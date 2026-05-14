#!/usr/bin/env python3
"""Tiny invariant checks for endpoint-collapsed internal branches."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_matrix(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        require(fields, f"Empty matrix header: {path}")
        gene_col = fields[0]
        rows = {row[gene_col]: row for row in reader}
        return fields[1:], rows


def parse_finite_numeric(value: str) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if not text or text.startswith("NA"):
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def load_crosswalk(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    rows = read_tsv(path)
    split_to_oracle: dict[str, str] = {}
    oracle_to_split: dict[str, str] = {}
    for row in rows:
        split = row.get("splitaligner_branch", "").strip()
        oracle = row.get("benchmark_unrooted_branch", "").strip()
        if not split or split == "-":
            continue
        require(oracle and oracle != "-", f"Missing benchmark_unrooted branch for {split}")
        split_to_oracle[split] = oracle
        oracle_to_split[oracle] = split
    return split_to_oracle, oracle_to_split


def classify_member_type(member: str) -> str:
    return "internal" if member.startswith("N_") else "terminal"


def normalize_member_set(text: str) -> tuple[str, ...]:
    return tuple(sorted([part for part in text.split("|") if part]))


def map_composite_column(column: str, split_to_oracle: dict[str, str]) -> tuple[str, ...]:
    members = []
    for part in column.split("|"):
        require(part in split_to_oracle, f"Composite column {column} contains unmapped branch {part}")
        members.append(split_to_oracle[part])
    return tuple(sorted(members))


def find_composite_column(
    fused_cols: list[str],
    fused_row: dict[str, str],
    split_to_oracle: dict[str, str],
    oracle_members: tuple[str, ...],
) -> tuple[str | None, float | None]:
    for col in fused_cols:
        if "|" not in col:
            continue
        if map_composite_column(col, split_to_oracle) != oracle_members:
            continue
        value = parse_finite_numeric(fused_row[col])
        if value is not None:
            return col, value
    return None, None


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scenario-root",
        required=True,
        help="Scenario root containing branch_label_map.tsv, benchmark_unrooted/, and splitaligner_perl/",
    )
    parser.add_argument(
        "--audit-dir",
        help="Optional directory for invariant outputs; defaults to <scenario-root>/../../audit/<scenario>",
    )
    args = parser.parse_args()

    scenario_root = Path(args.scenario_root)
    scenario_name = scenario_root.name
    benchmark_dir = scenario_root.parents[1]
    audit_dir = (
        Path(args.audit_dir)
        if args.audit_dir
        else benchmark_dir / "audit" / scenario_name
    )
    audit_dir.mkdir(parents=True, exist_ok=True)

    benchmark_unrooted = scenario_root / "benchmark_unrooted"
    splitaligner_perl = scenario_root / "splitaligner_perl"
    split_to_oracle, oracle_to_split = load_crosswalk(scenario_root / "branch_label_map.tsv")

    _, classified = read_matrix(splitaligner_perl / "benchmark_fix.fix.na_classified.txt")
    fused_cols, fused_matrix = read_matrix(splitaligner_perl / "benchmark.matrix_with_fuse.txt")
    oracle_groups = read_tsv(benchmark_unrooted / "oracle_fusion_groups.tsv")

    case_rows: list[dict[str, str]] = []
    warning_rows: list[dict[str, str]] = []

    # Case 1: endpoint-collapsed internal singleton must not remain numeric.
    singleton_gene = "main_step3"
    singleton_internal = "N_15"
    singleton_split = oracle_to_split[singleton_internal]
    singleton_value = classified[singleton_gene][singleton_split]
    singleton_pass = parse_finite_numeric(singleton_value) is None
    case_rows.append(
        {
            "case_id": "case1_internal_1k_singleton",
            "gene_id": singleton_gene,
            "oracle_members": singleton_internal,
            "splitaligner_columns": singleton_split,
            "expected": "non-numeric primitive state",
            "observed": singleton_value,
            "status": "PASS" if singleton_pass else "FAIL",
            "note": "Endpoint-collapsed internal singleton should not be emitted as primitive numeric branch.",
        }
    )

    # Case 2: endpoint-collapsed internal fused with terminal must emit fused pattern.
    fused_gene = "main_step2"
    fused_members = normalize_member_set("N_15|t4")
    fused_split_members = [oracle_to_split[m] for m in fused_members]
    fused_col, fused_len = find_composite_column(
        fused_cols, fused_matrix[fused_gene], split_to_oracle, fused_members
    )
    primitive_status_ok = all(
        classified[fused_gene][split_branch] == "NA_fuse"
        for split_branch in fused_split_members
    )
    fused_numeric_ok = fused_col is not None and fused_len is not None
    case_rows.append(
        {
            "case_id": "case2_internal_1k_fused_with_terminal",
            "gene_id": fused_gene,
            "oracle_members": "|".join(fused_members),
            "splitaligner_columns": "|".join(fused_split_members),
            "expected": "numeric fused column plus primitive NA_fuse",
            "observed": (
                f"composite={fused_col or '(missing)'};"
                f" composite_value={fused_matrix[fused_gene].get(fused_col, '') if fused_col else ''};"
                f" primitive_states={','.join(classified[fused_gene][b] for b in fused_split_members)}"
            ),
            "status": "PASS" if fused_numeric_ok and primitive_status_ok else "FAIL",
            "note": "Endpoint-collapsed internal fused with terminal should emit fused pattern and primitive NA_fuse states.",
        }
    )

    # Case 3: internal-only fusion groups are allowed as bookkeeping if documented.
    for row in oracle_groups:
        members = normalize_member_set(row["benchmark_unrooted_members"])
        if len(members) < 2:
            continue
        if any(classify_member_type(member) != "internal" for member in members):
            continue
        gene_id = row["gene_id"]
        composite_col, composite_val = find_composite_column(
            fused_cols, fused_matrix[gene_id], split_to_oracle, members
        )
        warning_rows.append(
            {
                "gene_id": gene_id,
                "step_id": row["step_id"],
                "merge_group_id": row["merge_group_id"],
                "benchmark_unrooted_members": "|".join(members),
                "splitaligner_composite_column": composite_col or "",
                "splitaligner_value": (
                    fused_matrix[gene_id].get(composite_col, "") if composite_col else ""
                ),
                "status": "PASS" if composite_col and composite_val is not None else "FAIL",
                "note": (
                    "Internal-only fused group retained as fused-path bookkeeping under "
                    "standard species-axis construction."
                ),
            }
        )

    overall_pass = all(row["status"] == "PASS" for row in case_rows) and all(
        row["status"] == "PASS" for row in warning_rows
    )

    write_tsv(
        audit_dir / "endpoint_collapse_cases.tsv",
        case_rows,
        [
            "case_id",
            "gene_id",
            "oracle_members",
            "splitaligner_columns",
            "expected",
            "observed",
            "status",
            "note",
        ],
    )
    write_tsv(
        audit_dir / "endpoint_collapse_internal_only_groups.tsv",
        warning_rows,
        [
            "gene_id",
            "step_id",
            "merge_group_id",
            "benchmark_unrooted_members",
            "splitaligner_composite_column",
            "splitaligner_value",
            "status",
            "note",
        ],
    )

    summary_lines = [
        f"Status: {'PASS' if overall_pass else 'FAIL'}",
        f"Scenario: {scenario_name}",
        f"Case 1 singleton check: {case_rows[0]['status']}",
        f"Case 2 fused-with-terminal check: {case_rows[1]['status']}",
        f"Internal-only fused groups checked: {len(warning_rows)}",
        f"Internal-only fused groups status: {'PASS' if all(row['status'] == 'PASS' for row in warning_rows) else 'FAIL'}",
        f"Case table: {audit_dir / 'endpoint_collapse_cases.tsv'}",
        f"Internal-only table: {audit_dir / 'endpoint_collapse_internal_only_groups.tsv'}",
    ]
    (audit_dir / "endpoint_collapse_summary.txt").write_text("\n".join(summary_lines) + "\n")

    if not overall_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
