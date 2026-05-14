#!/usr/bin/env python3
"""Compare SplitAligner packaged outputs against benchmark unrooted oracle."""

from __future__ import annotations

import argparse
import csv
import math
import tempfile
from pathlib import Path

TOL = 1e-8


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_tsv_with_fields(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return reader.fieldnames or [], rows


def read_matrix(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields:
            raise ValueError(f"Empty matrix header: {path}")
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


def status(value: str) -> str:
    if value is None:
        return "other"
    text = value.strip()
    if not text:
        return "other"
    if text.startswith("NA"):
        return text
    if parse_finite_numeric(text) is not None:
        return "observed"
    return "other"


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_crosswalk(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    mapping: dict[str, str] = {}
    for row in rows:
        b = row.get("splitaligner_branch") or row.get("SplitAligner") or ""
        unrooted = row.get("benchmark_unrooted_branch") or row.get("Benchmark_unrooted") or ""
        if not b or b == "-":
            continue
        require(b not in mapping, f"Duplicate SplitAligner branch in crosswalk: {b}")
        require(unrooted and unrooted != "-", f"Missing benchmark_unrooted branch for {b} in crosswalk")
        mapping[b] = unrooted
    return mapping


def sorted_keys(rows: dict[str, dict[str, str]]) -> list[str]:
    return sorted(rows.keys())


def summarize_primitive_integrity(
    oracle_cols: list[str],
    oracle_rows: dict[str, dict[str, str]],
    split_cols: list[str],
    split_rows: dict[str, dict[str, str]],
    branch_map: dict[str, str],
) -> dict[str, object]:
    oracle_col_set = set(oracle_cols)
    split_col_set = set(split_cols)
    crosswalk_split_cols = sorted(branch_map.keys())
    crosswalk_oracle_branches = sorted(set(branch_map.values()))

    unmapped_split_cols = sorted([col for col in split_cols if col not in branch_map])
    mapped_split_cols = sorted([col for col in split_cols if col in branch_map])
    missing_splitaligner_branches = sorted(set(crosswalk_split_cols) - split_col_set)
    missing_oracle_branches = sorted(
        {branch_map[col] for col in crosswalk_split_cols if branch_map[col] not in oracle_col_set}
    )
    uncovered_oracle_branches = sorted(oracle_col_set - set(crosswalk_oracle_branches))
    comparable_split_cols = [
        col for col in mapped_split_cols if branch_map[col] in oracle_col_set
    ]

    oracle_genes = set(oracle_rows)
    split_genes = set(split_rows)
    shared_genes = sorted(oracle_genes & split_genes)
    oracle_only_genes = sorted(oracle_genes - split_genes)
    splitaligner_only_genes = sorted(split_genes - oracle_genes)

    expected_compared_cells = len(shared_genes) * len(crosswalk_split_cols)

    return {
        "oracle_col_set": oracle_col_set,
        "split_col_set": split_col_set,
        "crosswalk_split_cols": crosswalk_split_cols,
        "crosswalk_oracle_branches": crosswalk_oracle_branches,
        "unmapped_split_cols": unmapped_split_cols,
        "mapped_split_cols": mapped_split_cols,
        "missing_splitaligner_branches": missing_splitaligner_branches,
        "missing_oracle_branches": missing_oracle_branches,
        "uncovered_oracle_branches": uncovered_oracle_branches,
        "comparable_split_cols": comparable_split_cols,
        "shared_genes": shared_genes,
        "oracle_only_genes": oracle_only_genes,
        "splitaligner_only_genes": splitaligner_only_genes,
        "expected_compared_cells": expected_compared_cells,
    }


def map_composite_column(column: str, branch_map: dict[str, str]) -> tuple[str, ...]:
    members = []
    for part in column.split("|"):
        require(part in branch_map, f"Composite column {column} contains unmapped SplitAligner branch {part}")
        members.append(branch_map[part])
    return tuple(sorted(members))


def compare_primitive_cells(
    oracle_dir: Path,
    perl_dir: Path,
    branch_map: dict[str, str],
    audit_dir: Path,
) -> dict[str, object]:
    oracle_cols, oracle = read_matrix(oracle_dir / "benchmark.table.txt")
    split_cols, split = read_matrix(perl_dir / "benchmark_fix.fix.na_classified.txt")
    integrity = summarize_primitive_integrity(oracle_cols, oracle, split_cols, split, branch_map)
    rows: list[dict[str, str]] = []

    for gene_id in integrity["shared_genes"]:
        oracle_row = oracle[gene_id]
        split_row = split[gene_id]
        for b_col in integrity["comparable_split_cols"]:
            branch = branch_map[b_col]
            expected_value = oracle_row[branch]
            observed_value = split_row[b_col]
            expected_status = status(expected_value)
            observed_status = status(observed_value)
            exp_num = parse_finite_numeric(expected_value)
            obs_num = parse_finite_numeric(observed_value)

            if exp_num is not None and obs_num is not None:
                exact = abs(exp_num - obs_num) < TOL
            else:
                exact = expected_value == observed_value

            if exact:
                artifact = "match"
            elif expected_status in {"NA_fuse", "NA_struct"} and observed_status == "observed":
                artifact = f"{expected_status}_numeric_leak"
            elif expected_status == "observed" and observed_status in {"NA_fuse", "NA_struct", "NA"}:
                artifact = "over_mask"
            else:
                artifact = "mismatch"

            rows.append(
                {
                    "gene_id": gene_id,
                    "splitaligner_branch_id": b_col,
                    "benchmark_unrooted_branch_id": branch,
                    "expected_value": expected_value,
                    "expected_status": expected_status,
                    "splitaligner_value": observed_value,
                    "splitaligner_status": observed_status,
                    "artifact_type": artifact,
                }
            )

    actual_compared_cells = len(rows)
    expected_cell_count = integrity["expected_compared_cells"]

    fieldnames = list(rows[0].keys()) if rows else [
        "gene_id",
        "splitaligner_branch_id",
        "benchmark_unrooted_branch_id",
        "expected_value",
        "expected_status",
        "splitaligner_value",
        "splitaligner_status",
        "artifact_type",
    ]
    write_tsv(audit_dir / "expected_vs_splitaligner.full.tsv", rows, fieldnames)

    diff_rows = [row for row in rows if row["artifact_type"] != "match"]
    write_tsv(audit_dir / "expected_vs_splitaligner.diff.tsv", diff_rows, fieldnames)
    write_tsv(audit_dir / "unexpected_mismatches.tsv", diff_rows, fieldnames)

    primitive_pass = (
        len(diff_rows) == 0
        and len(integrity["unmapped_split_cols"]) == 0
        and len(integrity["missing_splitaligner_branches"]) == 0
        and len(integrity["missing_oracle_branches"]) == 0
        and len(integrity["uncovered_oracle_branches"]) == 0
        and len(integrity["oracle_only_genes"]) == 0
        and len(integrity["splitaligner_only_genes"]) == 0
        and actual_compared_cells == expected_cell_count
    )

    return {
        "pass": primitive_pass,
        "expected_compared_cells": expected_cell_count,
        "actual_compared_cells": actual_compared_cells,
        "unexpected_mismatches": len(diff_rows),
        "shared_genes": len(integrity["shared_genes"]),
        "mapped_branches": len(integrity["crosswalk_split_cols"]),
        "unmapped_split_cols": integrity["unmapped_split_cols"],
        "missing_splitaligner_branches": integrity["missing_splitaligner_branches"],
        "missing_oracle_branches": integrity["missing_oracle_branches"],
        "uncovered_oracle_branches": integrity["uncovered_oracle_branches"],
        "oracle_only_genes": integrity["oracle_only_genes"],
        "splitaligner_only_genes": integrity["splitaligner_only_genes"],
        "full_table_path": str(audit_dir / "expected_vs_splitaligner.full.tsv"),
        "diff_table_path": str(audit_dir / "expected_vs_splitaligner.diff.tsv"),
    }


def compare_fusion_groups(
    oracle_dir: Path,
    perl_dir: Path,
    branch_map: dict[str, str],
    audit_dir: Path,
) -> dict[str, object]:
    oracle_path = oracle_dir / "oracle_fusion_groups.tsv"
    require(oracle_path.exists(), f"Missing oracle fusion group table: {oracle_path}")
    oracle_fields, oracle_groups = read_tsv_with_fields(oracle_path)

    fused_cols, fused_matrix = read_matrix(perl_dir / "benchmark.matrix_with_fuse.txt")
    composite_cols = [col for col in fused_cols if "|" in col]
    composite_members = {col: map_composite_column(col, branch_map) for col in composite_cols}

    member_col = next((key for key in oracle_fields if key.endswith("_members")), None)
    require(member_col is not None, "Could not find *_members column in oracle_fusion_groups.tsv")

    expected_by_gene: dict[str, list[dict[str, str]]] = {}
    expected_member_sets_by_gene: dict[str, set[tuple[str, ...]]] = {}
    for row in oracle_groups:
        expected_by_gene.setdefault(row["gene_id"], []).append(row)
        expected_member_sets_by_gene.setdefault(row["gene_id"], set()).add(
            tuple(sorted(row[member_col].split("|")))
        )

    full_rows: list[dict[str, str]] = []
    diff_rows: list[dict[str, str]] = []

    for gene_id, rows in expected_by_gene.items():
        require(gene_id in fused_matrix, f"Gene {gene_id} missing from benchmark.matrix_with_fuse.txt")
        fused_row = fused_matrix[gene_id]
        for row in rows:
            members = tuple(sorted(row[member_col].split("|")))
            expected_len = parse_finite_numeric(row["expected_fused_length"])
            candidate_cols = [col for col, mapped in composite_members.items() if mapped == members]
            numeric_candidates = [
                (col, parse_finite_numeric(fused_row[col]))
                for col in candidate_cols
                if parse_finite_numeric(fused_row[col]) is not None
            ]

            if not candidate_cols:
                artifact = "missing_fused_column"
                observed_col = ""
                observed_value = ""
            elif not numeric_candidates:
                artifact = "fused_column_not_numeric"
                observed_col = candidate_cols[0]
                observed_value = fused_row[observed_col]
            else:
                matched = None
                for col, val in numeric_candidates:
                    if expected_len is not None and val is not None and abs(val - expected_len) < TOL:
                        matched = (col, val)
                        break
                if matched is None:
                    observed_col = numeric_candidates[0][0]
                    observed_value = fused_row[observed_col]
                    artifact = "fused_length_mismatch"
                else:
                    observed_col = matched[0]
                    observed_value = fused_row[observed_col]
                    artifact = "match"

            out_row = {
                "gene_id": gene_id,
                "step_id": row["step_id"],
                "merge_group_id": row["merge_group_id"],
                member_col: "|".join(members),
                "group_size": row["group_size"],
                "expected_fused_length": row["expected_fused_length"],
                "splitaligner_composite_column": observed_col,
                "splitaligner_value": observed_value,
                "artifact_type": artifact,
            }
            full_rows.append(out_row)
            if artifact != "match":
                diff_rows.append(out_row)

    for gene_id, fused_row in fused_matrix.items():
        expected_member_sets = expected_member_sets_by_gene.get(gene_id, set())
        for col in composite_cols:
            numeric_value = parse_finite_numeric(fused_row[col])
            if numeric_value is None:
                continue
            mapped_members = composite_members[col]
            primitive_values = [parse_finite_numeric(fused_row.get(member, "")) for member in col.split("|")]
            if all(value is not None for value in primitive_values):
                synthetic_sum = sum(value for value in primitive_values if value is not None)
                artifact = (
                    "synthetic_composite_match"
                    if abs(numeric_value - synthetic_sum) < TOL
                    else "synthetic_composite_length_mismatch"
                )
                out_row = {
                    "gene_id": gene_id,
                    "step_id": "",
                    "merge_group_id": "",
                    member_col: "|".join(mapped_members),
                    "group_size": str(len(mapped_members)),
                    "expected_fused_length": "",
                    "splitaligner_composite_column": col,
                    "splitaligner_value": fused_row[col],
                    "artifact_type": artifact,
                }
                full_rows.append(out_row)
                if artifact != "synthetic_composite_match":
                    diff_rows.append(out_row)
                continue
            if mapped_members in expected_member_sets:
                continue
            out_row = {
                "gene_id": gene_id,
                "step_id": "",
                "merge_group_id": "",
                member_col: "|".join(mapped_members),
                "group_size": str(len(mapped_members)),
                "expected_fused_length": "",
                "splitaligner_composite_column": col,
                "splitaligner_value": fused_row[col],
                "artifact_type": "unexpected_active_fused_column",
            }
            full_rows.append(out_row)
            diff_rows.append(out_row)

    fieldnames = [
        "gene_id",
        "step_id",
        "merge_group_id",
        member_col,
        "group_size",
        "expected_fused_length",
        "splitaligner_composite_column",
        "splitaligner_value",
        "artifact_type",
    ]
    write_tsv(audit_dir / "fusion_groups.full.tsv", full_rows, fieldnames)
    write_tsv(audit_dir / "fusion_groups.diff.tsv", diff_rows, fieldnames)
    return {
        "pass": len(diff_rows) == 0,
        "rows_checked": len(full_rows),
        "unexpected_mismatches": len(diff_rows),
        "full_table_path": str(audit_dir / "fusion_groups.full.tsv"),
        "diff_table_path": str(audit_dir / "fusion_groups.diff.tsv"),
    }


def run_self_tests() -> None:
    require(status("NaN") == "other", "NaN must not be classified as observed")
    require(status("Inf") == "other", "Inf must not be classified as observed")
    require(status("-Inf") == "other", "-Inf must not be classified as observed")
    require(status("1.0") == "observed", "Finite numeric values must be classified as observed")
    require(status("NA_fuse") == "NA_fuse", "NA_fuse state detection must be preserved")

    with tempfile.TemporaryDirectory() as tmpdir_name:
        tmpdir = Path(tmpdir_name)
        oracle_dir = tmpdir / "oracle"
        perl_dir = tmpdir / "perl"
        audit_dir = tmpdir / "audit"
        oracle_dir.mkdir()
        perl_dir.mkdir()
        audit_dir.mkdir()

        (oracle_dir / "benchmark.table.txt").write_text(
            "gene_id\tN_12\nmid_step0\t0.5\n"
        )
        (perl_dir / "benchmark_fix.fix.na_classified.txt").write_text(
            "gene_id\tB14\tB99\nmid_step0\t0.5\t0.7\n"
        )
        branch_map = {"B14": "N_12"}
        primitive = compare_primitive_cells(oracle_dir, perl_dir, branch_map, audit_dir)
        require(not primitive["pass"], "Missing mapping branch must FAIL, not PASS")
        require(
            primitive["unmapped_split_cols"] == ["B99"],
            "Missing mapping branch should be reported explicitly",
        )

        (oracle_dir / "benchmark.table.txt").write_text(
            "gene_id\tU1\tU2\ng1\t0.1\t0.2\n"
        )
        (perl_dir / "benchmark_fix.fix.na_classified.txt").write_text(
            "gene_id\tB1\ng1\t0.1\n"
        )
        branch_map_missing = {"B1": "U1", "B2": "U2"}
        primitive = compare_primitive_cells(oracle_dir, perl_dir, branch_map_missing, audit_dir)
        require(
            not primitive["pass"],
            "Missing expected SplitAligner branch column must FAIL",
        )
        require(
            primitive["missing_splitaligner_branches"] == ["B2"],
            "Missing expected SplitAligner branch should be reported explicitly",
        )

        (oracle_dir / "oracle_fusion_groups.tsv").write_text(
            "gene_id\tstep_id\tmerge_group_id\texpected_fused_length\tgroup_size\tbenchmark_unrooted_members\n"
        )
        (perl_dir / "benchmark.matrix_with_fuse.txt").write_text(
            "gene_id\tB14\nmid_step0\t0.5\n"
        )
        (perl_dir / "benchmark_fix.fix.na_classified.txt").write_text(
            "gene_id\tB14\nmid_step0\t0.5\n"
        )
        fusion = compare_fusion_groups(oracle_dir, perl_dir, branch_map, audit_dir)
        require(fusion["pass"], "Empty oracle fusion table with no active numeric composite columns should PASS")
        require(fusion["rows_checked"] == 0, "Empty oracle fusion table should report zero checked fusion rows")

        (perl_dir / "benchmark.matrix_with_fuse.txt").write_text(
            "gene_id\tB14\tB99\tB14|B99\nmid_step0\t0.4\t0.5\t1.0\n"
        )
        (perl_dir / "benchmark_fix.fix.na_classified.txt").write_text(
            "gene_id\tB14\tB99\nmid_step0\t0.4\t0.5\n"
        )
        branch_map_extra = {"B14": "N_12", "B99": "N_99"}
        fusion = compare_fusion_groups(oracle_dir, perl_dir, branch_map_extra, audit_dir)
        require(
            not fusion["pass"],
            "Synthetic composite length mismatch should FAIL",
        )

        (perl_dir / "benchmark.matrix_with_fuse.txt").write_text(
            "gene_id\tB14\tB99\tB14|B99\nmid_step0\tNA\tNA\t0.9\n"
        )
        (perl_dir / "benchmark_fix.fix.na_classified.txt").write_text(
            "gene_id\tB14\tB99\nmid_step0\tNA\tNA\n"
        )
        fusion = compare_fusion_groups(oracle_dir, perl_dir, branch_map_extra, audit_dir)
        require(
            not fusion["pass"],
            "Unexpected active fused column should FAIL when primitive members are nonnumeric",
        )

        (oracle_dir / "oracle_fusion_groups.tsv").write_text(
            "gene_id\tstep_id\tmerge_group_id\texpected_fused_length\tgroup_size\tbenchmark_unrooted_members\n"
            "mid_step0\t0\tMG01\t0.6\t2\tU1|U2\n"
        )
        (perl_dir / "benchmark.matrix_with_fuse.txt").write_text(
            "gene_id\tB1|B2\tB1|B3\tB1\tB2\tB3\nmid_step0\t0.6\t0.4\t0.1\t0.5\t0.3\n"
        )
        (perl_dir / "benchmark_fix.fix.na_classified.txt").write_text(
            "gene_id\tB1\tB2\tB3\nmid_step0\t0.1\t0.2\t0.3\n"
        )
        branch_map_fusion = {"B1": "U1", "B2": "U2", "B3": "U3"}
        fusion = compare_fusion_groups(oracle_dir, perl_dir, branch_map_fusion, audit_dir)
        require(
            fusion["pass"],
            "Synthetic composite values with correct primitive sums should PASS",
        )

        (oracle_dir / "oracle_fusion_groups.tsv").write_text(
            "gene_id\tstep_id\tmerge_group_id\texpected_fused_length\tgroup_size\tbenchmark_unrooted_members\n"
        )
        (perl_dir / "benchmark.matrix_with_fuse.txt").write_text(
            "gene_id\tB1\tB3\tB1|B3\nmid_step0\t0.1\t0.3\t0.4\n"
        )
        (perl_dir / "benchmark_fix.fix.na_classified.txt").write_text(
            "gene_id\tB1\tB3\nmid_step0\t0.1\t0.3\n"
        )
        branch_map_synthetic = {"B1": "U1", "B3": "U3"}
        fusion = compare_fusion_groups(oracle_dir, perl_dir, branch_map_synthetic, audit_dir)
        require(
            fusion["pass"],
            "Synthetic composite with exact primitive sum should PASS",
        )

        (perl_dir / "benchmark.matrix_with_fuse.txt").write_text(
            "gene_id\tB1\tB3\tB1|B3\nmid_step0\t0.1\t0.3\t0.9\n"
        )
        fusion = compare_fusion_groups(oracle_dir, perl_dir, branch_map_synthetic, audit_dir)
        require(
            not fusion["pass"],
            "Synthetic composite with incorrect primitive sum should FAIL",
        )

        (perl_dir / "benchmark.matrix_with_fuse.txt").write_text(
            "gene_id\tB1\tB2\tB1|B2\nmid_step0\tNA\tNA\t0.6\n"
        )
        (perl_dir / "benchmark_fix.fix.na_classified.txt").write_text(
            "gene_id\tB1\tB2\nmid_step0\tNA_fuse\tNA_fuse\n"
        )
        (oracle_dir / "oracle_fusion_groups.tsv").write_text(
            "gene_id\tstep_id\tmerge_group_id\texpected_fused_length\tgroup_size\tbenchmark_unrooted_members\n"
            "mid_step0\t0\tMG01\t0.6\t2\tU1|U2\n"
        )
        branch_map_active = {"B1": "U1", "B2": "U2"}
        fusion = compare_fusion_groups(oracle_dir, perl_dir, branch_map_active, audit_dir)
        require(
            fusion["pass"],
            "Oracle-expected active fused coordinate with nonnumeric primitive members should PASS",
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--oracle-dir", help="benchmark_unrooted directory")
    parser.add_argument("--perl-dir", help="splitaligner_perl directory")
    parser.add_argument("--crosswalk", help="branch label map TSV")
    parser.add_argument("--audit-dir", help="output audit directory")
    parser.add_argument("--self-test", action="store_true", help="run built-in regression checks and exit")
    args = parser.parse_args()

    if args.self_test:
        run_self_tests()
        return

    require(args.oracle_dir is not None, "--oracle-dir is required unless --self-test is used")
    require(args.perl_dir is not None, "--perl-dir is required unless --self-test is used")
    require(args.crosswalk is not None, "--crosswalk is required unless --self-test is used")
    require(args.audit_dir is not None, "--audit-dir is required unless --self-test is used")

    oracle_dir = Path(args.oracle_dir)
    perl_dir = Path(args.perl_dir)
    crosswalk = Path(args.crosswalk)
    audit_dir = Path(args.audit_dir)
    audit_dir.mkdir(parents=True, exist_ok=True)

    branch_map = load_crosswalk(crosswalk)
    primitive = compare_primitive_cells(oracle_dir, perl_dir, branch_map, audit_dir)
    fusion = compare_fusion_groups(oracle_dir, perl_dir, branch_map, audit_dir)

    status_flag = "PASS" if primitive["pass"] and fusion["pass"] else "FAIL"
    summary_lines = [
        f"Status: {status_flag}",
        f"Primitive cell audit: {'PASS' if primitive['pass'] else 'FAIL'}",
        f"Fusion group audit: {'PASS' if fusion['pass'] else 'FAIL'}",
        f"Shared genes: {primitive['shared_genes']}",
        f"Mapped branches: {primitive['mapped_branches']}",
        f"Expected compared cells: {primitive['expected_compared_cells']}",
        f"Actual compared cells: {primitive['actual_compared_cells']}",
        f"Unexpected mismatches: {primitive['unexpected_mismatches']}",
        f"Unmapped SplitAligner branches: {','.join(primitive['unmapped_split_cols']) if primitive['unmapped_split_cols'] else '(none)'}",
        f"Missing expected SplitAligner branches: {','.join(primitive['missing_splitaligner_branches']) if primitive['missing_splitaligner_branches'] else '(none)'}",
        f"Missing oracle branches: {','.join(primitive['missing_oracle_branches']) if primitive['missing_oracle_branches'] else '(none)'}",
        f"Uncovered oracle branches: {','.join(primitive['uncovered_oracle_branches']) if primitive['uncovered_oracle_branches'] else '(none)'}",
        f"Oracle-only genes: {','.join(primitive['oracle_only_genes']) if primitive['oracle_only_genes'] else '(none)'}",
        f"SplitAligner-only genes: {','.join(primitive['splitaligner_only_genes']) if primitive['splitaligner_only_genes'] else '(none)'}",
        f"Full table path: {primitive['full_table_path']}",
        f"Diff table path: {primitive['diff_table_path']}",
        f"Fusion-group rows checked: {fusion['rows_checked']}",
        f"Fusion-group mismatches: {fusion['unexpected_mismatches']}",
        f"Fusion full table path: {fusion['full_table_path']}",
        f"Fusion diff table path: {fusion['diff_table_path']}",
        f"Crosswalk mapped branches: {len(branch_map)}",
    ]
    (audit_dir / "pass_fail_summary.txt").write_text("\n".join(summary_lines) + "\n")

    if status_flag != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
