"""Run the two-quartet example with native Windows or POSIX Python."""
from pathlib import Path
import subprocess
import sys


def main():
    if len(sys.argv) != 3:
        print("Usage: quickstart.py PATH_TO_SAF NEW_OUTPUT_DIRECTORY", file=sys.stderr)
        return 2
    example_dir = Path(__file__).resolve().parent
    return subprocess.call([
        sys.argv[1], "--reference", str(example_dir / "reference.nwk"),
        "--genes", str(example_dir / "genes.tsv"), "--scope", "full",
        "--role", "free", "--out", sys.argv[2],
    ])


if __name__ == "__main__":
    raise SystemExit(main())
