#!/bin/sh
set -eu
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 PATH_TO_SAF NEW_OUTPUT_DIRECTORY" >&2
    exit 2
fi
example_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "$1" --reference "$example_dir/reference.nwk" --genes "$example_dir/genes.tsv" \
    --scope full --role free --out "$2"
