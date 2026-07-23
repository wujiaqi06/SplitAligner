#!/usr/bin/env python3

import hashlib
import json
import pathlib
import sys


manifest_path = pathlib.Path(sys.argv[1])
state_path = pathlib.Path(sys.argv[2])
state_name = sys.argv[3]
output_path = pathlib.Path(sys.argv[4])

data = json.loads(manifest_path.read_text(encoding="utf-8"))
state_sha = hashlib.sha256(state_path.read_bytes()).hexdigest()
data["coordinate_state"]["filename"] = state_name
data["coordinate_state"]["sha256"] = state_sha
data["outputs"]["primitive_state"]["filename"] = state_name
data["outputs"]["primitive_state"]["sha256"] = state_sha
output_path.write_text(
    json.dumps(data, ensure_ascii=False, indent=3, sort_keys=True) + "\n",
    encoding="utf-8",
)
