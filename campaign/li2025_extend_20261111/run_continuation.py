#!/usr/bin/env python3
"""Use the frozen, hash-checked segment runner in this isolated directory."""
import hashlib
import json
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
manifest = json.loads((HERE / 'continuation_manifest.json').read_text())
for name, digest in manifest['dependencies'].items():
    if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
        raise SystemExit('Frozen dependency changed: ' + name)
import frozen_run_cp2k as runner
runner.BASE = HERE
runner.MANIFEST_CASES = {'continuation_manifest.json': ('interface_r1_seg3',)}
sys.exit(runner.main())
