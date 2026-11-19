#!/usr/bin/env python3
"""Isolated adapter: exact frozen runner with an explicit six-case whitelist."""
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
manifest = json.loads((HERE/'pilot_manifest.json').read_text())
for filename, expected in manifest['dependencies'].items():
    if hashlib.sha256((HERE/filename).read_bytes()).hexdigest() != expected:
        raise SystemExit('Frozen dependency hash mismatch: '+filename)
import frozen_run_cp2k as runner
runner.BASE = HERE
runner.MANIFEST_CASES = {'pilot_manifest.json': (
    'p_interface', 'p_interface_smoke', 'p_bulk', 'p_bulk_smoke',
    'r_interface_eq', 'r_interface_eq_smoke')}
sys.exit(runner.main())
