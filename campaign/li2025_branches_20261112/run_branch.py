#!/usr/bin/env python3
"""Run a fresh smoke and production pair through the frozen Li runner."""
import argparse
import json
from pathlib import Path
import frozen_run_cp2k as runner
from validate_stage import validate

BASE = Path(__file__).resolve().parent
BRANCHES = ('r_eq_metad', 'p_iface_v2', 'relax_v2')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=BRANCHES, required=True)
    parser.add_argument('--ntasks', type=int, default=64)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--workdir-root', type=Path, default=BASE)
    args = parser.parse_args()
    manifest = json.loads((BASE / 'manifest.json').read_text())
    for name, expected in manifest['dependencies'].items():
        if runner.sha256(BASE / name) != expected:
            raise RuntimeError('Dependency hash mismatch: ' + name)
    runner.BASE = BASE
    runner.MANIFEST_CASES = {'manifest.json': tuple(r['name'] for r in manifest['cases'])}
    root = args.workdir_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    cases = (args.case + '_smoke', args.case)
    if not args.dry_run and any((root / name).exists() for name in cases):
        raise RuntimeError('Existing branch workdirs; refuse duplicate workflow')
    for case in cases:
        work = root / case
        code = runner.run(case, str(work), args.ntasks, args.dry_run,
                          str(BASE / 'manifest.json'), None, str(root))
        if code:
            raise SystemExit(code)
        if not args.dry_run:
            record = next(r for r in manifest['cases'] if r['name'] == case)
            validate(work, record)
            print('STAGE_VALIDATED ' + case, flush=True)


if __name__ == '__main__':
    main()
