#!/usr/bin/env python3
"""Run exactly one selected candidate: fresh smoke, strict check, fresh 5 ps.

No submission action and no automatic loop over candidates.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import sys
import validate_stage as gate


def load_runner(base):
    path = Path(base) / 'frozen_run_cp2k.py'
    gate.demand(gate.sha(path) == gate.RUNNER_SHA, 'frozen runner changed')
    spec = importlib.util.spec_from_file_location('li_hb_frozen_runner', path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    return runner


def run_case(base, case, workdir_root, dry_run=False, ntasks=64):
    base, root = Path(base).resolve(), Path(workdir_root).resolve()
    gate.demand(ntasks == 64, 'this workflow is registered for 64 MPI tasks')
    manifest_path = base / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    gate.demand(manifest['schema'] == 'li2025_hbnetwork_workflow/v1', 'unregistered manifest schema')
    for name, digest in manifest['dependencies'].items():
        gate.demand(gate.sha(base / name) == digest, 'dependency hash mismatch: ' + name)
    records = {r['name']: r for r in manifest['cases']}
    gate.demand(case in manifest['candidate_names'] and case in records, 'explicit production candidate required')
    names = [case + '_smoke', case]
    gate.demand(all(n in records for n in names), 'missing smoke/production pair')
    smoke, production = [records[n] for n in names]
    gate.demand(smoke['steps'] == 20 and production['steps'] == 5000, 'stage duration changed')
    gate.demand(smoke['coordinate_file'] == production['coordinate_file']
                and smoke['seed'] == production['seed'] and smoke['source'] == production['source'],
                'smoke and production must start from the same frozen candidate/seed')
    runner = load_runner(base)
    runner.BASE = base
    runner.MANIFEST_CASES = {'manifest.json': tuple(records)}
    # Inspect both destinations before staging either one.
    for name in names:
        work = runner.resolve_workdir(root / name, root)
        runner.refuse_stale(work, name, records[name])
        gate.verify_input(base / 'inputs', records[name], base)
        if not dry_run:
            gate.demand(not work.exists(), 'existing stage directory; use a new isolated workdir root')
    root.mkdir(parents=True, exist_ok=True)
    for name in names:
        record = records[name]
        work = root / name
        code = runner.run(name, str(work), ntasks, dry_run, str(manifest_path), None, str(root))
        if code:
            return code
        if not dry_run:
            gate.validate(work, record, base, expected_job=os.environ.get('SLURM_JOB_ID'))
            print('STAGE_VALIDATED ' + name, flush=True)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--base', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--workdir-root', type=Path)
    parser.add_argument('--ntasks', type=int, default=int(os.environ.get('SLURM_NTASKS', '64')))
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    return run_case(args.base, args.case, args.workdir_root or args.base / 'runs', args.dry_run, args.ntasks)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, KeyError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(9)
