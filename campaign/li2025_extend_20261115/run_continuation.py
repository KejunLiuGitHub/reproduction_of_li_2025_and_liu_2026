#!/usr/bin/env python3
"""Run the byte-frozen runner, then require separate complete technical CHECKS."""
import argparse
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent


def main():
    manifest = json.loads((HERE / 'continuation_manifest.json').read_text())
    import hashlib
    for name, digest in manifest['dependencies'].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Frozen dependency changed: ' + name)
    import prepare
    import frozen_run_cp2k as runner
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', default='continuation_manifest.json')
    parser.add_argument('--case', default=prepare.CASE, choices=[prepare.CASE])
    parser.add_argument('--parent-dir', required=True, type=Path)
    parser.add_argument('--workdir', type=Path)
    parser.add_argument('--workdir-root', type=Path, default=HERE)
    parser.add_argument('--ntasks', type=int, default=int(os.environ.get('SLURM_NTASKS', '64')))
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if args.manifest != 'continuation_manifest.json':
        raise RuntimeError('Only the frozen manifest is allowed')
    record = manifest['cases'][0]
    work = args.workdir or HERE / args.case
    runner.BASE = HERE
    runner.MANIFEST_CASES = {'continuation_manifest.json': (prepare.CASE,)}
    # Refuse stale/unsafe destinations before any writes. Parent data are read only.
    work = runner.resolve_workdir(work, args.workdir_root)
    runner.refuse_stale(work, args.case, record)
    for name, digest in record['parent_evidence']['source_files_sha256'].items():
        if prepare.sha(args.parent_dir / name) != digest:
            raise RuntimeError('Validated parent changed: ' + name)
    result = runner.run(args.case, str(work), args.ntasks, args.dry_run,
                        args.manifest, args.parent_dir, args.workdir_root)
    if args.dry_run:
        return result
    # Preserve the original runner record even if this stronger gate rejects it.
    expected_inputs = {entry['target']: entry['sha256'] for entry in record['stage']
                       if entry['role'] not in ('parent_hills', 'parent_colvar')}
    try:
        prepare.validate_segment(work, args.case, 40000, 50000,
            os.environ.get('SLURM_JOB_ID'), HERE / 'inputs/interface_r1.xyz',
            args.parent_dir, expected_inputs, {}, work / 'CHECKS.json',
            record['technical_policy']['plumed_time_fs_per_step'])
    except (ValueError, OSError, KeyError) as error:
        print(str(error), file=sys.stderr)
        return result or 9
    return result


if __name__ == '__main__':
    sys.exit(main())
