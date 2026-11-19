#!/usr/bin/env python3
"""Package explicitly selected, stable Li296 frames for unbiased retention runs.

No HPC interaction, no submission, no coordinate modification, no implicit cases.
A new explicit output directory is mandatory; dry-run also exercises staging.
"""
import argparse
from contextlib import redirect_stdout
import io
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sys
import validate_stage as gate

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REFERENCE = ROOT / 'campaign/li2025_branches_20261112/inputs/p_iface_v2.inp'
RUNNER = ROOT / 'campaign/li2025_branches_20261112/frozen_run_cp2k.py'
ANALYZER = ROOT / 'campaign/li2025_branches_20261112/frozen_analyze_event_unbiased.py'
ALIASES = {'sha256': 'trajectory_sha256', 'source_CHECKS': 'source_checks',
           'source_CHECKS_sha256': 'source_checks_sha256', 'source_job': 'source_job_id'}


def normalize_selection(data):
    gate.demand(('cases' in data) != ('candidates' in data), 'provide exactly one cases/candidates list')
    items = data.get('cases', data.get('candidates'))
    gate.demand(isinstance(items, list) and 1 <= len(items) <= 4, 'explicit selection must contain 1 to 4 candidates')
    result = []
    for raw in items:
        item = dict(raw)
        for old, new in ALIASES.items():
            if old in item:
                gate.demand(new not in item or item[new] == item[old], 'conflicting candidate alias: ' + old)
                item[new] = item.pop(old)
        for name in ('name', 'trajectory', 'trajectory_sha256', 'frame_step', 'source_checks',
                     'source_checks_sha256', 'source_job_id', 'seed', 'source_prepared_under_bias'):
            gate.demand(name in item, 'candidate missing ' + name)
        gate.demand(re.fullmatch(r'[a-z][a-z0-9_]{0,55}', item['name']) and not item['name'].endswith('_smoke'), 'invalid candidate name')
        gate.demand(type(item['frame_step']) is int and item['frame_step'] >= 0, 'invalid frame step')
        gate.demand(type(item['source_prepared_under_bias']) is bool, 'explicit source bias history required')
        for name in ('trajectory_sha256', 'source_checks_sha256'):
            gate.demand(re.fullmatch(r'[0-9a-f]{64}', item[name]), 'invalid declared SHA: ' + name)
        gate.demand(type(item['seed']) is int and 0 < item['seed'] < 2147483647, 'invalid integer seed')
        result.append(item)
    gate.demand(len({x['name'] for x in result}) == len(result), 'duplicate candidate name')
    gate.demand(len({x['seed'] for x in result}) == len(result), 'candidate seeds must be distinct and explicitly selected')
    return result


def source_path(value):
    path = Path(value)
    path = path if path.is_absolute() else ROOT / path
    resolved = path.resolve()
    gate.demand('methanol-sulfonization' not in resolved.parts, 'forbidden source directory')
    return resolved


def science_signature(text):
    """Ignore only allowed initialization/bias differences when checking sources."""
    output, stack, skip = [], [], 0
    for raw in text.splitlines():
        line = raw.split('#', 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        if fields[0].upper() in ('&EXT_RESTART', '&FREE_ENERGY') and not skip:
            skip = 1
            continue
        if skip:
            if fields[0].startswith('&'):
                skip += -1 if fields[0].upper() == '&END' else 1
            continue
        if fields[0] in ('PROJECT', 'SEED', 'COORD_FILE_NAME', 'STEPS'):
            continue
        output.append(' '.join(fields))
    gate.demand(skip == 0, 'truncated source bias/restart section')
    return output


def verify_source(candidate, cache):
    trajectory = source_path(candidate['trajectory'])
    checks_path = source_path(candidate['source_checks'])
    key = (str(trajectory), candidate['trajectory_sha256'], str(checks_path), candidate['source_checks_sha256'])
    if key in cache:
        verified = cache[key]
    else:
        gate.demand(trajectory.parent == checks_path.parent, 'source trajectory/CHECKS must share a stable run directory')
        gate.demand(gate.sha(trajectory) == candidate['trajectory_sha256'], 'source trajectory SHA mismatch')
        gate.demand(gate.sha(checks_path) == candidate['source_checks_sha256'], 'source CHECKS SHA mismatch')
        official = json.loads(checks_path.read_text())
        checks = official.get('checks', {})
        gate.demand(official.get('status') == 'PASS' and checks and all((v.get('pass') is True if isinstance(v, dict) else v is True) for v in checks.values()), 'source CHECKS not all PASS')
        declared = official.get('files_sha256', {}).get(trajectory.name, official.get('trajectory_sha256'))
        gate.demand(declared == candidate['trajectory_sha256'], 'source CHECKS does not describe this trajectory')
        for name, digest in official.get('files_sha256', {}).items():
            gate.demand(Path(name).name == name and gate.sha(trajectory.parent / name) == digest, 'source CHECKS file hash mismatch: ' + name)
        run_path = trajectory.parent / 'run_record.json'
        run = json.loads(run_path.read_text())
        gate.demand(run['status'] == 'COMPLETED' and run['returncode'] == 0 and run['program_ended'] and run['abort_count'] == 0, 'source run is not normally completed')
        gate.demand(str(official['job_id']) == str(run['job_id']) and official['case'] == run['case'], 'source CHECKS/run identity mismatch')
        input_hashes = run.get('input_sha256', {})
        gate.demand(run['case'] + '.inp' in input_hashes, 'source run record lacks input SHA')
        bias_prefix_bytes = {}
        for name, digest in input_hashes.items():
            gate.demand(Path(name).name == name, 'unsafe source input name')
            path = trajectory.parent / name
            if name.startswith(('HILLS_', 'COLVAR_')):
                # These staged inputs are append-only during a metaD source run.
                running, size = hashlib.sha256(), 0
                matched = running.hexdigest() == digest
                with path.open('rb') as handle:
                    for line in handle:
                        if matched:
                            break
                        running.update(line)
                        size += len(line)
                        matched = running.hexdigest() == digest
                gate.demand(matched, 'source bias input prefix SHA mismatch: ' + name)
                bias_prefix_bytes[name] = size
            else:
                gate.demand(gate.sha(path) == digest, 'source run input SHA mismatch: ' + name)
        begin, end = run['md_step_first'], run['md_step_last']
        gate.demand(type(begin) is int and type(end) is int and 1 <= begin <= end, 'source MD range invalid')
        length = end - begin + 1
        log_path = trajectory.parent / (run['case'] + '.out')
        log = gate.read_log(log_path)
        gate.demand(log['steps'] == list(range(begin, end + 1)), 'source log MD coverage mismatch')
        gate.demand(log['scf_converged'] == length + 1 and log['scf_failed'] == 0 and log['scf_pending'] == 0 and log['scf_bins'] == [2] + [1] * (length - 1), 'source SCF coverage/failure mismatch')
        gate.demand(log['program_ended_count'] == 1 and log['abort_count'] == 0 and run['scf']['converged_runs'] == length + 1 and run['scf']['not_converged_runs'] == 0, 'source termination/SCF record mismatch')
        frames = gate.read_xyz(trajectory)
        first = 0 if begin == 1 else ((begin + 9) // 10) * 10
        gate.demand([f['step'] for f in frames] == list(range(first, end + 1, 10)), 'source XYZ full frame sequence mismatch')
        gate.demand(run['trajectory']['frames'] == len(frames) and run['trajectory']['first_step'] == first and run['trajectory']['last_step'] == end, 'source runner/trajectory mismatch')
        gate.demand(all(f['time'] is not None and f['energy'] is not None and abs(f['time'] - f['step']) < 0.0011 for f in frames), 'source XYZ time/energy mismatch')
        energy_path = trajectory.parent / run['energy']['file']
        energy = gate.read_energy(energy_path)
        energy_first = 0 if begin == 1 else begin
        gate.demand([r[0] for r in energy] == list(range(energy_first, end + 1)) and all(abs(r[0] - r[1]) < 2e-6 for r in energy), 'source energy full step coverage mismatch')
        gate.demand(all(abs(energy[f['step'] - energy_first][4] - f['energy']) < 2e-8 for f in frames), 'source XYZ/energy mismatch')
        project = trajectory.name.removesuffix('-pos-1.xyz')
        gate.demand(project != trajectory.name, 'source trajectory must retain its CP2K filename')
        restart_path = trajectory.parent / (project + '-1.restart')
        restart = gate.read_restart(restart_path, end, length, frames[-1]['coords'])
        source_input = trajectory.parent / (run['case'] + '.inp')
        gate.demand(science_signature(source_input.read_text()) == science_signature(REFERENCE.read_text()), 'source scientific model differs from registered Li template')
        source_biased = bool(re.search(r'(?mi)^\s*&FREE_ENERGY\b', source_input.read_text())) or bool(official.get('biased')) or bool(official.get('prepared_under_bias'))
        evidence = {'trajectory': str(trajectory), 'trajectory_sha256': gate.sha(trajectory),
            'source_checks': str(checks_path), 'source_checks_sha256': gate.sha(checks_path),
            'source_job_id': str(run['job_id']), 'source_case': run['case'],
            'source_range': [first, end], 'source_frames': len(frames),
            'source_was_biased_or_prepared_under_bias': source_biased,
            'independent_technical_recheck': 'PASS', 'restart': restart,
            'verified_input_sha256': input_hashes, 'bias_input_prefix_bytes': bias_prefix_bytes,
            'verified_files_sha256': {p.name: gate.sha(p) for p in (trajectory, checks_path, run_path, log_path, energy_path, restart_path, source_input)}}
        verified = {'frames': {f['step']: f for f in frames}, 'evidence': evidence}
        cache[key] = verified
    gate.demand(str(candidate['source_job_id']) == verified['evidence']['source_job_id'], 'candidate source job mismatch')
    gate.demand(not verified['evidence']['source_was_biased_or_prepared_under_bias'] or candidate['source_prepared_under_bias'], 'biased source history cannot be erased')
    gate.demand(candidate['frame_step'] in verified['frames'], 'selected frame missing from source')
    if 'frame_block_sha256' in candidate:
        gate.demand(hashlib.sha256(verified['frames'][candidate['frame_step']]['raw'].encode()).hexdigest() == candidate['frame_block_sha256'], 'selected frame block SHA mismatch')
    return verified['frames'][candidate['frame_step']], verified['evidence']


def build(selection_path, output_dir, dry_run=False):
    selection_path, output_dir = Path(selection_path).resolve(), Path(output_dir).resolve()
    gate.demand(not output_dir.exists(), 'output directory already exists; refuse overwrite')
    gate.demand(gate.sha(REFERENCE) == gate.REFERENCE_SHA, 'registered science reference changed')
    gate.demand(gate.sha(RUNNER) == gate.RUNNER_SHA, 'registered frozen runner changed')
    gate.demand(gate.sha(ANALYZER) == gate.ANALYZER_SHA, 'registered frozen distance analyzer changed')
    selection = json.loads(selection_path.read_text())
    candidates = normalize_selection(selection)
    cache, selected = {}, []
    # Verify all sources before writing a package.
    for candidate in candidates:
        frame, evidence = verify_source(candidate, cache)
        selected.append((candidate, frame, evidence))
    inputs = output_dir / 'inputs'
    inputs.mkdir(parents=True)
    shutil.copy2(REFERENCE, inputs / 'reference_unbiased.inp')
    shutil.copy2(RUNNER, output_dir / 'frozen_run_cp2k.py')
    shutil.copy2(ANALYZER, output_dir / 'frozen_analyze_event_unbiased.py')
    for name in ('run.py', 'validate_stage.py', 'build_campaign.py'):
        shutil.copy2(HERE / name, output_dir / name)
    shutil.copy2(selection_path, output_dir / 'selection.json')
    records = []
    for candidate, frame, evidence in selected:
        branch = candidate['name']
        xyz = inputs / (branch + '.xyz')
        xyz.write_text(frame['raw'])
        source = dict(evidence, frame_step=candidate['frame_step'], candidate_xyz_sha256=gate.sha(xyz),
                      source_prepared_under_bias=candidate['source_prepared_under_bias'])
        for suffix, steps in [('_smoke', 20), ('', 5000)]:
            name = branch + suffix
            inp = inputs / (name + '.inp')
            inp.write_text(gate.make_input(REFERENCE.read_text(), name, xyz.name, candidate['seed'], steps))
            record = {'name': name, 'branch': branch, 'unbiased': True,
                'prepared_under_bias': candidate['source_prepared_under_bias'],
                'seed': candidate['seed'], 'steps': steps, 'expected_frames': steps // 10 + 1,
                'cp2k_project': name, 'coordinate_file': xyz.name,
                'box_class': 'interface', 'box_angstrom': [15., 15., 30.], 'source': source,
                'selection_metadata': candidate, 'fes': {'status': 'not_applicable', 'note': 'Fully unbiased conditional-retention stage; recorded prior bias history retained.'},
                'cp2k': {'steps_initial_segment': steps, 'trajectory_every_fs': 10},
                'stage': [{'role': 'new', 'source': p.name, 'target': p.name, 'sha256': gate.sha(p)} for p in (xyz, inp)],
                'outputs': {'output': name + '.out', 'trajectory': name + '-pos-1.xyz', 'energy': name + '-1.ener', 'restart': name + '-1.restart'},
                'human_review': 'PENDING', 'structural_verdict': None}
            gate.verify_input(inputs, record, output_dir)
            records.append(record)
        script = ('#!/bin/bash\n#SBATCH --job-name=li25-' + branch + '\n#SBATCH --partition=compute\n#SBATCH --nodes=1\n#SBATCH --ntasks=64\n#SBATCH --mem=32G\n#SBATCH --time=48:00:00\nset -euo pipefail\n'
            'BASE="${LI_WORKFLOW_BASE:?Set the absolute uploaded package directory in a registered submission wrapper}"\n'
            '/public/home/temp/miniconda3/bin/python3 "$BASE/run.py" --base "$BASE" --case ' + branch + ' --ntasks "$SLURM_NTASKS"\n')
        (output_dir / ('submit_' + branch + '.slurm')).write_text(script)
    dependencies = {name: gate.sha(output_dir / name) for name in ('run.py', 'validate_stage.py', 'build_campaign.py', 'frozen_run_cp2k.py', 'frozen_analyze_event_unbiased.py', 'inputs/reference_unbiased.inp', 'selection.json')}
    dependencies.update({p.name: gate.sha(p) for p in output_dir.glob('submit_*.slurm')})
    manifest = {'schema': 'li2025_hbnetwork_workflow/v1', 'candidate_names': [c['name'] for c in candidates],
        'cases': records, 'dependencies': dependencies, 'selection_source': {'path': str(selection_path), 'sha256': gate.sha(selection_path)},
        'reference_input': {'path': str(REFERENCE), 'sha256': gate.REFERENCE_SHA}, 'runner_source': {'path': str(RUNNER), 'sha256': gate.RUNNER_SHA},
        'velocity_policy': 'Fresh MB velocities and zero counters for both smoke and production from the same candidate XYZ/seed; no restart or inherited velocity.',
        'submission_policy': 'No automatic submissions; explicit case selection and current resource review required. Resource request per case:64CPU/32G/48h.',
        'scope': '5 ps unbiased conditional retention of selected states after recorded history; no network causality, electronic bond, equilibrium spontaneity or FES verdict.',
        'human_review': 'PENDING', 'structural_verdict': None}
    gate.dump(output_dir / 'manifest.json', manifest)
    if dry_run:
        spec = __import__('importlib.util', fromlist=['util']).spec_from_file_location('li_hb_package_run', output_dir / 'run.py')
        module = __import__('importlib.util', fromlist=['util']).module_from_spec(spec)
        spec.loader.exec_module(module)
        logs = {}
        for name in manifest['candidate_names']:
            stream = io.StringIO()
            with redirect_stdout(stream):
                code = module.run_case(output_dir, name, output_dir / 'dryrun', True, 64)
            gate.demand(code == 0, 'staging dry-run failed')
            logs[name] = {'returncode': code, 'stdout': stream.getvalue()}
        gate.dump(output_dir / 'DRY_RUN.json', {'status': 'PASS', 'physical_simulation': False, 'cases': logs})
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', nargs='?', default='prepare', choices=['prepare'])
    parser.add_argument('--candidates', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    manifest = build(args.candidates, args.output_dir, args.dry_run)
    print(json.dumps({'status': 'PREPARED', 'candidates': manifest['candidate_names'], 'output_dir': str(args.output_dir), 'submitted': False}))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(9)
