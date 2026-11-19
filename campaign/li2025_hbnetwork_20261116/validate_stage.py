#!/usr/bin/env python3
"""Strict Li296 technical checks for fresh unbiased smoke/production stages.

This module produces no structural, bond, network-mechanism or FES verdict.
"""
import hashlib
import json
import math
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ORDER = ['H', 'H', 'O'] * 96 + ['C', 'O', 'H', 'H', 'S', 'O', 'O', 'O']
REFERENCE_SHA = '9536573c2158025949f53c2b03d31a0c5b2ee7dca41b2a08542aed4b53766319'
RUNNER_SHA = 'c3d65060daafbe1076304f092f72c2e40e738aea9cb6f4cc63d8263298397d61'
ANALYZER_SHA = '8529ce9c63f694fd93f86ed0fbac292aff4bf77a612b4069b2d45bb7c96c01cd'


def demand(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def numbers(tokens):
    values = [float(x.replace('D', 'E').replace('d', 'e')) for x in tokens]
    demand(all(math.isfinite(x) for x in values), 'nonfinite numeric value')
    return values


def read_xyz(path):
    frames = []
    with Path(path).open() as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            block = [header, handle.readline()] + [handle.readline() for _ in range(296)]
            demand(all(line and line.endswith('\n') for line in block), 'truncated XYZ frame/line')
            demand(header.strip() == '296', 'XYZ must contain 296 atoms')
            atoms = [line.split() for line in block[2:]]
            demand(all(len(a) == 4 for a in atoms) and [a[0] for a in atoms] == ORDER,
                   'XYZ atom order/width mismatch')
            coords = [numbers(a[1:]) for a in atoms]
            comment = block[1]
            step = re.search(r'\bi\s*=\s*(\d+)', comment)
            time = re.search(r'\btime\s*=\s*([^,\s]+)', comment)
            energy = re.search(r'\bE\s*=\s*([^,\s]+)', comment)
            frames.append({'step': int(step[1]) if step else None,
                'time': numbers([time[1]])[0] if time else None,
                'energy': numbers([energy[1]])[0] if energy else None,
                'coords': coords, 'raw': ''.join(block)})
    demand(bool(frames), 'empty XYZ')
    return frames


def read_energy(path):
    rows = []
    with Path(path).open() as handle:
        for line in handle:
            demand(line.endswith('\n'), 'truncated energy line')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            fields = line.split()
            demand(len(fields) == 7, 'energy row width mismatch')
            rows.append(numbers(fields))
    demand(bool(rows), 'empty energy table')
    return rows


def read_log(path):
    steps, bins, count, pending, failures, ends, aborts = [], [], 0, 0, 0, 0, 0
    with Path(path).open() as handle:
        for line in handle:
            demand(line.endswith('\n'), 'truncated output log line')
            if 'SCF run NOT converged' in line:
                failures += 1
            if 'SCF run converged in' in line:
                demand(re.search(r'SCF run converged in\s+\d+\s+steps', line), 'malformed SCF marker')
                count += 1
                pending += 1
            match = re.search(r'MD\|\s+Step number\s+(\d+)', line)
            if match:
                steps.append(int(match[1]))
                bins.append(pending)
                pending = 0
            ends += line.count('PROGRAM ENDED AT')
            aborts += int(bool(re.search(r'\bABORT\b', line)))
    return {'steps': steps, 'scf_bins': bins, 'scf_converged': count,
            'scf_pending': pending, 'scf_failed': failures,
            'program_ended_count': ends, 'abort_count': aborts}


def read_restart(path, expected_step, expected_steps, final_coords):
    """Parse balanced sections and all 296 finite positions/velocities/counters."""
    stack, blocks = [], {}
    with Path(path).open() as handle:
        for raw in handle:
            demand(raw.endswith('\n'), 'truncated restart line')
            line = raw.split('#', 1)[0].strip()
            if not line:
                continue
            fields = line.split()
            if line.startswith('&'):
                tag = fields[0][1:].upper()
                if tag == 'END':
                    demand(bool(stack) and (len(fields) == 1 or (len(fields) == 2 and fields[1].upper() == stack[-1])), 'unbalanced restart END')
                    stack.pop()
                else:
                    stack.append(tag)
                    key = tuple(stack)
                    demand(key not in blocks or 'KIND' in key, 'duplicate restart section')
                    blocks.setdefault(key, [])
            else:
                demand(bool(stack), 'restart data outside section')
                blocks[tuple(stack)].append(fields)
    demand(not stack, 'truncated restart/open section')
    coord = blocks.get(('FORCE_EVAL', 'SUBSYS', 'COORD'), [])
    velocity = blocks.get(('FORCE_EVAL', 'SUBSYS', 'VELOCITY'), [])
    demand(len(coord) == len(velocity) == 296, 'restart needs 296 coordinates/velocities')
    demand(all(len(a) == 4 for a in coord) and [a[0] for a in coord] == ORDER, 'restart atom order/row width')
    demand(all(len(v) == 3 for v in velocity), 'restart velocity row width')
    coords = [numbers(a[1:]) for a in coord]
    velocities = [numbers(v) for v in velocity]
    md_rows = blocks.get(('MOTION', 'MD'), [])
    md = {row[0].upper(): row[1:] for row in md_rows}
    demand(len(md) == len(md_rows), 'duplicate restart MD keyword')
    for key in ('STEP_START_VAL', 'TIME_START_VAL', 'TIMESTEP', 'STEPS', 'ECONS_START_VAL'):
        demand(len(md.get(key, [])) == 1, 'missing restart counter ' + key)
        numbers(md[key])
    demand(float(md['STEP_START_VAL'][0]) == expected_step, 'restart step mismatch')
    demand(abs(float(md['TIME_START_VAL'][0]) - expected_step) < 1e-5, 'restart time mismatch')
    demand(abs(float(md['TIMESTEP'][0]) - 1.0) < 1e-12 and int(md['STEPS'][0]) == expected_steps, 'restart timestep/STEPS mismatch')
    demand(md.get('ENSEMBLE') == ['NVT'], 'restart ensemble mismatch')
    maximum = max(abs(a - b) for x, y in zip(coords, final_coords) for a, b in zip(x, y))
    demand(maximum < 1e-7, 'restart coordinates differ from final XYZ')
    return {'step': expected_step, 'time_fs': float(md['TIME_START_VAL'][0]),
            'coordinates': len(coords), 'velocities': len(velocities),
            'finite_positions_and_velocities': True, 'max_final_frame_delta_A': maximum}


def make_input(reference, name, coordinate_file, seed, steps):
    demand(re.fullmatch(r'[a-z][a-z0-9_]{0,62}', name), 'unsafe case name')
    demand(Path(coordinate_file).name == coordinate_file and coordinate_file.endswith('.xyz'), 'unsafe coordinate filename')
    demand(type(seed) is int and 0 < seed < 2147483647, 'invalid integer seed')
    demand(steps in (20, 5000), 'only smoke20/production5000 are registered')
    body = reference
    for key, old, new in [('PROJECT', 'p_iface_v2', name), ('COORD_FILE_NAME', 'p_iface_v2.xyz', coordinate_file),
                          ('SEED', '2026111202', str(seed)), ('STEPS', '5000', str(steps))]:
        pattern = r'(?m)^(\s*' + key + r'\s+)' + re.escape(old) + r'(\s*)$'
        body, count = re.subn(pattern, lambda m: m[1] + new + m[2], body)
        demand(count == 1, 'reference keyword absent/ambiguous: ' + key)
    return body


def verify_input(work, record, base):
    ref = Path(base) / 'inputs/reference_unbiased.inp'
    demand(sha(ref) == REFERENCE_SHA, 'frozen reference input changed')
    expected = make_input(ref.read_text(), record['name'], record['coordinate_file'], record['seed'], record['steps'])
    actual = (Path(work) / (record['name'] + '.inp')).read_text()
    demand(actual == expected, 'input differs from frozen unbiased science/template')
    forbidden = r'(?mi)^\s*(?:&(?:FREE_ENERGY|METADYN|EXT_RESTART|CONSTRAINT|COLLECTIVE|VELOCITY)\b|PLUMED_INPUT_FILE\b|USE_PLUMED\b)'
    demand(not re.search(forbidden, actual), 'bias/constraint/restart/velocity section is forbidden')
    demand(record['unbiased'] is True, 'stage must be explicitly unbiased')
    return {'frozen_reference_sha256': REFERENCE_SHA, 'bias_and_constraints_absent': True,
            'restart_and_inherited_velocities_absent': True, 'seed': record['seed']}


def validate(work, record, base, expected_job=None):
    work, base = Path(work), Path(base)
    audit = {'schema': 'li2025_hbnetwork_technical_checks/v1', 'recorded_at': datetime.now(timezone.utc).isoformat(),
        'case': record['name'], 'status': 'FAIL', 'checks': {}, 'source': record.get('source'),
        'biased': False, 'prepared_under_bias': record['prepared_under_bias'],
        'smoke': record['steps'] == 20, 'structural_verdict': None, 'network_causal_verdict': None,
        'FES_verdict': None, 'human_review': 'PENDING', 'interpretation': 'Conditional retention after recorded source history; technical checks only.'}
    def check(name, function):
        value = function()
        audit['checks'][name] = {'pass': True, 'details': value}
        return value
    try:
        name, steps = record['name'], record['steps']
        demand(steps in (20, 5000), 'unsupported stage duration')
        run = json.loads((work / 'run_record.json').read_text())
        audit['job_id'] = run['job_id']
        demand(run['case'] == name and (expected_job is None or str(run['job_id']) == str(expected_job)), 'run identity mismatch')
        demand(run['status'] == 'COMPLETED' and run['returncode'] == 0 and run['program_ended'] and run['abort_count'] == 0, 'runner not normally completed')
        check('unbiased_frozen_input', lambda: verify_input(work, record, base))
        staged = {}
        for entry in record['stage']:
            actual = sha(work / entry['target'])
            demand(actual == entry['sha256'] == run['input_sha256'].get(entry['target']), 'staged/run-record input SHA mismatch')
            staged[entry['target']] = actual
        audit['checks']['input_hashes'] = {'pass': True, 'details': staged}
        demand(not any(list(work.glob(pattern)) for pattern in ('HILLS*', 'COLVAR*', '*.plumed', 'PLUMED.OUT')), 'unexpected bias artifacts')
        log = read_log(work / record['outputs']['output'])
        demand(log['steps'] == list(range(1, steps + 1)), 'MD full step coverage mismatch')
        demand(log['scf_converged'] == steps + 1 and log['scf_failed'] == 0 and log['scf_bins'] == [2] + [1] * (steps - 1) and log['scf_pending'] == 0, 'SCF per-step coverage/failure mismatch')
        demand(log['program_ended_count'] == 1 and log['abort_count'] == 0, 'log not uniquely normally terminated')
        demand(run['md_step_first'] == 1 and run['md_step_last'] == steps and run['scf']['converged_runs'] == steps + 1 and run['scf']['not_converged_runs'] == 0, 'run record disagrees with full log')
        audit['checks']['all_MD_SCF'] = {'pass': True, 'details': {'MD_steps': steps, 'SCF_converged': steps + 1, 'SCF_failed': 0}}
        trajectory = work / record['outputs']['trajectory']
        frames = read_xyz(trajectory)
        demand([f['step'] for f in frames] == list(range(0, steps + 1, 10)), 'XYZ missing/duplicate/out-of-order frame')
        demand(len(frames) == record['expected_frames'] and run['trajectory']['frames'] == len(frames) and run['trajectory']['first_step'] == 0 and run['trajectory']['last_step'] == steps, 'record/frame count mismatch')
        demand(all(f['time'] is not None and f['energy'] is not None and abs(f['time'] - f['step']) < 0.0011 for f in frames), 'XYZ time/energy missing or inconsistent')
        source = read_xyz(base / 'inputs' / record['coordinate_file'])
        demand(len(source) == 1, 'candidate source must contain one frame')
        maximum = max(abs(x - y) for a, b in zip(source[0]['coords'], frames[0]['coords']) for x, y in zip(a, b))
        demand(maximum < 1e-6, 'first frame differs from frozen candidate')
        audit['checks']['complete_XYZ_and_source'] = {'pass': True, 'details': {'frames': len(frames), 'first_step': 0, 'last_step': steps, 'atoms': 296, 'max_initial_delta_A': maximum}}
        energy = read_energy(work / record['outputs']['energy'])
        demand([r[0] for r in energy] == list(range(steps + 1)), 'energy full step coverage mismatch')
        demand(all(abs(r[0] - r[1]) < 2e-6 for r in energy), 'energy time mismatch')
        demand(all(abs(energy[f['step']][4] - f['energy']) < 2e-8 for f in frames), 'XYZ/energy potential mismatch')
        audit['checks']['complete_finite_energy'] = {'pass': True, 'details': {'rows': len(energy), 'first_step': 0, 'last_step': steps}}
        check('complete_restart', lambda: read_restart(work / record['outputs']['restart'], steps, steps, frames[-1]['coords']))
        audit['trajectory_sha256'] = sha(trajectory)
        audit['frames'], audit['first_step'], audit['last_step'] = len(frames), 0, steps
        audit['files_sha256'] = {p.name: sha(p) for p in [work / 'run_record.json'] + [work / n for n in record['outputs'].values()] + [work / x['target'] for x in record['stage']]}
        analyzer = base / 'frozen_analyze_event_unbiased.py'
        demand(sha(analyzer) == ANALYZER_SHA, 'frozen distance analyzer SHA mismatch')
        target = work / 'DISTANCE_ANALYSIS.json'
        subprocess.run([sys.executable, '-B', str(analyzer), '--trajectory', str(trajectory),
                        '--energy', str(work / record['outputs']['energy']), '--label', name,
                        '--box', '15,15,30', '--out', str(target)], check=True, capture_output=True, text=True)
        data = json.loads(target.read_text())
        data.update(source=record.get('source'), job_id=run['job_id'])
        data['interpretation'] = {'biased': False, 'unbiased': True,
            'prepared_under_bias': record['prepared_under_bias'], 'smoke': steps == 20,
            'human_review': 'PENDING', 'structural_verdict': None, 'network_causal_verdict': None,
            'description': 'Distance diagnostics for conditional retention after recorded source history; not spontaneous formation or a network causal comparison.'}
        dump(target, data)
        audit['distance_analysis'] = {'path': target.name, 'sha256': sha(target), 'analyzer_sha256': ANALYZER_SHA}
        audit['files_sha256'][target.name] = sha(target)
        audit['status'] = 'PASS'
    except (ValueError, OSError, KeyError, TypeError, IndexError, subprocess.CalledProcessError) as error:
        audit['error'] = str(error)
    dump(work / 'CHECKS.json', audit)
    if audit['status'] != 'PASS':
        raise ValueError('Stage technical gate FAIL: ' + audit['error'])
    return audit
