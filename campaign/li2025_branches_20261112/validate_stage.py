#!/usr/bin/env python3
"""Audit complete first-segment coverage and reuse the frozen distance analyzer."""
import json
import math
import re
import subprocess
import sys
from pathlib import Path
import frozen_run_cp2k as runner

BASE = Path(__file__).resolve().parent


def numeric_rows(path):
    if not path.exists():
        return []
    return [[float(x) for x in line.split()] for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]


def validate(work, record):
    case = record['name']
    end = record['steps']
    run = json.loads((work / 'run_record.json').read_text())
    outputs = record['outputs']
    trajectory = work / outputs['trajectory']
    lines = trajectory.read_text().splitlines()
    source = BASE / 'inputs' / (record['branch'] + '.xyz')
    source_atoms = [line.split() for line in source.read_text().splitlines()[2:]]
    steps = runner.trajectory_steps(trajectory)
    output = (work / outputs['output']).read_text()
    md_steps = [int(x) for x in re.findall(r'MD\|\s+Step number\s+(\d+)', output)]
    energy = numeric_rows(work / outputs['energy'])
    checks = {'runner_completed': run['status'] == 'COMPLETED' and run['returncode'] == 0,
              'normal_end': run['program_ended'] and run['abort_count'] == 0,
              'all_scf': run['scf']['not_converged_runs'] == 0 and run['scf']['converged_runs'] == end + 1,
              'md_coverage': md_steps == list(range(1, end + 1)),
              'trajectory_coverage': steps == list(range(0, end + 1, 10)),
              'complete_frames': len(lines) == record['expected_frames'] * 298,
              'energy_coverage': [int(row[0]) for row in energy] == list(range(end + 1)),
              'finite_energy': all(math.isfinite(x) for row in energy for x in row),
              'restart_end': runner.restart_md_start_step(work / outputs['restart']) == end}
    frame_valid = True
    for offset in range(0, len(lines), 298):
        block = lines[offset:offset + 298]
        if len(block) != 298 or block[0].strip() != '296':
            frame_valid = False
            break
        atoms = [line.split() for line in block[2:]]
        if any(len(a) != 4 or a[0] != b[0] or not all(math.isfinite(float(x)) for x in a[1:])
               for a, b in zip(atoms, source_atoms)):
            frame_valid = False
            break
    checks['frame_atom_order_and_finite_coordinates'] = frame_valid
    first = [line.split() for line in lines[2:298]]
    checks['first_frame_matches_source'] = len(first) == 296 and all(
        a[0] == b[0] and max(abs(float(x) - float(y)) for x, y in zip(a[1:], b[1:])) < 1e-6
        for a, b in zip(first, source_atoms))
    if not record['unbiased']:
        for kind, stride in [('hills', 30), ('colvar', 10)]:
            rows = numeric_rows(work / record['metad_files'][kind])
            expected = list(range(stride, end + 1, stride))
            checks[kind + '_coverage'] = len(rows) == len(expected) and all(
                abs(row[0] - time) < 0.002 for row, time in zip(rows, expected))
            checks[kind + '_finite'] = all(math.isfinite(x) for row in rows for x in row)
            if kind == 'hills':
                checks['hill_parameters'] = all(len(row) == 5 and abs(row[2] - 0.1) < 1e-8
                    and abs(row[3] - 0.03137545) < 1e-8 and row[4] == -1 for row in rows)
    audit = {'case': case, 'job_id': run['job_id'], 'status': 'PASS' if all(checks.values()) else 'FAIL',
             'checks': checks, 'source': record['source'], 'trajectory_sha256': runner.sha256(trajectory),
             'frames': len(steps), 'first_step': steps[0] if steps else None,
             'last_step': steps[-1] if steps else None, 'biased': not record['unbiased'],
             'prepared_under_bias': record['prepared_under_bias'], 'human_review': 'PENDING',
             'smoke': case.endswith('_smoke'), 'structural_verdict': None}
    (work / 'CHECKS.json').write_text(json.dumps(audit, indent=2) + '\n')
    if not all(checks.values()):
        raise RuntimeError('Stage checks failed: ' + ', '.join(k for k, v in checks.items() if not v))
    target = work / 'DISTANCE_ANALYSIS.json'
    subprocess.run([sys.executable, str(BASE / 'frozen_analyze_event_unbiased.py'),
                    '--trajectory', str(trajectory), '--energy', str(work / outputs['energy']),
                    '--label', case, '--box', '15,15,30', '--out', str(target)], check=True)
    data = json.loads(target.read_text())
    data['source'] = record['source']
    data['job_id'] = run['job_id']
    data['interpretation'] = {'biased': not record['unbiased'], 'prepared_under_bias': record['prepared_under_bias'],
                              'human_review': 'PENDING', 'structural_verdict': None,
                              'description': 'Distance diagnostics under recorded preparation and bias history',
                              'smoke': case.endswith('_smoke')}
    target.write_text(json.dumps(data, indent=2) + '\n')
    return audit
