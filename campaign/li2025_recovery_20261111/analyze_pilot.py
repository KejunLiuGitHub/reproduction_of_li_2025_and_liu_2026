#!/usr/bin/env python3
"""Reuse Li's distance analyzer; validate source/coverage and label public origins."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

base = Path(__file__).resolve().parent
case = sys.argv[1]
manifest = json.loads((base/'pilot_manifest.json').read_text())
record = next(r for r in manifest['cases'] if r['name'] == case)
work = base/case
run = json.loads((work/'run_record.json').read_text())
assert run['status'] == 'COMPLETED' and run['case'] == case
for filename, expected in manifest['dependencies'].items():
    assert hashlib.sha256((base/filename).read_bytes()).hexdigest() == expected
trajectory = work/record['outputs']['trajectory']
xyz = base/'inputs'/record['stage'][0]['source']
first = trajectory.read_text().splitlines()[:298]
reference = xyz.read_text().splitlines()
assert first[0].strip() == '296' and len(first) == len(reference) == 298
for a, b in zip(first[2:], reference[2:]):
    x, y = a.split(), b.split()
    assert x[0] == y[0] and max(abs(float(u)-float(v)) for u,v in zip(x[1:],y[1:])) < 1e-6
output = work/'PILOT_ANALYSIS.json'
subprocess.run([sys.executable, str(base/'frozen_analyze_event_unbiased.py'),
    '--trajectory', str(trajectory), '--energy', str(work/record['outputs']['energy']),
    '--label', case, '--box', ','.join(map(str,record['box_angstrom'])), '--out',str(output)], check=True)
data = json.loads(output.read_text())
assert data['frames']['count'] == record['expected_frames']
assert data['frames']['first_step'] == 0 and data['frames']['last_step'] == record['steps']
data['source'] = record['source']
data['job_id'] = run['job_id']
data['interpretation'] = dict(verdict=None, what_this_is='Distance diagnostics from public Li Figure 6b '+record['source']['label'],
    what_this_is_not=['not FES or rate', 'not confirmed product', 'not three-replica statistics',
                      'not proof of equilibrium or spontaneous formation; product-start tests only retention'],
    human_review='PENDING', smoke=case.endswith('_smoke'))
output.write_text(json.dumps(data,indent=2)+'\n')
