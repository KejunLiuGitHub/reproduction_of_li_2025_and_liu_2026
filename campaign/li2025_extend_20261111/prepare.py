#!/usr/bin/env python3
"""Freeze the completed Li interface 20 ps state into a separate 10 ps continuation."""
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / 'results/li2025_recovery_20261111/continuation_source'
OLD = ROOT / 'results/li2025_direct/interface_r1_seg2'
CASE = 'interface_r1_seg3'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
parent = json.loads((SOURCE / 'run_record.json').read_text())
assert parent['status'] == 'COMPLETED' and parent['program_ended'] and parent['returncode'] == 0
assert str(parent['job_id']) == '8565154' and parent['md_step_last'] == 20000
assert parent['scf']['not_converged_runs'] == 0
text = (SOURCE / 'interface_r1_seg2.inp').read_text()
assert sha(SOURCE / 'interface_r1_seg2.inp') == parent['input_sha256']['interface_r1_seg2.inp']
assert text.count('interface_r1_seg2_parent.restart') == 1
new = text.replace('interface_r1_seg2_parent.restart', CASE + '_parent.restart')
(HERE / 'inputs' / (CASE + '.inp')).write_text(new)
xyz = OLD / 'interface_r1.xyz'
assert sha(xyz) == parent['input_sha256']['interface_r1.xyz']
(HERE / 'inputs' / xyz.name).write_bytes(xyz.read_bytes())
assert sha(SOURCE / 'interface_r1.plumed') == parent['input_sha256']['interface_r1.plumed']
assert re.search(r'^RESTART$', (SOURCE / 'interface_r1.plumed').read_text(), re.M)
record = {
    'name': CASE, 'parent_case': 'interface_r1_seg2', 'parent_job_id': '8565154',
    'parent_workdir_remote': parent['workdir'], 'cp2k_project': 'interface_r1',
    'box_class': 'interface', 'box_angstrom': [15.0, 15.0, 30.0],
    'segment': {'start_step': 20000, 'steps': 10000, 'expected_first_step': 20001,
                'expected_end_step': 30000, 'parent_md_step_last': 20000,
                'segment_time_range_fs': [20000, 30000], 'cumulative_step_range': [1, 30000]},
    'outputs': {'output': CASE + '.out', 'trajectory': 'interface_r1-pos-1.xyz',
                'energy': 'interface_r1-1.ener', 'restart': 'interface_r1-1.restart'},
    'metad_files': {'hills': 'HILLS_interface_r1', 'colvar': 'COLVAR_interface_r1'},
    'metad_baseline': {k: parent['metad'][k] for k in
                      ('hills_rows', 'hills_last_time_fs', 'colvar_rows', 'colvar_last_time_fs')},
    'parent_evidence': {'run_record_sha256': sha(SOURCE / 'run_record.json'),
                        'local_source': str(SOURCE.relative_to(ROOT))},
    'stage': []}
for role, name, target in [
    ('new_input', CASE + '.inp', CASE + '.inp'),
    ('new_xyz', xyz.name, xyz.name),
    ('parent_restart', 'interface_r1-1_20000.restart', CASE + '_parent.restart'),
    ('parent_plumed', 'interface_r1.plumed', 'interface_r1.plumed'),
    ('parent_hills', 'HILLS_interface_r1', 'HILLS_interface_r1'),
    ('parent_colvar', 'COLVAR_interface_r1', 'COLVAR_interface_r1'),
    ('parent_run_record', 'run_record.json', 'parent_run_record.json')]:
    p = (HERE / 'inputs' if role.startswith('new') else SOURCE) / name
    record['stage'].append(dict(role=role, source=name, target=target, sha256=sha(p)))
manifest = {'schema': 'li2025_extend/v1', 'cases': [record],
    'dependencies': {'frozen_run_cp2k.py': sha(HERE / 'frozen_run_cp2k.py')},
    'scope': 'Ordinary biased metaD, same SO3^2- model and parameters; no FES claim.',
    'restart_policy': 'Same as seg2: POS/VEL/COUNTERS; CSVR parameters preserved, random state not restored; SCF_GUESS ATOMIC unchanged.'}
(HERE / 'continuation_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(CASE + ': 20001..30000; parent 8565154; HILLS/COLVAR retained')
