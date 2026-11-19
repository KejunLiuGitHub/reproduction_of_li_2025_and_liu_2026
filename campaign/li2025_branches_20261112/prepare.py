#!/usr/bin/env python3
"""Freeze three Li diagnostic branches from validated, unedited coordinates."""
import json
import shutil
from pathlib import Path
import frozen_prepare_event_unbiased as helper
import frozen_run_cp2k as runner

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
EVIDENCE = ROOT / 'results/li2025_autofollow/decision_20261112'
BOX = [15., 15., 30.]
CASES = [
    ('r_eq_metad', 'li2025_recovery_20261111/r_interface_eq', 5000, 10000, 2026111201),
    ('p_iface_v2', None, None, 5000, 2026111202),
    ('relax_v2', 'li2025_solvation_20261111/relax', 2000, 5000, 2026111202),
]


def main():
    inputs = BASE / 'inputs'
    reference = ROOT / 'campaign/li2025_direct/inputs/interface_r1.inp'
    plumed = ROOT / 'campaign/li2025_direct/inputs/interface_r1.plumed'
    records = []
    for name, parent, step, length, seed in CASES:
        xyz = inputs / (name + '.xyz')
        if parent:
            directory = EVIDENCE / parent
            record_path = directory / 'run_record.json'
            record = json.loads(record_path.read_text())
            assert record['status'] == 'COMPLETED' and record['program_ended']
            assert record['scf']['not_converged_runs'] == 0
            assert record['md_step_last'] == step and record['abort_count'] == 0
            assert record['scf']['converged_runs'] == step + 1
            if 'checks' in record:
                assert all(record['checks'].values())
            if 'unbiased' in record:
                assert record['unbiased']['step_check'] == 'PASS'
                assert record['unbiased']['trajectory_check'] == 'PASS'
            source = directory / (directory.name + '-pos-1.xyz')
            assert runner.trajectory_steps(source) == list(range(0, step + 1, 10))
            assert runner.count_frames(source) == step // 10 + 1
            analysis_name = 'PILOT_ANALYSIS.json' if name == 'r_eq_metad' else 'DISTANCE_ANALYSIS.json'
            analysis = json.loads((directory / analysis_name).read_text())
            assert analysis['trajectory_sha256'] == helper.sha256(source)
            frame = helper.find_event_frame(source, float(step))
            assert frame['step'] == step
            xyz.write_text('296\n' + frame['comment'] + '\n' + '\n'.join(frame['atom_lines']) + '\n')
            origin = {'trajectory': str(source.relative_to(ROOT)), 'trajectory_sha256': helper.sha256(source),
                      'frame_step': step, 'parent_job_id': record['job_id'],
                      'run_record': str(record_path.relative_to(ROOT)), 'run_record_sha256': helper.sha256(record_path)}
        else:
            source = ROOT / 'campaign/li2025_recovery_20261111/inputs/p_interface.xyz'
            shutil.copy2(source, xyz)
            origin = {'coordinate_file': str(source.relative_to(ROOT)), 'coordinate_sha256': helper.sha256(source),
                      'label': 'Public Figure 6b P_inter', 'workbook_range': 'Z2:AC297',
                      'workbook_sha256': 'f4e2bedc004f0abaf180f55e95bd2aa9f7ddc2aa86ef09a683fd2d5afe0a9992'}
        elements, coords = helper.parse_atoms(xyz.read_text().splitlines()[2:])
        geometry = helper.verify_structure(elements, coords, BOX)
        for suffix, steps in [('_smoke', 20), ('', length)]:
            case = name + suffix
            body, removed = helper.build_unbiased_input(reference.read_text(), seed, name,
                project=case, steps=steps, reference_project='interface_r1',
                reference_coord='interface_r1.xyz', expected_cell='ABC 15.0 15.0 30.0')
            files = [xyz]
            biased = name == 'r_eq_metad'
            if biased:
                block = '\n'.join(removed).replace('interface_r1.plumed', case + '.plumed') + '\n'
                body = body.replace('&MOTION\n', '&MOTION\n' + block)
                bias_file = inputs / (case + '.plumed')
                bias_file.write_text(plumed.read_text().replace('interface_r1', case))
                files.append(bias_file)
            inp = inputs / (case + '.inp')
            inp.write_text(body)
            files.append(inp)
            records.append({'name': case, 'branch': name, 'unbiased': not biased,
                'prepared_under_bias': name == 'relax_v2', 'seed': seed, 'steps': steps,
                'cp2k_project': case, 'expected_frames': steps // 10 + 1,
                'box_class': 'interface', 'box_angstrom': BOX, 'source': origin,
                'sc_distance_angstrom': geometry['sc_distance_angstrom'], 'geometry': geometry,
                'cp2k': {'steps_initial_segment': steps, 'trajectory_every_fs': 10},
                'stage': [{'role': 'new', 'source': p.name, 'target': p.name,
                           'sha256': helper.sha256(p), 'bytes': p.stat().st_size} for p in files],
                'outputs': {'output': case + '.out', 'trajectory': case + '-pos-1.xyz',
                            'energy': case + '-1.ener', 'restart': case + '-1.restart'},
                'metad_files': {'hills': 'HILLS_' + case, 'colvar': 'COLVAR_' + case}})
    dependencies = {p.name: helper.sha256(p) for p in BASE.glob('frozen_*.py')}
    dependencies.update({name: helper.sha256(BASE / name) for name in ['run_branch.py', 'validate_stage.py', 'run.slurm']})
    manifest = {'schema': 'li2025_branches/v1', 'cases': records, 'dependencies': dependencies,
                'reference_input_sha256': helper.sha256(reference), 'reference_plumed_sha256': helper.sha256(plumed),
                'velocity_policy': 'Fresh CP2K velocities with recorded GLOBAL SEED; positions unedited; no EXT_RESTART',
                'bias_history': 'r_eq_metad starts a separate new bias history after unbiased 5 ps R preparation',
                'interpretation': 'Conditional retention and initialization sensitivity; descriptive geometry; human review PENDING'}
    (BASE / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({r['name']: r['sc_distance_angstrom'] for r in records}, indent=2))


if __name__ == '__main__':
    main()
