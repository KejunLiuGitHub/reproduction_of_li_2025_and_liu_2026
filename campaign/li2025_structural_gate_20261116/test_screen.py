#!/usr/bin/env python3
"""Boundary tests; synthetic data are test fixtures only, never DFT evidence."""
import copy
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import screen

CFG = json.loads((screen.HERE / 'thresholds.json').read_text())
TEST_ROOT = screen.REPO / 'results/li2025_autofollow/20261116T234029+0800/structural_gate/tests'


def metric(t, sc=1.85, co=1.38, height=.30, clear=True):
    return {'time_fs': float(t), 'step': t, 'SC_A': sc, 'CO_A': co, 'carbon_pyramidal_height_A': height,
            'carbon_angle_min_deg': 100., 'carbon_angle_max_deg': 120., 'connectivity_guard_clear': clear}


class ScreenBoundaryTests(unittest.TestCase):
    def test_100fs_is_11_frames_not_10(self):
        ten = [metric(i) for i in range(0, 100, 10)]
        self.assertEqual(screen.continuous_segments(ten, [True] * 10, 10)[0]['span_fs'], 90)
        self.assertFalse(screen.classify(ten, CFG, {}, True, True)['threshold_sensitivity']['default']['terminal_geometry_candidate'])
        eleven = ten + [metric(100)]
        self.assertTrue(screen.classify(eleven, CFG, {}, True, True)['threshold_sensitivity']['default']['terminal_geometry_candidate'])

    def test_initial_short_remains_retention(self):
        rows = [metric(i) for i in range(0, 5010, 10)]
        report = screen.classify(rows, CFG, {'prepared_under_bias': True}, True, True)
        self.assertEqual(report['status'], 'DEFAULT_RETENTION_GEOMETRY_CANDIDATE')
        self.assertFalse(report['initial_reactant_window_confirmed'])

    def test_public_P_2p075_does_not_pass_default(self):
        r = metric(0, sc=2.075)
        self.assertFalse(all(screen.candidate_flags(r, CFG['default']).values()))
        self.assertTrue(all(screen.candidate_flags(r, CFG['loose']).values()))
        rows = [r] + [metric(i) for i in range(10, 5010, 10)]
        result = screen.classify(rows, CFG, {'public_P_label': True}, True, True)
        self.assertTrue(result['initial_context'].startswith('RETENTION'))
        self.assertNotIn('FORMATION', result['status'])

    def test_missing_t0_cannot_pass_full_coverage(self):
        rows = [metric(i) for i in range(2010, 7010, 10)]
        cov = screen.coverage(rows, CFG, 2000)
        self.assertFalse(cov['passed'])
        self.assertEqual(screen.classify(rows, CFG, {}, cov['passed'], True)['status'], 'UNDECIDED_INCOMPLETE_COVERAGE')

    def test_one_crossing_and_time_gap_break_segments(self):
        rows = [metric(i, sc=2.0 if i == 50 else 1.85) for i in range(0, 110, 10)]
        flags = [all(screen.candidate_flags(r, CFG['default']).values()) for r in rows]
        self.assertEqual([r['span_fs'] for r in screen.continuous_segments(rows, flags, 10)], [40, 40])
        rows = [metric(i) for i in (0, 10, 20, 40, 50, 60)]
        self.assertEqual([r['span_fs'] for r in screen.continuous_segments(rows, [True] * 6, 10)], [20, 20])

    def test_formation_window_is_disjoint_and_controlled_context_not_erased(self):
        rows = [metric(i, sc=3.5, co=1.23, height=.02) if i <= 100 else metric(i) for i in range(0, 5010, 10)]
        r = screen.classify(rows, CFG, {'controlled_initial_conditions': True}, True, True)
        self.assertEqual(r['status'], 'DEFAULT_FORMATION_GEOMETRY_CANDIDATE')
        self.assertGreater(r['threshold_sensitivity']['default']['terminal_sustained_segment']['start_time_fs'], 100)
        rows.pop(5)
        self.assertFalse(screen.classify(rows, CFG, {}, True, True)['initial_reactant_window_confirmed'])

    def test_connectivity_ambiguity_blocks_terminal_candidate(self):
        rows = [metric(i, clear=i != 1000) for i in range(0, 5010, 10)]
        r = screen.classify(rows, CFG, {}, True, True)
        self.assertEqual(r['status'], 'UNRESOLVED_CONNECTIVITY_AMBIGUITY')

    def test_threshold_sensitivity_is_explicit(self):
        rows = [metric(i, sc=1.98) for i in range(0, 5010, 10)]
        r = screen.classify(rows, CFG, {}, True, True)
        self.assertEqual(r['status'], 'DEFAULT_RETENTION_GEOMETRY_CANDIDATE')
        self.assertEqual(r['sensitivity_status'], 'UNRESOLVED_THRESHOLD_SENSITIVITY')

    def test_periodic_geometry_and_published_snapshot(self):
        parser = screen.load_parser()
        source = screen.REPO / 'results/li2025_zenodo/structure_compare/li2025_fig6b_P_iface.xyz'
        coords = np.asarray(parser.read_xyz(source)[0]['coords'])
        a = screen.geometry(coords, [15, 15, 30], CFG['connectivity_guard'])
        coords[screen.S] += [15, -30, 60]
        b = screen.geometry(coords, [15, 15, 30], CFG['connectivity_guard'])
        self.assertAlmostEqual(a['SC_A'], b['SC_A'], places=12)
        self.assertAlmostEqual(a['carbon_pyramidal_height_A'], b['carbon_pyramidal_height_A'], places=12)
        self.assertAlmostEqual(a['SC_A'], 2.075, delta=.002)
        self.assertFalse(all(screen.candidate_flags(a, CFG['default']).values()))

    def test_real_geometry_detects_parent_water_disruption(self):
        source = screen.REPO / 'results/li2025_zenodo/structure_compare/li2025_fig6b_P_iface.xyz'
        coords = np.asarray(screen.load_parser().read_xyz(source)[0]['coords'])
        coords[0] = coords[2] + [2., 0., 0.]
        report = screen.geometry(coords, [15, 15, 30], CFG['connectivity_guard'])
        self.assertFalse(report['connectivity_guard_clear'])
        self.assertTrue(any(i['kind'] == 'parent_water_OH_outside_guard' for i in report['connectivity_ambiguities']))

    def test_strict_parser_rejects_truncated_frame_and_wrong_order(self):
        source = screen.REPO / 'results/li2025_zenodo/structure_compare/li2025_fig6b_P_iface.xyz'
        raw = source.read_text()
        with tempfile.TemporaryDirectory(dir=TEST_ROOT) as tmp:
            p = Path(tmp) / 'bad.xyz'
            p.write_text(raw + '296\nincomplete\nH 0 0 0\n')
            with self.assertRaisesRegex(ValueError, 'truncated'):
                screen.load_parser().read_xyz(p)
            p.write_text(raw.replace('\nH ', '\nC ', 1))
            with self.assertRaisesRegex(ValueError, 'order'):
                screen.load_parser().read_xyz(p)

    def test_existing_output_cannot_be_overwritten(self):
        source = 'results/li2025_zenodo/structure_compare/li2025_fig6b_P_iface.xyz'
        with tempfile.TemporaryDirectory(dir=TEST_ROOT) as tmp:
            out = Path(tmp) / 'evidence'
            out.mkdir()
            keep = out / 'GEOMETRIC_SCREEN.json'
            keep.write_text('do not overwrite')
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                screen.analyze_case({'name': 'fixture', 'box_A': [15, 15, 30], 'files': {'trajectory': source}}, out)
            self.assertEqual(keep.read_text(), 'do not overwrite')

    def test_static_box_must_match_export_header(self):
        source = 'results/li2025_zenodo/structure_compare/li2025_fig6b_P_iface.xyz'
        with tempfile.TemporaryDirectory(dir=TEST_ROOT) as tmp:
            case = {'name': 'fixture', 'static': True, 'box_A': [17, 17, 47], 'files': {'trajectory': source}}
            with self.assertRaisesRegex(ValueError, 'box'):
                screen.analyze_case(case, Path(tmp) / 'bad_box')


class TechnicalBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=TEST_ROOT)
        self.root = Path(self.tmp.name)
        self.paths = {k: self.root / v for k, v in {'input': 'fixture.inp', 'source_xyz': 'source.xyz', 'trajectory': 'fixture-pos-1.xyz',
                                                   'run_record': 'run_record.json', 'checks': 'CHECKS.json', 'energy': 'fixture-1.ener', 'restart': 'fixture-1.restart'}.items()}
        for p in self.paths.values():
            p.write_text('fixture\n')
        self.paths['input'].write_text('CHARGE -2\nMULTIPLICITY 1\nTIMESTEP 1.0\nABC 15 15 30\n')
        self.record = {'case': 'fixture', 'job_id': '123', 'status': 'COMPLETED', 'program_ended': True, 'abort_count': 0, 'returncode': 0,
                       'md_steps_reported': 5000, 'md_step_first': 1, 'md_step_last': 5000,
                       'scf': {'converged_runs': 5001, 'not_converged_runs': 0},
                       'input_sha256': {'fixture.inp': screen.sha(self.paths['input'])}, 'input_cell_angstrom': [15., 15., 30.]}
        self.formal = {'case': 'fixture', 'job_id': '123', 'status': 'PASS', 'checks': {'normal': True}, 'trajectory_sha256': screen.sha(self.paths['trajectory'])}
        self.frames = [{'step': 0, 'time': 0., 'energy': -1., 'coords': [[0, 0, 0]] * 296},
                       {'step': 5000, 'time': 5000., 'energy': -1., 'coords': [[0, 0, 0]] * 296}]
        self.energy = [[i, float(i), 1., 298., -1., 0., 0.] for i in range(5001)]
        self.parser = SimpleNamespace(read_energy=lambda path: self.energy,
                                      read_restart=lambda *args: {'fixture': True},
                                      read_xyz=lambda path: [self.frames[0]])
        self.case = {'name': 'fixture', 'run_case': 'fixture', 'job_id': '123', 'start_step': 0, 'box_A': [15., 15., 30.],
                     'expected_source_sha256': screen.sha(self.paths['source_xyz'])}

    def tearDown(self):
        self.tmp.cleanup()

    def run_tech(self):
        screen.dump(self.paths['run_record'], self.record)
        screen.dump(self.paths['checks'], self.formal)
        return screen.technical_evidence(self.case, self.paths, self.frames, self.parser, {})

    def test_valid_fixture_certificate(self):
        self.assertTrue(self.run_tech()['passed'])

    def test_nonzero_or_missing_returncode_and_identity(self):
        self.record['returncode'] = 1
        self.assertFalse(self.run_tech()['checks']['terminal_certificate_bound_to_raw'])
        self.record.pop('returncode')
        self.assertFalse(self.run_tech()['checks']['terminal_certificate_bound_to_raw'])
        self.record['returncode'] = 0
        self.formal['job_id'] = '124'
        self.assertFalse(self.run_tech()['checks']['terminal_certificate_bound_to_raw'])

    def test_empty_formal_checks_not_pass(self):
        self.formal['checks'] = {}
        self.assertFalse(self.run_tech()['checks']['terminal_certificate_bound_to_raw'])
        self.record['trajectory_sha256'] = screen.sha(self.paths['trajectory'])
        self.formal = {}
        self.assertFalse(self.run_tech()['checks']['terminal_certificate_bound_to_raw'])

    def test_changed_input_hash_and_box(self):
        self.paths['input'].write_text(self.paths['input'].read_text() + '# changed\n')
        self.assertFalse(self.run_tech()['checks']['unbiased_input_and_recorded_sha'])
        self.record['input_sha256']['fixture.inp'] = screen.sha(self.paths['input'])
        self.case['box_A'] = [17., 17., 47.]
        self.assertFalse(self.run_tech()['checks']['unbiased_input_and_recorded_sha'])

    def test_fresh_energy_requires_t0_and_xyz_energy_agreement(self):
        self.energy.pop(0)
        self.assertFalse(self.run_tech()['checks']['complete_energy'])
        self.energy.insert(0, [0, 0., 1., 298., -1., 0., 0.])
        self.frames[-1]['energy'] = -2.
        self.assertFalse(self.run_tech()['checks']['complete_energy'])


if __name__ == '__main__':
    TEST_ROOT.mkdir(parents=True, exist_ok=True)
    unittest.main(verbosity=2)
