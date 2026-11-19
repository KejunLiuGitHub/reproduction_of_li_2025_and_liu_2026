#!/usr/bin/env python3
"""Local intervention tests; no physical simulation, HPC action or source edits.

Pure-Python geometry checks are independent of preparation/validation helpers.
Runtime corruption tests only modify isolated fixtures under the allowed outdir.
"""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import Mock, patch
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ARTIFACTS = ROOT / 'results/li2025_autofollow/20261116T234029+0800/intervention'
ALLOWED_OUT = ARTIFACTS / 'tests'
TEST_ROOT = ALLOWED_OUT
BOX = (15., 15., 30.)
ORDER = ['H', 'H', 'O'] * 96 + ['C', 'O', 'H', 'H', 'S', 'O', 'O', 'O']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def mic(a):
    return [x - b * round(x / b) for x, b in zip(a, BOX)]


def delta(a, b):
    return mic([x - y for x, y in zip(a, b)])


def matrix_vector(matrix, vector):
    return [dot(row, vector) for row in matrix]


def determinant(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def rodrigues(axis, angle_degrees):
    x, y, z = axis
    skew = [[0., -z, y], [z, 0., -x], [-y, x, 0.]]
    angle = math.radians(angle_degrees)
    c, s = math.cos(angle), math.sin(angle)
    return [[c * int(i == j) + (1 - c) * axis[i] * axis[j] + s * skew[i][j]
             for j in range(3)] for i in range(3)]


def xyz(path):
    lines = Path(path).read_text().splitlines()
    assert len(lines) == 298 and lines[0].strip() == '296'
    values = [line.split() for line in lines[2:]]
    assert [x[0] for x in values] == ORDER
    assert all(len(x) == 4 for x in values)
    coordinates = [[float(v) for v in row[1:]] for row in values]
    assert all(math.isfinite(v) for row in coordinates for v in row)
    return coordinates


def water_geometry(coords, oxygen, hydrogens):
    oh = [delta(coords[h - 1], coords[oxygen - 1]) for h in hydrogens]
    lengths = [norm(v) for v in oh]
    angle = math.degrees(math.acos(max(-1., min(1., dot(*oh) / (lengths[0] * lengths[1])))))
    return lengths + [angle, norm(delta(coords[hydrogens[0] - 1], coords[hydrogens[1] - 1]))]


def min_nonbonded(coords):
    excluded = set()
    for o in range(3, 289, 3):
        excluded.update((h - 1, o - 1) for h in (o - 2, o - 1))
    excluded.update((i - 1, j - 1) for i, j in [(289, 290), (289, 291), (289, 292), (293, 294), (293, 295), (293, 296)])
    best = (float('inf'), None)
    for i in range(296):
        for j in range(i + 1, 296):
            if (i, j) not in excluded:
                d = norm(delta(coords[i], coords[j]))
                if d < best[0]:
                    best = (d, [i + 1, j + 1])
    return best


def hb_to_Oc(coords, distance_limit, angle_limit):
    donors = []
    for oxygen in range(3, 289, 3):
        for hydrogen in (oxygen - 2, oxygen - 1):
            ho = delta(coords[oxygen - 1], coords[hydrogen - 1])
            ha = delta(coords[289], coords[hydrogen - 1])
            length = norm(ha)
            angle = math.degrees(math.acos(max(-1., min(1., dot(ho, ha) / (norm(ho) * length)))))
            if length <= distance_limit and angle >= angle_limit:
                donors.append(hydrogen)
    return donors


class IndependentGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.feasibility = read(ARTIFACTS / 'FEASIBILITY.json')
        cls.pairs = cls.feasibility['selected_pairs']
        assert len(cls.pairs) == 2
        cls.evidence = {'independent_math': 'Python standard library only; no preparation/gate geometry helper imported',
                        'feasibility_sha256': sha(ARTIFACTS / 'FEASIBILITY.json'), 'pairs': {}}
        cls.coords = {}
        for pair in cls.pairs:
            step = pair['source_step']
            cls.coords[step] = {name: xyz(info['path']) for name, info in pair['xyz'].items()}
            assert all(sha(info['path']) == info['sha256'] for info in pair['xyz'].values())
            cls.evidence['pairs'][str(step)] = {'source_xyz_sha256': pair['source_frame_sha256']}

    def test_rotation_matrices_are_proper_orthogonal(self):
        for pair in self.pairs:
            for mode in ('on', 'control'):
                with self.subTest(step=pair['source_step'], mode=mode):
                    matrix = pair[mode + '_rotation_matrix']
                    error = max(abs(dot(matrix[i], matrix[j]) - int(i == j)) for i in range(3) for j in range(3))
                    self.assertLess(error, 1e-12)
                    self.assertAlmostEqual(determinant(matrix), 1., delta=1e-12)
                    self.evidence['pairs'][str(pair['source_step'])][mode + '_matrix'] = {'max_RRT_error': error, 'determinant': determinant(matrix)}

    def test_axis_angles_and_inverse_rotation_agree(self):
        for pair in self.pairs:
            self.assertAlmostEqual(norm(pair['axis_unit']), 1., delta=1e-12)
            self.assertAlmostEqual(sum(pair['signed_angles_deg']), 0., delta=1e-12)
            for mode, angle in zip(('on', 'control'), pair['signed_angles_deg']):
                expected = rodrigues(pair['axis_unit'], angle)
                actual = pair[mode + '_rotation_matrix']
                self.assertLess(max(abs(a - b) for x, y in zip(expected, actual) for a, b in zip(x, y)), 1e-12)
            self.assertLess(max(abs(pair['on_rotation_matrix'][i][j] - pair['control_rotation_matrix'][j][i]) for i in range(3) for j in range(3)), 1e-12)

    def test_transformed_H_reconstruct_and_other_294_atoms_are_unchanged(self):
        for pair in self.pairs:
            step, oxygen, hydrogens = pair['source_step'], pair['water_O_1based'], pair['rotated_H_1based']
            source = self.coords[step]['source']
            self.assertEqual(hydrogens, [oxygen - 2, oxygen - 1])
            for mode in ('on', 'control'):
                target = self.coords[step][mode]
                for i in range(1, 297):
                    if i not in hydrogens:
                        self.assertEqual(source[i - 1], target[i - 1])
                    else:
                        moved = matrix_vector(pair[mode + '_rotation_matrix'], delta(source[i - 1], source[oxygen - 1]))
                        expected = [x + y for x, y in zip(source[oxygen - 1], moved)]
                        self.assertLess(max(abs(x - y) for x, y in zip(target[i - 1], expected)), 1e-9)

    def test_water_internal_geometry_preserved(self):
        for pair in self.pairs:
            step, oxygen, hydrogens = pair['source_step'], pair['water_O_1based'], pair['rotated_H_1based']
            before = water_geometry(self.coords[step]['source'], oxygen, hydrogens)
            for mode in ('on', 'control'):
                after = water_geometry(self.coords[step][mode], oxygen, hydrogens)
                self.assertLess(max(abs(x - y) for x, y in zip(before, after)), 1e-8)
                self.evidence['pairs'][str(step)][mode + '_OHH'] = after

    def test_equal_amplitude_H_displacement(self):
        for pair in self.pairs:
            coords = self.coords[pair['source_step']]
            for h in pair['rotated_H_1based']:
                lengths = [norm(delta(coords[mode][h - 1], coords['source'][h - 1])) for mode in ('on', 'control')]
                self.assertAlmostEqual(*lengths, delta=1e-9)

    def test_nonbonded_contact_floor(self):
        for pair in self.pairs:
            for mode in ('source', 'on', 'control'):
                minimum, atoms = min_nonbonded(self.coords[pair['source_step']][mode])
                self.assertGreaterEqual(minimum, 1.4)
                self.evidence['pairs'][str(pair['source_step'])][mode + '_minimum_nonbonded'] = {'angstrom': minimum, 'pair_1based': atoms}

    def test_three_threshold_HB_contrast(self):
        for pair in self.pairs:
            result = {}
            for criterion, limits in [('strict', (2.4, 160.)), ('standard', (2.5, 150.)), ('loose', (2.6, 140.))]:
                donors = {mode: hb_to_Oc(self.coords[pair['source_step']][mode], *limits) for mode in ('source', 'on', 'control')}
                self.assertEqual([len(donors[mode]) for mode in ('source', 'on', 'control')], [1, 2, 1])
                selected = set(pair['rotated_H_1based'])
                self.assertEqual([len(selected.intersection(donors[mode])) for mode in ('source', 'on', 'control')], [0, 1, 0])
                result[criterion] = donors
            self.evidence['pairs'][str(pair['source_step'])]['HB_Oc_donors_1based'] = result

    @classmethod
    def tearDownClass(cls):
        dump(TEST_ROOT / 'INDEPENDENT_GEOMETRY.json', cls.evidence)

PACKAGE = HERE / 'production_v1'


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def replace_coordinate(path, atom, coordinate):
    lines = Path(path).read_text().splitlines()
    symbol = lines[atom + 1].split()[0]
    lines[atom + 1] = symbol + ' ' + ' '.join('%.12f' % value for value in coordinate)
    Path(path).write_text('\n'.join(lines) + '\n')


class InterventionGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.package = PACKAGE.resolve()
        cls.manifest = read(cls.package / 'manifest.json')
        cls.original_hashes = {str(p.relative_to(cls.package)): sha(p) for p in cls.package.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
        load_module('frozen_validate_stage', cls.package / 'frozen_validate_stage.py')
        cls.gate = load_module('validate_stage', cls.package / 'validate_stage.py')
        cls.workflow = load_module('intervention_run_under_test', cls.package / 'run.py')
        cls.evidence = {}

    def setUp(self):
        self.directory = TEST_ROOT / self._testMethodName
        self.directory.mkdir()
        self.base = self.directory / 'package'
        shutil.copytree(self.package, self.base, ignore=shutil.ignore_patterns('dryrun', '__pycache__'))
        self.record = copy.deepcopy(self.manifest['cases'][0])

    def tearDown(self):
        # Keep compact evidence and portable dry-run package, discard mutation fixtures.
        if self._testMethodName != 'test_relocated_package_dryrun_without_workspace_reads':
            shutil.rmtree(self.base)

    def candidate_path(self):
        return self.base / 'inputs' / self.record['coordinate_file']

    def update_candidate_hashes(self):
        digest = sha(self.candidate_path())
        self.record['transformation']['transformed_xyz_sha256'] = digest
        for entry in self.record['stage']:
            if entry['target'] == self.record['coordinate_file']:
                entry['sha256'] = digest

    def reject(self, pattern):
        with self.assertRaisesRegex(ValueError, pattern) as context:
            self.gate.verify_transformation(self.record, self.base)
        dump(self.directory / 'REJECTION.json', {'status': 'EXPECTED_REJECTION', 'error': str(context.exception)})

    def test_all_eight_stage_inputs_and_four_transformations(self):
        result = {}
        self.assertEqual(len(self.manifest['cases']), 8)
        self.assertEqual(len(self.manifest['candidate_names']), 4)
        for record in self.manifest['cases']:
            detail = self.gate.verify_input(self.base / 'inputs', record, self.base)
            result[record['name']] = detail
            self.assertEqual(detail['controlled_transformation']['unrotated_atoms_exactly_unchanged'], 294)
            self.assertGreaterEqual(detail['controlled_transformation']['minimum_nonbonded_contact_A'], 1.4)
            self.assertTrue(record['unbiased'])
            self.assertTrue(record['controlled_initial_conditions'])
            self.assertTrue(record['prepared_under_bias'])
            self.assertEqual(record['human_review'], 'PENDING')
        dump(self.directory / 'STAGE_INPUT_CHECKS.json', result)

    def test_frozen_executor_and_technical_validator_unchanged(self):
        old = ROOT / 'campaign/li2025_hbnetwork_20261116/production_v1'
        for current, former in [('run.py', 'run.py'), ('frozen_validate_stage.py', 'validate_stage.py'),
                                ('frozen_run_cp2k.py', 'frozen_run_cp2k.py'), ('frozen_analyze_event_unbiased.py', 'frozen_analyze_event_unbiased.py')]:
            self.assertEqual(sha(self.base / current), sha(old / former))
        for name, expected in self.manifest['dependencies'].items():
            self.assertEqual(sha(self.base / name), expected)

    def test_pair_seed_protocol_and_fresh_smoke_production(self):
        records = self.manifest['cases']
        for step, seed in [(1150, 2026111603), (2200, 2026111604)]:
            pair = [r for r in records if r['pair_source_step'] == step]
            self.assertEqual(len(pair), 4)
            self.assertEqual({r['seed'] for r in pair}, {seed})
            for condition in ('oc_donor', 'reverse'):
                stages = [r for r in pair if r['condition'] == condition]
                self.assertEqual(sorted(r['steps'] for r in stages), [20, 5000])
                self.assertEqual(stages[0]['source'], stages[1]['source'])
                self.assertEqual(stages[0]['coordinate_file'], stages[1]['coordinate_file'])
                for stage in stages:
                    body = (self.base / 'inputs' / (stage['name'] + '.inp')).read_text()
                    self.assertNotRegex(body, r'(?i)EXT_RESTART|FREE_ENERGY|METADYN|PLUMED|CONSTRAINT|&VELOCITY')

    def test_bad_nonorthogonal_matrix(self):
        self.record['transformation']['rotation_matrix'][0][0] += .01
        self.reject('proper orthogonal rotation')

    def test_valid_rotation_with_wrong_signed_angle(self):
        self.record['transformation']['signed_angle_deg'] *= -1
        self.reject('matrix/axis/angle mismatch')

    def test_selected_H_coordinate_tamper_with_updated_declared_hash(self):
        atom = self.record['transformation']['rotated_H_1based'][0]
        coordinate = xyz(self.candidate_path())[atom - 1]
        coordinate[0] += .01
        replace_coordinate(self.candidate_path(), atom, coordinate)
        self.update_candidate_hashes()
        self.reject('do not reconstruct')

    def test_nonselected_atom_tamper_with_updated_declared_hash(self):
        coordinate = xyz(self.candidate_path())[0]
        coordinate[0] += .01
        replace_coordinate(self.candidate_path(), 1, coordinate)
        self.update_candidate_hashes()
        self.reject('unrotated atom changed')

    def test_source_frame_hash_mismatch(self):
        self.record['transformation']['source_frame_sha256'] = '0' * 64
        self.reject('source frame SHA mismatch')

    def test_audit_hash_mismatch(self):
        self.record['transformation']['source_audit_sha256'] = '0' * 64
        self.reject('source audit SHA mismatch')

    def test_forged_audit_with_failure_and_updated_declared_hash(self):
        path = self.base / 'SOURCE_AUDIT.json'
        audit = read(path)
        audit['checks']['all_MD_SCF_and_normal_end']['pass'] = False
        dump(path, audit)
        self.record['transformation']['source_audit_sha256'] = sha(path)
        self.reject('source audit not all PASS')

    def test_source_trajectory_audit_identity_mismatch(self):
        self.record['source']['trajectory_sha256'] = '0' * 64
        self.reject('source trajectory/audit identity mismatch')

    def test_source_frame_step_identity(self):
        self.record['source']['frame_step'] = 1160
        self.reject('source frame identity mismatch')

    def test_only_fixed_parent_water_H_pair_allowed(self):
        self.record['transformation']['rotated_H_1based'][0] -= 1
        self.reject('only a fixed-parent water H pair')

    def test_contact_floor_after_hash_and_reconstruction_gates(self):
        # Move the same unselected H in both source/candidate copies. This keeps
        # the submitted rotation exactly reconstructible while isolating the
        # nonbonded-contact check; it is explicitly a forged local test fixture.
        source = self.base / 'inputs' / self.record['transformation']['source_coordinate_file']
        position = xyz(source)[3]
        replace_coordinate(source, 1, position)
        replace_coordinate(self.candidate_path(), 1, position)
        self.record['transformation']['source_frame_sha256'] = sha(source)
        self.record['source']['original_frame_sha256'] = sha(source)
        self.update_candidate_hashes()
        self.reject('nonbonded contact below 1.4')

    def test_bias_input_rejected(self):
        inp = self.base / 'inputs' / (self.record['name'] + '.inp')
        inp.write_text(inp.read_text() + '\n&CONSTRAINT\n&END CONSTRAINT\n')
        next(x for x in self.record['stage'] if x['target'] == inp.name)['sha256'] = sha(inp)
        with self.assertRaisesRegex(ValueError, 'frozen unbiased science/template'):
            self.gate.verify_input(self.base / 'inputs', self.record, self.base)

    def test_transformation_failure_writes_FAIL_before_frozen_gate(self):
        work = self.directory / 'work'
        work.mkdir()
        self.record['transformation']['rotation_matrix'][0][0] += .01
        with patch.object(self.gate.frozen, 'validate') as frozen_gate:
            with self.assertRaisesRegex(ValueError, 'proper orthogonal rotation'):
                self.gate.validate(work, self.record, self.base)
        frozen_gate.assert_not_called()
        audit = read(work / 'CHECKS.json')
        self.assertEqual(audit['status'], 'FAIL')
        self.assertEqual(audit['human_review'], 'PENDING')
        self.assertIsNone(audit['structural_verdict'])

    def test_metadata_uses_long_SC_descriptive_label_and_refreezes_distance_sha(self):
        work = self.directory / 'work'
        work.mkdir()
        target = work / 'DISTANCE_ANALYSIS.json'
        old_descriptor = {'step': 0, 'criterion': 'd(S-C)>2.5 A for 10 frames'}
        dump(target, {'interpretation': {'description': 'old conditional retention', 'human_review': 'PENDING', 'structural_verdict': None},
                      'summary': {'sc': {'first_sustained_dissociation': old_descriptor}}})
        old_digest = sha(target)
        fake = {'status': 'PASS', 'checks': {}, 'files_sha256': {target.name: old_digest},
                'distance_analysis': {'path': target.name, 'sha256': old_digest}, 'human_review': 'PENDING', 'structural_verdict': None}
        with patch.object(self.gate.frozen, 'validate', return_value=fake):
            audit = self.gate.validate(work, self.record, self.base)
        data = read(target)
        self.assertNotIn('first_sustained_dissociation', data['summary']['sc'])
        self.assertEqual(data['summary']['sc']['first_sustained_SC_gt_2p5_A_10_frames_descriptive_only'], old_descriptor)
        self.assertIn('long-SC controlled initial water-H orientations', data['interpretation']['description'])
        self.assertTrue(data['interpretation']['controlled_initial_conditions'])
        self.assertFalse(data['interpretation']['source_segment_biased'])
        self.assertNotEqual(old_digest, sha(target))
        self.assertEqual(audit['distance_analysis']['sha256'], sha(target))
        self.assertEqual(audit['files_sha256'][target.name], sha(target))
        self.assertIsNone(audit['structural_verdict'])
        self.assertEqual(read(work / 'CHECKS.json')['status'], 'PASS')

    def test_smoke_validation_failure_never_launches_production(self):
        runner = Mock()
        runner.resolve_workdir.side_effect = lambda path, root: Path(path)
        runner.run.return_value = 0
        with patch.object(self.workflow, 'load_runner', return_value=runner), patch.object(self.workflow.gate, 'validate', side_effect=ValueError('synthetic smoke gate failure')):
            with self.assertRaisesRegex(ValueError, 'synthetic smoke gate failure'):
                self.workflow.run_case(self.base, 'r1150_oc_donor', self.directory / 'runtime')
        self.assertEqual([call.args[0] for call in runner.run.call_args_list], ['r1150_oc_donor_smoke'])

    def test_smoke_production_seed_mismatch_rejected_before_run(self):
        manifest = read(self.base / 'manifest.json')
        next(r for r in manifest['cases'] if r['name'] == 'r1150_oc_donor')['seed'] += 1
        dump(self.base / 'manifest.json', manifest)
        with patch.object(self.workflow, 'load_runner') as loader:
            with self.assertRaisesRegex(ValueError, 'same frozen candidate/seed'):
                self.workflow.run_case(self.base, 'r1150_oc_donor', self.directory / 'runtime')
        loader.assert_not_called()

    def test_relocated_package_dryrun_without_workspace_reads(self):
        import subprocess
        harness = r'''
import sys, json, os
from pathlib import Path
sys.dont_write_bytecode = True
base=Path(sys.argv[1]).resolve(); original_root=Path(sys.argv[2]).resolve()
reads=[]
def guard(event,args):
    if event == 'open' and isinstance(args[0],(str,bytes,os.PathLike)):
        path=Path(os.fsdecode(args[0])).resolve()
        if path.is_relative_to(original_root/'campaign') or path.is_relative_to(original_root/'results'):
            if not path.is_relative_to(base):
                raise RuntimeError('FORBIDDEN ORIGINAL WORKSPACE READ: '+str(path))
            reads.append(str(path.relative_to(base)))
sys.addaudithook(guard)
sys.path.insert(0,str(base))
import run
manifest=json.loads((base/'manifest.json').read_text())
for case in manifest['candidate_names']:
    code=run.run_case(base,case,base/'portable_dryrun',True,64)
    assert code==0
print('PORTABILITY_AUDIT='+json.dumps({'status':'PASS','cases':manifest['candidate_names'],'package_reads':sorted(set(reads)),'physical_simulation':False,'guard':'Any open under original workspace campaign/results outside relocated package raises'}))
'''
        process = subprocess.run([sys.executable, '-B', '-c', harness, str(self.base), str(ROOT)], cwd=self.directory,
                                 text=True, capture_output=True)
        (self.directory / 'stdout.txt').write_text(process.stdout)
        (self.directory / 'stderr.txt').write_text(process.stderr)
        self.assertEqual(process.returncode, 0, process.stderr)
        line = next(x for x in process.stdout.splitlines() if x.startswith('PORTABILITY_AUDIT='))
        audit = json.loads(line.split('=', 1)[1])
        self.assertEqual(audit['status'], 'PASS')
        self.assertEqual(len(audit['cases']), 4)
        self.assertFalse(audit['physical_simulation'])
        dump(self.directory / 'PORTABILITY_AUDIT.json', audit)

    @classmethod
    def tearDownClass(cls):
        after = {str(p.relative_to(cls.package)): sha(p) for p in cls.package.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
        dump(TEST_ROOT / 'PACKAGE_IMMUTABILITY.json', {'status': 'PASS' if after == cls.original_hashes else 'FAIL', 'before': cls.original_hashes, 'after': after})
        if after != cls.original_hashes:
            raise AssertionError('Original frozen package changed during tests')


class RecordingResult(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.records = []

    def addSuccess(self, test):
        super().addSuccess(test)
        self.records.append({'test': test.id(), 'status': 'PASS'})

    def addFailure(self, test, err):
        super().addFailure(test, err)
        self.records.append({'test': test.id(), 'status': 'FAIL', 'error': self._exc_info_to_string(err, test)})

    def addError(self, test, err):
        super().addError(test, err)
        self.records.append({'test': test.id(), 'status': 'ERROR', 'error': self._exc_info_to_string(err, test)})


def main():
    global TEST_ROOT, PACKAGE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=ALLOWED_OUT)
    parser.add_argument('--package', type=Path, default=PACKAGE)
    args = parser.parse_args()
    out = args.outdir.resolve()
    if not out.is_relative_to(ALLOWED_OUT.resolve()):
        parser.error('--outdir must remain within the delegated intervention/tests directory')
    out.mkdir(parents=True, exist_ok=True)
    TEST_ROOT = out / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    TEST_ROOT.mkdir()
    PACKAGE = args.package.resolve()
    frozen_sources = [ARTIFACTS / 'FEASIBILITY.json', ARTIFACTS / 'SOURCE_AUDIT.json', Path(__file__)]
    frozen_sources += [Path(v['path']) for p in read(ARTIFACTS / 'FEASIBILITY.json')['selected_pairs'] for v in p['xyz'].values()]
    before = {str(p): sha(p) for p in frozen_sources}
    suite = unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromTestCase(cls)
                               for cls in (IndependentGeometryTests, InterventionGateTests)])
    result = unittest.TextTestRunner(verbosity=2, resultclass=RecordingResult).run(suite)
    after = {str(p): sha(p) for p in frozen_sources}
    report = {'schema': 'li296_intervention_tests/v1', 'status': 'PASS' if result.wasSuccessful() and before == after else 'FAIL',
              'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
              'results': result.records, 'source_files_unchanged': before == after, 'source_sha256': before,
              'package_path': str(PACKAGE), 'package_manifest_sha256': sha(PACKAGE / 'manifest.json'),
              'package_code_sha256': {n: sha(PACKAGE / n) for n in ['run.py', 'validate_stage.py', 'frozen_validate_stage.py', 'build_campaign.py']},
              'physical_simulation': False, 'HPC_interaction': False, 'human_review': 'PENDING', 'structural_verdict': None,
              'limitations': ['No new CP2K smoke/production execution; wrapper sequencing and postprocessing tests use explicit mocks.', 'Unchanged 26-test Li technical-gate suite was not rerun; frozen gate and executor byte identity was checked.', 'Equal-angle orientation contrast is not an energy-matched or isolated-HB causal control.'],
              'artifact_directory': str(TEST_ROOT)}
    dump(TEST_ROOT / 'TEST_REPORT.json', report)
    dump(out / 'LATEST_TEST_REPORT.json', report)
    print('TEST_REPORT=' + str(TEST_ROOT / 'TEST_REPORT.json'))
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
