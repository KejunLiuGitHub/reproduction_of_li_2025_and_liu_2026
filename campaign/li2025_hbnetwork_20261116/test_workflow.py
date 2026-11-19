#!/usr/bin/env python3
"""Local technical workflow regression tests; never runs CP2K or submits jobs.

Golden data are the complete, previously verified 5 ps p_iface_v2 trajectory.
All mutations target isolated test fixtures under --outdir. The runner sequencing
checks use mocks; they are not new smoke or production simulations.
"""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import unittest
from unittest.mock import Mock, patch

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / 'results/li2025_autofollow/20261113T213847+0800/li2025_branches_20261112/p_iface_v2'
DEFAULT_OUT = ROOT / 'results/li2025_autofollow/20261116T064438+0800/workflow_audit/tests'
SOURCE_NAMES = ['CHECKS.json', 'run_record.json', 'p_iface_v2.out', 'p_iface_v2-pos-1.xyz',
                'p_iface_v2-1.ener', 'p_iface_v2-1.restart', 'p_iface_v2.inp', 'p_iface_v2.xyz']
TEST_ROOT = DEFAULT_OUT


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def read_json(path):
    return json.loads(Path(path).read_text())


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


gate = load('validate_stage', 'validate_stage.py')
builder = load('li_hb_test_builder', 'build_campaign.py')
workflow = load('li_hb_test_workflow', 'run.py')


def write_detached(path, content):
    path = Path(path)
    if path.is_symlink():
        path.unlink()
    path.write_text(content)


class WorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = TEST_ROOT
        cls.root.mkdir(parents=True, exist_ok=False)
        cls.source_hashes_before = {n: digest(SOURCE / n) for n in SOURCE_NAMES}
        cls.old_record = next(r for r in read_json(ROOT / 'campaign/li2025_branches_20261112/manifest.json')['cases'] if r['name'] == 'p_iface_v2')
        cls.selection = cls.root / 'valid_selection.json'
        dump(cls.selection, {'schema': 'li2025_hbnetwork_candidates/v1', 'candidates': [cls.candidate()]})
        cls.package = cls.root / 'valid_package'
        cls.manifest = builder.build(cls.selection, cls.package, dry_run=True)

    @classmethod
    def candidate(cls, source=SOURCE):
        return {'name': 'golden_candidate', 'trajectory': str(source / 'p_iface_v2-pos-1.xyz'),
                'trajectory_sha256': digest(source / 'p_iface_v2-pos-1.xyz'), 'frame_step': 5000,
                'source_checks': str(source / 'CHECKS.json'), 'source_checks_sha256': digest(source / 'CHECKS.json'),
                'source_job_id': '8569494', 'seed': 2026111601, 'source_prepared_under_bias': False,
                'purpose': 'technical fixture only', 'selection_rule': 'last completed golden frame'}

    def setUp(self):
        self.directory = self.root / self._testMethodName
        self.directory.mkdir()

    def tearDown(self):
        # Persist concise gate evidence; the test script recreates large mutated fixtures.
        for name in ('CHECKS.json', 'DISTANCE_ANALYSIS.json'):
            source = self.directory / 'work' / name
            if source.is_file():
                shutil.copy2(source, self.directory / name)
        for name in ('source', 'work', 'base'):
            target = self.directory / name
            if target.exists():
                shutil.rmtree(target)

    def selection_file(self, candidate, **extra):
        path = self.directory / 'selection.json'
        dump(path, dict(schema='li2025_hbnetwork_candidates/v1', candidates=[candidate], **extra))
        return path

    def source_copy(self):
        target = self.directory / 'source'
        target.mkdir()
        for name in SOURCE_NAMES:
            if name in ('p_iface_v2.out', 'p_iface_v2-1.ener', 'p_iface_v2-1.restart'):
                (target / name).symlink_to(SOURCE / name)
            else:
                shutil.copy2(SOURCE / name, target / name)
        return target

    def assert_rejected_build(self, candidate, regex):
        with self.assertRaisesRegex((ValueError, OSError), regex):
            builder.build(self.selection_file(candidate), self.directory / 'rejected_package')
        self.assertFalse((self.directory / 'rejected_package').exists())

    def golden(self):
        base = self.directory / 'base'
        work = self.directory / 'work'
        (base / 'inputs').mkdir(parents=True)
        work.mkdir()
        shutil.copy2(ROOT / 'campaign/li2025_branches_20261112/frozen_analyze_event_unbiased.py', base / 'frozen_analyze_event_unbiased.py')
        shutil.copy2(SOURCE / 'p_iface_v2.inp', base / 'inputs/reference_unbiased.inp')
        shutil.copy2(SOURCE / 'p_iface_v2.xyz', base / 'inputs/p_iface_v2.xyz')
        for name in SOURCE_NAMES:
            if name != 'CHECKS.json':
                (work / name).symlink_to(SOURCE / name)
        record = copy.deepcopy(self.old_record)
        record['coordinate_file'] = 'p_iface_v2.xyz'
        return base, work, record

    def rejected_gate(self, base, work, record, regex):
        with self.assertRaisesRegex(ValueError, regex):
            gate.validate(work, record, base, expected_job='8569494')
        result = read_json(work / 'CHECKS.json')
        self.assertEqual(result['status'], 'FAIL')
        self.assertIsNone(result['structural_verdict'])
        self.assertEqual(result['human_review'], 'PENDING')

    def test_builder_real_source_and_dryrun(self):
        manifest = self.manifest
        self.assertEqual(manifest['candidate_names'], ['golden_candidate'])
        self.assertEqual(digest(self.package / 'frozen_run_cp2k.py'), gate.RUNNER_SHA)
        self.assertEqual(digest(self.package / 'frozen_analyze_event_unbiased.py'), gate.ANALYZER_SHA)
        dryrun = read_json(self.package / 'DRY_RUN.json')
        self.assertEqual(dryrun['status'], 'PASS')
        self.assertFalse(dryrun['physical_simulation'])
        smoke, production = manifest['cases']
        self.assertEqual([smoke['steps'], production['steps']], [20, 5000])
        self.assertEqual([smoke['expected_frames'], production['expected_frames']], [3, 501])
        self.assertEqual(smoke['coordinate_file'], production['coordinate_file'])
        self.assertEqual(smoke['seed'], production['seed'])
        self.assertEqual(smoke['source'], production['source'])
        self.assertEqual(smoke['selection_metadata']['purpose'], 'technical fixture only')
        for record in manifest['cases']:
            name = record['name']
            staged = self.package / 'dryrun' / name
            self.assertTrue((staged / (name + '.inp')).is_file())
            body = (staged / (name + '.inp')).read_text()
            self.assertNotRegex(body, r'(?i)EXT_RESTART|FREE_ENERGY|METADYN|PLUMED|CONSTRAINT|&VELOCITY')
            self.assertEqual(digest(staged / record['coordinate_file']), record['source']['candidate_xyz_sha256'])
            self.assertTrue(record['unbiased'])
            self.assertFalse(any(staged.glob('*restart*')))
        for name, expected in manifest['dependencies'].items():
            self.assertEqual(digest(self.package / name), expected)

    def test_builder_refuses_overwrite(self):
        with self.assertRaisesRegex(ValueError, 'already exists'):
            builder.build(self.selection, self.package)

    def test_builder_wrong_source_sha(self):
        candidate = self.candidate()
        candidate['trajectory_sha256'] = '0' * 64
        self.assert_rejected_build(candidate, 'source trajectory SHA')

    def test_builder_missing_selected_frame(self):
        candidate = self.candidate()
        candidate['frame_step'] = 5001
        self.assert_rejected_build(candidate, 'selected frame missing')

    def test_builder_failed_checks(self):
        source = self.source_copy()
        checks = read_json(source / 'CHECKS.json')
        checks['status'] = 'FAIL'
        dump(source / 'CHECKS.json', checks)
        self.assert_rejected_build(self.candidate(source), 'CHECKS not all PASS')

    def test_builder_running_source(self):
        source = self.source_copy()
        run = read_json(source / 'run_record.json')
        run['status'] = 'RUNNING'
        dump(source / 'run_record.json', run)
        self.assert_rejected_build(self.candidate(source), 'not normally completed')

    def test_builder_forged_pass_missing_frame_with_updated_hash(self):
        source = self.source_copy()
        trajectory = source / 'p_iface_v2-pos-1.xyz'
        lines = trajectory.read_text().splitlines(keepends=True)
        del lines[298 * 50:298 * 51]
        trajectory.write_text(''.join(lines))
        checks = read_json(source / 'CHECKS.json')
        checks['trajectory_sha256'] = digest(trajectory)
        checks['frames'] = 500
        dump(source / 'CHECKS.json', checks)
        candidate = self.candidate(source)
        self.assertEqual(candidate['trajectory_sha256'], checks['trajectory_sha256'])
        self.assert_rejected_build(candidate, 'XYZ full frame sequence')

    def test_builder_source_input_hash_tamper(self):
        source = self.source_copy()
        inp = source / 'p_iface_v2.inp'
        inp.write_text(inp.read_text().replace('SEED 2026111202', 'SEED 2026111203'))
        self.assertNotEqual(digest(inp), read_json(source / 'run_record.json')['input_sha256'][inp.name])
        self.assert_rejected_build(self.candidate(source), '(?i)input.*hash|input.*SHA|SHA.*input')

    def test_builder_modern_checks_hash_mapping(self):
        source = self.source_copy()
        checks = read_json(source / 'CHECKS.json')
        checks.pop('trajectory_sha256')
        checks['files_sha256'] = {n: digest(source / n) for n in SOURCE_NAMES if n != 'CHECKS.json'}
        dump(source / 'CHECKS.json', checks)
        manifest = builder.build(self.selection_file(self.candidate(source)), self.directory / 'package')
        self.assertEqual(manifest['cases'][0]['source']['independent_technical_recheck'], 'PASS')

    def test_builder_bias_history_cannot_be_erased(self):
        source = self.source_copy()
        checks = read_json(source / 'CHECKS.json')
        checks['prepared_under_bias'] = True
        dump(source / 'CHECKS.json', checks)
        self.assert_rejected_build(self.candidate(source), 'biased source history cannot be erased')

    def test_selection_aliases_and_conflicts(self):
        canonical = self.candidate()
        alias = dict(canonical)
        for old, new in builder.ALIASES.items():
            alias[old] = alias.pop(new)
        self.assertEqual(builder.normalize_selection({'cases': [alias]}), [canonical])
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            builder.normalize_selection({'cases': [alias], 'candidates': [canonical]})
        alias['trajectory_sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'conflicting candidate alias'):
            builder.normalize_selection({'cases': [alias]})

    def test_validator_real_complete_5000_step_golden(self):
        base, work, record = self.golden()
        audit = gate.validate(work, record, base, expected_job='8569494')
        self.assertEqual(audit['status'], 'PASS')
        self.assertEqual([audit['first_step'], audit['last_step'], audit['frames']], [0, 5000, 501])
        self.assertEqual(audit['checks']['all_MD_SCF']['details']['SCF_converged'], 5001)
        self.assertEqual(audit['checks']['complete_restart']['details']['velocities'], 296)
        self.assertIsNone(audit['structural_verdict'])
        self.assertIsNone(audit['network_causal_verdict'])
        self.assertEqual(digest(work / 'DISTANCE_ANALYSIS.json'), audit['distance_analysis']['sha256'])
        self.assertEqual(audit['human_review'], 'PENDING')

    def test_validator_scf_failure(self):
        base, work, record = self.golden()
        log = work / record['outputs']['output']
        body, count = re.subn(r'SCF run converged in\s+\d+\s+steps', 'SCF run NOT converged', log.read_text(), count=1)
        self.assertEqual(count, 1)
        write_detached(log, body)
        self.rejected_gate(base, work, record, 'SCF per-step')

    def test_validator_missing_md_step(self):
        base, work, record = self.golden()
        log = work / record['outputs']['output']
        body, count = re.subn(r'^.*MD\|\s+Step number\s+2500[^\n]*\n', '', log.read_text(), count=1, flags=re.M)
        self.assertEqual(count, 1)
        write_detached(log, body)
        self.rejected_gate(base, work, record, 'MD full step coverage')

    def test_validator_atom_order(self):
        base, work, record = self.golden()
        trajectory = work / record['outputs']['trajectory']
        lines = trajectory.read_text().splitlines(keepends=True)
        fields = lines[2].split()
        fields[0] = 'O'
        lines[2] = ' '.join(fields) + '\n'
        write_detached(trajectory, ''.join(lines))
        self.rejected_gate(base, work, record, 'XYZ atom order')

    def test_validator_analyzer_failure_is_not_pass(self):
        base, work, record = self.golden()
        with patch.object(gate.subprocess, 'run', side_effect=subprocess.CalledProcessError(2, 'synthetic analyzer failure')):
            self.rejected_gate(base, work, record, 'synthetic analyzer failure')
        self.assertNotIn('distance_analysis', read_json(work / 'CHECKS.json'))

    def test_validator_missing_frame(self):
        base, work, record = self.golden()
        trajectory = work / record['outputs']['trajectory']
        lines = trajectory.read_text().splitlines(keepends=True)
        del lines[298 * 50:298 * 51]
        write_detached(trajectory, ''.join(lines))
        self.rejected_gate(base, work, record, 'XYZ missing')

    def test_validator_truncated_frame(self):
        base, work, record = self.golden()
        trajectory = work / record['outputs']['trajectory']
        write_detached(trajectory, trajectory.read_text()[:-10])
        self.rejected_gate(base, work, record, 'truncated XYZ')

    def test_validator_nonfinite_energy(self):
        base, work, record = self.golden()
        energy = work / record['outputs']['energy']
        lines = energy.read_text().splitlines(keepends=True)
        row = next(i for i, line in enumerate(lines) if line.strip() and not line.lstrip().startswith('#'))
        fields = lines[row].split()
        fields[4] = 'NaN'
        lines[row] = ' '.join(fields) + '\n'
        write_detached(energy, ''.join(lines))
        self.rejected_gate(base, work, record, 'nonfinite')

    def test_validator_restart_missing_velocity(self):
        base, work, record = self.golden()
        restart = work / record['outputs']['restart']
        lines = restart.read_text().splitlines(keepends=True)
        row = next(i for i, line in enumerate(lines) if line.strip() == '&VELOCITY')
        del lines[row + 1]
        write_detached(restart, ''.join(lines))
        self.rejected_gate(base, work, record, '296 coordinates/velocities')

    def test_validator_unbiased_input_rejects_restraint_even_with_matching_hashes(self):
        base, work, record = self.golden()
        inp = work / 'p_iface_v2.inp'
        write_detached(inp, inp.read_text() + '\n&CONSTRAINT\n&END CONSTRAINT\n')
        run = read_json(work / 'run_record.json')
        run['input_sha256'][inp.name] = digest(inp)
        write_detached(work / 'run_record.json', json.dumps(run) + '\n')
        next(x for x in record['stage'] if x['target'] == inp.name)['sha256'] = digest(inp)
        self.rejected_gate(base, work, record, 'frozen unbiased science/template')

    def test_validator_rejects_bias_artifact(self):
        base, work, record = self.golden()
        (work / 'HILLS_unexpected').write_text('#! FIELDS time sc\n')
        self.rejected_gate(base, work, record, 'unexpected bias artifacts')

    def mock_runner(self, returns):
        runner = Mock()
        runner.resolve_workdir.side_effect = lambda path, root: Path(path)
        runner.run.side_effect = returns
        return runner

    def test_smoke_nonzero_prevents_production(self):
        runner = self.mock_runner([7])
        with patch.object(workflow, 'load_runner', return_value=runner), patch.object(workflow.gate, 'validate') as validator:
            self.assertEqual(workflow.run_case(self.package, 'golden_candidate', self.directory / 'runs'), 7)
        self.assertEqual([x.args[0] for x in runner.run.call_args_list], ['golden_candidate_smoke'])
        validator.assert_not_called()

    def test_smoke_validation_failure_prevents_production(self):
        runner = self.mock_runner([0])
        with patch.object(workflow, 'load_runner', return_value=runner), patch.object(workflow.gate, 'validate', side_effect=ValueError('synthetic smoke failure')):
            with self.assertRaisesRegex(ValueError, 'synthetic smoke failure'):
                workflow.run_case(self.package, 'golden_candidate', self.directory / 'runs')
        self.assertEqual([x.args[0] for x in runner.run.call_args_list], ['golden_candidate_smoke'])

    def test_successful_sequencing_has_two_fresh_independent_inputs(self):
        runner = self.mock_runner([0, 0])
        with patch.object(workflow, 'load_runner', return_value=runner), patch.object(workflow.gate, 'validate', return_value={'status': 'PASS'}) as validator:
            self.assertEqual(workflow.run_case(self.package, 'golden_candidate', self.directory / 'runs'), 0)
        self.assertEqual([x.args[0] for x in runner.run.call_args_list], ['golden_candidate_smoke', 'golden_candidate'])
        self.assertEqual(validator.call_count, 2)
        self.assertNotEqual(runner.run.call_args_list[0].args[1], runner.run.call_args_list[1].args[1])
        for call in runner.run.call_args_list:
            self.assertFalse(call.args[3])
            self.assertIsNone(call.args[5])
        smoke, production = [x.args[1] for x in validator.call_args_list]
        self.assertEqual(smoke['source'], production['source'])
        self.assertEqual(smoke['coordinate_file'], production['coordinate_file'])
        self.assertEqual(smoke['seed'], production['seed'])

    def test_workflow_requires_explicit_production_case_and_64_tasks(self):
        with self.assertRaisesRegex(ValueError, 'explicit production candidate'):
            workflow.run_case(self.package, 'golden_candidate_smoke', self.directory / 'runs', True)
        with self.assertRaisesRegex(ValueError, '64 MPI tasks'):
            workflow.run_case(self.package, 'golden_candidate', self.directory / 'runs', True, ntasks=32)

    @classmethod
    def tearDownClass(cls):
        after = {n: digest(SOURCE / n) for n in SOURCE_NAMES}
        dump(cls.root / 'SOURCE_IMMUTABILITY.json', {'status': 'PASS' if after == cls.source_hashes_before else 'FAIL', 'before': cls.source_hashes_before, 'after': after})
        if after != cls.source_hashes_before:
            raise AssertionError('Original golden source was modified')


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
    global TEST_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.outdir.resolve()
    if not out.is_relative_to(DEFAULT_OUT.resolve()):
        parser.error('--outdir must remain inside the delegated workflow_audit/tests directory')
    out.mkdir(parents=True, exist_ok=True)
    TEST_ROOT = out / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    code_before = {name: digest(HERE / name) for name in ['build_campaign.py', 'run.py', 'validate_stage.py', 'test_workflow.py']}
    result = unittest.TextTestRunner(verbosity=2, resultclass=RecordingResult).run(unittest.defaultTestLoader.loadTestsFromTestCase(WorkflowTests))
    report = {'schema': 'li2025_hbnetwork_workflow_tests/v1', 'status': 'PASS' if result.wasSuccessful() else 'FAIL',
              'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
              'results': result.records, 'physical_simulation': False, 'HPC_interaction': False,
              'golden_source': str(SOURCE), 'golden_job': '8569494',
              'limitations': ['No complete local 20-step raw smoke fixture was available; smoke staging/sequencing are tested without a new physical simulation.', 'Mocked workflow control-flow tests are not a smoke or production PASS.'],
              'tested_code_sha256': {name: digest(HERE / name) for name in ['build_campaign.py', 'run.py', 'validate_stage.py', 'test_workflow.py']},
              'artifact_directory': str(TEST_ROOT), 'code_sha256_before': code_before}
    report['code_unchanged_during_run'] = code_before == report['tested_code_sha256']
    if not report['code_unchanged_during_run']:
        report['status'] = 'FAIL'
    if TEST_ROOT.exists():
        dump(TEST_ROOT / 'TEST_REPORT.json', report)
    dump(out / 'LATEST_TEST_REPORT.json', report)
    print('TEST_REPORT=' + str(TEST_ROOT / 'TEST_REPORT.json'))
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
