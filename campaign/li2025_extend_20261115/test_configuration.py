#!/usr/bin/env python3
"""Target tests for the seg5 configuration and new retrieval/CHECKS provenance.

The inherited content gates are byte-identical to the tested seg4 code. This
test exercises new adapters using isolated altered copies, never parent files.
"""
import ast
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import prepare as gate


def main():
    base, evidence = gate.HERE, gate.DEFAULT_AUDIT.parent
    old = base.parent / 'li2025_extend_20261113'
    manifest = json.loads((base / 'continuation_manifest.json').read_text())
    record = manifest['cases'][0]
    old_text = (old / 'prepare.py').read_text()
    new_text = (base / 'prepare.py').read_text()
    checks = {
        'frozen_runner_identical': (old / 'frozen_run_cp2k.py').read_bytes() == (base / 'frozen_run_cp2k.py').read_bytes(),
        'all_inherited_integrity_helpers_identical': old_text[old_text.index('def sha('):old_text.index('def prepare(')]
            == new_text[new_text.index('def sha('):new_text.index('def retrieval_hashes(')],
        'scientific_input_only_changes_parent_restart_filename':
            (base / 'inputs/interface_r1_seg5.inp').read_bytes()
            == (gate.DEFAULT_SOURCE / 'interface_r1_seg4.inp').read_bytes().replace(
                b'interface_r1_seg4_parent.restart', b'interface_r1_seg5_parent.restart'),
        'manifest_identity': record['name'] == 'interface_r1_seg5'
            and record['parent_case'] == 'interface_r1_seg4' and record['parent_job_id'] == '8571154',
        'segment_range': record['segment']['start_step'] == 40000
            and record['segment']['expected_first_step'] == 40001
            and record['segment']['expected_end_step'] == 50000
            and record['segment']['steps'] == 10000,
        'bias_rows': record['metad_baseline']['hills_rows'] == 1333
            and record['technical_policy']['hills_new_rows'] == 333
            and record['metad_baseline']['colvar_rows'] == 4000
            and record['technical_policy']['colvar_new_rows'] == 1000,
        'stop_rule': 'STOP at cumulative 50 ps' in manifest['stop_rule'],
        'no_inherited_submission': not (base / 'SUBMISSION.json').exists(),
        'dependency_hashes': all(gate.sha(base / name) == digest for name, digest in manifest['dependencies'].items()),
    }
    for path in base.glob('*.py'):
        ast.parse(path.read_text())
    subprocess.run(['bash', '-n', str(base / 'run.slurm')], check=True)
    checks['python_and_shell_syntax'] = True
    slurm = (base / 'run.slurm').read_text()
    checks['slurm_config'] = all(value in slurm for value in (
        '--partition=compute', '--ntasks=64', '--mem=32G', '--time=48:00:00',
        'saa-vib-v1/li2025_extend_20261115\n',
        'saa-vib-v1/li2025_extend_20261113/interface_r1_seg4\n', '--case interface_r1_seg5'))
    formal = json.loads((gate.DEFAULT_SOURCE / 'CHECKS.json').read_text())
    required = set(formal['files_sha256']) | {'CHECKS.json'}
    remote = gate.retrieval_hashes(gate.DEFAULT_SOURCE_SHA, required)
    metadata = json.loads(gate.DEFAULT_SOURCE_SHA.read_text())
    negative = []
    with tempfile.TemporaryDirectory(prefix='configuration-tests-', dir=evidence) as temporary:
        fixtures = Path(temporary)
        for label in ('live_metadata', 'bad_remote_line', 'missing_checkpoint',
                      'formal_failed', 'formal_wrong_job', 'formal_bad_source_hash'):
            fixture = fixtures / label
            fixture.mkdir()
            badmeta = json.loads(json.dumps(metadata))
            badformal = json.loads(json.dumps(formal))
            expected = ''
            if label == 'live_metadata':
                badmeta['classification'] = 'LIVE_SNAPSHOT'
                expected = 'source is not stable terminal'
            elif label == 'bad_remote_line':
                badmeta['remote_local_sha'][0]['remote_digest_line'] += 1
                expected = 'raw remote SHA evidence disagreement'
            elif label == 'missing_checkpoint':
                badmeta['remote_local_sha'] = [item for item in badmeta['remote_local_sha']
                    if not item['file'].endswith('interface_r1-1_40000.restart')]
                expected = 'remote source SHA index misses required evidence'
            elif label == 'formal_failed':
                badformal['status'] = 'FAIL'
                expected = 'parent CHECKS is not all PASS'
            elif label == 'formal_wrong_job':
                badformal['job_id'] = '8568711'
                expected = 'parent CHECKS identity/range mismatch'
            else:
                badformal['files_sha256']['interface_r1-pos-1.xyz'] = '0' * 64
                expected = 'parent CHECKS/remote SHA mismatch'
            caught = ''
            try:
                if label.startswith('formal_'):
                    gate.dump(fixture / 'CHECKS.json', badformal)
                    # Simulate valid transport of a semantically inconsistent CHECKS.
                    expected_remote = dict(remote, **{'CHECKS.json': gate.sha(fixture / 'CHECKS.json')})
                    for name in formal['files_sha256']:
                        (fixture / name).symlink_to(gate.DEFAULT_SOURCE / name)
                    gate.verify_parent_checks(fixture, expected_remote, required)
                else:
                    gate.dump(fixture / 'RETRIEVAL_METADATA.json', badmeta)
                    for name in ('HASHES.json', 'REMOTE_EVIDENCE.txt'):
                        (fixture / name).symlink_to(gate.DEFAULT_SOURCE_SHA.parent / name)
                    gate.retrieval_hashes(fixture / 'RETRIEVAL_METADATA.json', required)
            except ValueError as error:
                caught = str(error)
            negative.append({'case': label, 'expected_error': expected, 'actual_error': caught,
                             'pass': expected in caught})
        dry = subprocess.run([sys.executable, str(base / 'run_continuation.py'),
            '--parent-dir', str(gate.DEFAULT_SOURCE), '--workdir-root', str(evidence),
            '--workdir', str(evidence / 'staging_smoke'), '--dry-run'], check=True, capture_output=True, text=True)
        (evidence / 'DRY_RUN.json').write_text(dry.stdout)
        staged = json.loads(dry.stdout)
        checks['dry_run'] = staged['dry_run'] and staged['segment'] == record['segment']
        checks['staged_hashes'] = all(gate.sha(evidence / 'staging_smoke' / entry['target']) == entry['sha256']
                                      for entry in record['stage'])
    checks['new_provenance_negatives'] = all(case['pass'] for case in negative)
    result = {'recorded_at': datetime.now(timezone.utc).isoformat(),
        'status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks,
        'negative_cases': negative, 'manifest_sha256': gate.sha(base / 'continuation_manifest.json'),
        'parent_checks_sha256': gate.sha(gate.DEFAULT_AUDIT), 'dry_run_sha256': gate.sha(evidence / 'DRY_RUN.json'),
        'parent_source_files': len(remote), 'parent_files_modified': False,
        'prior_content_gate_negative_tests': '20261113T213847+0800/continuation_validation/NEGATIVE_TESTS.json; 8 passed, unchanged helpers not retested'}
    gate.dump(evidence / 'IMPLEMENTATION_CHECKS.json', result)
    print(json.dumps({'status': result['status'], 'checks': checks, 'manifest_sha256': result['manifest_sha256']}, indent=2))
    gate.demand(result['status'] == 'PASS', 'configuration/provenance tests failed')


if __name__ == '__main__':
    main()
