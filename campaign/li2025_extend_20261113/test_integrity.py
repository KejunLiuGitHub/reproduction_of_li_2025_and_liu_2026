#!/usr/bin/env python3
"""Exercise technical refusal gates against corrupted local copies of real seg3.

Writes only to the supplied evidence directory. Source files are never modified.
Parser-isolation cases deliberately omit output SHA expectations so rejection
must come from the tested content gate; a separate case proves the hash gate.
"""
import argparse
import json
from pathlib import Path
import tempfile

import prepare as gate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=gate.DEFAULT_SOURCE)
    parser.add_argument('--outdir', type=Path, default=gate.DEFAULT_AUDIT.parent)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    source = args.source.resolve()
    prior = json.loads((gate.HERE.parent / 'li2025_extend_20261111/continuation_manifest.json').read_text())['cases'][0]
    inputs = {entry['target']: entry['sha256'] for entry in prior['stage']
              if entry['role'] not in ('parent_hills', 'parent_colvar')}
    time_scale = gate.rows(gate.DEFAULT_HISTORY / 'HILLS_interface_r1', 5)[0][0] / 30
    cases = [
        ('scf_failure', 'interface_r1_seg3.out',
         lambda s: s.replace('SCF run converged in', 'SCF run NOT converged in', 1),
         'SCF failure/count mismatch'),
        ('truncated_xyz', 'interface_r1-pos-1.xyz',
         lambda s: '\n'.join(s.splitlines()[:-1]) + '\n', 'truncated or malformed XYZ atom row'),
        ('missing_interior_frame', 'interface_r1-pos-1.xyz',
         lambda s: '\n'.join(s.splitlines()[:2980] + s.splitlines()[3278:]) + '\n',
         'trajectory missing/duplicate/out-of-order frame'),
        ('changed_bias_prefix', 'HILLS_interface_r1',
         lambda s: s.replace('3.890372470471124', '3.790372470471124', 1),
         'parent byte prefix modified'),
        ('truncated_restart', 'interface_r1-1_30000.restart',
         lambda s: '\n'.join(s.splitlines()[:-1]) + '\n', 'truncated restart: open sections'),
        ('missing_energy_step', 'interface_r1-1.ener',
         lambda s: '\n'.join(s.splitlines()[:50] + s.splitlines()[51:]) + '\n',
         'energy step coverage'),
        ('nonfinite_velocity', 'interface_r1-1_30000.restart',
         lambda s: s[:s.index('&VELOCITY')] + s[s.index('&VELOCITY'):].replace(
             s[s.index('&VELOCITY'):].splitlines()[1].split()[0], 'NaN', 1),
         'nonfinite numeric data'),
        ('wrong_source_hash', None, None, 'source SHA256 mismatch'),
    ]
    outcomes = []
    with tempfile.TemporaryDirectory(prefix='negative-fixtures-', dir=args.outdir) as temporary:
        root = Path(temporary)
        for label, filename, transform, expected_error in cases:
            work = root / label
            work.mkdir()
            for path in source.iterdir():
                if path.is_file():
                    (work / path.name).symlink_to(path)
            if filename:
                original = (source / filename).read_text()
                mutated = transform(original)
                gate.demand(mutated != original, 'test did not change content: ' + label)
                (work / filename).unlink()  # Never write through a source symlink.
                (work / filename).write_text(mutated)
                del original, mutated
            source_hashes = {'interface_r1.xyz': '0' * 64} if label == 'wrong_source_hash' else {}
            report = args.outdir / ('NEGATIVE_' + label + '.json')
            rejected = False
            error = ''
            try:
                gate.validate_segment(work, 'interface_r1_seg3', 20000, 30000, '8568711',
                    source / 'interface_r1.xyz', gate.DEFAULT_HISTORY, inputs, source_hashes,
                    report, time_scale)
            except ValueError as failure:
                error = str(failure)
                rejected = expected_error in error
            outcomes.append({'case': label, 'expected_error': expected_error, 'actual_error': error,
                             'rejected_by_expected_gate': rejected, 'evidence': str(report)})
            print(label + ': ' + ('PASS' if rejected else 'FAIL'), flush=True)
    result = {'status': 'PASS' if all(item['rejected_by_expected_gate'] for item in outcomes) else 'FAIL',
              'cases': outcomes, 'source_files_modified': False,
              'scope': 'Negative technical integrity checks on altered copies; no simulation or scientific verdict.'}
    gate.dump(args.outdir / 'NEGATIVE_TESTS.json', result)
    gate.demand(result['status'] == 'PASS', 'negative integrity test failed')


if __name__ == '__main__':
    main()
