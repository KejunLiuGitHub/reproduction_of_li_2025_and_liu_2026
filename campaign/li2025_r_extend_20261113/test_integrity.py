#!/usr/bin/env python3
"""Exercise technical refusal gates against corrupted local copies of real initial R segment.

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
    prior = next(case for case in json.loads((gate.HERE.parent / 'li2025_branches_20261112/manifest.json').read_text())['cases'] if case['name'] == 'r_eq_metad')
    inputs = {entry['target']: entry['sha256'] for entry in prior['stage']
              if entry['role'] not in ('parent_hills', 'parent_colvar')}
    time_scale = gate.rows(gate.DEFAULT_HISTORY / 'HILLS_r_eq_metad', 5)[0][0] / 30
    cases = [
        ('scf_failure', 'r_eq_metad.out',
         lambda s: s.replace('SCF run converged in', 'SCF run NOT converged in', 1),
         'SCF failure/count mismatch'),
        ('truncated_xyz', 'r_eq_metad-pos-1.xyz',
         lambda s: '\n'.join(s.splitlines()[:-1]) + '\n', 'truncated or malformed XYZ atom row'),
        ('missing_interior_frame', 'r_eq_metad-pos-1.xyz',
         lambda s: '\n'.join(s.splitlines()[:2980] + s.splitlines()[3278:]) + '\n',
         'trajectory missing/duplicate/out-of-order frame'),
        ('truncated_restart', 'r_eq_metad-1_10000.restart',
         lambda s: '\n'.join(s.splitlines()[:-1]) + '\n', 'truncated restart: open sections'),
        ('missing_energy_step', 'r_eq_metad-1.ener',
         lambda s: '\n'.join(s.splitlines()[:50] + s.splitlines()[51:]) + '\n',
         'energy step coverage'),
        ('nonfinite_velocity', 'r_eq_metad-1_10000.restart',
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
            source_hashes = {'r_eq_metad.xyz': '0' * 64} if label == 'wrong_source_hash' else {}
            report = args.outdir / ('NEGATIVE_' + label + '.json')
            rejected = False
            error = ''
            try:
                gate.validate_segment(work, 'r_eq_metad', 0, 10000, '8569493',
                    source / 'r_eq_metad.xyz', gate.DEFAULT_HISTORY, inputs, source_hashes,
                    report, time_scale)
            except ValueError as failure:
                error = str(failure)
                rejected = expected_error in error
            outcomes.append({'case': label, 'expected_error': expected_error, 'actual_error': error,
                             'rejected_by_expected_gate': rejected, 'evidence': str(report)})
            print(label + ': ' + ('PASS' if rejected else 'FAIL'), flush=True)
    # The continued bias history must preserve its parent byte prefix exactly.
    # Test this primitive directly: the genuine initial segment has no earlier bias prefix.
    bad_bias = args.outdir / 'changed_prefix_HILLS_fixture'
    original = (source / 'HILLS_r_eq_metad').read_text()
    bad_bias.write_text('# modified parent prefix\n' + original)
    prefix_report = args.outdir / 'NEGATIVE_changed_bias_prefix.json'
    error = ''
    try:
        gate.bias_check(bad_bias, 'hills', 10000, time_scale, source / 'HILLS_r_eq_metad')
    except ValueError as failure:
        error = str(failure)
    rejected = 'parent byte prefix modified' in error
    item = {'case': 'changed_bias_prefix', 'expected_error': 'parent byte prefix modified',
            'actual_error': error, 'rejected_by_expected_gate': rejected,
            'evidence': str(prefix_report), 'fixture_scope': 'direct continuation prefix gate'}
    gate.dump(prefix_report, item)
    outcomes.append(item)
    bad_bias.unlink()
    print('changed_bias_prefix: ' + ('PASS' if rejected else 'FAIL'), flush=True)
    result = {'status': 'PASS' if all(item['rejected_by_expected_gate'] for item in outcomes) else 'FAIL',
              'cases': outcomes, 'source_files_modified': False,
              'scope': 'Negative technical integrity checks on altered copies; no simulation or scientific verdict.'}
    gate.dump(args.outdir / 'NEGATIVE_TESTS.json', result)
    gate.demand(result['status'] == 'PASS', 'negative integrity test failed')


if __name__ == '__main__':
    main()
