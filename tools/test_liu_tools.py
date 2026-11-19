#!/usr/bin/env python3
"""Exercise the Liu replay in an isolated copy using one released trajectory."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PortableLiuTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='liu-replay-test-')
        cls.repo = Path(cls.temp.name) / 'repository'
        cls.repo.mkdir()
        paths = ['tools/_common.py', 'tools/reproduce_liu_case.py', 'metadata/LIU2026_SOURCE_FILES.json']
        paths += ['campaign/repl_liu2026/iface/' + name for name in
                  ['analyze_cs_coupling.py', 'analyze_unbiased.py', 'analyze_iface.py', 'iface_rc_common.py']]
        for rel in paths:
            target = cls.repo / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / rel, target)
        catalog = json.loads((ROOT / 'data/LIU2026_DATASETS.json').read_text())
        cls.real_case = dict(catalog['cases'][CASE])
        for key in ('trajectory', 'energy', 'input', 'coordinates'):
            if cls.real_case.get(key):
                rel = cls.real_case[key]
                target = cls.repo / rel
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(ROOT / rel, target)
        cls.xyz = cls.repo / cls.real_case['trajectory']
        lines = cls.xyz.read_text().splitlines(keepends=True)
        if len(lines) < 792 or lines[0].strip() != '394':
            raise RuntimeError('The test case needs at least two complete Liu394 frames')
        cls.blocks = [''.join(lines[i:i + 396]) for i in range(0, len(lines), 396) if lines[i].strip()]
        cls.initial_hash = digest(cls.xyz)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def invoke(self, case, blocks=None, expect=0, energy_text=None, out=None):
        test_id = self.id().split('.')[-1]
        case = dict(case)
        case['id'] = test_id
        if blocks is not None:
            xyz = self.repo / ('data/' + test_id + '.xyz')
            xyz.write_text(''.join(blocks))
            case['trajectory'] = str(xyz.relative_to(self.repo))
            case['sha256'] = {'trajectory': digest(xyz)}
        if energy_text is not None:
            energy = self.repo / ('data/' + test_id + '.ener')
            energy.write_text(energy_text)
            case['energy'] = str(energy.relative_to(self.repo))
            case.setdefault('sha256', {})['energy'] = digest(energy)
        (self.repo / 'data/LIU2026_DATASETS.json').write_text(json.dumps({'cases': {test_id: case}}))
        out = out or Path(self.temp.name) / test_id
        proc = subprocess.run([sys.executable, str(self.repo / 'tools/reproduce_liu_case.py'),
                               '--case', test_id, '--out', str(out)], cwd=self.temp.name,
                              text=True, capture_output=True)
        self.assertEqual(proc.returncode, expect, proc.stdout + proc.stderr)
        self.assertEqual(digest(self.xyz), self.initial_hash)
        return out, proc

    def fixture(self):
        return {'atom_count': 394, 'box_angstrom': [17, 17, 47], 'biased': False,
                'history': 'Synthetic parser test based on released coordinates.', 'status': 'Test fixture'}

    def first_step(self):
        import re
        return int(re.search(r'\bi\s*=\s*(\d+)', self.blocks[0].splitlines()[1]).group(1))

    def test_full_trajectory_and_independent_distances(self):
        out, _ = self.invoke(self.real_case)
        summary = json.loads((out / 'SUMMARY.json').read_text())
        with (out / 'FRAME_METRICS.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        by_step = {}
        import re
        for block in self.blocks:
            lines = block.splitlines()
            step = int(re.search(r'\bi\s*=\s*(\d+)', lines[1]).group(1))
            by_step[step] = [tuple(map(float, line.split()[1:4])) for line in lines[2:]]
        self.assertEqual(len(rows), summary['coverage']['frames'])
        for row in rows:
            coordinates = by_step[int(row['step'])]
            for field, a, b in [('SC_A', 388, 384), ('CO_A', 384, 385)]:
                delta = [coordinates[a][i] - coordinates[b][i] for i in range(3)]
                expected = math.sqrt(sum(((x + length / 2) % length - length / 2) ** 2
                                         for x, length in zip(delta, [17, 17, 47])))
                self.assertAlmostEqual(float(row[field]), expected, places=12)
        self.assertTrue(json.loads((out / 'PROVENANCE.json').read_text())['source_bytes_unchanged'])
        self.assertNotIn('verdict', summary)

    def test_existing_output_refused(self):
        out = Path(self.temp.name) / 'existing'
        out.mkdir()
        marker = out / 'keep.txt'
        marker.write_text('Preserve this file.')
        self.invoke(self.real_case, out=out, expect=2)
        self.assertEqual(marker.read_text(), 'Preserve this file.')

    def test_checksum_mismatch(self):
        case = dict(self.real_case)
        case['sha256'] = {'trajectory': '0' * 64}
        _, proc = self.invoke(case, expect=2)
        self.assertIn('checksum mismatch', proc.stderr)

    def test_full_atom_identity(self):
        lines = self.blocks[0].splitlines(keepends=True)
        lines[2] = lines[2].replace('O', 'N', 1)
        _, proc = self.invoke(self.fixture(), [''.join(lines)], expect=2)
        self.assertIn('Atom order/composition', proc.stderr)

    def test_nonfinite_coordinate(self):
        lines = self.blocks[0].splitlines(keepends=True)
        tokens = lines[2].split()
        tokens[1] = 'nan'
        lines[2] = ' '.join(tokens) + '\n'
        _, proc = self.invoke(self.fixture(), [''.join(lines)], expect=2)
        self.assertIn('Non-finite', proc.stderr)

    def test_truncated_frame(self):
        _, proc = self.invoke(self.fixture(), ['\n'.join(self.blocks[0].splitlines()[:-1]) + '\n'], expect=2)
        self.assertIn('truncated frame', proc.stderr)

    def test_duplicate_endpoint_recorded(self):
        out, _ = self.invoke(self.fixture(), self.blocks[:2] + [self.blocks[1]])
        coverage = json.loads((out / 'SUMMARY.json').read_text())['coverage']
        self.assertEqual(coverage['frames'], 2)
        self.assertEqual(coverage['deduplication']['duplicate_frames'], 1)
        self.assertTrue(coverage['raw_frame_checks']['duplicate_occurrences'][0]['identical_coordinates'])

    def test_different_duplicate_last_occurrence(self):
        lines = self.blocks[1].splitlines(keepends=True)
        tokens = lines[386].split()
        tokens[1] = str(float(tokens[1]) + 0.25)
        lines[386] = ' '.join(tokens) + '\n'
        out, _ = self.invoke(self.fixture(), self.blocks[:2] + [''.join(lines)])
        summary = json.loads((out / 'SUMMARY.json').read_text())
        self.assertEqual(summary['coverage']['frames'], 2)
        self.assertFalse(summary['coverage']['raw_frame_checks']['duplicate_occurrences'][0]['identical_coordinates'])

    def test_mixed_stages_preserved(self):
        case = self.fixture()
        first = self.first_step()
        case.update(biased='mixed', stages=[{'first_step': first, 'last_step': first, 'biased': True, 'role': 'A'},
                                         {'first_step': first + 1, 'last_step': first + 10, 'biased': False, 'role': 'B'}])
        out, _ = self.invoke(case, self.blocks[:2])
        summary = json.loads((out / 'SUMMARY.json').read_text())
        self.assertEqual([stage['biased'] for stage in summary['stages']], [True, False])
        self.assertEqual(summary['biased'], 'mixed')

    def test_ambiguous_stage_refused(self):
        case = self.fixture()
        first = self.first_step()
        case.update(biased='mixed', stages=[{'first_step': first, 'last_step': first + 10, 'biased': True},
                                         {'first_step': first, 'last_step': first + 10, 'biased': False}])
        _, proc = self.invoke(case, self.blocks[:2], expect=2)
        self.assertIn('Exactly one bias stage', proc.stderr)

    def test_temperature_time_mismatch_not_guessed(self):
        first = self.first_step()
        out, _ = self.invoke(self.fixture(), self.blocks[:1], energy_text=f'{first} 999999 0.4 300.0\n')
        summary = json.loads((out / 'SUMMARY.json').read_text())
        self.assertEqual(summary['metrics']['temperature_K']['n'], 0)
        self.assertEqual(len(summary['coverage']['temperature_time_mismatches']), 1)

    def test_no_implicit_mediator_map(self):
        out, _ = self.invoke(self.fixture(), self.blocks[:1])
        summary = json.loads((out / 'SUMMARY.json').read_text())
        self.assertIsNone(summary['mediator_map_1based'])
        self.assertEqual(summary['metrics']['d2_water_OH_A']['n'], 0)
        self.assertEqual(summary['metrics']['O386_H393_A']['n'], 1)

    def test_unresolved_GDS_retains_frame(self):
        lines = self.blocks[0].splitlines(keepends=True)
        for i in range(2, 386):
            tokens = lines[i].split()
            tokens[3] = '0.0'
            lines[i] = ' '.join(tokens) + '\n'
        out, _ = self.invoke(self.fixture(), [''.join(lines)])
        summary = json.loads((out / 'SUMMARY.json').read_text())
        self.assertEqual(summary['coverage']['frames'], 1)
        self.assertEqual(summary['coverage']['gds_resolved_frames'], 0)
        self.assertEqual(summary['metrics']['SC_A']['n'], 1)

    def test_wrong_cell_refused(self):
        case = self.real_case.copy()
        case['box_angstrom'] = [15, 15, 30]
        _, proc = self.invoke(case, expect=2)
        self.assertIn('17 x 17 x 47', proc.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True, help='A released trajectory with at least two saved frames')
    args, remainder = parser.parse_known_args()
    CASE = args.case
    unittest.main(argv=[sys.argv[0]] + remainder)
