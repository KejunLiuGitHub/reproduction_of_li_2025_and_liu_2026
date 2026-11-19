#!/usr/bin/env python3
"""Describe a published Liu394 trajectory using the original parsing functions."""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
from _common import ROOT, fresh_directory, sha256, source_path, write_json

CORE = ROOT / 'campaign/repl_liu2026/iface'
sys.path.insert(0, str(CORE))
import analyze_iface as iface
import analyze_unbiased as ub

EXPECTED_ELEMENTS = ['O', 'H', 'H'] * 128 + ['C', 'O', 'H', 'H', 'S', 'O', 'O', 'O', 'H', 'K']
CORE_NAMES = ('analyze_cs_coupling.py', 'analyze_unbiased.py', 'analyze_iface.py', 'iface_rc_common.py')
FIXED_INDICES = {'C': 385, 'O_carbonyl': 386, 'S': 389, 'O_hydroxyl': 390, 'H_hydroxyl': 393, 'K': 394}


def describe(values):
    values = [value for value in values if value is not None]
    if not values:
        return {'n': 0, 'first': None, 'last': None, 'min': None, 'max': None, 'mean': None, 'median': None}
    return {'n': len(values), 'first': values[0], 'last': values[-1], 'min': min(values),
            'max': max(values), 'mean': statistics.fmean(values), 'median': statistics.median(values)}


def read_case(case_id, manifest):
    data = json.loads(manifest.read_text())
    if case_id not in data['cases']:
        raise ValueError('Unknown Liu case: ' + case_id)
    case = data['cases'][case_id]
    if case.get('id', case_id) != case_id:
        raise ValueError('Case id does not match its dictionary key')
    for field in ('trajectory', 'atom_count', 'box_angstrom', 'biased', 'history', 'status'):
        if field not in case:
            raise ValueError('Missing case field: ' + field)
    if case['atom_count'] != 394 or case['box_angstrom'] != [17, 17, 47]:
        raise ValueError('The original density function requires 394 atoms and the 17 x 17 x 47 A cell')
    if type(case['biased']) is not bool and case['biased'] != 'mixed':
        raise ValueError('biased must be a boolean or the string mixed')
    return case


def stage_for_step(case, step):
    stages = case.get('stages', case.get('bias_segments', []))
    if not stages:
        if case['biased'] == 'mixed':
            raise ValueError('A mixed-bias case must declare non-overlapping stages')
        return {'biased': case['biased'], 'role': 'trajectory'}
    for stage in stages:
        if (type(stage.get('first_step')) is not int or type(stage.get('last_step')) is not int
                or stage['first_step'] > stage['last_step'] or type(stage.get('biased')) is not bool):
            raise ValueError('Each stage requires integer bounds and a boolean biased field')
    matches = [stage for stage in stages if stage['first_step'] <= step <= stage['last_step']]
    if len(matches) != 1:
        raise ValueError('Exactly one bias stage must cover saved step ' + str(step))
    stage = matches[0]
    if type(case['biased']) is bool and stage['biased'] != case['biased']:
        raise ValueError('Stage bias contradicts the case bias flag')
    return stage


def check_raw_frames(path):
    """Check every frame, including any duplicate that the core parser discards."""
    frames = iface.read_frames(path)
    lines = path.read_text().splitlines()
    previous = {}
    duplicate_details = []
    raw_steps = []
    for frame in frames:
        if [atom[0] for atom in frame.atoms] != EXPECTED_ELEMENTS:
            raise ValueError('Atom order/composition mismatch at XYZ line ' + str(frame.start_line))
        if not all(math.isfinite(value) for atom in frame.atoms for value in atom[1:]):
            raise ValueError('Non-finite coordinate at XYZ line ' + str(frame.start_line))
        step, timestamp = ub.parse_comment(lines[frame.start_line], path, frame.start_line + 1)
        raw_steps.append(step)
        if step in previous:
            prior_time, prior_atoms = previous[step]
            duplicate_details.append({'step': step, 'same_time': prior_time == timestamp,
                                      'identical_coordinates': prior_atoms == frame.atoms,
                                      'kept_occurrence': 'last in file'})
        previous[step] = (timestamp, frame.atoms)
    return {'raw_frames': len(frames), 'duplicate_occurrences': duplicate_details,
            'raw_steps_nondecreasing': all(b >= a for a, b in zip(raw_steps, raw_steps[1:])),
            'atom_count': 394, 'all_atoms_finite': True, 'full_element_order_matches': True}


def mediator_map(case):
    """Use explicit atom identities only; no nearest-water identity inference."""
    mapping = case.get('atom_map')
    if mapping is None:
        return None
    if not isinstance(mapping, dict):
        raise ValueError('atom_map must be a dictionary of 1-based atom indices')
    keys = ('water_O', 'water_H1', 'water_H2')
    present = [mapping.get(key) for key in keys]
    if all(value is None for value in present):
        return None
    if not all(type(value) is int and 1 <= value <= 384 for value in present):
        raise ValueError('A mediator map requires all three water atom indices')
    oxygen, hydrogen1, hydrogen2 = present
    if oxygen % 3 != 1 or {hydrogen1, hydrogen2} != {oxygen + 1, oxygen + 2}:
        raise ValueError('Mediator identities must be one original O,H,H water triplet')
    return dict(zip(keys, present))


def measure_frame(record, box, stage, mediator, temperature):
    atoms = record.atoms

    def distance(i, j):
        return ub.periodic_distance(atoms[i - 1][1:4], atoms[j - 1][1:4], box)

    nearest_h = min((distance(386, index + 1), index + 1)
                    for index, atom in enumerate(atoms) if atom[0] == 'H')
    row = {'step': record.step, 'time_fs': record.time, 'stage': stage.get('role', 'trajectory'),
           'biased': stage['biased'], 'SC_A': record.sc, 'CO_A': distance(385, 386),
           'O390_H393_A': record.oh, 'O386_H393_A': distance(386, 393),
           'O386_nearest_H_A': nearest_h[0], 'nearest_H_index_1based': nearest_h[1],
           'd2_water_OH_A': None, 'd3_waterH_O386_A': None, 'd5_waterO_H393_A': None,
           'temperature_K': temperature}
    if mediator:
        row.update(d2_water_OH_A=distance(mediator['water_O'], mediator['water_H1']),
                   d3_waterH_O386_A=distance(mediator['water_H1'], 386),
                   d5_waterO_H393_A=distance(mediator['water_O'], 393))
    for key in ('gds_lower_A', 'gds_upper_A', 'slab_thickness_A', 'slab_center_A'):
        row[key] = None
    for name in ('C385', 'S389', 'O386', 'K394'):
        row[name + '_signed_GDS_A'] = None
        row[name + '_on_vapour_side'] = None
    row['gds_error'] = ''
    try:
        metrics = iface.frame_metrics(atoms).result
    except RuntimeError as error:
        row['gds_error'] = str(error)
    else:
        row.update(gds_lower_A=metrics['z_lower_A'], gds_upper_A=metrics['z_upper_A'],
                   slab_thickness_A=metrics['thickness_A'], slab_center_A=metrics['center_A'])
        for name, value in metrics['solutes'].items():
            row[name + '_signed_GDS_A'] = value['signed_distance_to_gds']
            row[name + '_on_vapour_side'] = value['on_vapour_side']
    return row


def summarize_rows(rows):
    keys = ('SC_A', 'CO_A', 'O390_H393_A', 'O386_H393_A', 'O386_nearest_H_A',
            'd2_water_OH_A', 'd3_waterH_O386_A', 'd5_waterO_H393_A', 'temperature_K')
    summary = {key: describe([row[key] for row in rows]) for key in keys}
    summary['signed_GDS'] = {}
    for name in ('C385', 'S389', 'O386', 'K394'):
        observed = [row[name + '_on_vapour_side'] for row in rows
                    if row[name + '_on_vapour_side'] is not None]
        summary['signed_GDS'][name] = {
            'distance_A': describe([row[name + '_signed_GDS_A'] for row in rows]),
            'vapour_side_frames': sum(observed), 'resolved_frames': len(observed),
            'vapour_side_fraction': sum(observed) / len(observed) if observed else None}
    return summary


def run(case_id, out, manifest):
    case = read_case(case_id, manifest)
    sources = {'dataset_manifest': manifest.resolve(), 'trajectory': source_path(case['trajectory'])}
    for key in ('energy', 'input', 'coordinates'):
        if case.get(key):
            sources[key] = source_path(case[key])
    for name in CORE_NAMES:
        sources['code/' + name] = CORE / name
    sources['code/reproduce_liu_case.py'] = Path(__file__).resolve()
    sources['code/_common.py'] = ROOT / 'tools/_common.py'
    sources['code_source_manifest'] = ROOT / 'metadata/LIU2026_SOURCE_FILES.json'
    before = {key: sha256(path) for key, path in sources.items()}
    expected = case.get('sha256', {})
    if not isinstance(expected, dict) or not expected.get('trajectory'):
        raise ValueError('Case sha256 must specify the trajectory checksum')
    for key, value in expected.items():
        if key in before and before[key] != value:
            raise ValueError('Dataset checksum mismatch for ' + key)
    for item in json.loads(sources['code_source_manifest'].read_text())['files']:
        if sha256(source_path(item['path'])) != item['sha256']:
            raise ValueError('Original core source checksum mismatch: ' + item['path'])
    out = fresh_directory(out)
    raw_check = check_raw_frames(sources['trajectory'])
    records, deduplication = ub.load_trajectory(sources['trajectory'], case['box_angstrom'])
    original_frames = len(records)
    records = [record for record in records
               if case.get('first_step', records[0].step) <= record.step <= case.get('last_step', records[-1].step)]
    if not records:
        raise ValueError('The case step window contains no frames')
    mediator = mediator_map(case)
    energy_records, energy_info = (ub.read_energy(sources['energy']) if 'energy' in sources else ([], None))
    by_step = {record.step: record for record in energy_records}
    rows = []
    temperature_mismatch = []
    for record in records:
        energy = by_step.get(record.step)
        temperature = None
        if energy:
            if abs(energy.time - record.time) <= 1e-6:
                temperature = energy.temp
            else:
                temperature_mismatch.append({'step': record.step, 'XYZ_time_fs': record.time,
                                             'energy_time_fs': energy.time})
        rows.append(measure_frame(record, case['box_angstrom'], stage_for_step(case, record.step), mediator, temperature))
    intervals = [{'previous_step': a.step, 'step': b.step, 'delta_step': b.step - a.step,
                  'delta_time_fs': b.time - a.time} for a, b in zip(records, records[1:])]
    coverage = {'frames': len(records), 'first_step': records[0].step, 'last_step': records[-1].step,
                'first_time_fs': records[0].time, 'last_time_fs': records[-1].time,
                'observed_span_fs': records[-1].time - records[0].time,
                'distinct_step_intervals': sorted({item['delta_step'] for item in intervals}),
                'distinct_time_intervals_fs': sorted({round(item['delta_time_fs'], 9) for item in intervals}),
                'intervals_other_than_10fs': [item for item in intervals if item['delta_step'] != 10 or abs(item['delta_time_fs'] - 10) > 1e-6],
                'source_unique_frames': original_frames, 'frames_excluded_by_case_window': original_frames - len(records),
                'deduplication': deduplication, 'raw_frame_checks': raw_check,
                'energy': energy_info, 'temperature_matched_frames': sum(row['temperature_K'] is not None for row in rows),
                'temperature_time_mismatches': temperature_mismatch,
                'gds_resolved_frames': sum(not row['gds_error'] for row in rows),
                'gds_unresolved': [{'step': row['step'], 'reason': row['gds_error']} for row in rows if row['gds_error']],
                'note': 'Coverage is the observed saved-frame span, not an inferred simulation duration.'}
    groups = {}
    for row in rows:
        key = (row['stage'], row['biased'])
        groups.setdefault(key, []).append(row)
    stages = [{'role': role, 'biased': biased, 'frames': len(group),
               'first_step': group[0]['step'], 'last_step': group[-1]['step'],
               'observed_span_fs': group[-1]['time_fs'] - group[0]['time_fs'],
               'metrics': summarize_rows(group)} for (role, biased), group in groups.items()]
    summary = {'schema': 'liu394-descriptive-replay-v1', 'case_id': case_id,
               'label': case.get('label', case_id), 'job_ids': case.get('job_ids', [case['job_id']] if case.get('job_id') else []),
               'status_as_recorded': case['status'], 'biased': case['biased'], 'history': case['history'],
               'prepared_under_bias': case.get('prepared_under_bias'),
               'controlled_initial_conditions': case.get('controlled_initial_conditions'),
               'source_segment_biased': case.get('source_segment_biased'),
               'atom_count': 394, 'box_angstrom': case['box_angstrom'],
               'atom_indices_1based': FIXED_INDICES, 'mediator_map_1based': mediator,
               'coverage': coverage, 'metrics': summarize_rows(rows), 'stages': stages,
               'methods': {'distance': 'Orthorhombic minimum image on all three axes; original analyze_unbiased.periodic_distance.',
                           'nearest_H': 'Nearest of all 259 hydrogen atoms to O386; a geometric descriptor with the atom index retained.',
                           'GDS': 'Original analyze_iface.frame_metrics: 0.5 A water-oxygen density bins, three-bin moving average, half maximum. Positive is vapour; negative is inside the water slab. An operational surface estimate.',
                           'temperature': 'Original read_energy; exact step and matching time, with no interpolation.',
                           'water_coordination': 'Not calculated: the four archived core files do not define a reusable water-coordination metric.',
                           'scope': 'Saved-frame distances and positions; no bond, reaction, rate or reproduction verdict is generated.'}}
    write_json(out / 'CASE_MANIFEST.json', case)
    write_json(out / 'SUMMARY.json', summary)
    with (out / 'FRAME_METRICS.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    after = {key: sha256(path) for key, path in sources.items()}
    if after != before:
        raise ValueError('A source file changed during replay')
    provenance = {'command': ['python3', 'tools/reproduce_liu_case.py', '--case', case_id, '--out', str(out)],
                  'python_version': sys.version, 'sources': {
                      key: {'path': str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                            'sha256_before': before[key], 'sha256_after': after[key]}
                      for key, path in sources.items()}, 'source_bytes_unchanged': True,
                  'original_core_commit': json.loads(sources['code_source_manifest'].read_text())['original_commit']}
    write_json(out / 'PROVENANCE.json', provenance)
    write_json(out / 'OUTPUT_SHA256.json', {path.name: sha256(path) for path in sorted(out.iterdir()) if path.is_file()})
    print(json.dumps({'case': case_id, 'frames': len(rows), 'observed_span_fs': coverage['observed_span_fs'],
                      'SC_mean_A': summary['metrics']['SC_A']['mean'], 'out': str(out)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True, help='Key in data/LIU2026_DATASETS.json')
    parser.add_argument('--out', type=Path, required=True, help='A new output directory')
    args = parser.parse_args()
    try:
        run(args.case, args.out, ROOT / 'data/LIU2026_DATASETS.json')
    except (ValueError, RuntimeError, OSError, KeyError) as error:
        parser.exit(2, str(error) + '\n')


if __name__ == '__main__':
    main()
