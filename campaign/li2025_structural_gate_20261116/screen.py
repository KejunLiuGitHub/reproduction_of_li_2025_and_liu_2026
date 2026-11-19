#!/usr/bin/env python3
"""Retrospective Li296 geometric screening; never an electronic bond verdict.

Python >=3.9 + numpy. Source files are read only. All parser dependencies are
hash checked. CLI consumes an explicit source manifest; no auto-discovery.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PARSER = REPO / 'campaign/li2025_hbnetwork_20261116/production_v1/validate_stage.py'
PARSER_SHA = '3df8173e5b37eb2d1aba2215026710ae1e868d60de3d29090fed8f00bfcc282a'
ORDER = ['H', 'H', 'O'] * 96 + ['C', 'O', 'H', 'H', 'S', 'O', 'O', 'O']
C, OC, S = 288, 289, 292
WATER_O = np.arange(2, 288, 3)
WATER_H = np.array([i for i in range(288) if i % 3 != 2])
ALL_H = np.array([i for i, e in enumerate(ORDER) if e == 'H'])
ALL_O = np.array([i for i, e in enumerate(ORDER) if e == 'O'])


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load_parser():
    if sha(PARSER) != PARSER_SHA:
        raise ValueError('Frozen technical parser SHA mismatch')
    spec = importlib.util.spec_from_file_location('li296_frozen_parser', PARSER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def geometry(coords, box, guard):
    """Same MIC/pyramidal-height definition as registered analyze_solvation.

    Fixed-parent checks deliberately do not infer molecular identities anew.
    Unexpected short distances are ambiguity flags, not newly assigned bonds.
    """
    x = np.asarray(coords, dtype=float)
    box = np.asarray(box, dtype=float)
    d = x[:, None, :] - x[None, :, :]
    d -= box * np.rint(d / box)
    r = np.sqrt(np.einsum('ijk,ijk->ij', d, d))
    v = -d[C, [OC, 290, 291, S]]
    normal = np.cross(v[1] - v[0], v[2] - v[0])
    norm = np.linalg.norm(normal)
    if norm < 1.e-12:
        raise ValueError('degenerate Oc/H/H plane at carbon')
    height = float(abs(np.dot(v[0], normal)) / norm)
    lengths = np.linalg.norm(v, axis=1)
    if min(lengths) < 1.e-12:
        raise ValueError('zero carbon-neighbour distance')
    angles = [float(np.degrees(np.arccos(np.clip(np.dot(v[i], v[j]) /
              (lengths[i] * lengths[j]), -1, 1)))) for i in range(4) for j in range(i + 1, 4)]
    issues = []

    def flag(name, pairs, lower=None, upper=None):
        for a, b in pairs:
            value = float(r[a, b])
            bad = (lower is not None and value < lower) or (upper is not None and value > upper)
            if bad:
                issues.append({'kind': name, 'atoms_1based': [int(a + 1), int(b + 1)], 'distance_A': value})

    parent_pairs = [(int(h), int(h - h % 3 + 2)) for h in WATER_H]
    flag('parent_water_OH_outside_guard', parent_pairs, guard['water_OH_min'], guard['water_OH_max'])
    flag('parent_CH_outside_guard', [(C, 290), (C, 291)], guard['CH_min'], guard['CH_max'])
    flag('parent_SO_outside_guard', [(S, i) for i in (293, 294, 295)], guard['SO_min'], guard['SO_max'])
    flag('parent_CO_outside_guard', [(C, OC)], guard['CO_min'], guard['CO_max'])

    def short_contacts(name, aa, bb, threshold, allowed=frozenset(), triangle=False):
        pairs = np.argwhere(r[np.ix_(aa, bb)] < threshold)
        for ii, jj in pairs:
            a, b = int(aa[ii]), int(bb[jj])
            if (triangle and a >= b) or (a, b) in allowed:
                continue
            issues.append({'kind': name, 'atoms_1based': [a + 1, b + 1], 'distance_A': float(r[a, b])})

    short_contacts('unexpected_OH_close_contact', ALL_H, ALL_O, guard['unexpected_OH_max_exclusive'], frozenset(parent_pairs))
    short_contacts('unexpected_CO_close_contact', [C], ALL_O, guard['unexpected_CO_max_exclusive'], frozenset([(C, OC)]))
    short_contacts('unexpected_SO_close_contact', [S], ALL_O, guard['unexpected_SO_max_exclusive'], frozenset([(S, i) for i in (293, 294, 295)]))
    short_contacts('unexpected_CH_close_contact', [C], ALL_H, guard['unexpected_CH_max_exclusive'], frozenset([(C, 290), (C, 291)]))
    short_contacts('unexpected_SH_close_contact', [S], ALL_H, guard['unexpected_SH_max_exclusive'])
    short_contacts('OO_collision', ALL_O, ALL_O, guard['OO_collision_max_exclusive'], triangle=True)
    short_contacts('HH_collision', ALL_H, ALL_H, guard['HH_collision_max_exclusive'], triangle=True)
    return {'SC_A': float(r[C, S]), 'CO_A': float(r[C, OC]),
            'CH_A': [float(r[C, j]) for j in (290, 291)],
            'SO_A': [float(r[S, j]) for j in (293, 294, 295)],
            'carbon_pyramidal_height_A': height, 'carbon_six_angles_deg': angles,
            'carbon_angle_min_deg': min(angles), 'carbon_angle_max_deg': max(angles),
            'parent_water_OH_max_A': max(float(r[a, b]) for a, b in parent_pairs),
            'connectivity_guard_clear': not issues, 'connectivity_ambiguities': issues}


def candidate_flags(m, cut):
    return {
        'SC': m['SC_A'] < cut['SC_max_exclusive'],
        'CO': cut['CO_min'] <= m['CO_A'] <= cut['CO_max'],
        'pyramidal_height': cut['pyramidal_height_min'] <= m['carbon_pyramidal_height_A'] <= cut['pyramidal_height_max'],
        'four_neighbour_angles': cut['C_angle_min'] <= m['carbon_angle_min_deg'] and m['carbon_angle_max_deg'] <= cut['C_angle_max'],
        'connectivity': m['connectivity_guard_clear']}


def continuous_segments(rows, flags, stride):
    """Split on every false frame or non-exact gap. No interpolation/dedup."""
    segments, start = [], None
    for i, yes in enumerate(flags):
        joined = i == 0 or abs(rows[i]['time_fs'] - rows[i - 1]['time_fs'] - stride) < 1.e-6
        if start is not None and (not yes or not joined):
            segments.append((start, i - 1))
            start = None
        if yes and start is None:
            start = i
    if start is not None:
        segments.append((start, len(rows) - 1))
    return [{'first_index': a, 'last_index': b, 'frames': b - a + 1,
             'start_time_fs': rows[a]['time_fs'], 'end_time_fs': rows[b]['time_fs'],
             'span_fs': rows[b]['time_fs'] - rows[a]['time_fs'],
             'left_censored': a == 0, 'right_censored': b == len(rows) - 1}
            for a, b in segments]


def coverage(rows, cfg, start_step, static=False):
    if static:
        return {'passed': False, 'status': 'NOT_APPLICABLE_STATIC_SNAPSHOT', 'frames': len(rows)}
    steps = [r['step'] for r in rows]
    times = [r['time_fs'] for r in rows]
    expected = list(range(start_step, start_step + int(cfg['duration_fs']) + 1, int(cfg['output_stride_fs'])))
    checks = {'full_step_sequence': steps == expected,
              'time_matches_step_1fs': all(abs(t - s) < 1.e-5 for t, s in zip(times, steps)),
              'observed_span_at_least_5ps': times[-1] - times[0] >= cfg['duration_fs'] - 1.e-6,
              'no_duplicate_or_irregular_gap': all(abs(b - a - cfg['output_stride_fs']) < 1.e-6 for a, b in zip(times, times[1:]))}
    return {'passed': all(checks.values()), 'status': 'PASS' if all(checks.values()) else 'UNDECIDED_INCOMPLETE_COVERAGE',
            'checks': checks, 'frames': len(rows), 'first_step': steps[0], 'last_step': steps[-1],
            'observed_span_fs': times[-1] - times[0], 'required_start_step': start_step}


def classify(rows, cfg, context, covered, technical, static=False):
    """Classification is derived once from default; sensitivity never replaces it."""
    first = rows[0]
    react = cfg['initial_reactant']
    if context.get('public_P_label') or first['SC_A'] <= react['SC_min_exclusive']:
        kind = 'RETENTION_CONTEXT_INITIAL_SHORT_OR_PUBLIC_P'
        initial_ok = False
    else:
        end = first['time_fs'] + cfg['initial_reactant_span_fs']
        initial = [r for r in rows if r['time_fs'] <= end + 1.e-6]
        initial_ok = (abs(initial[-1]['time_fs'] - end) < 1.e-6 and
                      all(r['SC_A'] > react['SC_min_exclusive'] and r['CO_A'] <= react['CO_max'] and
                          r['carbon_pyramidal_height_A'] <= react['pyramidal_height_max'] and
                          r['connectivity_guard_clear'] for r in initial) and
                      all(abs(b['time_fs'] - a['time_fs'] - cfg['output_stride_fs']) < 1.e-6 for a, b in zip(initial, initial[1:])))
        kind = 'FORMATION_CONTEXT_CONFIRMED_REACTANT_WINDOW' if initial_ok else 'UNRESOLVED_INITIAL_STATE'
    schemes = {}
    for label in ('strict', 'default', 'loose'):
        flags = [all(candidate_flags(r, cfg[label]).values()) for r in rows]
        segments = continuous_segments(rows, flags, cfg['output_stride_fs'])
        sustained = [s for s in segments if s['span_fs'] >= cfg['sustained_span_fs'] - 1.e-6]
        # Formation candidate window must not share the reactant window endpoint.
        eligible = [s for s in sustained if not initial_ok or s['start_time_fs'] > first['time_fs'] + cfg['initial_reactant_span_fs'] + 1.e-6]
        terminal = next((s for s in eligible if s['right_censored']), None)
        schemes[label] = {'candidate_frames': sum(flags), 'fraction': sum(flags) / len(flags),
                          'segments': segments, 'sustained_segments': sustained,
                          'eligible_sustained_segments': eligible, 'terminal_sustained_segment': terminal,
                          'terminal_geometry_candidate': terminal is not None,
                          'criterion_failure_counts': {key: sum(not candidate_flags(r, cfg[label])[key] for r in rows)
                                                       for key in candidate_flags(rows[0], cfg[label])}}
    d = schemes['default']
    ambiguities = sum(not r['connectivity_guard_clear'] for r in rows)
    if static:
        status = 'STATIC_GEOMETRY_CANDIDATE' if all(candidate_flags(first, cfg['default']).values()) else 'STATIC_GEOMETRY_OUTSIDE_DEFAULT_SCREEN'
    elif not covered:
        status = 'UNDECIDED_INCOMPLETE_COVERAGE'
    elif not technical:
        status = 'UNRESOLVED_TECHNICAL_PROVENANCE'
    elif ambiguities:
        status = 'UNRESOLVED_CONNECTIVITY_AMBIGUITY'
    elif kind == 'UNRESOLVED_INITIAL_STATE':
        status = 'UNRESOLVED_INITIAL_STATE'
    elif d['terminal_geometry_candidate']:
        status = 'DEFAULT_RETENTION_GEOMETRY_CANDIDATE' if kind.startswith('RETENTION') else 'DEFAULT_FORMATION_GEOMETRY_CANDIDATE'
    else:
        status = 'UNRESOLVED_NO_TERMINAL_GEOMETRIC_SEGMENT'
    return {'status': status, 'initial_context': kind, 'initial_reactant_window_confirmed': initial_ok,
            'connectivity_ambiguous_frames': ambiguities, 'threshold_sensitivity': schemes,
            'terminal_candidate_same_under_all_thresholds': len({s['terminal_geometry_candidate'] for s in schemes.values()}) == 1,
            'sensitivity_status': 'CONSISTENT_TERMINAL_SCREEN' if len({s['terminal_geometry_candidate'] for s in schemes.values()}) == 1 else 'UNRESOLVED_THRESHOLD_SENSITIVITY',
            'makes_negative_reaction_claim': False,
            'interpretation': 'Only geometric screening. No candidate is not proof of no reaction; threshold crossings are not molecular bond breaking.'}


def technical_evidence(case, paths, frames, parser, seen):
    details, checks = {}, {}

    def record_check(name, fn):
        try:
            value = fn()
            checks[name] = True
            details[name] = value
        except (ValueError, KeyError, AssertionError, OSError) as exc:
            checks[name] = False
            details[name] = str(exc)

    def need(condition, msg):
        if not condition:
            raise ValueError(msg)

    record = json.loads(paths['run_record'].read_text()) if 'run_record' in paths else {}
    formal = json.loads(paths['checks'].read_text()) if 'checks' in paths else {}
    distance_certificate = json.loads(paths['distance_certificate'].read_text()) if 'distance_certificate' in paths else {}
    nstep = 5000
    start = case.get('start_step', 0)

    def input_check():
        text = '\n'.join(line.split('#')[0] for line in paths['input'].read_text().splitlines())
        need(re.search(r'(?mi)^\s*CHARGE\s+-2\s*$', text), 'explicit charge -2 absent')
        need(re.search(r'(?mi)^\s*MULTIPLICITY\s+1\s*$', text), 'explicit singlet absent')
        need(not re.search(r'(?mi)^\s*(?:&(?:FREE_ENERGY|METADYN|CONSTRAINT|COLLECTIVE)\b|PLUMED_INPUT_FILE\b|USE_PLUMED\b)', text), 'bias/constraint keyword present')
        need(re.search(r'(?mi)^\s*TIMESTEP\s+1(?:\.0*)?\s*$', text), '1 fs timestep absent')
        cells = re.findall(r'(?mi)^\s*ABC\s+([^\n]+)', text)
        need(len(cells) == 1, 'input ABC absent/ambiguous')
        cell = [float(n) for n in cells[0].split()]
        need(len(cell) == 3 and all(abs(a - b) < 1.e-8 for a, b in zip(cell, case['box_A'])), 'declared box differs from input ABC')
        need(not re.search(r'(?mi)^\s*(?:ALPHA_BETA_GAMMA|A|B|C)\s+[-+.0-9]', text), 'non-orthorhombic/alternative cell unsupported')
        need(record.get('input_cell_angstrom', record.get('box_angstrom', cell)) == cell, 'record cell differs from input')
        for file, expected in record.get('input_sha256', {}).items():
            p = paths['input'].parent / file
            need(p.is_file() and sha(p) == expected, 'recorded input hash mismatch: ' + file)
            seen[p] = sha(p)
        need(record.get('input_sha256', {}).get(paths['input'].name) == sha(paths['input']), 'input not bound by run_record SHA')
        if re.search(r'(?mi)^\s*&EXT_RESTART\b', text):
            need('source_restart' in paths, 'EXT_RESTART source not supplied')
            ref = re.search(r'(?mi)^\s*RESTART_FILE_NAME\s+(\S+)', text)
            need(ref and paths['source_restart'].name == ref[1], 'restart source file mismatch')
        return {'charge': -2, 'multiplicity': 1, 'unbiased': True, 'sha256': sha(paths['input'])}
    record_check('unbiased_input_and_recorded_sha', input_check)

    def terminal_record():
        need(record.get('status') == 'COMPLETED' and record.get('program_ended') is True and record.get('abort_count') == 0 and record.get('returncode') == 0, 'run_record not normal completed/zero returncode')
        need(record.get('case') == case['run_case'] and str(record.get('job_id')) == str(case['job_id']), 'case/job identity mismatch')
        need(record.get('md_steps_reported') == nstep and record.get('md_step_first') == start + 1 and record.get('md_step_last') == start + nstep, 'run_record MD count/endpoint mismatch')
        need(record.get('scf', {}).get('converged_runs') == nstep + 1 and record.get('scf', {}).get('not_converged_runs') == 0, 'run_record SCF mismatch')
        bound = formal.get('trajectory_sha256', record.get('trajectory_sha256', distance_certificate.get('trajectory_sha256')))
        need(bound == sha(paths['trajectory']), 'trajectory not bound to formal certificate SHA')
        if 'checks' in paths:
            need(formal.get('status') == 'PASS', 'formal technical certificate did not pass')
            need(bool(formal.get('checks')), 'formal checks empty')
            need(formal.get('case') == record.get('case') and str(formal.get('job_id')) == str(record.get('job_id')), 'formal/run_record case/job mismatch')
            need(all(v.get('pass') is True if isinstance(v, dict) else v is True for v in formal.get('checks', {}).values()), 'formal checks contain failure')
            for file, expected in formal.get('files_sha256', {}).items():
                p = paths['trajectory'].parent / file
                need(p.is_file() and sha(p) == expected, 'formal internal file SHA mismatch/missing: ' + file)
                seen[p] = sha(p)
        return {'job_id': record.get('job_id'), 'SCF_converged': nstep + 1, 'SCF_failed': 0,
                'trajectory_bound_to_certificate': True, 'formal_checks_present': 'checks' in paths}
    record_check('terminal_certificate_bound_to_raw', terminal_record)

    def energy_check():
        rows = parser.read_energy(paths['energy'])
        steps = [int(r[0]) for r in rows]
        allowed = [list(range(start, start + nstep + 1))]
        if start != 0:
            allowed.append(list(range(start + 1, start + nstep + 1)))
        need(steps in allowed, 'energy incomplete sequence (fresh stages require t0)')
        need(all(abs(r[1] - r[0]) < 1.e-5 for r in rows), 'energy time inconsistent')
        by_step = {int(row[0]): row for row in rows}
        need(all(f['energy'] is not None and f['step'] in by_step and abs(f['energy'] - by_step[f['step']][4]) < 1.e-8 for f in frames), 'XYZ energy disagrees with energy table')
        return {'rows': len(rows), 'first_step': steps[0], 'last_step': steps[-1], 'finite': True}
    record_check('complete_energy', energy_check)
    record_check('final_restart', lambda: parser.read_restart(paths['restart'], start + nstep, nstep, frames[-1]['coords']))

    def source_check():
        if frames[0]['step'] != start:
            raise ValueError('t0 missing: source cannot be substituted for a missing observed frame')
        source = parser.read_xyz(paths['source_xyz'])
        need(len(source) == 1, 'source must be one frame')
        delta = float(np.max(np.abs(np.asarray(source[0]['coords']) - np.asarray(frames[0]['coords']))))
        need(delta < 1.e-7, 'first frame/source mismatch')
        expected = case.get('expected_source_sha256')
        need(expected and sha(paths['source_xyz']) == expected, 'registered source SHA missing/mismatch')
        if 'parent_trajectory' in paths:
            need(sha(paths['parent_trajectory']) == case['expected_parent_trajectory_sha256'], 'parent trajectory SHA mismatch')
            parent = [f for f in parser.read_xyz(paths['parent_trajectory']) if f['step'] == case['parent_step']]
            need(len(parent) == 1, 'parent step absent/duplicated')
            need(float(np.max(np.abs(np.asarray(parent[0]['coords']) - np.asarray(source[0]['coords'])))) < 1.e-7, 'parent source coordinates mismatch')
        return {'source_sha256': expected, 'max_delta_A': delta}
    record_check('first_frame_and_source', source_check)

    if 'log' in paths:
        def raw_log():
            r = parser.read_log(paths['log'])
            need(r['steps'] == list(range(start + 1, start + nstep + 1)), 'raw log MD sequence')
            need(r['scf_bins'] == [2] + [1] * (nstep - 1) and r['scf_converged'] == nstep + 1 and r['scf_failed'] == 0 and r['scf_pending'] == 0, 'raw SCF sequence')
            need(r['program_ended_count'] == 1 and r['abort_count'] == 0, 'raw log abnormal end')
            return {k: v for k, v in r.items() if k not in ('steps', 'scf_bins')}
        record_check('raw_log_recheck', raw_log)
        level = 'RAW_LOG_PLUS_BOUND_FORMAL_CERTIFICATE'
    else:
        level = 'BOUND_FORMAL_CERTIFICATE_LOG_NOT_AVAILABLE_LOCALLY'
        details['raw_log_limit'] = 'No raw log reread is claimed; terminal MD/SCF evidence comes from the hashed run_record/certificate.'
    return {'passed': all(checks.values()), 'checks': checks, 'details': details, 'evidence_level': level,
            'failed_items': [k for k, v in checks.items() if not v]}


def analyze_case(case, out, thresholds=HERE / 'thresholds.json'):
    cfg = json.loads(Path(thresholds).read_text())
    parser = load_parser()
    box = case.get('box_A')
    if not isinstance(box, list) or len(box) != 3 or not all(isinstance(x, (int, float)) and math.isfinite(x) and x > 0 for x in box):
        raise ValueError('box_A requires three finite positive orthorhombic lengths')
    paths = {key: (REPO / value).resolve() for key, value in case['files'].items()}
    out = Path(out).resolve()
    if any(p == out or out in p.parents for p in paths.values()):
        raise ValueError('source inside output directory is forbidden')
    if out.exists() and any(out.iterdir()):
        raise ValueError('output directory is not empty: evidence overwrite refused')
    out.mkdir(parents=True, exist_ok=True)
    seen = {p: sha(p) for p in paths.values()}
    seen[PARSER] = sha(PARSER)
    seen[Path(thresholds).resolve()] = sha(thresholds)
    seen[Path(__file__).resolve()] = sha(__file__)
    for key, expected in case.get('expected_file_sha256', {}).items():
        if sha(paths[key]) != expected:
            raise ValueError('registered file SHA mismatch: ' + key)
    frames = parser.read_xyz(paths['trajectory'])
    static = case.get('static', False)
    if static and len(frames) != 1:
        raise ValueError('static mode requires exactly one frame')
    if static:
        header_box = re.search(r'\bbox=([0-9.,]+)', frames[0]['raw'].splitlines()[1])
        if not header_box or len(header_box[1].split(',')) != 3 or any(abs(a - float(b)) > 1.e-7 for a, b in zip(box, header_box[1].split(','))):
            raise ValueError('static snapshot box is absent or differs from exported header')
    rows = []
    for i, frame in enumerate(frames):
        if not static and (frame['step'] is None or frame['time'] is None):
            raise ValueError('trajectory step/time absent')
        metrics = geometry(frame['coords'], case['box_A'], cfg['connectivity_guard'])
        rows.append({'index': i, 'step': frame['step'], 'time_fs': frame['time'] if not static else 0.0, **metrics})
    cov = coverage(rows, cfg, case.get('start_step', 0), static)
    tech = {'passed': False, 'status': 'NOT_APPLICABLE_STATIC_REFERENCE'} if static else technical_evidence(case, paths, frames, parser, seen)
    classification = classify(rows, cfg, case.get('context', {}), cov['passed'], tech['passed'], static)
    after = {p: sha(p) for p in seen}
    unchanged = seen == after
    if not unchanged:
        classification['status'] = 'UNRESOLVED_SOURCE_CHANGED_DURING_ANALYSIS'
    for label in ('strict', 'default', 'loose'):
        for row in rows:
            row[label + '_criteria'] = candidate_flags(row, cfg[label])
            row[label + '_candidate'] = all(row[label + '_criteria'].values())
    report = {'schema': 'li296_geometric_screen/v1', 'recorded_at': datetime.now(timezone.utc).isoformat(),
              'case': case['name'], 'retrospective': True, 'case_manifest': case,
              'thresholds': cfg, 'thresholds_sha256': sha(thresholds),
              'composition': {'formula': 'C1H194O100S1', 'water_count': 96, 'solute': 'HCHO + SO3^2-', 'atom_order_verified': True},
              'coverage': cov, 'technical_evidence': tech, 'geometry_screen': classification,
              'source_before_after_sha256': [{'path': str(p.relative_to(REPO)) if p.is_relative_to(REPO) else str(p),
                                             'before': seen[p], 'after': after[p]} for p in seen],
              'all_sources_unchanged': unchanged,
              'geometry_summary': {key: {'first': rows[0][key], 'last': rows[-1][key], 'mean': sum(r[key] for r in rows) / len(rows),
                                        'min': min(r[key] for r in rows), 'max': max(r[key] for r in rows)}
                                   for key in ('SC_A', 'CO_A', 'carbon_pyramidal_height_A', 'carbon_angle_min_deg', 'carbon_angle_max_deg')},
              'structural_verdict': None, 'electronic_bond_verdict': None, 'FES_verdict': None,
              'vision_model_check': 'PENDING', 'human_review': 'PENDING',
              'limits': ['Geometric candidate is not electronic bond proof or an approved product.',
                         '10 fs saved frames cannot rule out sub-frame excursions; persistence is observed geometric span, not kinetic lifetime.',
                         'No signed-GDS gate here: interface promotion/localization requires separate registered solvation evidence.',
                         'Initial/preparation history is preserved; conditional retention is not equilibrium spontaneous formation.',
                         'No product/formation conclusion is released before separate visual screening and explicit human approval.']}
    dump(out / 'GEOMETRIC_SCREEN.json', report)
    dump(out / 'FRAME_GEOMETRY.json', rows)
    fields = ['index', 'step', 'time_fs', 'SC_A', 'CO_A', 'carbon_pyramidal_height_A', 'carbon_angle_min_deg', 'carbon_angle_max_deg', 'connectivity_guard_clear', 'strict_candidate', 'default_candidate', 'loose_candidate']
    with (out / 'FRAME_GEOMETRY.tsv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True, help='Explicit single-case JSON; paths relative to repository')
    parser.add_argument('--out', type=Path, required=True, help='New output directory, separate from raw sources')
    args = parser.parse_args()
    case = json.loads(args.manifest.read_text())
    try:
        report = analyze_case(case, args.out)
    except (ValueError, OSError, KeyError) as exc:
        if args.out.exists() and any(args.out.iterdir()):
            print('Refusing to overwrite existing output; error: ' + str(exc))
            raise SystemExit(2)
        args.out.mkdir(parents=True, exist_ok=True)
        dump(args.out / 'ERROR.json', {'status': 'UNRESOLVED_INPUT_ERROR', 'error': str(exc), 'human_review': 'PENDING', 'structural_verdict': None})
        raise SystemExit(2)
    print(case['name'] + ': ' + report['geometry_screen']['status'])


if __name__ == '__main__':
    main()
