#!/usr/bin/env python3
"""Reproducible retrospective suite, read-only sources, new evidence directory."""
import argparse
import json
from pathlib import Path
import screen


def relative(p):
    return str(Path(p).resolve().relative_to(screen.REPO))


def cases():
    result = []
    old = screen.REPO / 'results/li2025_autofollow/20261113T213847+0800'
    natural = screen.REPO / 'results/li2025_autofollow/20261116T222349+0800/retrieval/li2025_hbnetwork_20261116/production_v1/runs'
    for name, folder, manifest in [
        ('natural_original40', natural / 'natural_original40', 'campaign/li2025_hbnetwork_20261116/production_v1/manifest.json'),
        ('natural_r19830', natural / 'natural_r19830', 'campaign/li2025_hbnetwork_20261116/production_v1/manifest.json'),
        ('p_iface_v2', old / 'li2025_branches_20261112/p_iface_v2', 'campaign/li2025_branches_20261112/manifest.json'),
        ('relax_v2', old / 'li2025_branches_20261112/relax_v2', 'campaign/li2025_branches_20261112/manifest.json')]:
        source_manifest = screen.REPO / manifest
        registration = next(c for c in json.loads(source_manifest.read_text())['cases'] if c['name'] == name)
        record = json.loads((folder / 'run_record.json').read_text())
        case = {'name': name, 'run_case': name, 'job_id': record['job_id'], 'box_A': [15., 15., 30.], 'start_step': 0,
                'context': {'prepared_under_bias': name != 'p_iface_v2', 'public_P_label': name == 'p_iface_v2',
                            'controlled_initial_conditions': False,
                            'description': 'Public-P-labelled initial configuration, new velocities' if name == 'p_iface_v2' else 'Conditional retention after recorded biased preparation; not equilibrium spontaneous formation'},
                'source_history': registration['source'],
                'files': {'trajectory': relative(folder / (name + '-pos-1.xyz')), 'input': relative(folder / (name + '.inp')),
                          'source_xyz': relative(folder / (name + '.xyz')), 'energy': relative(folder / (name + '-1.ener')),
                          'restart': relative(folder / (name + '-1.restart')), 'run_record': relative(folder / 'run_record.json'),
                          'checks': relative(folder / 'CHECKS.json'), 'source_manifest': manifest},
                'expected_source_sha256': next(s['sha256'] for s in registration['stage'] if s['target'] == name + '.xyz'),
                'expected_file_sha256': {'input': next(s['sha256'] for s in registration['stage'] if s['target'] == name + '.inp')}}
        if (folder / (name + '.out')).exists():
            case['files']['log'] = relative(folder / (name + '.out'))
        src = registration['source']
        if 'trajectory' in src:
            case['files']['parent_trajectory'] = relative(screen.REPO / src['trajectory'])
            case['expected_parent_trajectory_sha256'] = src['trajectory_sha256']
            case['parent_step'] = src['frame_step']
        result.append(case)
    release = old / 'li2025_solvation_20261111/release'
    result.append({'name': 'legacy_release', 'run_case': 'release', 'job_id': '8568949', 'box_A': [15., 15., 30.], 'start_step': 2000,
                   'context': {'prepared_under_bias': True, 'public_P_label': False, 'controlled_initial_conditions': False,
                               'description': '5 ps release following 2 ps constrained preparation; saved XYZ begins at 2010. No synthetic t0 is inserted. event10410.xyz is overridden by EXT_RESTART.'},
                   'files': {k: relative(release / v) for k, v in {'trajectory': 'release-pos-1.xyz', 'input': 'release.inp',
                             'source_xyz': 'event10410.xyz', 'source_restart': 'relax_parent.restart', 'source_provenance': 'parent_provenance.json',
                             'energy': 'release-1.ener', 'restart': 'release-1.restart', 'run_record': 'run_record.json'}.items()},
                   'expected_file_sha256': {'input': '24688c668f967a29a855c184ea1795b48643a316e07e1acb90d36bfa4bdf47c4',
                                            'source_restart': '534b7879ae8f2d8551079e49afb6e22725c8a23d293f4dfac5d871ea7941160a'}})
    rc = screen.REPO / 'results/li2025_autofollow/decision_20261112/li2025_recovery_20261111/r_interface_eq'
    result.append({'name': 'public_RC_5ps', 'run_case': 'r_interface_eq', 'job_id': '8568513', 'box_A': [15., 15., 30.], 'start_step': 0,
                   'context': {'prepared_under_bias': False, 'public_P_label': False, 'controlled_initial_conditions': False,
                               'description': 'Public Fig6b R_interface, unbiased 5 ps; single trajectory, no negative reaction conclusion'},
                   'files': {k: relative(rc / v) for k, v in {'trajectory': 'r_interface_eq-pos-1.xyz', 'input': 'r_interface_eq.inp',
                             'source_xyz': 'r_interface_eq.xyz', 'energy': 'r_interface_eq-1.ener', 'restart': 'r_interface_eq-1.restart',
                             'run_record': 'run_record.json', 'distance_certificate': 'PILOT_ANALYSIS.json', 'log': 'r_interface_eq.out'}.items()},
                   'expected_source_sha256': 'b44766d837cae154f7281df37f17a179d6ec607b47ef0cebae578759f9c580b5',
                   'expected_file_sha256': {'input': '41bbc6f1d569180f1cc841da3685e392fda3bd9b5d5b0e0da985a141f1d713e9'}})
    workbook = screen.REPO / 'results/li2025_zenodo/source-data.xlsx'
    for label in ('P_iface', 'P_bulk', 'RC_iface', 'RC_bulk'):
        p = screen.REPO / ('results/li2025_zenodo/structure_compare/li2025_fig6b_' + label + '.xyz')
        box = [14.2296] * 3 if label.endswith('bulk') else [15., 15., 30.]
        result.append({'name': 'public_snapshot_' + label, 'static': True, 'box_A': box,
                       'context': {'public_P_label': label.startswith('P_'), 'published_label': label,
                                   'description': 'Static Fig6b coordinate label, not a verified TS/product calculation or time trajectory'},
                       'files': {'trajectory': relative(p), 'source_workbook': relative(workbook),
                                 'export_metadata': 'results/li2025_zenodo/structure_compare/LI2025_STRUCTURE_COMPARISON.json'},
                       'expected_file_sha256': {'source_workbook': 'f4e2bedc004f0abaf180f55e95bd2aa9f7ddc2aa86ef09a683fd2d5afe0a9992'}})
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        raise SystemExit('Refusing nonempty suite output')
    args.out.mkdir(parents=True, exist_ok=True)
    selected = cases()
    screen.dump(args.out / 'SUITE_MANIFEST.json', {'retrospective': True, 'threshold_sha256': screen.sha(screen.HERE / 'thresholds.json'), 'cases': selected})
    summaries = []
    for case in selected:
        screen.dump(args.out / (case['name'] + '_manifest.json'), case)
        report = screen.analyze_case(case, args.out / case['name'])
        g = report['geometry_screen']
        brief = {'case': case['name'], 'status': g['status'], 'technical_pass': report['technical_evidence']['passed'],
                 'coverage': report['coverage'], 'context': g['initial_context'], 'sensitivity_status': g['sensitivity_status'],
                 'ambiguous_frames': g['connectivity_ambiguous_frames'], 'sources_unchanged': report['all_sources_unchanged'],
                 'by_threshold': {k: {'frames': v['candidate_frames'], 'sustained_segments': len(v['sustained_segments']),
                                      'longest_span_fs': max([s['span_fs'] for s in v['segments']] or [0]),
                                      'terminal_span_fs': v['terminal_sustained_segment']['span_fs'] if v['terminal_sustained_segment'] else None}
                                  for k, v in g['threshold_sensitivity'].items()},
                 'structural_verdict': None, 'human_review': 'PENDING'}
        summaries.append(brief)
        print(case['name'], brief['status'], brief['by_threshold']['default'], flush=True)
    screen.dump(args.out / 'SUMMARY.json', {'schema': 'li296_geometry_screen_suite/v1', 'retrospective': True, 'cases': summaries,
                                           'human_review': 'PENDING', 'structural_verdict': None})


if __name__ == '__main__':
    main()
