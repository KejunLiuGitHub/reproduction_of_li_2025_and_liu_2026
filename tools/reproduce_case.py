#!/usr/bin/env python3
"""Replay the registered Li296 distance, solvation and HB-network analyses."""
import argparse
import gzip
import json
import math
from pathlib import Path
import re
import statistics
import subprocess
import sys
import time

sys.dont_write_bytecode = True
from _common import ROOT, fresh_directory, load_case, load_module, sha256, source_path, write_json

DISTANCE = ROOT / 'campaign/li2025_direct/analyze_event_unbiased.py'
SOLVATION = ROOT / 'campaign/li2025_recovery_20261111/analyze_solvation.py'
NETWORK = ROOT / 'campaign/li2025_hbnetwork_20261116/analyze_network.py'
PARSER = ROOT / 'campaign/li2025_hbnetwork_20261116/production_v1/validate_stage.py'


def complete_prefixes(xyz, energy, out, through_step=None):
    """Separate an unfinished final frame/line; all complete records are retained."""
    parser = load_module(PARSER, 'li296_public_parsers')
    xyz_target, energy_target = out / 'COMPLETE.xyz', out / 'COMPLETE.ener'
    xyz_tail = b''
    with xyz.open('rb') as source, xyz_target.open('wb') as target:
        while first := source.readline():
            block = [first] + [source.readline() for _ in range(297)]
            if not all(line and line.endswith(b'\n') for line in block):
                xyz_tail = b''.join(block) + source.read()
                break
            if first.strip() != b'296':
                raise ValueError('Unexpected atom count in complete XYZ block')
            target.writelines(block)
    energy_tail = b''
    with energy.open('rb') as source, energy_target.open('wb') as target:
        for line in source:
            if not line.endswith(b'\n'):
                energy_tail = line + source.read()
                break
            target.write(line)
    (out / 'XYZ_EXCLUDED_TAIL.bin').write_bytes(xyz_tail)
    (out / 'ENERGY_EXCLUDED_TAIL.bin').write_bytes(energy_tail)
    frames = parser.read_xyz(xyz_target)
    energies = parser.read_energy(energy_target)
    original_coverage = {'XYZ_frames': len(frames), 'XYZ_last_step': frames[-1]['step'],
                         'energy_rows': len(energies), 'energy_last_step': int(energies[-1][0])}
    steps = [frame['step'] for frame in frames]
    if any(step is None for step in steps) or any(b-a != 10 for a, b in zip(steps, steps[1:])):
        raise ValueError('XYZ steps must be consecutive 10 fs samples')
    if any(frame['time'] is None or abs(frame['time']-frame['step']) > .0011 for frame in frames):
        raise ValueError('XYZ time does not match the 1 fs timestep')
    if any(b[0]-a[0] != 1 for a, b in zip(energies, energies[1:])):
        raise ValueError('Energy rows are not consecutive')
    if any(int(row[0]) != row[0] or abs(row[0]-row[1]) > 2e-6 for row in energies):
        raise ValueError('Energy step/time mismatch')
    by_step = {int(row[0]): row for row in energies}
    shared = [frame for frame in frames if frame['step'] in by_step]
    if not shared:
        raise ValueError('XYZ and energy have no shared steps')
    for frame in shared:
        if frame['energy'] is not None and abs(frame['energy']-by_step[frame['step']][4]) > 2e-8:
            raise ValueError('XYZ comment and energy table disagree at step ' + str(frame['step']))
    if through_step is not None:
        if through_step < 0:
            raise ValueError('--through-step must be nonnegative')
        frames = [frame for frame in frames if frame['step'] <= through_step]
        energies = [row for row in energies if row[0] <= through_step]
        shared = [frame for frame in shared if frame['step'] <= through_step]
        if not frames or not energies or not shared:
            raise ValueError('--through-step selects no shared complete observations')
        xyz_target.write_text(''.join(frame['raw'] for frame in frames))
        kept = []
        for line in energy_target.read_bytes().splitlines(keepends=True):
            if not line.strip() or line.lstrip().startswith(b'#') or float(line.split()[0]) <= through_step:
                kept.append(line)
        energy_target.write_bytes(b''.join(kept))
        steps = [frame['step'] for frame in frames]
    coverage = {'XYZ': {'frames': len(frames), 'first_step': steps[0], 'last_step': steps[-1], 'excluded_tail_bytes': len(xyz_tail)},
                'energy': {'rows': len(energies), 'first_step': int(energies[0][0]), 'last_step': int(energies[-1][0]), 'excluded_tail_bytes': len(energy_tail)},
                'shared_XYZ_energy_first_step': shared[0]['step'], 'shared_XYZ_energy_last_step': shared[-1]['step'],
                'log_coverage': 'Not inferred by this descriptive replay; use the release run record for MD/SCF coverage.',
                'selection': {'through_step_requested': through_step, 'original_complete_coverage': original_coverage,
                              'frames_omitted_after_requested_step': original_coverage['XYZ_frames']-len(frames),
                              'selected_XYZ_sha256': sha256(xyz_target), 'selected_energy_sha256': sha256(energy_target)}}
    return xyz_target, energy_target, coverage, energies


def stats(values):
    return {'mean': statistics.fmean(values), 'min': min(values), 'max': max(values), 'first': values[0], 'last': values[-1]}


def reproduce(case_id, out, datasets=None, network=True, through_step=None):
    case, manifest_path = load_case(case_id, datasets)
    xyz, energy = source_path(case['trajectory']), source_path(case['energy'])
    for key, path in [('trajectory', xyz), ('energy', energy)]:
        if key in case.get('sha256', {}) and sha256(path) != case['sha256'][key]:
            raise ValueError('Dataset SHA mismatch for ' + key)
    input_sources = []
    if case.get('input'):
        inp = source_path(case['input'])
        cells = re.findall(r'(?mi)^\s*ABC\s+([0-9.eEdD+\-]+)\s+([0-9.eEdD+\-]+)\s+([0-9.eEdD+\-]+)', inp.read_text())
        if len(cells) != 1 or any(abs(float(v.replace('D','E').replace('d','e'))-length) > 1e-8 for v, length in zip(cells[0], case['box_angstrom'])):
            raise ValueError('Dataset cell differs from the archived CP2K input')
        input_sources.append(inp)
    code = [Path(__file__).resolve(), ROOT / 'tools/_common.py', DISTANCE, SOLVATION, PARSER]
    if network:
        code += [NETWORK, NETWORK.with_name('NETWORK_DESIGN.md')]
    source_hashes = {str(p): sha256(p) for p in [manifest_path, xyz, energy, *input_sources, *code]}
    out = fresh_directory(out)
    write_json(out / 'SOURCES_BEFORE.json', source_hashes)
    complete_xyz, complete_energy, coverage, energies = complete_prefixes(xyz, energy, out, through_step)
    box = ','.join(str(v) for v in case['box_angstrom'])
    commands = [[sys.executable, '-B', str(DISTANCE), '--trajectory', str(complete_xyz), '--energy', str(complete_energy),
                 '--label', str(case['label']), '--box', box, '--out', str(out / 'DISTANCE_ANALYSIS_RAW.json')],
                [sys.executable, '-B', str(SOLVATION), str(complete_xyz), '--out', str(out / 'SOLVATION.json')]]
    if case['environment'] == 'bulk':
        commands[1].append('--bulk')
    if network:
        commands.append([sys.executable, '-B', str(NETWORK), '--xyz', str(complete_xyz), '--solvation', str(out / 'SOLVATION.json'),
                         '--outdir', str(out / 'network'), '--label', str(case['label']), '--job-id', str(case['job_id']) if case['job_id'] is not None else 'archive:' + case_id,
                         '--history', str(case['history']), '--biased', str(case['biased']).lower()])
    execution = {'cwd': str(ROOT), 'commands': commands, 'case_id': case_id, 'results': []}
    write_json(out / 'COMMANDS.json', execution)
    for index, command in enumerate(commands):
        start = time.monotonic()
        run = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
        (out / f'command_{index}.stdout.txt').write_text(run.stdout)
        (out / f'command_{index}.stderr.txt').write_text(run.stderr)
        execution['results'].append({'index': index, 'returncode': run.returncode, 'elapsed_seconds': time.monotonic()-start})
        write_json(out / 'COMMANDS.json', execution)
        if run.returncode:
            raise ValueError(f'Analysis command {index} failed; see its stderr file')
    distance = json.loads((out / 'DISTANCE_ANALYSIS_RAW.json').read_text())
    normalized = dict(distance)
    normalized['summary'] = {**distance['summary'], 'sc': dict(distance['summary']['sc'])}
    normalized['summary']['sc']['first_SC_gt2p5_A_10_sample_span'] = normalized['summary']['sc'].pop('first_sustained_dissociation', None)
    normalized['interpretation'] = {'biased': case['biased'], 'history': case['history'],
                                    'description': 'Periodic distance observations over the supplied complete frames.'}
    normalized['raw_output'] = {'path': 'DISTANCE_ANALYSIS_RAW.json', 'sha256': sha256(out / 'DISTANCE_ANALYSIS_RAW.json')}
    write_json(out / 'DISTANCE_ANALYSIS.json', normalized)
    solvation = json.loads((out / 'SOLVATION.json').read_text())
    rows = solvation['rows']
    if len(rows) != coverage['XYZ']['frames'] or len(rows) != len(distance['per_frame']):
        raise ValueError('Analyzer frame coverage differs')
    for sol, dist in zip(rows, distance['per_frame']):
        if sol['step'] != dist['step'] or abs(sol['SC']-dist['sc']) > 1e-9 or abs(sol['CO']-dist['co']) > 1e-9:
            raise ValueError('Distance and solvation per-frame geometry differ')
    summary = {'schema': 'li296_public_replay/v1', 'case': case_id, 'label': case['label'], 'job_id': str(case['job_id']) if case['job_id'] is not None else None,
               'capture_status': case['status'], 'biased': case['biased'], 'history': case['history'],
               'box_angstrom': case['box_angstrom'], 'environment': case['environment'], 'coverage': coverage,
               'SC_A': stats([r['SC'] for r in rows]), 'CO_A': stats([r['CO'] for r in rows]),
               'temperature_K': stats([r[3] for r in energies if rows[0]['step'] <= r[0] <= rows[-1]['step']]),
               'hydrogen_bonds': {criterion: {role: stats([r['hb_'+role+suffix] for r in rows]) for role in ('Oc', 'Os1', 'Os2', 'Os3')}
                                  for criterion, suffix in [('strict', '_strict'), ('standard', ''), ('loose', '_loose')]},
               'water_O_within_3p5_A': {role: stats([r['water_O_3p5_'+role] for r in rows]) for role in ('Oc', 'Os1', 'Os2', 'Os3')},
               'carbon_pyramidal_height_A': stats([r['carbon_pyramidal_height'] for r in rows]),
               'signed_GDS_A': {role: {**stats([r['gds_'+role] for r in rows]), 'vapour_fraction': sum(r['gds_'+role] > 0 for r in rows)/len(rows)} for role in ('C','S','Oc')} if case['environment'] == 'interface' else None,
               'normal_cosines': {key: stats([r[key] for r in rows]) for key in ('CO_normal_cos','CS_normal_cos')} if case['environment'] == 'interface' else None,
               'network_computed': network, 'description': 'Complete-frame geometrical observations; status and preparation history are retained from the selected dataset.'}
    if network:
        with gzip.open(out / 'network/NETWORK.json.gz', 'rt') as handle:
            net = json.load(handle)
        if [r['registered_geometry'] for r in net['rows']] != rows or net['biased'] != case['biased']:
            raise ValueError('Network geometry or bias metadata differs')
        summary['network'] = net['summary']
        summary['target_water204_to_Oc_counts'] = {threshold: {str(h): sum([204,h,290] in row['networks'][threshold]['edges_1based_D_H_A'] for row in net['rows']) for h in (202,203)} for threshold in net['criteria']}
    # Raw analyzer files remain unchanged. This public record carries the actual
    # case history instead of adopting the generic legacy distance interpretation.
    summary['distance_field_note'] = 'Legacy first_sustained_dissociation records SC>2.5 A for 10 sampled frames; it is a distance-span statistic.'
    write_json(out / 'SUMMARY.json', summary)
    after = {p: sha256(Path(p)) for p in source_hashes}
    if after != source_hashes:
        raise ValueError('A source changed during replay')
    write_json(out / 'SOURCE_SHA256.json', {p: {'before': h, 'after': after[p]} for p, h in source_hashes.items()})
    write_json(out / 'OUTPUT_SHA256.json', {str(p.relative_to(out)): sha256(p) for p in out.rglob('*') if p.is_file() and p.name != 'OUTPUT_SHA256.json'})
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--datasets', type=Path, help='Alternative dataset manifest; default data/DATASETS.json')
    parser.add_argument('--no-network', action='store_true')
    parser.add_argument('--through-step', type=int, help='Use all complete frames/energy rows at or before this step; default all available')
    args = parser.parse_args()
    result = reproduce(args.case, args.out, args.datasets, not args.no_network, args.through_step)
    print(json.dumps({'case': result['case'], 'frames': result['coverage']['XYZ']['frames'], 'SC_A': result['SC_A'], 'output': str(args.out.resolve())}))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(2)
