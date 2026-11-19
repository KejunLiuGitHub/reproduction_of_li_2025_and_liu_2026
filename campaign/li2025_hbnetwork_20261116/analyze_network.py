#!/usr/bin/env python3
"""Li296 local directed hydrogen-bond graphs; descriptive, no structural verdict."""
import argparse
from collections import Counter, defaultdict, deque
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
SOL = ROOT / 'campaign/li2025_recovery_20261111/analyze_solvation.py'
SOL_SHA = '2194d0b4c19c43f4cddae901e8eb71a069fb18167ddf6874fb0a2adff6c7d2f9'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
if sha(SOL) != SOL_SHA:
    raise ValueError('Registered solvation source changed; review before use')
spec = importlib.util.spec_from_file_location('li296_solvation', SOL)
sol = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sol)
THRESHOLDS = {'strict': (2.4, 160), 'standard': (2.5, 150), 'loose': (2.6, 140)}
ACCEPTORS = np.r_[sol.WATER_O, [sol.ROLES[k] for k in ('Oc', 'Os1', 'Os2', 'Os3')]]
OWNERS = sol.WATER_H - sol.WATER_H % 3 + 2
WATERS = set((sol.WATER_O + 1).tolist())
ROLES = {k: v + 1 for k, v in sol.ROLES.items() if k.startswith('O')}


def bond_geometry(x, box):
    """All donor-H/acceptor-O pairs, fixed parent ownership, periodic vectors."""
    owner = sol.mic(x[OWNERS] - x[sol.WATER_H], box)
    vectors = sol.mic(x[ACCEPTORS][None, :, :] - x[sol.WATER_H][:, None, :], box)
    distance = np.linalg.norm(vectors, axis=2)
    lengths = np.linalg.norm(owner, axis=1)
    denom = distance * lengths[:, None]
    cosine = np.divide(np.sum(vectors * owner[:, None, :], axis=2), denom,
                       out=np.ones_like(distance), where=denom > 0)
    angles = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
    valid = (ACCEPTORS[None, :] != OWNERS[:, None]) & (denom > 0)
    return distance, angles, valid, lengths


def edges_from_geometry(geometry, threshold):
    distance, angles, valid, _ = geometry
    cut, angle = threshold
    ii, jj = np.where(valid & (distance <= cut) & (angles >= angle))
    return [(int(OWNERS[i] + 1), int(sol.WATER_H[i] + 1), int(ACCEPTORS[j] + 1),
             float(distance[i, j]), float(angles[i, j])) for i, j in zip(ii, jj)]


def shortest_path(adjacency, sources, targets):
    queue = deque((x, 0) for x in sorted(sources))
    seen = set(sources)
    while queue:
        node, depth = queue.popleft()
        if node in targets:
            return depth
        for nxt in adjacency.get(node, ()):
            if nxt not in seen:
                seen.add(nxt)
                queue.append((nxt, depth + 1))
    return None


def local_graph(edges):
    donors = {k: sorted({e[0] for e in edges if e[2] == atom}) for k, atom in ROLES.items()}
    first = set().union(*(set(v) for v in donors.values()))
    water_edges = [e for e in edges if e[2] in WATERS]
    adjacency = defaultdict(set)
    for donor, _, acceptor, *_ in water_edges:
        adjacency[donor].add(acceptor)
        adjacency[acceptor].add(donor)
    second = set().union(*(adjacency[w] for w in first)) - first if first else set()
    nodes = first | second
    local = [e for e in edges if e[0] in nodes and (e[2] in nodes or e[2] not in WATERS)]
    directed, undirected = defaultdict(set), defaultdict(set)
    for d, _, a, *_ in local:
        if a in WATERS:
            directed[d].add(a)
            undirected[d].add(a)
            undirected[a].add(d)
    outgoing, incoming = Counter(e[0] for e in edges), Counter(e[2] for e in edges)
    families = {str(w): f'D{outgoing[w]}A{incoming[w]}' for w in sorted(nodes)}
    oc, os = set(donors['Oc']), set().union(*(set(donors[k]) for k in ('Os1', 'Os2', 'Os3')))
    path = shortest_path(undirected, oc, os)
    counts = {k: sum(e[2] == atom for e in local) for k, atom in ROLES.items()}
    return {'donor_waters_1based': donors, 'first_shell': sorted(first), 'second_shell': sorted(second),
            'edges_1based_D_H_A': [list(e[:3]) for e in local],
            'water_degree_families_global': families, 'family_histogram': dict(sorted(Counter(families.values()).items())),
            'solute_HB_counts': counts, 'local_water_edges': sum(e[2] in WATERS for e in local),
            'min_water_path_undirected': path,
            'min_water_path_Oc_donor_waters_to_Os_donor_waters': shortest_path(directed, oc, os),
            'min_water_path_Os_donor_waters_to_Oc_donor_waters': shortest_path(directed, os, oc),
            'class': 'path_le2' if path is not None and path <= 2 else 'path_gt2_or_absent',
            'edge_geometry': {'H_acceptor_A': stats([e[3] for e in local]),
                              'O_H_acceptor_deg': stats([e[4] for e in local])}}


def stats(values):
    if not values:
        return None
    return {'n': len(values), 'mean': float(np.mean(values)), 'min': float(min(values)),
            'median': float(np.median(values)), 'max': float(max(values))}


def episodes(indices, steps, gap=0):
    """Endpoint-span episodes; gap counts allowed absent sampled frames."""
    groups = []
    for i in indices:
        if not groups or i - groups[-1][-1] > gap + 1:
            groups.append([i])
        else:
            groups[-1].append(i)
    return [{'first_step': steps[g[0]], 'last_step': steps[g[-1]],
             'present_frames': len(g), 'absent_frames_within_span': g[-1] - g[0] + 1 - len(g),
             'endpoint_span_fs': steps[g[-1]] - steps[g[0]],
             'left_censored': g[0] == 0, 'right_censored': g[-1] == len(steps) - 1} for g in groups]


def lifetime_description(graphs, steps):
    occurrences = defaultdict(list)
    for i, graph in enumerate(graphs):
        for edge in graph['edges_1based_D_H_A']:
            occurrences[tuple(edge)].append(i)
    table = []
    for edge, indices in sorted(occurrences.items()):
        table.append({'edge_1based_D_H_A': edge, 'present_frames': len(indices),
                      'continuous': episodes(indices, steps), 'intermittent_gap_1_frame': episodes(indices, steps, 1)})
    return {'edge_count': len(table), 'edges': table,
            'note': 'Local membership and geometric HB presence jointly define observation; not physical rates or reaction lifetimes. 10 fs resolution; gap=1 allows one missing sampled frame; endpoint censoring retained.'}


def summarize(rows, criterion, clean=False):
    chosen = [r for r in rows if not clean or not r['parent_OH_warning_H_1based']]
    graphs = [r['networks'][criterion] for r in chosen]
    return {'frames': len(chosen), 'classes': dict(Counter(g['class'] for g in graphs)),
            'mean_local_water_edges': float(np.mean([g['local_water_edges'] for g in graphs])) if graphs else None,
            'mean_first_shell_waters': float(np.mean([len(g['first_shell']) for g in graphs])) if graphs else None,
            'mean_second_shell_waters': float(np.mean([len(g['second_shell']) for g in graphs])) if graphs else None,
            'HB_counts': {k: stats([g['solute_HB_counts'][k] for g in graphs]) for k in ROLES}}


def analyze(xyz, reference, label, job_id, history, biased):
    ref = json.loads(Path(reference).read_text())
    if sha(xyz) != ref['sha256']:
        raise ValueError('XYZ SHA does not match registered all-frame solvation source')
    old = ref['rows']
    rows, steps = [], []
    box = np.array(ref['box'])
    for index, (step, x) in enumerate(sol.frames(xyz)):
        if index >= len(old) or step != old[index]['step']:
            raise ValueError('Registered row/XYZ step mismatch')
        geometry = bond_geometry(x, box)
        networks = {}
        for name, threshold in THRESHOLDS.items():
            edges = edges_from_geometry(geometry, threshold)
            graph = local_graph(edges)
            suffix = '' if name == 'standard' else '_' + name
            for role, atom in ROLES.items():
                if graph['solute_HB_counts'][role] != old[index]['hb_' + role + suffix]:
                    raise ValueError(f'Solvation HB mismatch: {step}/{name}/{role}')
                if name == 'standard' and sorted(e[1] for e in edges if e[2] == atom) != old[index]['donor_H_1based'][role]:
                    raise ValueError('Registered donor identity mismatch')
            networks[name] = graph
        warnings = (sol.WATER_H[geometry[3] > 1.3] + 1).tolist()
        if len(warnings) != old[index]['water_parent_assignment_warning']:
            raise ValueError('Registered OH warning mismatch')
        rows.append({'step': step, 'registered_geometry': old[index], 'parent_OH_warning_H_1based': warnings,
                     'networks': networks})
        steps.append(step)
    if not rows or len(rows) != len(old) or any(b - a != 10 for a, b in zip(steps, steps[1:])):
        raise ValueError('Incomplete reference coverage or nonuniform 10 fs samples')
    starts = [i for i in range(1, len(rows)) if old[i]['SC'] < 2 <= old[i-1]['SC']]
    digest_sources = [Path(xyz), Path(reference), SOL, Path(__file__), Path(__file__).with_name('NETWORK_DESIGN.md')]
    return {'schema': 'li296_local_hbond_network/v1', 'label': label, 'source_job': job_id,
            'history': history, 'biased': biased, 'source_xyz': str(Path(xyz).resolve()),
            'source_sha256': {str(p.resolve()): sha(p) for p in digest_sources},
            'coverage': {'frames': len(rows), 'first_step': steps[0], 'last_step': steps[-1], 'sample_interval_fs': 10,
                         'registered_direct_HB_and_donor_ids_match_all_frames': True},
            'criteria': THRESHOLDS, 'short_SC_access_onset_steps': [steps[i] for i in starts],
            'segment_start_already_short_left_censored': old[0]['SC'] < 2,
            'warning_steps': [r['step'] for r in rows if r['parent_OH_warning_H_1based']],
            'summary': {name: {'all_frames': summarize(rows, name), 'nonwarning_frames': summarize(rows, name, True),
                              'short_SC_access_onsets': summarize([rows[i] for i in starts], name)} for name in THRESHOLDS},
            'lifetimes': {name: lifetime_description([r['networks'][name] for r in rows], steps) for name in THRESHOLDS},
            'rows': rows, 'structural_verdict': None, 'FES_verdict': None, 'human_review': 'PENDING',
            'interpretation': 'Descriptive local HB geometry/topology. Correlated frames; inherited history retained; no electronic bond, proton transfer, spontaneity, equilibrium, causal or rate verdict.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xyz', type=Path, required=True)
    p.add_argument('--solvation', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--label', required=True)
    p.add_argument('--job-id', required=True)
    p.add_argument('--history', required=True)
    p.add_argument('--biased', choices=['true', 'false'], required=True)
    a = p.parse_args()
    result = analyze(a.xyz, a.solvation, a.label, a.job_id, a.history, a.biased == 'true')
    result['replay_argv'] = [sys.executable, '-B', str(Path(__file__).resolve()), *sys.argv[1:]]
    a.outdir.mkdir(parents=True, exist_ok=True)
    full = a.outdir / 'NETWORK.json.gz'
    with full.open('wb') as handle:
        with gzip.GzipFile(fileobj=handle, mode='wb', mtime=0, filename='') as compressed:
            compressed.write(json.dumps(result, separators=(',', ':'), allow_nan=False).encode())
    small = {k: v for k, v in result.items() if k not in ('rows', 'lifetimes')}
    small.update(full_result=str(full.resolve()), full_result_sha256=sha(full))
    (a.outdir / 'SUMMARY.json').write_text(json.dumps(small, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'label': a.label, 'coverage': result['coverage'], 'summary': result['summary']['standard']}))


if __name__ == '__main__':
    main()
