#!/usr/bin/env python3
"""Freeze strictly matched candidates or explicitly unmatched natural retention pilots."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT / 'results/li2025_autofollow/20261115T100111+0800'
CASES = {
    'original_30_40ps': ('8571154', 'li2025_extend_20261113/interface_r1_seg4', 'interface_r1', 2026111601),
    'r_10_20ps': ('8571156', 'li2025_r_extend_20261113/r_eq_metad_seg2', 'r_eq_metad', 2026111602),
}
CALIPERS = {'SC': .05, 'CO': .03, 'CO_normal_cos': .10, 'CS_normal_cos': .10,
            'gds_C': .30, 'gds_Oc': .30, 'gds_S': .50, 'OCS_angle': 5.,
            'carbon_pyramidal_height': .03, 'temperature_K': 10.}
EXACT = ['hb_Oc', 'water_O_3p5_Oc']
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def compare(a, b):
    deltas = {k: abs(a[k] - b[k]) for k in CALIPERS}
    failed = [k for k, diff in deltas.items() if diff > CALIPERS[k]]
    failed += [k for k in EXACT if a[k] != b[k]]
    return {'differences': deltas, 'exact_pairs': {k: [a[k], b[k]] for k in EXACT},
            'failed_fields': failed, 'pass': not failed}


def eligible_onsets(network, temperature):
    starts = set(network['short_SC_access_onset_steps'])
    return [dict(r['registered_geometry'], temperature_K=temperature[r['step']],
                 network_class=r['networks']['standard']['class']) for r in network['rows']
            if r['step'] in starts and not r['parent_OH_warning_H_1based']]


def matches(rows):
    pairs, failures, tested = [], Counter(), 0
    for i, a in enumerate(rows):
        for b in rows[i+1:]:
            if a['network_class'] == b['network_class'] or b['step'] - a['step'] < 1000:
                continue
            tested += 1
            comparison = compare(a,b)
            failures.update(comparison['failed_fields'])
            if comparison['pass']:
                pairs.append({'steps': [a['step'],b['step']], 'classes': [a['network_class'],b['network_class']],
                              'matching': comparison})
    return {'eligible_onsets': len(rows), 'class_counts': dict(Counter(r['network_class'] for r in rows)),
            'tested_opposite_class_pairs_separated_ge1ps': tested,
            'failed_field_counts_nonexclusive': dict(failures), 'strict_matches': pairs,
            'independence': '1 ps separation only rejects adjacent samples; it does not establish independence.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,required=True)
    a=parser.parse_args()
    hashes, loaded, reports, selected = {}, {}, {}, []
    for label,(job,rel,project,seed) in CASES.items():
        full=a.results/label/'NETWORK.json.gz'
        with gzip.open(full,'rt') as f: network=json.load(f)
        work=PARENT/rel;xyz=work/(project+'-pos-1.xyz');checks_path=work/'CHECKS.json'
        checks=json.loads(checks_path.read_text())
        if checks['status']!='PASS' or checks['job_id']!=job or sha(xyz)!=checks['files_sha256'][xyz.name]:
            raise ValueError('Formal source CHECKS/trajectory mismatch')
        if sha(xyz)!=network['source_sha256'][network['source_xyz']]:
            raise ValueError('Graph source does not match formal trajectory')
        energy=work/(project+'-1.ener')
        if sha(energy)!=checks['files_sha256'][energy.name]:
            raise ValueError('Formal energy source mismatch')
        temperature={int(r[0]):float(r[3]) for r in np.loadtxt(energy)}
        onsets=eligible_onsets(network,temperature)
        reports[label]=matches(onsets)
        loaded[label]=(network,temperature,xyz,checks_path)
        for p in [full,xyz,checks_path,energy]: hashes[str(p.resolve())]=sha(p)
    group_available=all(reports[k]['strict_matches'] for k in CASES)
    # A four-case mechanism group additionally needs justified independence.
    # Temporal separation is deliberately not upgraded to that claim here.
    mode='TWO_UNMATCHED_NATURAL_RETENTION_PILOTS'
    reason=('Strict geometric pairs exist but independent network preparations have not been established; use pilots only.'
            if group_available else 'The two distinct histories do not both supply a strictly matched opposite-network pair; no four-case mechanism group.')
    for label,(job,rel,project,seed) in CASES.items():
        network,temperature,xyz,checks_path=loaded[label]
        row=next(r for r in reversed(network['rows']) if r['registered_geometry']['SC']<2 and not r['parent_OH_warning_H_1based'])
        step=row['step'];geometry=dict(row['registered_geometry'],temperature_K=temperature[step])
        index=next(i for i,r in enumerate(network['rows']) if r['step']==step)
        lines=xyz.read_bytes().splitlines(keepends=True);frame=b''.join(lines[index*298:(index+1)*298])
        graph=row['networks']['standard']
        selected.append({'name':'natural_original40' if label.startswith('original') else 'natural_r19830',
            'trajectory':str(xyz.resolve()),'sha256':sha(xyz),'frame_step':step,'frame_index_0based':index,
            'frame_block_sha256':hashlib.sha256(frame).hexdigest(),'source_CHECKS':str(checks_path.resolve()),
            'source_CHECKS_sha256':sha(checks_path),'source_job':job,'seed':seed,'source_prepared_under_bias':True,
            'source_biased':True,'production_biased':False,'production_steps':5000,
            'selection_rule':'Latest frame of this formal complete segment with SC<2 A and no fixed-parent OH warning. No future duration/retention/network favorable outcome is used.',
            'bias_history':network['history'],'purpose':'Independent-history natural-configuration unbiased retention pilot; unmatched, not a network causal group or spontaneous formation test.',
            'descriptor':{k:geometry[k] for k in ['step',*CALIPERS,*EXACT,'hb_Oc_strict','hb_Oc_loose','max_parent_water_OH']},
            'network_descriptor':{k:graph[k] for k in ['class','donor_waters_1based','first_shell','second_shell','local_water_edges','min_water_path_undirected','family_histogram']},
            'network_result':str((a.results/label/'NETWORK.json.gz').resolve()),
            'human_review':'PENDING','structural_verdict':None})
    difference=compare(selected[0]['descriptor'],selected[1]['descriptor'])
    for p in [Path(__file__),Path(__file__).with_name('NETWORK_DESIGN.md'),Path(__file__).with_name('analyze_network.py')]:
        hashes[str(p.resolve())]=sha(p)
    report={'schema':'li296_frozen_natural_candidates/v1','status':'FROZEN_CANDIDATE_SELECTION_NOT_SUBMISSION',
            'mode':mode,'reason':reason,'cases':selected,'strict_group_available_before_independence_check':group_available,
            'matching':reports,'calipers':CALIPERS,'exact_matching_fields':EXACT,'pilot_pair_matching':difference,
            'independence':'Two distinct prepared histories, one configuration per pilot; neither n=2 per condition nor proof of equilibrium independence. Same-coordinate velocity replicas are not new ICs.',
            'limitations':['No network causal comparison from these unmatched pilots.','Ordinary metaD preparation retained; production>=5ps unbiased tests retention only.',
                           'Negative GDS or zero vapour is not interface failure.','Original chain stops at50ps for review; no automatic further extension.'],
            'source_sha256':hashes,'replay_argv':[sys.executable,'-B',str(Path(__file__).resolve()),*sys.argv[1:]],
            'human_review':'PENDING','structural_verdict':None,'FES_verdict':None}
    (a.results/'CANDIDATES.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'mode':mode,'reason':reason,'matching':reports,'cases':[{'name':c['name'],'step':c['frame_step'],'descriptor':c['descriptor']} for c in selected]},indent=2))


if __name__=='__main__':
    main()
