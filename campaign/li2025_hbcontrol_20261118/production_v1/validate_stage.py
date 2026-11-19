#!/usr/bin/env python3
"""Bind unrotated R controls to audited exact frames and the frozen Li gate.

Local source provenance is verified at packaging time; the uploaded package
contains the immutable audit and exact source frames, not fabricated CHECKS.
No chemical, electronic-bond, causality or equilibrium verdict is produced.
"""
import json
import math
from pathlib import Path
import shutil
import frozen_validate_stage as frozen
from frozen_validate_stage import demand, sha, dump, read_xyz, make_input, REFERENCE_SHA, RUNNER_SHA, ANALYZER_SHA

FROZEN_VALIDATOR_SHA = '3df8173e5b37eb2d1aba2215026710ae1e868d60de3d29090fed8f00bfcc282a'
SOURCE_AUDIT_SHA = 'c8078b39f550bc643bec50656615b5e9ffdb88fa9d0985f08cde13b19dc0ca07'
SOURCE_TRAJECTORY_SHA = '80eedcaaa7de00f6f0b04cee30162a6e0ce834d25d01f2c8fea8b507d8219d4d'
APPROVED = {
    'r1150_original': (1150, 2026111603, '5e1b93a126e25373541c1a71f0c11c58987aeb1ef19af65d7c8ddd79c91c7b61'),
    'r2200_original': (2200, 2026111604, 'd21984232c1a00e1baf2087b0cd1f61bfe87370a662e1ec18b1071b06580f26b'),
}
BOX = (15., 15., 30.)


def relative(a, b):
    return [x-y-L*round((x-y)/L) for x, y, L in zip(a, b, BOX)]


def norm(v):
    return math.sqrt(sum(x*x for x in v))


def verify_source(record, base):
    base = Path(base)
    branch = record['branch']
    demand(branch in APPROVED, 'only the two registered original controls are allowed')
    step, seed, digest = APPROVED[branch]
    demand(record['name'] in (branch, branch+'_smoke'), 'unregistered stage name')
    demand(record['steps'] == (20 if record['name'].endswith('_smoke') else 5000), 'registered stage duration changed')
    demand(record['seed'] == seed and record['pair_source_step'] == step, 'paired source/seed changed')
    demand(record['controlled_initial_conditions'] is True and record['prepared_under_bias'] is True,
           'controlled initial conditions and public preparation history must be retained')
    demand(record['source_segment_biased'] is False and record['condition'] == 'original_unrotated', 'original-control flags changed')
    demand('transformation' not in record, 'original control must not carry a rotation transform')
    source = record['source']
    demand(source['source_segment_biased'] is False and source['prepared_under_bias'] is True
           and source['controlled_initial_conditions'] is True, 'source provenance flags changed')
    demand(str(source['source_job_id']) == '8568513' and source['frame_step'] == step,
           'wrong source job/frame')
    demand(source['public_snapshot_origin'] == 'Figure6b metaD R_inter snapshot', 'public preparation origin changed')
    demand(source['trajectory_sha256'] == SOURCE_TRAJECTORY_SHA and source['original_frame_sha256'] == digest,
           'wrong source trajectory/frame fingerprint')
    demand('methanol-sulfonization' not in Path(source['trajectory']).parts, 'forbidden source directory')
    demand(sha(base/'frozen_validate_stage.py') == FROZEN_VALIDATOR_SHA, 'frozen validator changed')
    audit_path = base/'SOURCE_AUDIT.json'
    demand(sha(audit_path) == SOURCE_AUDIT_SHA == source['source_audit_sha256'], 'source audit changed')
    audit = json.loads(audit_path.read_text())
    demand(audit['status'] == 'PASS' and str(audit['source_job']) == '8568513'
           and all(v['pass'] is True for v in audit['checks'].values()), 'source audit not all PASS')
    demand(audit['checks']['complete_XYZ']['details']['trajectory']['sha256'] == SOURCE_TRAJECTORY_SHA,
           'audit does not describe registered source')
    coordinate = record['coordinate_file']
    demand(coordinate == branch+'.xyz', 'unregistered coordinate filename')
    candidate = base/'inputs'/coordinate
    original = base/'inputs'/('source_frame_%d.xyz' % step)
    demand(sha(original) == sha(candidate) == digest and original.read_bytes() == candidate.read_bytes(),
           'original-control coordinates must equal the exact audited frame')
    frames = read_xyz(candidate)
    demand(len(frames) == 1 and frames[0]['step'] == step, 'wrong source frame identity')
    x = frames[0]['coords']
    bonds = {tuple(sorted((3*i+h, 3*i+2))) for i in range(96) for h in (0,1)}
    bonds |= {(288, j) for j in (289,290,291)} | {(292, j) for j in (293,294,295)}
    minimum, pair = min((norm(relative(x[i], x[j])), [i+1, j+1])
                        for i in range(296) for j in range(i+1,296) if (i,j) not in bonds)
    demand(minimum >= 1.4, 'source has a nonbonded contact below 1.4 A')
    return {'source_audit_sha256': SOURCE_AUDIT_SHA, 'source_trajectory_sha256': SOURCE_TRAJECTORY_SHA,
            'source_frame_sha256': digest, 'source_job_id': '8568513', 'source_frame_step': step,
            'unchanged_atoms': 296, 'coordinate_change': 'none; byte-identical original source frame',
            'seed_matches_existing_orientation_pair': seed, 'minimum_nonbonded_contact_A': minimum,
            'minimum_contact_pair': pair, 'SC_A': norm(relative(x[292], x[288])),
            'CO_A': norm(relative(x[289], x[288])), 'source_segment_biased': False,
            'prepared_under_bias': True, 'controlled_initial_conditions': True}


def verify_input(work, record, base):
    result = frozen.verify_input(work, record, base)
    result['original_source'] = verify_source(record, base)
    return result


def validate(work, record, base, expected_job=None):
    try:
        provenance = verify_source(record, base)
        verify_input(work, record, base)
    except (ValueError, OSError, KeyError, TypeError) as error:
        dump(Path(work)/'CHECKS.json', {'status':'FAIL', 'case':record['name'], 'error':str(error),
             'human_review':'PENDING', 'structural_verdict':None})
        raise
    audit = frozen.validate(work, record, base, expected_job)
    target = Path(work)/'DISTANCE_ANALYSIS.json'
    raw = Path(work)/'DISTANCE_ANALYSIS_RAW.json'
    shutil.copyfile(target, raw)
    data = json.loads(target.read_text())
    data['frozen_legacy_output'] = {'path':raw.name, 'sha256':sha(raw),
        'scope':'Unchanged frozen technical-validator distance output; normalized interpretation is below.'}
    data['interpretation'].update(controlled_initial_conditions=True, source_segment_biased=False,
        public_snapshot_origin=record['source']['public_snapshot_origin'],
        description='Long-SC original-coordinate control with fresh velocities matched by seed/protocol to the prior orientation pair. Distances do not imply dissociation, electronic bond formation, equilibrium spontaneity or isolated HB causality.')
    sc = data['summary']['sc']
    sc['first_sustained_SC_gt_2p5_A_10_frames_descriptive_only'] = sc.pop('first_sustained_dissociation', None)
    dump(target, data)
    audit['distance_analysis']['sha256'] = sha(target)
    audit['files_sha256'][target.name] = sha(target)
    audit['files_sha256'][raw.name] = sha(raw)
    audit['checks']['original_source_coordinates_and_pair_seed'] = {'pass':True, 'details':provenance}
    audit['controlled_initial_conditions'] = True
    audit['source_segment_biased'] = False
    audit['interpretation'] = 'Fully unbiased long-SC original-coordinate control after public metaD preparation history; technical completeness only.'
    dump(Path(work)/'CHECKS.json', audit)
    return audit
