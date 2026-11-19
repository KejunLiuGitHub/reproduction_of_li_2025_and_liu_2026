#!/usr/bin/env python3
"""Build exactly two original-coordinate third arms; no HPC or submission."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys
import validate_stage as gate

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p/'AGENTS.md').exists())
OLD = ROOT/'campaign/li2025_hbintervention_20261116/production_v1'
REMOTE = '/public/home/temp/projects/sulfoacetic-vib/saa-vib-v1/li2025_hbcontrol_20261118/production_v1'


def verify_local_sources():
    audit_path = OLD/'SOURCE_AUDIT.json'
    gate.demand(gate.sha(audit_path) == gate.SOURCE_AUDIT_SHA, 'pinned source audit changed')
    audit = json.loads(audit_path.read_text())
    gate.demand(audit['status'] == 'PASS' and all(v['pass'] is True for v in audit['checks'].values()), 'source audit failed')
    verified = {}
    for path, digest in audit['sources_sha256'].items():
        gate.demand('methanol-sulfonization' not in Path(path).parts, 'forbidden source')
        gate.demand(gate.sha(path) == digest, 'audited local source changed: '+path)
        verified[path] = digest
    source = Path(audit['checks']['complete_XYZ']['details']['trajectory']['path'])
    gate.demand(gate.sha(source) == gate.SOURCE_TRAJECTORY_SHA, 'full source trajectory changed')
    frames = gate.read_xyz(source)
    gate.demand(len(frames) == 501 and [f['step'] for f in frames] == list(range(0,5001,10)), 'source coverage changed')
    selected = {}
    old_manifest = json.loads((OLD/'manifest.json').read_text())
    for branch, (step, seed, digest) in gate.APPROVED.items():
        original = OLD/'inputs'/('source_frame_%d.xyz' % step)
        matches = [f for f in frames if f['step'] == step]
        gate.demand(len(matches) == 1 and original.read_text() == matches[0]['raw']
                    and gate.sha(original) == digest, 'source frame is not an exact unique block')
        pair = [r for r in old_manifest['cases'] if r['name'] in ('r%d_oc_donor'%step, 'r%d_reverse'%step)]
        gate.demand(len(pair) == 2 and all(r['seed'] == seed for r in pair), 'seed differs from original pair')
        selected[branch] = {'step':step, 'seed':seed, 'frame_sha256':digest,
                            'pair_cases':[r['name'] for r in pair]}
    return audit, verified, selected


def build(out, evidence):
    out, evidence = Path(out).resolve(), Path(evidence).resolve()
    gate.demand(not out.exists(), 'output exists; refuse overwrite')
    gate.demand(not evidence.exists(), 'evidence directory exists; refuse overwrite')
    audit, source_hashes, selected = verify_local_sources()
    evidence.mkdir(parents=True)
    inputs = out/'inputs'
    inputs.mkdir(parents=True)
    for name in ('run.py','validate_stage.py','frozen_validate_stage.py','build_campaign.py'):
        shutil.copy2(HERE/name, out/name)
    for name, digest in [('frozen_run_cp2k.py',gate.RUNNER_SHA), ('frozen_analyze_event_unbiased.py',gate.ANALYZER_SHA)]:
        gate.demand(gate.sha(OLD/name) == digest, 'frozen dependency changed')
        shutil.copy2(OLD/name, out/name)
    ref = OLD/'inputs/reference_unbiased.inp'
    gate.demand(gate.sha(ref) == gate.REFERENCE_SHA, 'science reference changed')
    shutil.copy2(ref, inputs/'reference_unbiased.inp')
    shutil.copy2(OLD/'SOURCE_AUDIT.json', out/'SOURCE_AUDIT.json')
    old_manifest = json.loads((OLD/'manifest.json').read_text())
    records = []
    geometry = {}
    for branch, choice in selected.items():
        step, seed = choice['step'], choice['seed']
        original = OLD/'inputs'/('source_frame_%d.xyz'%step)
        shutil.copy2(original, inputs/original.name)
        xyz = inputs/(branch+'.xyz')
        shutil.copy2(original, xyz)
        template = next(r for r in old_manifest['cases'] if r['name'] == 'r%d_oc_donor'%step)
        for suffix, steps in [('_smoke',20),('',5000)]:
            name = branch+suffix
            inp = inputs/(name+'.inp')
            inp.write_text(gate.make_input(ref.read_text(), name, xyz.name, seed, steps))
            record = copy.deepcopy(template)
            record.pop('transformation')
            record.update(name=name, branch=branch, condition='original_unrotated', seed=seed, steps=steps,
                          expected_frames=steps//10+1, cp2k_project=name, coordinate_file=xyz.name,
                          source_segment_biased=False, cp2k={'steps_initial_segment':steps,'trajectory_every_fs':10})
            record['source']['source_audit_sha256'] = gate.SOURCE_AUDIT_SHA
            record['source']['coordinate_change'] = 'none; all 296 coordinates byte-identical to audited frame'
            record['stage'] = [{'role':'new','source':p.name,'target':p.name,'sha256':gate.sha(p)} for p in (xyz,inp)]
            record['outputs'] = {'output':name+'.out','trajectory':name+'-pos-1.xyz',
                                 'energy':name+'-1.ener','restart':name+'-1.restart'}
            record['fes'] = {'status':'not_applicable','note':'Original-coordinate long-SC control; fully unbiased dynamics with preparation history retained.'}
            geometry[name] = gate.verify_input(inputs, record, out)
            records.append(record)
        script = ('#!/bin/bash\n#SBATCH --job-name=li25-'+branch+'\n#SBATCH --partition=compute\n#SBATCH --nodes=1\n'
                  '#SBATCH --ntasks=64\n#SBATCH --mem=32G\n#SBATCH --time=48:00:00\nset -euo pipefail\n'
                  'BASE="${LI_WORKFLOW_BASE:?Use the frozen launch wrapper with an absolute package directory}"\n'
                  '/public/home/temp/miniconda3/bin/python3 "$BASE/run.py" --base "$BASE" --case '+branch+' --ntasks "$SLURM_NTASKS"\n')
        (out/('submit_'+branch+'.slurm')).write_text(script)
        wrapper = ('#!/bin/bash\n#SBATCH --job-name=li25-'+branch+'\n#SBATCH --partition=compute\n#SBATCH --nodes=1\n'
                   '#SBATCH --ntasks=64\n#SBATCH --mem=32G\n#SBATCH --time=48:00:00\nset -euo pipefail\n'
                   'export LI_WORKFLOW_BASE="'+REMOTE+'"\nsource "$LI_WORKFLOW_BASE/submit_'+branch+'.slurm"\n')
        (out/('launch_'+branch+'.slurm')).write_text(wrapper)
    selection = {'schema':'li296_original_third_arm_selection/v1','source_job':'8568513', 'selected':selected,
                 'coordinate_policy':'Exact original R1150/2200 frames; no rotation/translation/optimization.',
                 'selection_reason':'Fill the absent original-coordinate arm at each already-tested source frame and same seed; no search for favorable outcomes.',
                 'prepared_under_bias':True,'controlled_initial_conditions':True,'source_segment_biased':False,
                 'n_per_condition':1,'independent_equilibrium_replicates':False}
    gate.dump(out/'SELECTION.json',selection)
    deps = {str(p.relative_to(out)):gate.sha(p) for p in out.rglob('*') if p.is_file()}
    manifest = {'schema':'li2025_hbnetwork_workflow/v1','campaign_type':'li296_original_coordinate_third_arm',
                'candidate_names':list(selected),'cases':records,'dependencies':deps,
                'prepared_under_bias':True,'controlled_initial_conditions':True,'source_segment_biased':False,
                'public_snapshot_origin':'Figure6b metaD R_inter snapshot',
                'velocity_policy':'Fresh Maxwell-Boltzmann generation; same seed/protocol as each prior donor/reverse pair. Independent smoke and fresh production, no velocity/restart inheritance.',
                'comparison':'Within each source frame compare original/donor/reverse over the same full5ps and fixed0–0.1,(0.1,1],(1,5]ps windows; all outcomes retained.',
                'scope':'Two original-coordinate controls, each n=1, correlated source history; no equilibrium spontaneity, independent replication or isolated HB causality.',
                'resource_per_case':{'partition':'compute','ntasks':64,'mem':'32G','walltime':'48:00:00'},
                'stop':'Each20-step smoke must pass before fresh5000-step production; stop at5ps without extension, extra seeds, restraints or new candidates.',
                'submission':'Not submitted by builder. Requires frozen commit, current queue/ledger/resource audit and one receipt per case; skip any recorded/live/completed duplicate.',
                'human_review':'PENDING','structural_verdict':None}
    gate.dump(out/'manifest.json',manifest)
    dryruns = []
    for branch in selected:
        command = [sys.executable,'-B',str(out/'run.py'),'--base',str(out),'--case',branch,
                   '--workdir-root',str(evidence/('dryrun_'+branch)),'--dry-run']
        result = subprocess.run(command,capture_output=True,text=True)
        dryruns.append({'case':branch,'argv':command,'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr})
        gate.demand(result.returncode == 0,'isolated dry-run failed: '+result.stderr)
    for path,digest in source_hashes.items():
        gate.demand(gate.sha(path)==digest,'source changed during build: '+path)
    gate.dump(evidence/'BUILD_VALIDATION.json',{'status':'PASS','physical_simulation':False,'submitted':False,
        'local_source_files_verified_before_after':source_hashes,'source_audit_sha256':gate.SOURCE_AUDIT_SHA,
        'source_schema_note':'Original R source has PILOT_ANALYSIS, not production CHECKS. The previously pinned full technical SOURCE_AUDIT is retained unchanged; no replacement certificate fabricated.',
        'input_and_original_geometry':geometry,'dryruns':dryruns,'manifest_sha256':gate.sha(out/'manifest.json')})
    return manifest


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--evidence',type=Path,required=True)
    args=p.parse_args()
    result=build(args.out,args.evidence)
    print(json.dumps({'status':'LOCAL_READY_NOT_SUBMITTED','cases':result['candidate_names'],'job_ids':[]}))
