#!/usr/bin/env python3
"""Freeze the approved two-frame, +/- water-H intervention into four workflows."""
import argparse,copy,hashlib,io,json,shutil,subprocess,sys
from pathlib import Path
from contextlib import redirect_stdout
import validate_stage as gate
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
REF=ROOT/'campaign/li2025_branches_20261112/inputs/p_iface_v2.inp'
LEGACY=ROOT/'campaign/li2025_hbnetwork_20261116/production_v1'

def build(evidence,out,dry_run=False):
 evidence,out=Path(evidence).resolve(),Path(out).resolve();gate.demand(not out.exists(),'output exists; refuse overwrite')
 feasibility=evidence/'FEASIBILITY.json';audit=evidence/'SOURCE_AUDIT.json';f=json.loads(feasibility.read_text());a=json.loads(audit.read_text())
 gate.demand(f['status']=='FEASIBLE_FOR_DESIGN_REVIEW' and [p['source_step'] for p in f['selected_pairs']]==[1150,2200],'only the approved frames may be built')
 gate.demand(a['status']=='PASS' and all(v['pass'] for v in a['checks'].values()),'source audit failed')
 source=Path(a['checks']['complete_XYZ']['details']['trajectory']['path']);gate.demand(gate.sha(source)==a['checks']['complete_XYZ']['details']['trajectory']['sha256'],'full source trajectory changed')
 blocks=source.read_text().splitlines(keepends=True);parsed=gate.read_xyz(source);gate.demand(len(parsed)==501 and [q['step'] for q in parsed]==list(range(0,5001,10)), 'source frame coverage changed');source_frames={int(parsed[i]['step']):''.join(blocks[i*298:(i+1)*298]) for i in range(501)}
 for n,h in f['source_sha256'].items():gate.demand(gate.sha(n)==h,'feasibility source/code changed: '+n)
 gate.demand(gate.sha(REF)==gate.REFERENCE_SHA,'reference science changed')
 inputs=out/'inputs';inputs.mkdir(parents=True)
 for name in ('run.py','validate_stage.py','frozen_validate_stage.py','build_campaign.py','prepare_intervention.py'):
  shutil.copy2(HERE/name,out/name)
 for name,expected in [('frozen_run_cp2k.py',gate.RUNNER_SHA),('frozen_analyze_event_unbiased.py',gate.ANALYZER_SHA)]:
  gate.demand(gate.sha(LEGACY/name)==expected,'frozen executor/analyzer changed');shutil.copy2(LEGACY/name,out/name)
 shutil.copy2(REF,inputs/'reference_unbiased.inp');shutil.copy2(audit,out/'SOURCE_AUDIT.json')
 selection={'schema':'li296_controlled_orientation_selection/v1','selection_rules':f['selection_rules'],'counts':{k:f[k] for k in ('frames_screened','frames_after_filters','rotation_attempts','feasible_frames')},'full_feasibility_evidence':{'path':str(feasibility),'sha256':gate.sha(feasibility)},'selected_pairs':copy.deepcopy(f['selected_pairs'])}
 for selected in selection['selected_pairs']:
  for threshold,edges in selected['network_edges_and_differences'].items():
   for key in ('source_edges','on_edges','control_edges'):edges.pop(key)
  selected['full_directed_edge_arrays']='See full_feasibility_evidence, selected_pairs matching source_step; complete original/on/reverse edges retained there.'
 gate.dump(out/'SELECTION.json',selection)
 records=[];names=[]
 for p in f['selected_pairs']:
  step=p['source_step'];seed={1150:2026111603,2200:2026111604}[step];source_name='source_frame_%d.xyz'%step
  original=Path(p['xyz']['source']['path']);gate.demand(original.read_text()==source_frames[step] and gate.sha(original)==p['source_frame_sha256'],'selected raw frame differs from audited full source')
  shutil.copy2(original,inputs/source_name)
  for role,coordinate,rotation,theta in [('oc_donor','on',p['on_rotation_matrix'],p['signed_angles_deg'][0]),('reverse','control',p['control_rotation_matrix'],p['signed_angles_deg'][1])]:
   branch='r%d_%s'%(step,role);names.append(branch);xyz=inputs/(branch+'.xyz');shutil.copy2(p['xyz'][coordinate]['path'],xyz)
   gate.demand(gate.sha(xyz)==p['xyz'][coordinate]['sha256'],'selected transformed coordinates changed')
   transform={'source_audit_sha256':gate.sha(audit),'source_coordinate_file':source_name,'source_frame_sha256':p['source_frame_sha256'],'transformed_xyz_sha256':gate.sha(xyz),'water_O_1based':p['water_O_1based'],'rotated_H_1based':p['rotated_H_1based'],'axis_unit':p['axis_unit'],'rotation_matrix':rotation,'signed_angle_deg':theta}
   source_record={'trajectory':str(source),'trajectory_sha256':gate.sha(source),'source_job_id':'8568513','frame_step':step,'source_segment_biased':False,'public_snapshot_origin':'Figure6b metaD R_inter snapshot','prepared_under_bias':True,'controlled_initial_conditions':True,'original_frame_sha256':p['source_frame_sha256']}
   for suffix,steps in [('_smoke',20),('',5000)]:
    name=branch+suffix;inp=inputs/(name+'.inp');inp.write_text(gate.make_input(REF.read_text(),name,xyz.name,seed,steps))
    record={'name':name,'branch':branch,'unbiased':True,'prepared_under_bias':True,'controlled_initial_conditions':True,'condition':role,'pair_source_step':step,'seed':seed,'steps':steps,'expected_frames':steps//10+1,'cp2k_project':name,'coordinate_file':xyz.name,'box_class':'interface','box_angstrom':[15.,15.,30.],'source':source_record,'transformation':transform,'cp2k':{'steps_initial_segment':steps,'trajectory_every_fs':10},'stage':[{'role':'new','source':q.name,'target':q.name,'sha256':gate.sha(q)} for q in [xyz,inp]],'outputs':{'output':name+'.out','trajectory':name+'-pos-1.xyz','energy':name+'-1.ener','restart':name+'-1.restart'},'fes':{'status':'not_applicable','note':'Controlled initial orientation, fully unbiased dynamics; no FES verdict.'},'human_review':'PENDING','structural_verdict':None}
    gate.verify_input(inputs,record,out);records.append(record)
   script='#!/bin/bash\n#SBATCH --job-name=li25-'+branch+'\n#SBATCH --partition=compute\n#SBATCH --nodes=1\n#SBATCH --ntasks=64\n#SBATCH --mem=32G\n#SBATCH --time=48:00:00\nset -euo pipefail\nBASE="${LI_WORKFLOW_BASE:?Set absolute uploaded package directory in registered wrapper}"\n/public/home/temp/miniconda3/bin/python3 "$BASE/run.py" --base "$BASE" --case '+branch+' --ntasks "$SLURM_NTASKS"\n'
   (out/('submit_'+branch+'.slurm')).write_text(script)
 deps={str(p.relative_to(out)):gate.sha(p) for p in out.rglob('*') if p.is_file()}
 manifest={'schema':'li2025_hbnetwork_workflow/v1','campaign_type':'li296_controlled_water_hydrogen_orientation_intervention','candidate_names':names,'cases':records,'dependencies':deps,'selection':'Approved deterministic first feasible frames1150/2200fs, same R history, no future trajectory outcome used.','velocity_policy':'Fresh MB sampling with identical protocol and same seed within each pair; independent smoke and fresh production from frozen transformed XYZ; no restart inheritance. Exact evolving velocity equality is not asserted.','pair_seeds':{'1150':2026111603,'2200':2026111604},'scope':'Two correlated source-frame pairs, n=1 per condition; water-H orientations controlled, all dynamics unbiased. Equal-angle is not equal-energy or an isolated Oc-HB intervention.','source_segment_biased':False,'public_snapshot_origin':'Figure6b metaD R_inter snapshot','prepared_under_bias':True,'controlled_initial_conditions':True,'resource_per_case':{'partition':'compute','ntasks':64,'mem':'32G','walltime':'48:00:00'},'submission':'None by this builder. Root registers frozen wrapper and explicit job receipts. Three64CPU jobs may run and fourth waits under200CPU quota.','human_review':'PENDING','structural_verdict':None}
 gate.dump(out/'manifest.json',manifest)
 if dry_run:
  result=[]
  for name in names:
   command=[sys.executable,'-B',str(out/'run.py'),'--base',str(out),'--case',name,'--workdir-root',str(out/'dryrun'),'--dry-run'];r=subprocess.run(command,text=True,capture_output=True)
   result.append({'case':name,'argv':command,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr});gate.demand(r.returncode==0,'dry-run failed: '+r.stderr)
  gate.dump(out/'DRY_RUN.json',{'status':'PASS','physical_simulation':False,'cases':result})
 return manifest
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--dry-run',action='store_true');a=p.parse_args();m=build(a.evidence,a.out,a.dry_run);print(json.dumps({'status':'PREPARED','cases':m['candidate_names'],'submitted':False}))
