#!/usr/bin/env python3
"""Add registered rigid-transform provenance to the unchanged Li technical gate."""
import json,math
from pathlib import Path
import frozen_validate_stage as frozen
from frozen_validate_stage import demand,sha,dump,read_xyz,make_input,REFERENCE_SHA,RUNNER_SHA,ANALYZER_SHA
FROZEN_VALIDATOR_SHA='3df8173e5b37eb2d1aba2215026710ae1e868d60de3d29090fed8f00bfcc282a'
BOX=(15.,15.,30.)
def mic(v):return [x-L*round(x/L) for x,L in zip(v,BOX)]
def norm(v):return math.sqrt(sum(x*x for x in v))
def relative(a,b):return mic([x-y for x,y in zip(a,b)])
def angle(a,b):return math.degrees(math.acos(max(-1.,min(1.,sum(x*y for x,y in zip(a,b))/(norm(a)*norm(b))))))
def safe(base,name):
 demand(Path(name).name==name,'unsafe provenance filename');return Path(base)/'inputs'/name

def verify_transformation(record,base):
 base=Path(base);t=record['transformation'];source=record['source']
 demand(record['controlled_initial_conditions'] is True and record['prepared_under_bias'] is True,'controlled/history flags required')
 demand(source['source_segment_biased'] is False and source['public_snapshot_origin']=='Figure6b metaD R_inter snapshot','source history mismatch')
 demand(sha(base/'frozen_validate_stage.py')==FROZEN_VALIDATOR_SHA,'frozen validator changed')
 auditpath=base/'SOURCE_AUDIT.json';demand(sha(auditpath)==t['source_audit_sha256'],'source audit SHA mismatch')
 audit=json.loads(auditpath.read_text());demand(audit['status']=='PASS' and str(audit['source_job'])=='8568513' and all(v['pass'] is True for v in audit['checks'].values()),'source audit not all PASS')
 original=safe(base,t['source_coordinate_file']);candidate=safe(base,record['coordinate_file'])
 demand(sha(original)==t['source_frame_sha256'],'source frame SHA mismatch')
 demand(sha(candidate)==t['transformed_xyz_sha256'],'transformed XYZ SHA mismatch')
 demand(source['trajectory_sha256']==audit['checks']['complete_XYZ']['details']['trajectory']['sha256'],'source trajectory/audit identity mismatch')
 old=read_xyz(original);new=read_xyz(candidate);demand(len(old)==len(new)==1 and old[0]['step']==source['frame_step'],'source frame identity mismatch')
 x,y=old[0]['coords'],new[0]['coords'];o=t['water_O_1based']-1;hs=[h-1 for h in t['rotated_H_1based']]
 demand(o in range(2,288,3) and hs==[o-2,o-1],'only a fixed-parent water H pair may rotate')
 R=t['rotation_matrix'];demand(len(R)==3 and all(len(row)==3 for row in R) and all(math.isfinite(v) for row in R for v in row),'invalid rotation matrix')
 orth=max(abs(sum(R[k][i]*R[k][j] for k in range(3))-(i==j)) for i in range(3) for j in range(3))
 determinant=R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1])-R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0])+R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0])
 demand(orth<1e-12 and abs(determinant-1)<1e-12,'matrix is not a proper orthogonal rotation')
 axis=t['axis_unit'];theta=math.radians(t['signed_angle_deg']);demand(len(axis)==3 and abs(norm(axis)-1)<1e-12 and 10<=abs(t['signed_angle_deg'])<=60,'registered axis/angle invalid')
 a,b,c=axis;K=[[0,-c,b],[c,0,-a],[-b,a,0]]
 expected=[[math.cos(theta)*(i==j)+(1-math.cos(theta))*axis[i]*axis[j]+math.sin(theta)*K[i][j] for j in range(3)] for i in range(3)]
 demand(max(abs(R[i][j]-expected[i][j]) for i in range(3) for j in range(3))<1e-12,'matrix/axis/angle mismatch')
 maximum=0.
 for i in range(296):
  if i not in hs:demand(x[i]==y[i],'unrotated atom changed')
  else:
   v=relative(x[i],x[o]);reconstructed=[x[o][j]+sum(R[j][k]*v[k] for k in range(3)) for j in range(3)]
   maximum=max(maximum,max(abs(a-b) for a,b in zip(reconstructed,y[i])))
 demand(maximum<1e-9,'submitted coordinates do not reconstruct from source rotation')
 before=[relative(x[h],x[o]) for h in hs];after=[relative(y[h],y[o]) for h in hs]
 demand(max(abs(norm(a)-norm(b)) for a,b in zip(before,after))<1e-9 and abs(angle(*before)-angle(*after))<1e-7,'water internal geometry changed')
 bonds={tuple(sorted((3*i+h,3*i+2))) for i in range(96) for h in [0,1]}|{(288,j) for j in [289,290,291]}|{(292,j) for j in [293,294,295]}
 minimum,pair=min((norm(relative(y[i],y[j])),[i+1,j+1]) for i in range(296) for j in range(i+1,296) if (i,j) not in bonds)
 demand(minimum>=1.4,'nonbonded contact below 1.4 A')
 return {'source_audit_sha256':sha(auditpath),'source_frame_sha256':sha(original),'transformed_xyz_sha256':sha(candidate),'source_frame_step':source['frame_step'],'rotated_H_1based':t['rotated_H_1based'],'unrotated_atoms_exactly_unchanged':294,'max_reconstruction_delta_A':maximum,'minimum_nonbonded_contact_A':minimum,'minimum_contact_pair':pair,'rotation_matrix':R,'signed_angle_deg':t['signed_angle_deg'],'controlled_initial_conditions':True,'source_segment_biased':False,'prepared_under_bias':True,'source_public_snapshot_origin':source['public_snapshot_origin']}

def verify_input(work,record,base):
 result=frozen.verify_input(work,record,base);result['controlled_transformation']=verify_transformation(record,base);return result

def validate(work,record,base,expected_job=None):
 # Run this before the frozen validator can emit PASS. The original first-MD
 # frame equality check remains and uses the transformed frozen input XYZ.
 try:transform=verify_transformation(record,base)
 except (ValueError,OSError,KeyError,TypeError) as error:
  dump(Path(work)/'CHECKS.json',{'status':'FAIL','case':record['name'],'error':str(error),'human_review':'PENDING','structural_verdict':None});raise
 audit=frozen.validate(work,record,base,expected_job)
 target=Path(work)/'DISTANCE_ANALYSIS.json';data=json.loads(target.read_text())
 data['interpretation'].update(controlled_initial_conditions=True,source_segment_biased=False,public_snapshot_origin=record['source']['public_snapshot_origin'],description='Distance diagnostics from long-SC controlled initial water-H orientations; not retention of a preselected short-SC state, chemical formation, equilibrium spontaneity, or a pure Oc-HB causal intervention.')
 sc=data['summary']['sc'];sc['first_sustained_SC_gt_2p5_A_10_frames_descriptive_only']=sc.pop('first_sustained_dissociation',None)
 dump(target,data);audit['distance_analysis']['sha256']=sha(target);audit['files_sha256'][target.name]=sha(target)
 audit['checks']['controlled_source_rotation']={'pass':True,'details':transform}
 audit['controlled_initial_conditions']=True;audit['interpretation']='Fully unbiased dynamics from controlled water-H orientations; inherited public metaD snapshot history retained. Technical checks do not imply chemical formation, equilibrium spontaneity, or a pure HB causal effect.'
 dump(Path(work)/'CHECKS.json',audit)
 return audit
