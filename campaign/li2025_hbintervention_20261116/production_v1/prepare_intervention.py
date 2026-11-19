#!/usr/bin/env python3
"""Local deterministic feasibility search: real R frames, matched +/- water rotations."""
import argparse, hashlib, importlib.util, json, math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'results/li2025_autofollow/decision_20261112/li2025_recovery_20261111/r_interface_eq/r_interface_eq-pos-1.xyz'
SOLVATION=ROOT/'results/li2025_autofollow/decision_20261112/R_eq_solvation.json'
NET=ROOT/'campaign/li2025_hbnetwork_20261116/analyze_network.py'
spec=importlib.util.spec_from_file_location('li_network_reference',NET);net=importlib.util.module_from_spec(spec);spec.loader.exec_module(net)
sol=net.sol;BOX=np.array([15.,15.,30.]);ORDER=sol.ORDER
RULES={'source_step_min':1000,'SC_min_A':3.0,'CO_interval_A':[1.1,1.3],'C_Oc_signed_GDS_interval_A':[-3.0,.5], 'S_signed_GDS_interval_A':[-5.0,.5],'parent_OH_warning_max':0,'nonbonded_min_A':1.4,'water_O_to_Oc_interval_A':[2.65,3.3],'rotation_abs_angle_deg_interval':[10.,60.],'pair_frame_separation_min_fs':1000,'frame_order':'ascending step; first feasible, then first feasible at least1000fs later','water_order':'ascending minimal absolute angle, O index, selected H index','contrast':'selected water contributes an Oc donor in strict/standard/loose for +theta and none even loose for -theta; it is not an Oc donor in the original even loose','scope':'Controlled initial-condition pilot; same source history, separated frames are correlated, not independent equilibrium samples.'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def mic(v):return sol.mic(v,BOX)
BONDS=[(3*i+h,3*i+2) for i in range(96) for h in [0,1]]+[(288,j) for j in [289,290,291]]+[(292,j) for j in [293,294,295]]
MASK=np.triu(np.ones((296,296),bool),1)
for i,j in BONDS:MASK[min(i,j),max(i,j)]=False

def contact(x):
 d=np.linalg.norm(mic(x[:,None,:]-x[None,:,:]),axis=2);d[~MASK]=np.inf;i,j=np.unravel_index(np.argmin(d),d.shape)
 return {'minimum_A':float(d[i,j]),'pair_1based':[int(i+1),int(j+1)],'elements':[ORDER[i],ORDER[j]],'exclusions':'Only fixed parent water O-H, HCHO C-O/C-H and SO3 S-O covalent bonds; all other pairs including within-molecule nonbonded pairs checked, MIC.'}
def rotation(axis,theta):
 axis=np.asarray(axis,float);axis=axis/np.linalg.norm(axis);a,b,c=axis;K=np.array([[0,-c,b],[c,0,-a],[-b,a,0]])
 return np.eye(3)*math.cos(theta)+(1-math.cos(theta))*np.outer(axis,axis)+math.sin(theta)*K

def rotate_water(x,o,R):
 y=x.copy()
 for h in [o-2,o-1]:y[h]=x[o]+R@mic(x[h]-x[o])
 return y

def water_geometry(x,o):
 vectors=[mic(x[h]-x[o]) for h in [o-2,o-1]]
 return {'OH_lengths_A':[float(np.linalg.norm(v)) for v in vectors],'HOH_angle_deg':sol.angle(*vectors),'HH_distance_A':float(np.linalg.norm(mic(x[o-2]-x[o-1])))}
def graph_edges(x):return {name:net.edges_from_geometry(net.bond_geometry(x,BOX),threshold) for name,threshold in net.THRESHOLDS.items()}
def donor(edges,o):return [e for e in edges if e[0]==o+1 and e[2]==290]
def write_xyz(path,block,x):
 lines=block.splitlines();old=np.array([[float(v) for v in line.split()[1:]] for line in lines[2:]])
 for i in range(296):
  if not np.array_equal(x[i],old[i]):lines[i+2]=' %s % .12f % .12f % .12f'%(ORDER[i],*x[i])
 Path(path).write_text('\n'.join(lines)+'\n')
def frame_blocks(path):
 lines=path.read_text().splitlines(keepends=True)
 if len(lines)%298:raise ValueError('Incomplete source XYZ')
 return [''.join(lines[i:i+298]) for i in range(0,len(lines),298)]
def scan(out):
 out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
 if (out/'FEASIBILITY.json').exists():raise ValueError('Refuse overwrite existing feasibility')
 reference=json.loads(SOLVATION.read_text());assert sha(SOURCE)==reference['sha256']
 frames=list(sol.frames(SOURCE));blocks=frame_blocks(SOURCE);assert len(frames)==501==len(reference['rows'])
 all_frames=[];attempts=[];feasible=[]
 for (step,x),r,block in zip(frames,reference['rows'],blocks):
  assert step==r['step'];c=contact(x);reasons=[]
  if step<1000:reasons.append('before_registered_1ps_offset')
  if r['SC']<3:reasons.append('SC_below_long_distance_filter')
  if not 1.1<=r['CO']<=1.3:reasons.append('CO_outside_filter')
  if any(not -3<=r[k]<=.5 for k in ['gds_C','gds_Oc']):reasons.append('C_or_Oc_GDS_outside_filter')
  if not -5<=r['gds_S']<=.5:reasons.append('S_GDS_outside_filter')
  if r['water_parent_assignment_warning']:reasons.append('parent_OH_warning')
  if c['minimum_A']<1.4:reasons.append('source_contact_below_1p4')
  record={'step':step,'SC':r['SC'],'CO':r['CO'],'gds_C':r['gds_C'],'gds_Oc':r['gds_Oc'],'gds_S':r['gds_S'],'CO_normal_cos':r['CO_normal_cos'],'CS_normal_cos':r['CS_normal_cos'],'hb_Oc':r['hb_Oc'],'contact':c,'rejections':reasons,'feasible_rotations':0};all_frames.append(record)
  if reasons:continue
  original=graph_edges(x);options=[]
  for o in sol.WATER_O:
   o=int(o);target=mic(x[289]-x[o]);length=float(np.linalg.norm(target))
   if not 2.65<=length<=3.3 or donor(original['loose'],o):continue
   for h in [o-2,o-1]:
    u=mic(x[h]-x[o]);u=u/np.linalg.norm(u);v=target/length;theta=math.acos(float(np.clip(u@v,-1,1)));deg=math.degrees(theta)
    if not 10<=deg<=60:continue
    axis=np.cross(u,v);axis/=np.linalg.norm(axis);R=rotation(axis,theta);on=rotate_water(x,o,R);off=rotate_water(x,o,R.T);cp,cm=contact(on),contact(off)
    a={'step':step,'water_O_1based':o+1,'aligned_H_1based':h+1,'theta_deg':deg,'on_min_contact_A':cp['minimum_A'],'control_min_contact_A':cm['minimum_A'],'pass':False,'rejections':[]}
    if min(cp['minimum_A'],cm['minimum_A'])<1.4:a['rejections'].append('rotated_contact_below_1p4')
    gp=graph_edges(on);gm=graph_edges(off)
    if any(not donor(gp[k],o) for k in net.THRESHOLDS) or donor(gm['loose'],o):a['rejections'].append('insufficient_three_threshold_contrast')
    a['pass']=not a['rejections'];attempts.append(a)
    if a['pass']:options.append((deg,o,h,axis,R,on,off,original,gp,gm,block,r,cp,cm))
  record['feasible_rotations']=len(options)
  if options:feasible.append((step,sorted(options,key=lambda z:z[:3])[0]))
 selected=[]
 for step,option in feasible:
  if not selected or step-selected[-1][0]>=1000:selected.append((step,option))
  if len(selected)==2:break
 pairs=[]
 for number,(step,opt) in enumerate(selected,1):
  deg,o,h,axis,R,on,off,original,gp,gm,block,r,cp,cm=opt;base=out/('frame_%05d'%step);base.mkdir()
  (base/'source.xyz').write_text(block);write_xyz(base/'on.xyz',block,on);write_xyz(base/'control.xyz',block,off)
  x=dict(frames)[step];before=water_geometry(x,o);after=[water_geometry(y,o) for y in (on,off)]
  for a in after:
   assert max(abs(a['OH_lengths_A'][i]-before['OH_lengths_A'][i]) for i in range(2))<1e-12
   assert abs(a['HOH_angle_deg']-before['HOH_angle_deg'])<1e-10
  unchanged=[i for i in range(296) if i not in (o-2,o-1)];assert np.array_equal(on[unchanged],x[unchanged]) and np.array_equal(off[unchanged],x[unchanged])
  changes={}
  for key in net.THRESHOLDS:
   ids=lambda es:set(tuple(e[:3]) for e in es);b,p,m=map(ids,[original[key],gp[key],gm[key]])
   changes[key]={'source_edges':original[key],'on_edges':gp[key],'control_edges':gm[key],'on_added':sorted(p-b),'on_removed':sorted(b-p),'control_added':sorted(m-b),'control_removed':sorted(b-m),'on_vs_control_only_on':sorted(p-m),'on_vs_control_only_control':sorted(m-p)}
  pairs.append({'pair':number,'source_step':step,'source_frame_sha256':sha(base/'source.xyz'),'source_directory':str(SOURCE.parent),'source_job_id':'8568513','source_segment_biased':False,'public_snapshot_origin':'Figure6b metaD R_inter snapshot','prepared_under_bias':True,'controlled_initial_conditions':True,'water_O_1based':o+1,'rotated_H_1based':[o-1,o],'aligned_H_1based':h+1,'axis_unit':axis.tolist(),'on_rotation_matrix':R.tolist(),'control_rotation_matrix':R.T.tolist(),'signed_angles_deg':[deg,-deg],'translation_A':[0,0,0],'water_geometry_before':before,'water_geometry_after':after,'contact_on':cp,'contact_control':cm,'source_descriptors':r,'on_descriptors':sol.describe(on,BOX),'control_descriptors':sol.describe(off,BOX),'heavyatom_and_all_other_atoms_exactly_unchanged':True,'network_edges_and_differences':changes,'xyz':{k:{'path':str(base/(k+'.xyz')),'sha256':sha(base/(k+'.xyz'))} for k in ('source','on','control')},'human_review':'PENDING','structural_verdict':None})
 result={'schema':'li296_hb_intervention_feasibility/v1','status':'FEASIBLE_FOR_DESIGN_REVIEW' if len(pairs)==2 else 'INSUFFICIENT_FEASIBLE_PAIRS','selection_rules':RULES,'source_sha256':{str(p):sha(p) for p in (SOURCE,SOLVATION,NET,Path(__file__))},'frames_screened':len(all_frames),'frames_after_filters':sum(not x['rejections'] for x in all_frames),'rotation_attempts':len(attempts),'feasible_frames':len(feasible),'all_frame_candidates':all_frames,'all_rotation_attempts':attempts,'selected_pairs':pairs,'limitations':['Same R history; separation does not establish equilibrium independence.','Equal-angle rotations do not equalize potential energy or isolate one HB edge; all changes retained.','No CP2K smoke/energy or chemical approval supplied here; dynamics would be fully unbiased from controlled ICs.'],'human_review':'PENDING','structural_verdict':None}
 dump(out/'FEASIBILITY.json',result);print(json.dumps({k:result[k] for k in ('status','frames_screened','frames_after_filters','rotation_attempts','feasible_frames')}));print([(p['source_step'],p['water_O_1based'],p['signed_angles_deg'][0],p['contact_on']['minimum_A'],p['contact_control']['minimum_A']) for p in pairs])
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();scan(a.out)
