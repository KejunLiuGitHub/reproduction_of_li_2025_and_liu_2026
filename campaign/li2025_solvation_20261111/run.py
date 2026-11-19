#!/usr/bin/env python3
"""Run smoke -> restrained 2ps -> unbiased 5ps with strict coverage and SCF checks."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
import frozen_run_cp2k as cp2k

BASE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def write(path,data): path.write_text(json.dumps(data,indent=2)+'\n')

def validate(case,work,rc,elapsed,manifest):
    spec=manifest['cases'][case]; text=(work/(case+'.out')).read_text()
    parsed=cp2k.parse_output(text)
    trajectory=work/(case+'-pos-1.xyz')
    steps=cp2k.trajectory_steps(trajectory)
    expected=list(range(spec['start_step']+(10 if case.startswith('release') else 0),spec['expected_end_step']+1,10))
    checks={'program_ended':parsed['program_ended'],'returncode':rc==0,
        'all_scf_converged':parsed['scf']['not_converged_runs']==0 and parsed['scf']['converged_runs']>0,
        'md_steps':parsed['md_step_first']==spec['start_step']+1 and parsed['md_step_last']==spec['expected_end_step'],
        'trajectory_steps':steps==expected and cp2k.count_frames(trajectory)==len(expected),
        'energy_step':cp2k.energy_last_step(work/(case+'-1.ener'))==spec['expected_end_step'],
        'restart_step':cp2k.restart_md_start_step(work/(case+'-1.restart'))==spec['expected_end_step']}
    if spec['biased']:
        rows=[line.split() for line in (work/'COLVAR_hold').read_text().splitlines() if line.strip() and not line.startswith('#')]
        checks['restraint_output']=len(rows)>=spec['steps']//10 and all(len(r)>=3 for r in rows)
        checks['restraint_formula']=all(abs(float(r[2])-50*(float(r[1])-1.95)**2)<0.002 for r in rows)
    result=dict(case=case,job_id=os.environ.get('SLURM_JOB_ID'),status='COMPLETED' if all(checks.values()) else 'FAILED',
        biased=spec['biased'],prepared_under_bias=case.startswith('release'),checks=checks,wall_seconds=elapsed,
        input_sha256={p.name:sha(p) for p in work.iterdir() if p.suffix in ('.inp','.plumed')},
        trajectory_sha256=sha(trajectory),**parsed)
    write(work/'run_record.json',result)
    if not all(checks.values()): raise RuntimeError('Stage failed checks: '+json.dumps(checks))
    return result

def run_stage(case,manifest,ntasks):
    work=BASE/case
    if work.exists(): raise RuntimeError('Refuse to overwrite existing stage '+str(work))
    work.mkdir()
    for filename in (case+'.inp','event10410.xyz'):
        shutil.copy2(BASE/'inputs'/filename,work/filename)
    if not case.startswith('release'): shutil.copy2(BASE/'inputs/hold_sc.plumed',work/'hold_sc.plumed')
    else:
        parent=BASE/('smoke' if case=='release_smoke' else 'relax'); record=json.loads((parent/'run_record.json').read_text())
        assert record['status']=='COMPLETED' and all(record['checks'].values())
        restart=parent/(parent.name+'-1.restart')
        assert cp2k.restart_md_start_step(restart)==manifest['cases'][case]['start_step']
        shutil.copy2(restart,work/'relax_parent.restart')
        body=(work/(case+'.inp')).read_text().upper()
        assert not any(x in body for x in ('PLUMED','FREE_ENERGY','CONSTRAINT','RESTRAINT'))
        write(work/'parent_provenance.json',{'restart_sha256':sha(work/'relax_parent.restart'),
             'run_record_sha256':sha(parent/'run_record.json'),'source':str(parent),'start_step':manifest['cases'][case]['start_step']})
    env=cp2k.cp2k_environment();cmd=cp2k.command_for(case,ntasks,env)
    start=time.time()
    with (work/'runner.log').open('w') as log:
        rc=subprocess.call(cmd,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT)
    validate(case,work,rc,time.time()-start,manifest)
    out=work/'DISTANCE_ANALYSIS.json'
    subprocess.run([sys.executable,str(BASE/'frozen_analyze_event_unbiased.py'),'--trajectory',str(work/(case+'-pos-1.xyz')),
        '--energy',str(work/(case+'-1.ener')),'--label',case,'--box','15,15,30','--out',str(out)],check=True)
    data=json.loads(out.read_text());data['source']=manifest['source']
    data['job_id']=os.environ.get('SLURM_JOB_ID')
    data['interpretation']={'verdict':None,'biased':not case.startswith('release'),'prepared_under_bias':case.startswith('release'),
        'human_review':'PENDING','smoke':'smoke' in case,'what_this_is':manifest['intervention']['interpretation'],
        'stage_time_ps':(manifest['cases'][case]['expected_end_step']-manifest['cases'][case]['start_step'])/1000}
    write(out,data)
    print('STAGE_COMPLETED',case,flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--dry-run',action='store_true');p.add_argument('--ntasks',type=int,default=64)
    args=p.parse_args();manifest=json.loads((BASE/'manifest.json').read_text())
    for name,digest in {**manifest['files'],**manifest['dependencies']}.items():
        assert sha(BASE/name)==digest, 'Hash mismatch: '+name
    if args.dry_run:
        print(json.dumps({'hash_check':'PASS','sequence':['smoke','release_smoke','relax','release'],'ntasks':args.ntasks}));return
    if any((BASE/c).exists() for c in manifest['cases']): raise RuntimeError('Existing stage directories; refuse duplicate workflow')
    for case in ('smoke','release_smoke','relax','release'): run_stage(case,manifest,args.ntasks)

if __name__=='__main__': main()
