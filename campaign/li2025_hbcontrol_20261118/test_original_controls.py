#!/usr/bin/env python3
"""Exercise source/seed/science mutation refusal on isolated package copies."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import validate_stage as gate


def run(package, out):
    package, out = Path(package).resolve(), Path(out).resolve()
    gate.demand(not out.exists(), 'test output exists; refuse overwrite')
    out.mkdir(parents=True)
    manifest=json.loads((package/'manifest.json').read_text())
    originals={r['name']:r for r in manifest['cases']}
    results=[]
    for case in manifest['candidate_names']:
        record=copy.deepcopy(originals[case])
        details=gate.verify_input(package/'inputs',record,package)
        results.append({'test':case+'_unchanged_source_pass','pass':True,'details':details})
    def refusal(name, mutate):
        base=out/name
        shutil.copytree(package,base)
        record=copy.deepcopy(originals['r1150_original'])
        mutate(base,record)
        try:
            gate.verify_input(base/'inputs',record,base)
        except (ValueError,KeyError,TypeError,OSError) as error:
            results.append({'test':name,'pass':True,'rejection':str(error)})
        else:
            raise AssertionError('mutant accepted: '+name)
    def seed(base,record):
        record['seed']+=1
        (base/'inputs/r1150_original.inp').write_text(gate.make_input(
            (base/'inputs/reference_unbiased.inp').read_text(),record['name'],record['coordinate_file'],record['seed'],record['steps']))
    refusal('different_seed_even_with_consistent_input',seed)
    def coordinate(base,record):
        for name in ('r1150_original.xyz','source_frame_1150.xyz'):
            p=base/'inputs'/name
            lines=p.read_text().splitlines();a=lines[2].split();a[1]=str(float(a[1])+0.01);lines[2]=' '.join(a)
            p.write_text('\n'.join(lines)+'\n')
        record['source']['original_frame_sha256']=gate.sha(base/'inputs/r1150_original.xyz')
    refusal('changed_source_and_candidate_with_updated_digest',coordinate)
    def audit(base,record):
        p=base/'SOURCE_AUDIT.json';d=json.loads(p.read_text());d['source_job']='1';gate.dump(p,d)
        record['source']['source_audit_sha256']=gate.sha(p)
    refusal('changed_audit_with_updated_digest',audit)
    def cutoff(base,record):
        p=base/'inputs/r1150_original.inp';text=p.read_text();gate.demand('CUTOFF 300' in text,'test cutoff absent')
        p.write_text(text.replace('CUTOFF 300','CUTOFF 301',1))
    refusal('scientific_cutoff_change',cutoff)
    def bias(base,record):
        p=base/'inputs/r1150_original.inp';p.write_text(p.read_text()+'\n&FREE_ENERGY\n&END FREE_ENERGY\n')
    refusal('added_bias_section',bias)
    refusal('erased_bias_history',lambda base,r:r.update(prepared_under_bias=False))
    refusal('identity_rotation_not_original',lambda base,r:r.update(transformation={'rotation_matrix':[[1,0,0],[0,1,0],[0,0,1]]}))
    refusal('wrong_frame',lambda base,r:r['source'].update(frame_step=2200))
    refusal('forbidden_source_path',lambda base,r:r['source'].update(trajectory='/tmp/methanol-sulfonization/source.xyz'))
    gate.dump(out/'TEST_REPORT.json',{'status':'PASS','tests':results,'count':len(results),
        'scope':'Actual unchanged package positive cases and deliberately corrupted isolated copies; no CP2K/HPC/production execution.'})
    return len(results)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();print(json.dumps({'status':'PASS','tests':run(a.package,a.out)}))
