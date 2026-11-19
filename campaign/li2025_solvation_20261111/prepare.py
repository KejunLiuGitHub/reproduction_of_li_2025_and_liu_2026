#!/usr/bin/env python3
"""Controlled solvent/solute relaxation from the unchanged Li event10410 frame."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
reference=ROOT/'campaign/li2025_direct/inputs/interface_event_unbiased_s2.inp'
text=reference.read_text()
assert not any(x in text for x in ('EXT_RESTART','FREE_ENERGY','CONSTRAINT','PLUMED'))
seed=int(re.search(r'^\s*SEED\s+(\d+)',text,re.M)[1])
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
reference_xyz=ROOT/'campaign/li2025_direct/inputs/interface_r1_event10410.xyz'
assert sha(HERE/'inputs/event10410.xyz')==sha(reference_xyz)
records={}
for case,steps in [('smoke',20),('release_smoke',20),('relax',2000),('release',5000)]:
    inp=re.sub(r'(?m)^(\s*PROJECT)\s+\S+',r'\1 '+case,text)
    inp=inp.replace('interface_r1_event10410.xyz','event10410.xyz')
    inp=re.sub(r'(?m)^(\s*STEPS)\s+\d+',r'\g<1> '+str(steps),inp)
    is_release=case.startswith('release')
    if not is_release:
        inp=inp.replace('&MOTION\n','&MOTION\n  &FREE_ENERGY\n    &METADYN\n      USE_PLUMED .TRUE.\n      PLUMED_INPUT_FILE ./hold_sc.plumed\n    &END METADYN\n  &END FREE_ENERGY\n')
    else:
        inp=inp.replace('&END GLOBAL\n','&END GLOBAL\n&EXT_RESTART\n  RESTART_FILE_NAME relax_parent.restart\n  RESTART_DEFAULT .FALSE.\n  RESTART_POS .TRUE.\n  RESTART_VEL .TRUE.\n  RESTART_COUNTERS .TRUE.\n&END EXT_RESTART\n',1)
    path=HERE/'inputs'/(case+'.inp');path.write_text(inp)
    start=2000 if case=='release' else 20 if case=='release_smoke' else 0
    records[case]={'steps':steps,'start_step':start,
                   'expected_end_step':start+steps,
                   'biased':not is_release,'input_sha256':sha(path)}
plumed=HERE/'inputs/hold_sc.plumed'
plumed.write_text('UNITS LENGTH=A ENERGY=kcal/mol TIME=fs\nd_sc: DISTANCE ATOMS=293,289\nhold: RESTRAINT ARG=d_sc AT=1.95 KAPPA=100.0\nPRINT ARG=d_sc,hold.bias STRIDE=10 FILE=COLVAR_hold\n')
manifest={'schema':'li2025_solvation/v1','cases':records,'seed':seed,
 'source':{'input':str(reference.relative_to(ROOT)),'input_sha256':sha(reference),
           'xyz':str(reference_xyz.relative_to(ROOT)),'xyz_sha256':sha(reference_xyz),
           'coordinate_change':'none; no atom moved, added or removed'},
 'files':{str(p.relative_to(HERE)):sha(p) for p in sorted((HERE/'inputs').glob('*'))},
 'dependencies':{p:sha(HERE/p) for p in ['frozen_run_cp2k.py','frozen_analyze_event_unbiased.py']},
 'intervention':{'stage_A':'2 ps S-C harmonic restraint; solute orientation and all water coordinates free',
  'AT_A':1.95,'KAPPA_kcal_mol_A2':100.,'definition':'V=0.5*KAPPA*(SC-AT)^2',
  'parameter_origin':'project diagnostic choice; target near public P 3ps mean SC=1.96 A; not Li original protocol',
  'stage_B':'5 ps fully unbiased, inheriting A positions/velocities/counters; no PLUMED or geometric constraint',
  'thermostat':'CSVR parameters retained; random internal thermostat state not restored',
  'interpretation':'guided preparation then unbiased retention; not spontaneous reactant-to-product formation; no FES'},
 'controls':'existing direct-release s2 has identical initial positions and seed; this intervention allows both solute and water relaxation, not a solvent-only causal isolation',
 'smoke':'20 restrained steps +20 release steps in separate directories; production starts fresh at step0, not from smoke',
 'plumed_reference':'https://www.plumed.org/doc-v2.9/user-doc/html/_r_e_s_t_r_a_i_n_t.html'}
(HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'cases':records,'seed':seed}))
