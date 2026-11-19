#!/usr/bin/env python3
"""Check periodic invariance and both surface signs against actual Li input geometry."""
import json
from pathlib import Path
import numpy as np
import analyze_solvation as a
ROOT=Path(__file__).resolve().parents[2]
box=np.array([15.,15.,30.])
event=ROOT/'campaign/li2025_direct/inputs/interface_r1_event10410.xyz'
public=ROOT/'campaign/li2025_recovery_20261111/inputs/p_interface.xyz'
_,x=next(a.frames(event)); ref=a.describe(x,box)
keys=[k for k,v in ref.items() if isinstance(v,(int,float))]
for y in [(x+np.array([3.1,-17.2,12.7]))%box,x*np.array([1.,1.,-1.])]:
    got=a.describe(y,box)
    assert all(abs(ref[k]-got[k])<1e-8 for k in keys)
ds,_,g=a.surface(x,box)
for z,sign in [(g['lower']-1,1),(g['upper']+1,1),((g['lower']+g['upper'])/2,-1)]:
    trial=x.copy();trial[288,2]=(z+g['origin'])%30
    assert a.surface(trial,box)[0][288]*sign>0
_,p=next(a.frames(public)); pub=a.describe(p,box)
assert ref['hb_Oc']==2 and pub['hb_Oc']==3
assert abs(ref['SC']-1.8081340822266159)<1e-10
assert abs(pub['SC']-2.074714787994484)<1e-10
assert ref['water_O_3p5_Oc']==2 and pub['water_O_3p5_Oc']==4
report={'periodic_translation_wrapping':'PASS','slab_z_reflection':'PASS',
 'both_interfaces_signed_positive_outside':'PASS','independent_distance_crosscheck':'PASS',
 'initial_hbond_crosscheck':'PASS',
 'gds_sigma_sensitivity':{label:{str(sig):a.surface(xyz,box,sig)[0][[288,289,292]].tolist()
    for sig in (.6,.8,1.)} for label,xyz in [('event',x),('public',p)]}}
out=ROOT/'results/li2025_solvation_20261111/GEOMETRY_CHECKS.json'
out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n')
print('PASS: '+str(out))
