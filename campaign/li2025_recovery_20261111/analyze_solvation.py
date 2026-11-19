#!/usr/bin/env python3
"""Li Figure 6b distance, hydration and signed interface diagnostics; no bond verdict.

Strict fixed Li atom order, minimum-image geometry, all frames retained. GDS is
an operational half-density plane: center the slab in its largest water-O gap,
Gaussian-smooth O density (sigma 0.8 A), use half its maximum and the outer
crossings. This is not a molecular intrinsic surface or proof of equilibrium.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path
import numpy as np

ORDER = ['H', 'H', 'O'] * 96 + ['C', 'O', 'H', 'H', 'S', 'O', 'O', 'O']
ROLES = {'C': 288, 'Oc': 289, 'S': 292, 'Os1': 293, 'Os2': 294, 'Os3': 295}
WATER_O = np.arange(2, 288, 3)
WATER_H = np.array([i for i in range(288) if i % 3 != 2])

def frames(path):
    with open(path) as f:
        while line := f.readline():
            if not line.strip():
                continue
            if int(line) != 296:
                raise ValueError('Expected 296 atoms')
            comment = f.readline()
            rows = [f.readline().split() for _ in range(296)]
            if any(len(r) != 4 for r in rows) or [r[0] for r in rows] != ORDER:
                raise ValueError('Incomplete XYZ or wrong atom order')
            coords = np.array([[float(x) for x in r[1:]] for r in rows])
            if not np.isfinite(coords).all():
                raise ValueError('Nonfinite coordinates')
            m = re.search(r'\bi\s*=\s*(\d+)', comment)
            yield (int(m[1]) if m else 0), coords

def mic(v, box):
    return v - box * np.rint(v / box)

def angle(u, v):
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v) /
                   (np.linalg.norm(u) * np.linalg.norm(v)), -1, 1))))

def surface(x, box, sigma=0.8):
    zs = np.sort(x[WATER_O, 2] % box[2])
    gaps = np.diff(np.r_[zs, zs[0] + box[2]])
    k = np.argmax(gaps)
    origin = (zs[k] + gaps[k] / 2) % box[2]
    z = (x[:, 2] - origin) % box[2]
    grid = np.arange(0.05, box[2], 0.1)
    delta = grid[:, None] - z[WATER_O]
    delta -= box[2] * np.rint(delta / box[2])
    density = np.exp(-0.5 * (delta / sigma)**2).sum(axis=1)
    density /= np.sqrt(2*np.pi)*sigma*box[0]*box[1]
    threshold = density.max()/2
    selected = np.flatnonzero(density >= threshold)
    a, b = selected[0], selected[-1]
    if a == 0 or b == len(grid)-1:
        raise ValueError('No two resolved surfaces')
    lo = float(np.interp(threshold, density[a-1:a+1], grid[a-1:a+1]))
    hi = float(np.interp(threshold, density[b:b+2][::-1], grid[b:b+2][::-1]))
    ds = np.maximum(lo-z, z-hi)  # inside negative; outside positive on BOTH sides
    normal = 1.0 if z[288] >= (lo+hi)/2 else -1.0
    return ds, normal, {'lower': lo, 'upper': hi, 'origin': float(origin),
                        'thickness': hi-lo, 'peak_density': float(density.max())}

def describe(x, box, interface=True):
    dist = lambda i,j: float(np.linalg.norm(mic(x[j]-x[i],box)))
    result = {'SC': dist(288,292), 'CO': dist(288,289),
              'OCS_angle': angle(mic(x[289]-x[288],box),mic(x[292]-x[288],box))}
    plane = mic(x[[289,290,291]]-x[288],box)
    n = np.cross(plane[1]-plane[0],plane[2]-plane[0])
    result['carbon_pyramidal_height'] = float(abs(np.dot(plane[0],n))/np.linalg.norm(n))
    donors = {}
    h_to_owner = mic(x[WATER_H-WATER_H%3+2]-x[WATER_H],box)
    result['max_parent_water_OH'] = float(np.linalg.norm(h_to_owner,axis=1).max())
    result['water_parent_assignment_warning'] = int(np.sum(np.linalg.norm(h_to_owner,axis=1)>1.3))
    for label, idx in ROLES.items():
        if not label.startswith('O'):
            continue
        h_to_acc = mic(x[idx]-x[WATER_H],box)
        distances = np.linalg.norm(h_to_acc,axis=1)
        cosine = np.sum(h_to_acc*h_to_owner,axis=1)/(distances*np.linalg.norm(h_to_owner,axis=1))
        angles = np.degrees(np.arccos(np.clip(cosine,-1,1)))
        for suffix,cut,deg in [('',2.5,150),('_strict',2.4,160),('_loose',2.6,140)]:
            result['hb_'+label+suffix] = int(np.sum((distances<=cut)&(angles>=deg)))
        keep = (distances<=2.5)&(angles>=150)
        donors[label] = [int(i+1) for i in WATER_H[keep]]
        result['nearest_H_'+label] = float(distances.min())
        result['water_O_3p5_'+label] = int(np.sum(np.linalg.norm(mic(x[WATER_O]-x[idx],box),axis=1)<=3.5))
    oc_w = {(i-1)//3 for i in donors['Oc']}
    os_w = {(i-1)//3 for label in ('Os1','Os2','Os3') for i in donors[label]}
    result['shared_donor_waters'] = len(oc_w & os_w)
    result['donor_H_1based'] = donors
    if interface:
        signed, normal, gds = surface(x,box)
        for label,idx in ROLES.items():
            result['gds_'+label] = float(signed[idx])
        for label, idx in [('CO',289),('CS',292)]:
            v = mic(x[idx]-x[288],box)
            result[label+'_normal_cos'] = float(normal*v[2]/np.linalg.norm(v))
        result['surface'] = gds
    return result

def aggregate(rows):
    return {k: {'mean':float(np.mean(v)), 'median':float(np.median(v)),
                'p10':float(np.percentile(v,10)), 'p90':float(np.percentile(v,90)),
                **({'vapour_fraction':float(np.mean(np.array(v)>0))} if k.startswith('gds_') else {})}
            for k in rows[0] if isinstance(rows[0][k],(int,float)) and k!='step'
            for v in [[r[k] for r in rows]]}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('xyz',type=Path); p.add_argument('--out',type=Path,required=True)
    p.add_argument('--bulk',action='store_true'); args=p.parse_args()
    box=np.array([14.2296]*3 if args.bulk else [15.,15.,30.])
    rows=[dict(step=step,**describe(x,box,not args.bulk)) for step,x in frames(args.xyz)]
    if not rows: raise ValueError('Empty trajectory')
    if len(rows)>1 and any(b['step']<=a['step'] for a,b in zip(rows,rows[1:])):
        raise ValueError('Non-monotonic or duplicate steps')
    start=rows[0]['step']
    data={'source':str(args.xyz.resolve()),'sha256':hashlib.sha256(args.xyz.read_bytes()).hexdigest(),
          'box':box.tolist(),'frames':len(rows),'initial':rows[0],'all':aggregate(rows),
          'first_300fs':aggregate([r for r in rows if r['step']-start<=300]),
          'after_1000fs':aggregate([r for r in rows if r['step']-start>=1000]) if rows[-1]['step']-start>=1000 else None,
          'criteria':{'hbond':'H...O<=2.5 A, Owater-H...O>=150deg; strict 2.4/160; loose 2.6/140',
                      'surface':'periodic slab centered at largest O gap; 0.1 A grid, Gaussian sigma=0.8 A, half maximum; diagnostic estimate'},
          'interpretation':'descriptive correlated-frame geometry only; no causal or structural verdict',
          'human_review':'PENDING','rows':rows}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'frames':len(rows),'out':str(args.out),'initial':rows[0]}))

if __name__=='__main__': main()
