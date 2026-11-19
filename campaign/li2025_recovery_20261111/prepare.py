#!/usr/bin/env python3
"""Freeze public Figure 6b snapshots and reuse the audited Li runner, unchanged."""
import hashlib
import json
import shutil
import sys
from pathlib import Path
from openpyxl import load_workbook

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LEGACY = HERE.parent / 'li2025_direct'
sys.path.insert(0, str(LEGACY))
import prepare_li2025_direct as source
import prepare_event_unbiased as unbiased

CASES = [('p_interface', 26, 'P_inter', 'interface', 3000, 2026111101),
         ('p_bulk', 11, 'P_bulk', 'bulk', 3000, 2026111102),
         ('r_interface_eq', 16, 'R_inter', 'interface', 5000, 2026111103)]

def main():
    assert source.md5(source.SOURCE) == source.SOURCE_MD5
    ws = load_workbook(source.SOURCE, read_only=True, data_only=True)[source.SHEET]
    inputs = HERE / 'inputs'
    inputs.mkdir(exist_ok=True)
    records = []
    for name, col, label, env, steps, seed in CASES:
        assert ws.cell(1, col).value == label, (label, ws.cell(1, col).value)
        elements, coords = source.read_block(ws, col)
        box = (15., 15., 30.) if env == 'interface' else (14.2296,) * 3
        metrics = source.verify_structure(elements, coords, box, label)
        text = source.xyz_text(elements, coords, 'Li2025 Figure 6b %s; q=-2; original coordinates' % label)
        xyz = inputs / (name + '.xyz')
        xyz.write_text(text)
        # Independent cross-check against the previously archived extraction.
        old = ROOT / 'results/li2025_zenodo/structure_compare' / ('li2025_fig6b_%s_%s.xyz' % ('RC' if label.startswith('R') else 'P', 'iface' if env == 'interface' else 'bulk'))
        old_atoms = [line.split() for line in old.read_text().splitlines()[2:]]
        assert [a[0] for a in old_atoms] == elements
        assert max(abs(float(a[j+1])-c[j]) for a,c in zip(old_atoms,coords) for j in range(3)) < 1e-12
        assert max(abs(float('%.6f' % v)-v) for c in coords for v in c) < 1e-12
        ref = (LEGACY / 'inputs' / (env + '_r1.inp')).read_text()
        for suffix, nsteps in [('', steps), ('_smoke', 20)]:
            case = name + suffix
            inp = inputs / (case + '.inp')
            body, removed = unbiased.build_unbiased_input(ref, seed, name, project=case,
                steps=nsteps, reference_project=env+'_r1', reference_coord=env+'_r1.xyz')
            inp.write_text(body)
            records.append(dict(name=case, unbiased=True, cp2k_project=case,
                box_class=env, box_angstrom=list(box), seed=seed, steps=nsteps,
                expected_frames=nsteps//10+1, sc_distance_angstrom=metrics['sc_distance_angstrom'],
                source=dict(sheet=source.SHEET, column=col, rows=[2,297], label=label),
                geometry=metrics, stage=[dict(role='new', source=p.name, target=p.name,
                    sha256=source.sha256(p), bytes=p.stat().st_size) for p in [xyz,inp]],
                outputs=dict(output=case+'.out', trajectory=case+'-pos-1.xyz',
                    energy=case+'-1.ener', restart=case+'-1.restart')))
    dependencies = {}
    for filename in ['run_cp2k.py', 'analyze_event_unbiased.py']:
        shutil.copy2(LEGACY / filename, HERE / ('frozen_'+filename))
        dependencies['frozen_'+filename] = source.sha256(LEGACY / filename)
    manifest = dict(schema='li2025_recovery/v1', cases=records, dependencies=dependencies,
        source=dict(workbook=str(source.SOURCE.relative_to(ROOT)), md5=source.SOURCE_MD5,
                    sha256=source.sha256(source.SOURCE), doi=source.PAPER_DOI, zenodo=source.ZENODO_DOI),
        policy='First diagnostic replicas only; not spontaneous formation, FES, or confirmed product. R equilibration requires later stationarity/interface review. Smoke never enters scientific statistics.')
    (HERE/'pilot_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({r['name']:r['geometry']['sc_distance_angstrom'] for r in records},indent=2))

if __name__ == '__main__':
    main()
