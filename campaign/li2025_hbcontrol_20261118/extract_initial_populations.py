#!/usr/bin/env python3
"""Extract Li296 CP2K populations with explicit MD/energy/source binding.

This is a population-partition diagnostic, with no charge-transfer, bond-order,
structural or free-energy verdict. Inputs are read only; output must be new.
"""
import argparse
import hashlib
import importlib.util
import io
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PARSER = ROOT / 'campaign/li2025_hbnetwork_20261116/production_v1/validate_stage.py'
PARSER_SHA = '3df8173e5b37eb2d1aba2215026710ae1e868d60de3d29090fed8f00bfcc282a'
ORDER = ['H', 'H', 'O'] * 96 + ['C', 'O', 'H', 'H', 'S', 'O', 'O', 'O']
VALENCE = {'H': 1, 'O': 6, 'C': 4, 'S': 6}
KINDS = {'H': 1, 'O': 2, 'C': 3, 'S': 4}
GROUPS = {'HCHO_atoms_289_292': list(range(289, 293)),
          'SO3_atoms_293_296': list(range(293, 297)),
          'water_202_203_204': [202, 203, 204],
          'all_water_atoms_1_288': list(range(1, 289))}


def demand(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def number(token):
    value = float(token.replace('D', 'E').replace('d', 'e'))
    demand(math.isfinite(value), 'nonfinite population/energy')
    return value


def dump(path, value):
    path = Path(path)
    demand(not path.exists(), 'Refuse existing output: ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def read_population(lines, scheme, start_line):
    atoms = []
    total = None
    precision = 6 if scheme == 'mulliken' else 3
    rounding = 0.5 * 10 ** (-precision)
    for line_number, line in lines:
        fields = line.split()
        if fields and fields[0].isdigit():
            width = 5 if scheme == 'mulliken' else 6
            demand(len(fields) == width, 'population row width mismatch')
            index = int(fields[0])
            demand(index == len(atoms) + 1 and index <= 296, 'population atom index missing/duplicate')
            element = fields[1]
            demand(element == ORDER[index - 1] and int(fields[2]) == KINDS[element], 'population element/kind mismatch')
            values = [number(x) for x in fields[3:]]
            population, charge = values[-2:]
            reference = VALENCE[element] if scheme == 'mulliken' else values[0]
            demand(abs(reference - VALENCE[element]) < 1.e-10, 'unexpected valence reference')
            demand(abs(reference - population - charge) <= 2 * rounding + 1.e-10, 'population/charge arithmetic mismatch')
            atoms.append({'atom_1based': index, 'element': element, 'kind': int(fields[2]),
                          'reference_valence_e': reference, 'population_e': population, 'net_charge_e': charge})
        if re.search(r'^\s*#?\s*Total\s+[Cc]harge\s+', line):
            values = [number(x) for x in re.sub(r'^\s*#?\s*Total\s+[Cc]harge\s+', '', line).split()]
            demand(len(values) == (2 if scheme == 'mulliken' else 1), 'population total width mismatch')
            total = values[-1]
            break
    demand(len(atoms) == 296 and total is not None, 'incomplete 296-atom population table')
    summed = sum(a['net_charge_e'] for a in atoms)
    tolerance = (len(atoms) + 1) * rounding + 1.e-10
    demand(abs(summed - total) <= tolerance, 'atomic charge sum differs from printed total beyond rounding')
    if scheme == 'mulliken':
        demand(abs(total + 2) <= rounding + 1.e-10, 'Mulliken total does not match physical charge -2')
        demand(abs(sum(a['population_e'] for a in atoms) - values[0]) <= tolerance, 'Mulliken population sum mismatch')
    return {'scheme': scheme, 'units': 'elementary charge e; populations in electrons',
            'printed_decimals': precision, 'printed_total_charge_e': total,
            'sum_printed_atomic_charges_e': summed, 'sum_minus_printed_total_e': summed - total,
            'rounding_sum_tolerance_e': tolerance, 'physical_system_charge_e': -2,
            'printed_total_minus_physical_charge_e': total + 2,
            'total_charge_consistency': 'MATCHES_PHYSICAL_CHARGE' if abs(total + 2) <= rounding + 1.e-10 else 'PARTITION_TOTAL_RESIDUAL_RETAINED',
            'normalization_applied': False, 'line_start': start_line, 'line_end': line_number,
            'atoms': atoms,
            'group_rounding_bound_e': {name: len(indices) * rounding for name, indices in GROUPS.items()},
            'group_net_charge_e': {name: sum(atoms[i - 1]['net_charge_e'] for i in indices) for name, indices in GROUPS.items()}}


def parse_log(handle, wanted):
    """Read only through required MD block; explicitly bind t0 at MD_INI.

    CP2K prints t0 SCF, MD_INI, step1 SCF, then MD step1. Counting the two
    SCFs before MD1 as one state would associate the wrong charges with t0.
    """
    wanted = set(wanted)
    demand(wanted and all(type(n) is int and n >= 0 for n in wanted), 'invalid requested MD steps')
    end_step = max(1, max(wanted))
    tables, pending, active = {}, None, None
    result, trace = {}, []
    scf_pending = total_scf = energies = 0
    last_assigned = -1
    lines = enumerate(handle, start=1)
    for line_number, line in lines:
        demand(line.endswith('\n'), 'truncated output line in analyzed prefix')
        demand('SCF run NOT converged' not in line and not re.search(r'\bABORT\b', line), 'SCF failure/abort in analyzed prefix')
        if 'SCF run converged in' in line:
            scf_pending += 1
            total_scf += 1
        for title, scheme in [('Mulliken Population Analysis', 'mulliken'), ('Hirshfeld Charges', 'hirshfeld')]:
            if title in line:
                demand(scheme not in tables and pending is None, 'duplicate/unbound population table')
                tables[scheme] = read_population(lines, scheme, line_number)
        match = re.search(r'ENERGY\| Total FORCE_EVAL\s+\( QS \) energy \[a.u.\]:\s*(\S+)', line)
        if match:
            demand(set(tables) == {'mulliken', 'hirshfeld'} and scf_pending == 1, 'force evaluation lacks unique converged SCF/two population tables')
            demand(pending is None and active is None, 'energy block lacks preceding MD assignment')
            pending = {'energy_hartree': number(match[1]), 'energy_line': line_number,
                       'populations': tables, 'scf_ordinal_1based': total_scf}
            tables, scf_pending = {}, 0
            energies += 1
        step_match = re.search(r'MD\|\s+Step number\s+(\d+)', line)
        if 'MD_INI| MD initialization' in line or step_match:
            step = int(step_match[1]) if step_match else 0
            demand(pending is not None and active is None, 'MD assignment lacks unique preceding force evaluation')
            demand(step == last_assigned + 1 and pending['scf_ordinal_1based'] == step + 1 and energies == step + 1, 'nonsequential/ambiguous SCF-to-MD mapping')
            active, pending = pending, None
            active.update(step=step, md_marker_line=line_number, md_time_fs=0.0 if step == 0 else None)
            last_assigned = step
        time_match = re.search(r'^\s*MD\|\s+Time \[fs\]\s+(\S+)', line)
        if time_match:
            demand(active is not None, 'MD time without assigned population block')
            active['md_time_fs'] = number(time_match[1])
        energy_match = re.search(r'^\s*(MD_INI|MD)\|\s+Potential energy \[hartree\]\s+(\S+)', line)
        if energy_match:
            demand(active is not None, 'MD potential without assigned population block')
            demand(active['md_time_fs'] is not None and abs(active['md_time_fs'] - active['step']) < 1.e-6, 'MD time/step mismatch')
            md_energy = number(energy_match[2])
            demand(abs(md_energy - active['energy_hartree']) < 1.e-7, 'MD potential and FORCE_EVAL energy mismatch')
            active['md_potential_energy_hartree'] = md_energy
            active['md_potential_line'] = line_number
            if active['step'] in wanted:
                result[active['step']] = active
            if active['step'] <= 1 or active['step'] in wanted:
                trace.append({k: active[k] for k in ['step', 'scf_ordinal_1based', 'energy_line', 'md_marker_line', 'md_potential_line']})
            complete_step = active['step']
            active = None
            if complete_step == end_step:
                demand(set(result) == wanted, 'requested MD block not found')
                return result, {'scope': 'validated log prefix through MD potential of step ' + str(end_step),
                                'last_validated_step': end_step, 'SCF_evaluations': total_scf,
                                'mapping_trace': trace, 'SCFs_before_MD1': 2,
                                'whole_file_completion_revalidated': False}
    raise ValueError('requested complete MD/population block absent')


def frozen_parser():
    demand(sha(PARSER) == PARSER_SHA, 'registered read-only parser SHA mismatch')
    spec = importlib.util.spec_from_file_location('li296_source_parser', PARSER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def extract(manifest_path):
    manifest = json.loads(Path(manifest_path).read_text())
    cases = manifest['cases']
    parser = frozen_parser()
    source_paths = {Path(manifest_path).resolve(), Path(__file__).resolve(), PARSER}
    for case in cases:
        source_paths.update(Path(p).resolve() for p in case['files'].values())
    before = {str(p): sha(p) for p in source_paths}
    for case in cases:
        demand(set(case['expected_file_sha256']) == set(case['files']), 'all input paths need declared expected SHA')
        for key, path in case['files'].items():
            demand(before[str(Path(path).resolve())] == case['expected_file_sha256'][key], 'registered source SHA mismatch: ' + key)
    parsed_logs = {}
    for log in sorted({c['files']['log'] for c in cases}):
        wanted = [c['step'] for c in cases if c['files']['log'] == log]
        with Path(log).open() as handle:
            parsed_logs[log] = parse_log(handle, wanted)
    results = []
    for case in cases:
        paths = {key: Path(value) for key, value in case['files'].items()}
        step = case['step']
        blocks, mapping = parsed_logs[case['files']['log']]
        block = blocks[step]
        energies = parser.read_energy(paths['energy'])
        demand([int(row[0]) for row in energies] == list(range(len(energies))), 'energy table must contain fresh sequential steps from 0')
        demand(all(abs(row[0] - row[1]) < 1.e-6 for row in energies), 'energy time/step inconsistent')
        demand(step < len(energies), 'target energy row missing')
        frames = [frame for frame in parser.read_xyz(paths['trajectory']) if frame['step'] == step]
        demand(len(frames) == 1 and frames[0]['energy'] is not None, 'target XYZ missing/duplicated/without E')
        frame = frames[0]
        demand(abs(frame['time'] - step) < 1.e-6, 'XYZ time/step mismatch')
        demand(abs(energies[step][4] - block['energy_hartree']) < 1.e-8 and abs(frame['energy'] - block['energy_hartree']) < 1.e-8, 'population FORCE_EVAL energy differs from XYZ/ener')
        source = parser.read_xyz(paths['source_xyz'])
        demand(len(source) == 1, 'initial source XYZ must have one frame')
        delta = max(abs(a-b) for x,y in zip(source[0]['coords'],frame['coords']) for a,b in zip(x,y))
        demand(delta < 1.e-7, 'target XYZ differs from frozen source coordinates')
        results.append({'name': case['name'], 'job_id': case['job_id'], 'source_history': case['source_history'],
                        'files': case['files'], 'step': step, 'max_source_coordinate_delta_A': delta,
                        'energy_binding': {'force_eval_hartree': block['energy_hartree'], 'ener_hartree': energies[step][4], 'XYZ_hartree': frame['energy']},
                        'log_mapping': mapping, 'population_block': block})
    after = {str(p): sha(p) for p in source_paths}
    demand(before == after, 'source changed during extraction')
    by_name = {r['name']:r for r in results}
    comparisons = []
    for pair in manifest.get('comparisons', []):
        reference, changed = by_name[pair['reference']], by_name[pair['changed']]
        delta_energy = changed['energy_binding']['force_eval_hartree'] - reference['energy_binding']['force_eval_hartree']
        comparisons.append({'reference':reference['name'],'changed':changed['name'], 'delta_potential_hartree':delta_energy,
                            'delta_potential_kcal_mol':delta_energy * 627.509474,
                            'delta_group_net_charge_e':{scheme:{group:changed['population_block']['populations'][scheme]['group_net_charge_e'][group] - reference['population_block']['populations'][scheme]['group_net_charge_e'][group] for group in GROUPS} for scheme in ['mulliken','hirshfeld']},
                            'delta_group_rounding_bound_e': {scheme: {group: reference['population_block']['populations'][scheme]['group_rounding_bound_e'][group] + changed['population_block']['populations'][scheme]['group_rounding_bound_e'][group] for group in GROUPS} for scheme in ['mulliken', 'hirshfeld']},
                            'interpretation':'Different-geometry population-partition and potential-energy differences, with no charge-transfer, bond-order or barrier interpretation'})
    return {'schema':'li296_initial_population_diagnostic/v1','recorded_at':datetime.now(timezone.utc).isoformat(),
            'cases':results,'comparisons':comparisons,'source_sha256':[{'path':p,'before':h,'after':after[p]} for p,h in sorted(before.items())],
            'all_sources_unchanged':True,'physical_charge_e':-2,'population_scheme_caveat':'Mulliken and Hirshfeld are distinct partitions with different printed precision. Retain Hirshfeld total residual and rounded atomic sums; no normalization or electronic charge-transfer inference.',
            'scope':'Requested initial states only; full trajectory numerical completion is certified separately.',
            'structural_verdict':None,'charge_transfer_verdict':None,'bond_order_verdict':None,'FES_verdict':None,'human_review':'PENDING',
            'replay_argv':[sys.executable,'-B',str(Path(__file__).resolve()),*sys.argv[1:]]}


def self_test():
    def population(scheme):
        header = 'Mulliken Population Analysis\n' if scheme == 'mulliken' else 'Hirshfeld Charges\n'
        rows = []
        for index, element in enumerate(ORDER,1):
            charge = -2.0 if index == 293 else 0.0
            values = [VALENCE[element]-charge,charge] if scheme == 'mulliken' else [VALENCE[element],VALENCE[element]-charge,charge]
            rows.append(f'{index} {element} {KINDS[element]} ' + ' '.join(f'{v:.6f}' for v in values) + '\n')
        tail = '# Total charge 806.000000 -2.000000\n' if scheme == 'mulliken' else 'Total Charge -2.000\n'
        return header + ''.join(rows) + tail
    def block(step):
        energy = -10.0-step
        prefix = 'SCF run converged in 8 steps\n' + population('mulliken') + population('hirshfeld')
        prefix += f'ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: {energy}\n'
        marker = 'MD_INI| MD initialization\n' if step == 0 else f'MD| Step number {step}\nMD| Time [fs] {step}.0\n'
        return prefix + marker + ('MD_INI' if step == 0 else 'MD') + f'| Potential energy [hartree] {energy}\n'
    valid = ''.join(block(step) for step in range(3))
    records = []
    values, mapping = parse_log(io.StringIO(valid),[0,1,2])
    demand([values[n]['energy_hartree'] for n in [0,1,2]] == [-10.,-11.,-12.] and mapping['SCFs_before_MD1'] == 2, 'mapping test failed')
    records.append({'test':'initial_and_two_MD_blocks_map_to_distinct_energies','pass':True})
    residual_fixture = valid.replace('Total Charge -2.000', 'Total Charge -1.977', 1)
    residual_values, _ = parse_log(io.StringIO(residual_fixture), [0])
    residual_table = residual_values[0]['populations']['hirshfeld']
    demand(abs(residual_table['printed_total_minus_physical_charge_e'] - .023) < 1.e-10 and not residual_table['normalization_applied'], 'Hirshfeld residual was concealed')
    records.append({'test': 'Hirshfeld_partition_residual_retained_without_normalization', 'pass': True})
    malformed = {
        'missing_MD_INI':valid.replace('MD_INI| MD initialization\n','',1),
        'duplicate_MD_step':valid.replace('MD| Step number 2','MD| Step number 1'),
        'wrong_atom_identity':valid.replace('289 C 3','289 O 3',1),
        'nonfinite_charge':valid.replace('1 H 1 1.000000 0.000000','1 H 1 nan 0.000000',1),
        'wrong_total_charge':valid.replace('# Total charge 806.000000 -2.000000','# Total charge 806.000000 -1.000000',1),
        'missing_converged_SCF':valid.replace('SCF run converged in 8 steps\n','',1),
        'MD_energy_disagrees':valid.replace('MD| Potential energy [hartree] -11.0','MD| Potential energy [hartree] -12.0'),
        'truncated_target_block':valid[:valid.rfind('MD| Potential energy')],
    }
    for name,text in malformed.items():
        try:
            parse_log(io.StringIO(text),[0,1,2])
        except ValueError as exc:
            records.append({'test':name,'pass':True,'rejected_with':str(exc)})
        else:
            raise AssertionError('malformed fixture accepted: '+name)
    return {'schema':'li296_initial_population_parser_tests/v1','tests':records,'all_pass':True,'parser_sha256':sha(__file__)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    demand(bool(args.manifest) != args.self_test, 'choose manifest or self-test')
    result = self_test() if args.self_test else extract(args.manifest)
    dump(args.out,result)
    print(json.dumps({'out':str(args.out),'cases':len(result.get('cases',[])),'tests':len(result.get('tests',[])),'sources_unchanged':result.get('all_sources_unchanged')}))


if __name__ == '__main__':
    main()
