#!/usr/bin/env python3
"""Freeze one independent 40→50 ps continuation after technical integrity gates.

These preparation/run checks make no structural, kinetic or FES verdict.
The frozen runner and all old production directories remain byte-for-byte intact.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CASE = 'interface_r1_seg5'
PROJECT = 'interface_r1'
FROZEN_SHA = 'c3d65060daafbe1076304f092f72c2e40e738aea9cb6f4cc63d8263298397d61'
DEFAULT_SOURCE = ROOT / 'results/li2025_autofollow/20261115T100111+0800/li2025_extend_20261113/interface_r1_seg4'
DEFAULT_HISTORY = ROOT / 'results/li2025_autofollow/20261113T213847+0800/li2025_extend_20261111/interface_r1_seg3'
DEFAULT_AUDIT = ROOT / 'results/li2025_autofollow/20261115T100111+0800/continuation_validation/PARENT_CHECKS.json'
DEFAULT_SOURCE_SHA = DEFAULT_SOURCE.parents[1] / 'RETRIEVAL_METADATA.json'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def dump(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def demand(condition, label):
    if not condition:
        raise ValueError(label)


def numbers(tokens):
    values = [float(token.replace('D', 'E').replace('d', 'e')) for token in tokens]
    demand(all(math.isfinite(value) for value in values), 'nonfinite numeric data')
    return values


def rows(path, width):
    data = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split()
        demand(len(fields) == width, '%s: incorrect column count' % path)
        data.append(numbers(fields))
    demand(bool(data), '%s: empty table' % path)
    return data


def read_xyz(path, expected_symbols=None):
    frames = []
    with Path(path).open() as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            demand(header.strip() == '296', 'XYZ atom count must be 296')
            comment = handle.readline()
            demand(bool(comment), 'truncated XYZ comment')
            atoms = []
            for _ in range(296):
                fields = handle.readline().split()
                demand(len(fields) == 4, 'truncated or malformed XYZ atom row')
                atoms.append((fields[0], numbers(fields[1:])))
            symbols = [atom[0] for atom in atoms]
            demand(expected_symbols is None or symbols == expected_symbols,
                   'XYZ atom order differs from frozen initial XYZ')
            step = re.search(r'\bi\s*=\s*(\d+)', comment)
            time = re.search(r'\btime\s*=\s*([^,\s]+)', comment)
            energy = re.search(r'\bE\s*=\s*([^,\s]+)', comment)
            frames.append({'step': int(step[1]) if step else None,
                           'time': numbers([time[1]])[0] if time else None,
                           'energy': numbers([energy[1]])[0] if energy else None,
                           'atoms': atoms})
    demand(bool(frames), 'empty XYZ')
    return frames


def read_restart(path, symbols, expected_step, final_atoms=None):
    """Check section balance and all restarted POS/VEL/MD counters explicitly.

    CP2K 9.1 prints unused native-METADYN sentinel NaN values even when using
    PLUMED. They are not restarted by POS/VEL/COUNTERS and are not interpreted.
    """
    stack, blocks = [], {}
    for raw in Path(path).read_text().splitlines():
        line = raw.split('#', 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        if line.startswith('&'):
            tag = fields[0][1:].upper()
            if tag == 'END':
                demand(bool(stack), 'restart unmatched END')
                demand(len(fields) <= 2 and (len(fields) == 1 or fields[1].upper() == stack[-1]),
                       'restart mismatched END')
                stack.pop()
            else:
                stack.append(tag)
                key = tuple(stack)
                demand(key not in blocks or 'KIND' in key, 'duplicate restart section: ' + '/'.join(key))
                blocks.setdefault(key, [])
            continue
        demand(bool(stack), 'restart data outside a section')
        blocks[tuple(stack)].append(fields)
    demand(not stack, 'truncated restart: open sections')
    coord = blocks.get(('FORCE_EVAL', 'SUBSYS', 'COORD'), [])
    vel = blocks.get(('FORCE_EVAL', 'SUBSYS', 'VELOCITY'), [])
    demand(len(coord) == len(vel) == 296, 'restart needs 296 coordinates and velocities')
    demand(all(len(row) == 4 for row in coord) and [row[0] for row in coord] == symbols,
           'restart coordinate atom order')
    coords = [numbers(row[1:]) for row in coord]
    demand(all(len(row) == 3 for row in vel), 'restart velocity row width')
    velocities = [numbers(row) for row in vel]
    md_rows = blocks.get(('MOTION', 'MD'), [])
    md = {row[0].upper(): row[1:] for row in md_rows}
    demand(len(md) == len(md_rows), 'duplicate restart MD keyword')
    for key in ('STEP_START_VAL', 'TIME_START_VAL', 'TIMESTEP', 'STEPS', 'ECONS_START_VAL'):
        demand(len(md.get(key, [])) == 1, 'restart missing/ambiguous ' + key)
        numbers(md[key])
    demand(float(md['STEP_START_VAL'][0]) == expected_step, 'restart step mismatch')
    demand(abs(float(md['TIME_START_VAL'][0]) - expected_step) < 1e-5,
           'restart time counter mismatch')
    demand(abs(float(md['TIMESTEP'][0]) - 1.) < 1e-12 and int(md['STEPS'][0]) == 10000,
           'restart timestep/segment length changed')
    demand(md.get('ENSEMBLE') == ['NVT'], 'restart ensemble changed')
    maximum = None
    if final_atoms is not None:
        maximum = max(abs(x - y) for values, atom in zip(coords, final_atoms)
                      for x, y in zip(values, atom[1]))
        demand(maximum < 1e-7, 'restart coordinates disagree with final trajectory frame')
    return {'step': expected_step, 'time_fs': float(md['TIME_START_VAL'][0]),
            'coordinate_count': len(coords), 'velocity_count': len(velocities),
            'max_final_frame_coordinate_difference_angstrom': maximum,
            'sha256': sha(path)}


def bias_check(path, kind, end, time_per_step, prefix=None):
    width, stride = (5, 30) if kind == 'hills' else (3, 10)
    data = rows(path, width)
    expected = list(range(stride, end + 1, stride))
    demand(len(data) == len(expected), kind + ': incomplete cumulative row count')
    demand(all(abs(row[0] - step * time_per_step) < 2e-6 for row, step in zip(data, expected)),
           kind + ': time sequence mismatch')
    if kind == 'hills':
        demand(all(abs(row[2] - 0.1) < 1e-12 and abs(row[3] - 0.03137545) < 1e-12
                   and row[4] == -1 for row in data), 'HILLS parameters changed')
    prefix_rows = 0
    if prefix is not None:
        prior = rows(prefix, width)
        demand(Path(path).read_bytes().startswith(Path(prefix).read_bytes()),
               kind + ': parent byte prefix modified')
        demand(data[:len(prior)] == prior, kind + ': parent numeric prefix modified')
        prefix_rows = len(prior)
    return {'rows': len(data), 'last_time_fs': data[-1][0],
            'parent_rows': prefix_rows, 'new_rows': len(data) - prefix_rows,
            'sha256': sha(path), 'byte_prefix_preserved': prefix is not None}


def validate_segment(work, case, start, end, expected_job, reference_xyz, history,
                     expected_input_hashes, expected_source_hashes, audit_out, time_per_step):
    """One fail-closed technical preparation/run gate; persist failure evidence."""
    work, history = Path(work), Path(history)
    audit = {'schema': 'li2025_continuation_technical_checks/v1',
             'recorded_at': datetime.now(timezone.utc).isoformat(), 'case': case,
             'job_id': str(expected_job), 'source_directory': str(work),
             'step_range': [start + 1, end], 'status': 'FAIL', 'checks': {},
             'structural_verdict': None, 'fes_verdict': None, 'biased': True,
             'interpretation': 'Technical completeness only; ordinary biased metaD.'}
    checked_files = set()

    def check(label, action):
        try:
            result = action()
            audit['checks'][label] = {'pass': True, 'details': result}
            return result
        except (ValueError, KeyError, OSError, TypeError) as error:
            audit['checks'][label] = {'pass': False, 'error': str(error)}
            raise

    def hashes(expected):
        observed = {}
        for name, digest in expected.items():
            observed[name] = sha(work / name)
            checked_files.add(name)
            demand(observed[name] == digest, 'source SHA256 mismatch: ' + name)
        return observed

    try:
        check('frozen_source_hashes', lambda: hashes(expected_source_hashes))
        run = json.loads((work / 'run_record.json').read_text())
        checked_files.add('run_record.json')
        def terminal():
            demand(run['case'] == case and str(run['job_id']) == str(expected_job), 'run identity mismatch')
            demand(run['status'] == 'COMPLETED' and run['returncode'] == 0
                   and run['program_ended'] and run['abort_count'] == 0, 'runner not normally completed')
            demand(run['md_step_first'] == start + 1 and run['md_step_last'] == end,
                   'run record step mismatch')
            return {'returncode': run['returncode'], 'program_ended': run['program_ended']}
        check('runner_terminal', terminal)
        def inputs():
            observed = hashes(expected_input_hashes)
            demand(all(run['input_sha256'].get(name) == digest for name, digest in observed.items()),
                   'run_record staged hash mismatch')
            return observed
        check('input_provenance', inputs)
        source = read_xyz(reference_xyz)
        demand(len(source) == 1, 'initial XYZ must have one frame')
        symbols = [atom[0] for atom in source[0]['atoms']]
        out_name = case + '.out'
        output = (work / out_name).read_text()
        checked_files.add(out_name)
        def log():
            md_matches = list(re.finditer(r'MD\|\s+Step number\s+(\d+)', output))
            steps = [int(match[1]) for match in md_matches]
            demand(steps == list(range(start + 1, end + 1)), 'MD log full step coverage')
            converged = list(re.finditer(r'SCF run converged in\s+(\d+)\s+steps', output))
            failures = len(re.findall(r'SCF run NOT converged', output))
            demand(failures == 0 and len(converged) == end - start + 1, 'SCF failure/count mismatch')
            # CP2K writes one initial SCF, then one SCF before each MD step report.
            bins, ci, previous = [], 0, 0
            for match in md_matches:
                while ci < len(converged) and converged[ci].start() < match.start():
                    ci += 1
                bins.append(ci - previous)
                previous = ci
            demand(bins == [2] + [1] * (end - start - 1) and ci == len(converged),
                   'SCF coverage per MD step mismatch')
            demand(output.count('PROGRAM ENDED AT') == 1 and not re.search(r'\bABORT\b', output),
                   'log lacks unique normal termination or has ABORT')
            demand(run['scf']['converged_runs'] == len(converged)
                   and run['scf']['not_converged_runs'] == failures, 'record SCF disagrees with actual log')
            return {'md_steps': len(steps), 'scf_converged': len(converged),
                    'scf_failed': failures, 'initial_scf': 1, 'per_md_step_scf': 1}
        check('all_md_and_scf_log_coverage', log)
        trajectory_name = PROJECT + '-pos-1.xyz'
        frames = check('complete_296_atom_frames', lambda: read_xyz(work / trajectory_name, symbols))
        # Avoid embedding all coordinates in the audit.
        audit['checks']['complete_296_atom_frames']['details'] = {'frames': len(frames), 'atoms_per_frame': 296}
        checked_files.add(trajectory_name)
        def trajectory():
            expected = list(range(start + 10, end + 1, 10))
            demand([frame['step'] for frame in frames] == expected, 'trajectory missing/duplicate/out-of-order frame')
            demand(all(frame['time'] is not None and frame['energy'] is not None
                       and abs(frame['time'] - frame['step']) < 0.0011 for frame in frames),
                   'trajectory time/energy comment missing or inconsistent')
            demand(run['trajectory']['frames'] == len(frames)
                   and run['trajectory']['first_step'] == expected[0]
                   and run['trajectory']['last_step'] == end, 'record trajectory mismatch')
            return {'frames': len(frames), 'first_step': expected[0], 'last_step': end, 'stride': 10}
        check('trajectory_sequence', trajectory)
        energy_name = PROJECT + '-1.ener'
        energy = rows(work / energy_name, 7)
        checked_files.add(energy_name)
        def energies():
            demand([row[0] for row in energy] == list(range(start + 1, end + 1)), 'energy step coverage')
            demand(all(abs(row[1] - row[0]) < 2e-6 for row in energy), 'energy time coverage')
            demand(all(abs(energy[frame['step'] - start - 1][4] - frame['energy']) < 2e-8
                       for frame in frames), 'XYZ/energy potential disagreement')
            return {'rows': len(energy), 'first_step': int(energy[0][0]), 'last_step': int(energy[-1][0])}
        check('energy_all_steps_and_finite', energies)
        for kind, name in [('hills', 'HILLS_' + PROJECT), ('colvar', 'COLVAR_' + PROJECT)]:
            result = check(kind + '_full_history', lambda k=kind, n=name: bias_check(
                work / n, k, end, time_per_step, history / n))
            stride = 30 if kind == 'hills' else 10
            demand(result['parent_rows'] == start // stride
                   and result['new_rows'] == end // stride - start // stride,
                   kind + ': wrong parent/new row count')
            demand(run['metad'][kind + '_rows'] == result['rows']
                   and abs(run['metad'][kind + '_last_time_fs'] - result['last_time_fs']) < 1e-8,
                   kind + ': run record differs from actual table')
            checked_files.add(name)
        for name in (PROJECT + '-1.restart', PROJECT + '-1_%d.restart' % end):
            check('checkpoint_' + name, lambda n=name: read_restart(work / n, symbols, end, frames[-1]['atoms']))
            checked_files.add(name)
        def matching_restart_payload():
            paths = [work / (PROJECT + '-1.restart'), work / (PROJECT + '-1_%d.restart' % end)]
            payloads = ['\n'.join(line for line in path.read_text().splitlines()
                                  if not line.lstrip().startswith('#')) for path in paths]
            demand(payloads[0] == payloads[1], 'numbered/final checkpoint state payloads differ')
            return {'same_noncomment_payload': True,
                    'note': 'CP2K writes distinct millisecond timestamps; both full-file SHA values retained.'}
        check('final_numbered_checkpoint_consistency', matching_restart_payload)
        audit['status'] = 'PASS'
    except (ValueError, KeyError, OSError, TypeError) as error:
        audit['error'] = str(error)
    audit['files_sha256'] = {name: sha(work / name) for name in sorted(checked_files) if (work / name).is_file()}
    dump(audit_out, audit)
    if audit['status'] != 'PASS':
        raise ValueError('Technical gate FAIL: ' + audit['error'])
    return audit


def retrieval_hashes(source_sha, required):
    """Resolve this stable segment against local hashes and raw remote SHA lines."""
    source_sha = Path(source_sha)
    metadata = json.loads(source_sha.read_text())
    local_index = json.loads((source_sha.parent / 'HASHES.json').read_text())
    prefix = 'li2025_extend_20261113/interface_r1_seg4/'
    demand(metadata['classification'] == 'STABLE_TERMINAL', 'source is not stable terminal')
    entries = [entry for entry in metadata['remote_local_sha'] if entry['file'].startswith(prefix)]
    demand(entries and all(entry['match'] is True for entry in entries), 'remote/local source check failed')
    hashes = {}
    for entry in entries:
        name = entry['file'][len(prefix):]
        demand(Path(name).name == name and name not in hashes, 'duplicate/unsafe remote SHA path')
        digest = entry['remote_sha256']
        demand(re.fullmatch('[0-9a-f]{64}', digest) is not None, 'malformed remote SHA')
        demand(local_index[entry['file']]['sha256'] == digest, 'local/remote SHA index disagreement')
        evidence_name = entry['remote_digest_evidence']
        demand(Path(evidence_name).name == evidence_name, 'unsafe remote evidence path')
        evidence_path = source_sha.parent / evidence_name
        demand(sha(evidence_path) == local_index[evidence_name]['sha256'], 'remote evidence SHA mismatch')
        lines = evidence_path.read_text().splitlines()
        line_number = entry['remote_digest_line']
        demand(1 <= line_number <= len(lines), 'remote SHA line outside evidence')
        demand(lines[line_number - 1].split() == [digest, entry['file']]
               and entry['remote_digest_lookup_key'] == entry['file'], 'raw remote SHA evidence disagreement')
        hashes[name] = digest
    demand(required <= set(hashes), 'remote source SHA index misses required evidence')
    return hashes


def verify_parent_checks(source, remote_hashes, required):
    """Require the saved parent PASS to describe exactly the retrieved files."""
    source = Path(source)
    path = source / 'CHECKS.json'
    demand(sha(path) == remote_hashes['CHECKS.json'], 'parent CHECKS source SHA mismatch')
    record = json.loads(path.read_text())
    demand(record['schema'] == 'li2025_continuation_technical_checks/v1'
           and record['case'] == 'interface_r1_seg4' and str(record['job_id']) == '8571154'
           and record['step_range'] == [30001, 40000], 'parent CHECKS identity/range mismatch')
    demand(record['status'] == 'PASS' and record['checks']
           and all(item['pass'] is True for item in record['checks'].values()), 'parent CHECKS is not all PASS')
    hashes = record['files_sha256']
    demand(required - {'CHECKS.json'} <= set(hashes), 'parent CHECKS omits required source hashes')
    for name, digest in hashes.items():
        demand(Path(name).name == name and remote_hashes.get(name) == digest,
               'parent CHECKS/remote SHA mismatch: ' + name)
        demand(sha(source / name) == digest, 'parent CHECKS/local source SHA mismatch: ' + name)
    return {'status': 'PASS', 'case': record['case'], 'job_id': str(record['job_id']),
            'recorded_at': record['recorded_at'], 'sha256': sha(path),
            'all_check_items_pass': True, 'source_files_sha256_verified': hashes}


def prepare(source, history, audit_out, validate_only=False, source_sha=DEFAULT_SOURCE_SHA):
    source, history = Path(source).resolve(), Path(history).resolve()
    old_manifest = json.loads((HERE.parent / 'li2025_extend_20261113/continuation_manifest.json').read_text())
    old = old_manifest['cases'][0]
    expected_inputs = {entry['target']: entry['sha256'] for entry in old['stage']
                       if entry['role'] not in ('parent_hills', 'parent_colvar')}
    # Prefix files must themselves be exactly the historical files frozen for seg4.
    for entry in old['stage']:
        if entry['role'] in ('parent_hills', 'parent_colvar'):
            demand(sha(history / entry['source']) == entry['sha256'], 'historical bias SHA mismatch')
    initial_xyz = source / (PROJECT + '.xyz')
    time_scale = rows(history / ('HILLS_' + PROJECT), 5)[0][0] / 30.
    # Verify historical cumulative sequences before using them as trusted prefixes.
    for kind, name in [('hills', 'HILLS_' + PROJECT), ('colvar', 'COLVAR_' + PROJECT)]:
        bias_check(history / name, kind, 30000, time_scale)
    required = set(expected_inputs) | {'CHECKS.json', 'run_record.json', 'interface_r1_seg4.out',
        PROJECT + '-pos-1.xyz', PROJECT + '-1.ener', PROJECT + '-1.restart',
        PROJECT + '-1_40000.restart', 'HILLS_' + PROJECT, 'COLVAR_' + PROJECT}
    remote_hashes = retrieval_hashes(source_sha, required)
    formal = verify_parent_checks(source, remote_hashes, required)
    audit = validate_segment(source, 'interface_r1_seg4', 30000, 40000, '8571154',
                             initial_xyz, history, expected_inputs, remote_hashes, audit_out, time_scale)
    audit['parent_formal_checks'] = formal
    audit['retrieval_evidence_sha256'] = {name: sha(Path(source_sha).parent / name)
        for name in ('RETRIEVAL_METADATA.json', 'HASHES.json', 'REMOTE_EVIDENCE.txt')}
    dump(audit_out, audit)
    if validate_only:
        return audit
    demand(sha(HERE / 'frozen_run_cp2k.py') == FROZEN_SHA, 'frozen runner differs from seg4')
    old_input = (source / 'interface_r1_seg4.inp').read_text()
    demand(old_input.count('interface_r1_seg4_parent.restart') == 1, 'parent restart reference missing/ambiguous')
    new_input = old_input.replace('interface_r1_seg4_parent.restart', CASE + '_parent.restart')
    (HERE / 'inputs').mkdir(exist_ok=True)
    (HERE / 'inputs' / (CASE + '.inp')).write_text(new_input)
    (HERE / 'inputs' / initial_xyz.name).write_bytes(initial_xyz.read_bytes())
    run = json.loads((source / 'run_record.json').read_text())
    record = {'name': CASE, 'parent_case': 'interface_r1_seg4', 'parent_job_id': '8571154',
              'parent_workdir_remote': run['workdir'], 'cp2k_project': PROJECT,
              'box_class': 'interface', 'box_angstrom': [15., 15., 30.],
              'segment': {'start_step': 40000, 'steps': 10000, 'expected_first_step': 40001,
                          'expected_end_step': 50000, 'parent_md_step_last': 40000,
                          'segment_time_range_fs': [40000, 50000], 'cumulative_step_range': [1, 50000]},
              'outputs': {'output': CASE + '.out', 'trajectory': PROJECT + '-pos-1.xyz',
                          'energy': PROJECT + '-1.ener', 'restart': PROJECT + '-1.restart'},
              'metad_files': {'hills': 'HILLS_' + PROJECT, 'colvar': 'COLVAR_' + PROJECT},
              'metad_baseline': {key: run['metad'][key] for key in
                                 ('hills_rows', 'hills_last_time_fs', 'colvar_rows', 'colvar_last_time_fs')},
              'parent_evidence': {'local_source': str(source.relative_to(ROOT)),
                                  'remote_sha_check_path': str(Path(source_sha).resolve().relative_to(ROOT)),
                                  'remote_sha_check_sha256': sha(source_sha),
                                  'technical_checks_path': str(Path(audit_out).resolve().relative_to(ROOT)),
                                  'technical_checks_sha256': sha(audit_out),
                                  'source_files_sha256': audit['files_sha256'],
                                  'formal_checks_sha256': formal['sha256'],
                                  'retrieval_evidence_sha256': audit['retrieval_evidence_sha256'],
                                  'technical_checks_status': 'PASS'},
              'technical_policy': {'plumed_time_fs_per_step': time_scale,
                                   'plumed_time_tolerance_fs': 2e-6,
                                   'trajectory_first_step': 40010, 'trajectory_frames': 1000,
                                   'md_steps': 10000, 'energy_rows': 10000, 'scf_converged': 10001,
                                   'scf_failed': 0, 'hills_new_rows': 333, 'colvar_new_rows': 1000},
              'stage': []}
    for role, name, target in [
            ('new_input', CASE + '.inp', CASE + '.inp'),
            ('new_xyz', initial_xyz.name, initial_xyz.name),
            ('parent_restart', PROJECT + '-1_40000.restart', CASE + '_parent.restart'),
            ('parent_plumed', PROJECT + '.plumed', PROJECT + '.plumed'),
            ('parent_hills', 'HILLS_' + PROJECT, 'HILLS_' + PROJECT),
            ('parent_colvar', 'COLVAR_' + PROJECT, 'COLVAR_' + PROJECT),
            ('parent_run_record', 'run_record.json', 'parent_run_record.json')]:
        path = (HERE / 'inputs' if role.startswith('new') else source) / name
        record['stage'].append({'role': role, 'source': name, 'target': target, 'sha256': sha(path)})
    manifest = {'schema': 'li2025_extend/v2', 'cases': [record],
                'dependencies': {name: sha(HERE / name) for name in
                                 ('frozen_run_cp2k.py', 'prepare.py', 'run_continuation.py', 'run.slurm')},
                'scope': 'Ordinary biased metaD; unchanged science; no structural/FES verdict.',
                'stop_rule': 'STOP at cumulative 50 ps; scientific review required before any further continuation.',
                'restart_policy': 'Same as seg2/seg3/seg4: POS/VEL/COUNTERS only; CSVR parameters preserved, random state not restored; SCF_GUESS ATOMIC unchanged.'}
    dump(HERE / 'continuation_manifest.json', manifest)
    print(CASE + ': parent technical PASS; manifest frozen for 40001..50000')
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--history', type=Path, default=DEFAULT_HISTORY)
    parser.add_argument('--audit-out', type=Path, default=DEFAULT_AUDIT)
    parser.add_argument('--source-sha', type=Path, default=DEFAULT_SOURCE_SHA)
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    prepare(args.source, args.history, args.audit_out, args.validate_only, args.source_sha)


if __name__ == '__main__':
    main()
