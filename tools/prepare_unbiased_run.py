#!/usr/bin/env python3
"""Prepare a fresh portable input for either unrotated original-coordinate control."""
import argparse
import difflib
import json
from pathlib import Path
import re
import shutil
import sys

sys.dont_write_bytecode = True
from _common import ROOT, fresh_directory, load_case, sha256, source_path, write_json

PACKAGE = ROOT / 'campaign/li2025_hbcontrol_20261118/production_v1'
ALLOWED_CASES = {'r1150_original', 'r2200_original'}
DATA_FILES = {'BASIS_SET_FILE_NAME': 'BASIS_MOLOPT', 'POTENTIAL_FILE_NAME': 'GTH_POTENTIALS', 'PARAMETER_FILE_NAME': 'dftd3.dat'}


def portable_text(original, project, data_dir):
    """Change only four approved runtime fields, preserving every other line."""
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,62}', project):
        raise ValueError('PROJECT must be a short alphanumeric/underscore name')
    if any(char.isspace() for char in str(data_dir)):
        raise ValueError('Use a CP2K data directory whose path contains no whitespace')
    replacements = {**{key: str(data_dir / name) for key, name in DATA_FILES.items()}, 'PROJECT': project}
    counts = {key: 0 for key in replacements}
    lines = []
    for line in original.splitlines(keepends=True):
        match = re.match(r'^(\s*)([A-Z_]+)\s+([^!#\r\n]+?)(\s*(?:[!#].*)?)(\r?\n)?$', line)
        if match and match[2] in replacements:
            key = match[2]
            counts[key] += 1
            lines.append(match[1] + key + ' ' + replacements[key] + match[4] + (match[5] or ''))
        else:
            lines.append(line)
    if any(value != 1 for value in counts.values()):
        raise ValueError('Expected exactly one occurrence of each portable input field: ' + str(counts))
    rendered = ''.join(lines)
    # A line-wise allowlist catches accidental changes to seeds, steps, science,
    # coordinates, velocities or bias options, including future renderer edits.
    old, new = original.splitlines(), rendered.splitlines()
    if len(old) != len(new):
        raise ValueError('Unexpected input line-count change')
    for a, b in zip(old, new):
        if a != b and a.strip().split()[0] not in replacements:
            raise ValueError('Nonportable scientific input field changed')
    return rendered


def prepare(case_id, out, data_dir, project=None, datasets=None):
    if case_id not in ALLOWED_CASES:
        raise ValueError('Supported cases: ' + ', '.join(sorted(ALLOWED_CASES)))
    case, datasets_path = load_case(case_id, datasets)
    if case['biased'] or case['environment'] != 'interface':
        raise ValueError('Original controls must use the documented unbiased interface model')
    original_input = source_path(case['input'])
    coordinates = source_path(case['coordinates'])
    frozen_input = PACKAGE / 'inputs' / (case_id + '.inp')
    frozen_coordinates = PACKAGE / 'inputs' / (case_id + '.xyz')
    manifest_path = PACKAGE / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    record = next(row for row in manifest['cases'] if row['name'] == case_id)
    expected = {row['target']: row['sha256'] for row in record['stage']}
    if sha256(original_input) != sha256(frozen_input) or sha256(original_input) != expected[case_id + '.inp']:
        raise ValueError('Selected original input differs from the frozen production input')
    if sha256(coordinates) != sha256(frozen_coordinates) or sha256(coordinates) != expected[case_id + '.xyz']:
        raise ValueError('Selected coordinates differ from the frozen original frame')
    original = original_input.read_text()
    forbidden = r'(?mi)^\s*(?:&(?:FREE_ENERGY|METADYN|EXT_RESTART|CONSTRAINT|COLLECTIVE|VELOCITY)\b|PLUMED_INPUT_FILE\b|USE_PLUMED\b)'
    if re.search(forbidden, original):
        raise ValueError('Bias, constraint, restart or inherited velocity section found')
    data_dir = Path(data_dir).resolve()
    data_files = [data_dir / name for name in DATA_FILES.values()]
    if not all(path.is_file() for path in data_files):
        raise ValueError('CP2K data directory must contain BASIS_MOLOPT, GTH_POTENTIALS and dftd3.dat')
    coordinate_names = re.findall(r'(?mi)^\s*COORD_FILE_NAME\s+(\S+)', original)
    if coordinate_names != [case_id + '.xyz']:
        raise ValueError('Unexpected coordinate filename in frozen original input')
    sources = [datasets_path, original_input, coordinates, frozen_input, frozen_coordinates, manifest_path, *data_files,
               Path(__file__).resolve(), ROOT / 'tools/_common.py']
    before = {str(path): sha256(path) for path in sources}
    project = project or case_id + '_portable'
    rendered = portable_text(original, project, data_dir)
    out = fresh_directory(out)
    input_target = out / (case_id + '.inp')
    input_target.write_text(rendered)
    shutil.copyfile(coordinates, out / coordinate_names[0])
    diff = ''.join(difflib.unified_diff(original.splitlines(keepends=True), rendered.splitlines(keepends=True),
                                        fromfile='original/' + original_input.name, tofile='portable/' + input_target.name))
    (out / 'INPUT_PATH_CHANGES.diff').write_text(diff)
    after = {path: sha256(Path(path)) for path in before}
    if before != after:
        raise ValueError('A preparation source changed')
    result = {'schema': 'li296_portable_input/v1', 'case': case_id, 'project': project,
              'input': input_target.name, 'coordinates': coordinate_names[0], 'steps': record['steps'], 'seed': record['seed'],
              'box_angstrom': case['box_angstrom'], 'biased': False, 'history': case['history'],
              'source_segment_biased': record['source_segment_biased'], 'prepared_under_bias': record['prepared_under_bias'],
              'controlled_initial_conditions': record['controlled_initial_conditions'],
              'allowed_changed_fields': [*DATA_FILES, 'PROJECT'], 'scientific_input_lines_unchanged': True,
              'source_sha256': {path: {'before': digest, 'after': after[path]} for path, digest in before.items()},
              'output_sha256': {input_target.name: sha256(input_target), coordinate_names[0]: sha256(out / coordinate_names[0]),
                                'INPUT_PATH_CHANGES.diff': sha256(out / 'INPUT_PATH_CHANGES.diff')},
              'runtime': {'recorded_engine': 'CP2K 9.1', 'mpi_tasks': 64, 'omp_threads': 1, 'data_directory': str(data_dir)},
              'run_status': 'PREPARED_NOT_RUN', 'original_case_status': case['status'],
              'note': 'Fresh input copy with runtime paths and output PROJECT changed. No historical run certificate is copied.'}
    write_json(out / 'PREPARATION.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True, choices=sorted(ALLOWED_CASES))
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--cp2k-data-dir', type=Path, required=True)
    parser.add_argument('--project')
    parser.add_argument('--datasets', type=Path)
    args = parser.parse_args()
    result = prepare(args.case, args.out, args.cp2k_data_dir, args.project, args.datasets)
    print(json.dumps({'case': result['case'], 'steps': result['steps'], 'seed': result['seed'], 'output': str(args.out.resolve()), 'status': result['run_status']}))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError, StopIteration) as error:
        print(str(error), file=sys.stderr)
        sys.exit(2)
