"""Path, provenance and JSON helpers for the public replay commands."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def source_path(value, root=ROOT):
    path = Path(value)
    if path.is_absolute():
        raise ValueError('Dataset paths must be relative to the repository')
    result = (root / path).resolve()
    if not result.is_relative_to(root.resolve()) or not result.is_file():
        raise ValueError('Missing source or path outside repository: ' + str(value))
    return result


def load_case(case_id, datasets=None):
    manifest = Path(datasets) if datasets else ROOT / 'data/DATASETS.json'
    data = json.loads(manifest.read_text())
    cases = data['cases']
    if not isinstance(cases, dict) or case_id not in cases:
        raise ValueError('Unknown case: ' + case_id)
    case = dict(cases[case_id])
    if case.get('id', case_id) != case_id:
        raise ValueError('Case id does not match dictionary key')
    for key in ('label', 'job_id', 'trajectory', 'energy', 'box_angstrom', 'biased', 'history', 'status'):
        if key not in case:
            raise ValueError('Missing case field: ' + key)
    if not isinstance(case['biased'], bool):
        raise ValueError('biased must be a JSON boolean')
    box = case['box_angstrom']
    if not isinstance(box, list) or len(box) != 3 or not all(isinstance(v, (int, float)) and math.isfinite(v) and v > 0 for v in box):
        raise ValueError('box_angstrom must contain three positive finite lengths')
    if all(abs(a-b) < 1e-8 for a, b in zip(box, (15., 15., 30.))):
        case['environment'] = 'interface'
    elif all(abs(a-14.2296) < 1e-8 for a in box):
        case['environment'] = 'bulk'
    else:
        raise ValueError('The registered solvation CLI supports only the two documented Li cells')
    return case, manifest.resolve()


def fresh_directory(path):
    path = Path(path).resolve()
    if path.exists():
        raise ValueError('Output already exists; choose a new directory: ' + str(path))
    path.mkdir(parents=True, exist_ok=False)
    return path


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
