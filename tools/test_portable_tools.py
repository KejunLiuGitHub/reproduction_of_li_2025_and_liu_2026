#!/usr/bin/env python3
"""Isolated integration and boundary tests; no CP2K or network calls."""
import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_tests(report):
    cases = json.loads((ROOT / 'data/DATASETS.json').read_text())['cases']
    selected = ['natural_original40', 'r1150_original', 'r2200_original', 'p_bulk_smoke', 'r_eq_metad_smoke']
    results = []
    with tempfile.TemporaryDirectory(prefix='li296-public-test-') as temp:
        root = Path(temp) / 'repository'
        shutil.copytree(ROOT, root, ignore=shutil.ignore_patterns('.git', '.venv', '__pycache__', 'data', 'evidence'))
        (root / 'data').mkdir()
        fixture = {'cases': {name: cases[name] for name in selected}}
        (root / 'data/DATASETS.json').write_text(json.dumps(fixture))
        copied = {}
        for case in fixture['cases'].values():
            for key in ('trajectory', 'energy', 'input', 'coordinates'):
                if not case.get(key):
                    continue
                path = Path(case[key]); target = root / path
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(ROOT / path, target)
                copied[str(path)] = sha(target)
            directory = Path(case['trajectory']).parent
            for name in ('DISTANCE_ANALYSIS.json', 'run_record.json'):
                if (ROOT / directory / name).is_file():
                    shutil.copyfile(ROOT / directory / name, root / directory / name)
        def command(argv, succeeds=True):
            result = subprocess.run([sys.executable, '-B', *argv], cwd=root, text=True, capture_output=True)
            if (result.returncode == 0) != succeeds:
                raise AssertionError(str(argv) + '\n' + result.stdout + '\n' + result.stderr)
            return result
        output = root / 'replay'
        for name in ('natural_original40', 'r1150_original', 'p_bulk_smoke', 'r_eq_metad_smoke'):
            dest = output / name
            command(['tools/reproduce_case.py', '--case', name, '--out', str(dest)])
            summary = json.loads((dest / 'SUMMARY.json').read_text())
            expected = cases[name]['coverage']['frames']
            assert summary['coverage']['XYZ']['frames'] == expected
            assert summary['biased'] == cases[name]['biased']
            assert summary['environment'] == ('bulk' if name == 'p_bulk_smoke' else 'interface')
            assert (summary['signed_GDS_A'] is None) == (name == 'p_bulk_smoke')
            assert all(sha(dest / path) == digest for path, digest in json.loads((dest / 'OUTPUT_SHA256.json').read_text()).items())
            results.append({'test': 'isolated_full_replay_' + name, 'pass': True, 'frames': expected,
                            'SC_A': summary['SC_A'], 'CO_A': summary['CO_A'], 'summary_sha256': sha(dest / 'SUMMARY.json')})
            archived = root / Path(cases[name]['trajectory']).parent / 'DISTANCE_ANALYSIS.json'
            if archived.exists():
                prior = json.loads(archived.read_text())
                new = json.loads((dest / 'DISTANCE_ANALYSIS_RAW.json').read_text())
                assert prior['per_frame'] == new['per_frame']
                results.append({'test': 'archived_distance_all_frames_' + name, 'pass': True, 'reference_sha256': sha(archived)})
        original = json.loads((output / 'r1150_original/SOLVATION.json').read_text())['rows']
        first = original[:258]
        assert first[-1]['step'] == 2570
        assert abs(sum(r['SC'] for r in first)/len(first)-3.6770440693864246) < 1e-12
        assert abs(sum(r['CO'] for r in first)/len(first)-1.2379804748654546) < 1e-12
        results.append({'test': 'original_prior_258_frame_scalar_reference', 'pass': True,
                        'reference': '20261119T210251+0800 analysis/SUMMARY.json SHA7491a7a170ba26841446a707834c23dbaf41246ebb75b4994132bf1d5522800b; constants carried into public test, no private file read'})
        command(['tools/reproduce_case.py','--case','r1150_original','--through-step','100','--out',str(root/'first100')])
        bounded=json.loads((root/'first100/SUMMARY.json').read_text())
        assert bounded['coverage']['XYZ']['frames']==11 and bounded['coverage']['XYZ']['last_step']==100
        assert bounded['coverage']['energy']['rows']==101 and bounded['coverage']['energy']['last_step']==100
        results.append({'test':'through_step100_has11frames101energies','pass':True})
        command(['tools/reproduce_case.py','--case','r1150_original','--out',str(output/'r1150_original')], False)
        results.append({'test':'refuse_existing_output','pass':True})
        def altered_manifest(name, mutate):
            data=copy.deepcopy(fixture);mutate(data['cases'][name]);path=root/(name+'_altered.json');path.write_text(json.dumps(data));return path
        path=altered_manifest('r1150_original',lambda c:c.update(trajectory='../outside.xyz'))
        command(['tools/reproduce_case.py','--case','r1150_original','--datasets',str(path),'--out',str(root/'bad_path')],False)
        results.append({'test':'reject_source_path_escape','pass':True})
        path=altered_manifest('r1150_original',lambda c:c.update(box_angstrom=[14.2296,14.2296,14.2296]))
        command(['tools/reproduce_case.py','--case','r1150_original','--datasets',str(path),'--out',str(root/'bad_box')],False)
        results.append({'test':'reject_manifest_cell_different_from_input','pass':True})
        path=altered_manifest('r1150_original',lambda c:c['sha256'].update(trajectory='0'*64))
        command(['tools/reproduce_case.py','--case','r1150_original','--datasets',str(path),'--out',str(root/'bad_hash')],False)
        results.append({'test':'reject_declared_SHA_mismatch','pass':True})
        data=copy.deepcopy(fixture);case=data['cases']['p_bulk_smoke']
        for key,suffix in [('trajectory',b'296\nunfinished frame\n'),('energy',b' 21 21.0')]:
            path=root/'data'/('partial_'+key)
            path.write_bytes((root/case[key]).read_bytes()+suffix)
            case[key]=str(path.relative_to(root));case['sha256'][key]=sha(path)
        data['cases']['p_bulk_smoke']['job_id']=None
        path=root/'partial.json';path.write_text(json.dumps(data))
        command(['tools/reproduce_case.py','--case','p_bulk_smoke','--datasets',str(path),'--out',str(root/'partial')])
        partial=json.loads((root/'partial/SUMMARY.json').read_text())
        assert partial['coverage']['XYZ']['frames']==3 and partial['coverage']['XYZ']['excluded_tail_bytes']>0
        assert partial['coverage']['energy']['excluded_tail_bytes']>0 and partial['job_id'] is None
        with gzip.open(root/'partial/network/NETWORK.json.gz','rt') as handle:assert json.load(handle)['source_job']=='archive:p_bulk_smoke'
        results.append({'test':'partial_final_records_and_unknown_job','pass':True})
        data_dir=root/'cp2k_data_for_rendering_test';data_dir.mkdir()
        for name in ('BASIS_MOLOPT','GTH_POTENTIALS','dftd3.dat'):
            (data_dir/name).write_text('Rendering-test placeholder only; never used for CP2K.\n')
        for name in ('r1150_original','r2200_original'):
            dest=root/'prepared'/name
            command(['tools/prepare_unbiased_run.py','--case',name,'--cp2k-data-dir',str(data_dir),'--out',str(dest)])
            prep=json.loads((dest/'PREPARATION.json').read_text())
            assert prep['steps']==5000 and prep['seed']==(2026111603 if name=='r1150_original' else 2026111604)
            assert prep['scientific_input_lines_unchanged'] and prep['run_status']=='PREPARED_NOT_RUN'
            assert not (dest/'CHECKS.json').exists() and not (dest/'run_record.json').exists()
            assert sha(dest/prep['coordinates'])==sha(root/cases[name]['coordinates'])
            assert all(sha(dest/n)==h for n,h in prep['output_sha256'].items())
            results.append({'test':'portable_input_only_runtime_changes_'+name,'pass':True})
        source_input=root/cases['r1150_original']['input'];saved=source_input.read_bytes()
        source_input.write_bytes(saved.replace(b'CUTOFF 300',b'CUTOFF 301'))
        command(['tools/prepare_unbiased_run.py','--case','r1150_original','--cp2k-data-dir',str(data_dir),'--out',str(root/'wrong_input')],False)
        source_input.write_bytes(saved)
        results.append({'test':'reject_changed_scientific_input','pass':True})
        assert all(sha(root/path)==digest for path,digest in copied.items())
        results.append({'test':'all_copied_sources_unchanged','pass':True,'files':len(copied)})
        report_data={'status':'PASS','tests':results,'count':len(results),'isolated_copy':True,'outside_repository_scientific_data_reads':False,
                     'CP2K_executed':False,'test_data_files_note':'CP2K parameter placeholders are used solely to test input rendering; no dynamics or engine validation is performed.',
                     'source_dataset_manifest_sha256':sha(ROOT/'data/DATASETS.json'),'python':sys.version}
        if report:
            Path(report).parent.mkdir(parents=True,exist_ok=True)
            Path(report).write_text(json.dumps(report_data,indent=2)+'\n')
        print(json.dumps({'status':'PASS','count':len(results),'replayed_cases':4,'report':str(report) if report else None}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path)
    run_tests(parser.parse_args().report)
