#!/usr/bin/env python3
"""Download, unpack and verify selected Li 2025 or Liu 2026 data."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import tarfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def destination(root, name):
    posix = PurePosixPath(name)
    if posix.is_absolute() or '..' in posix.parts or '\\' in name:
        raise ValueError(f'Unsafe archive path: {name}')
    path = root.joinpath(*posix.parts)
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f'Path escapes repository: {name}')
    return path


def extract_checked(archive, records, root=ROOT):
    expected = {r['path']: r for r in records}
    seen = set()
    with tarfile.open(archive, mode='r|gz') as tar:
        for member in tar:
            if not member.isfile() or member.name not in expected or member.name in seen:
                raise ValueError(f'Unexpected archive member: {member.name}')
            rec = expected[member.name]
            if member.size != rec['bytes']:
                raise ValueError(f'Size mismatch: {member.name}')
            dest = destination(root, member.name)
            stream = tar.extractfile(member)
            if dest.exists():
                if not dest.is_file() or digest(dest) != rec['sha256']:
                    raise ValueError(f'Existing file differs; refusing overwrite: {member.name}')
                seen.add(member.name)
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            tmp = dest.with_name(dest.name + '.part')
            h = hashlib.sha256()
            if tmp.exists():
                raise ValueError(f'Existing partial file; refusing overwrite: {tmp}')
            try:
                with tmp.open('xb') as f:
                    for block in iter(lambda: stream.read(1024 * 1024), b''):
                        h.update(block)
                        f.write(block)
                if h.hexdigest() != rec['sha256']:
                    raise ValueError(f'Checksum mismatch: {member.name}')
                os.replace(tmp, dest)
            except Exception:
                tmp.unlink(missing_ok=True)
                raise
            seen.add(member.name)
    if seen != set(expected):
        raise ValueError(f'Missing archive members: {sorted(set(expected) - seen)}')
    return len(seen)


def verify(records, root=ROOT):
    bad = []
    for rec in records:
        p = destination(root, rec['path'])
        if not p.is_file() or p.stat().st_size != rec['bytes'] or digest(p) != rec['sha256']:
            bad.append(rec['path'])
    if bad:
        raise ValueError('Missing or changed data files: ' + ', '.join(bad[:12]))
    return len(records)


def fetch_one(asset, archive_argument=None, verify_only=False):
    manifest = destination(ROOT, asset.get('manifest', 'metadata/DATA_FILES.json'))
    records = json.loads(manifest.read_text())['files']
    label = asset.get('dataset', 'li2025')
    if verify_only:
        print(f"{label}: verified {verify(records)} files.")
        return
    archive = archive_argument.resolve() if archive_argument else ROOT / 'downloads' / asset['filename']
    if not archive.exists():
        if archive_argument:
            raise FileNotFoundError(archive)
        archive.parent.mkdir(parents=True, exist_ok=True)
        tmp = archive.with_suffix(archive.suffix + '.part')
        if tmp.exists():
            raise ValueError(f'Existing partial download; refusing overwrite: {tmp}')
        request = urllib.request.Request(asset['url'], headers={'User-Agent': 'interfacial-cs-reproduction/1'})
        try:
            with urllib.request.urlopen(request, timeout=120) as source, tmp.open('xb') as target:
                for block in iter(lambda: source.read(1024 * 1024), b''):
                    target.write(block)
            os.replace(tmp, archive)
        except Exception:
            tmp.unlink(missing_ok=True)
            raise
    if archive.stat().st_size != asset['bytes'] or digest(archive) != asset['sha256']:
        raise ValueError(f'{label}: archive size or SHA-256 mismatch')
    n = extract_checked(archive, records)
    verify(records)
    print(f'{label}: extracted and verified {n} files from {archive.name}.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', choices=('li2025', 'liu2026', 'all'), default='li2025',
                        help='Data collection to retrieve (default: li2025)')
    parser.add_argument('--archive', type=Path, help='Use a local archive for one selected dataset')
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    assets = json.loads((ROOT / 'metadata/RELEASE_ASSETS.json').read_text())['assets']
    selected = [a for a in assets if args.dataset == 'all' or a.get('dataset', 'li2025') == args.dataset]
    if not selected:
        raise ValueError(f'No archive registered for {args.dataset}')
    if args.archive and len(selected) != 1:
        parser.error('--archive requires one dataset, not all')
    for asset in selected:
        fetch_one(asset, args.archive, args.verify_only)


if __name__ == '__main__':
    main()
