#!/usr/bin/env python3
"""Download, unpack and verify the selected Li 2025 data archive."""
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, help='Use a previously downloaded archive')
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    records = json.loads((ROOT / 'metadata/DATA_FILES.json').read_text())['files']
    if args.verify_only:
        print(f'Verified {verify(records)} files.')
        return
    assets = json.loads((ROOT / 'metadata/RELEASE_ASSETS.json').read_text())['assets']
    if len(assets) != 1:
        raise ValueError('This snapshot expects one data archive')
    asset = assets[0]
    archive = args.archive.resolve() if args.archive else ROOT / 'downloads' / asset['filename']
    if not archive.exists():
        if args.archive:
            raise FileNotFoundError(archive)
        archive.parent.mkdir(parents=True, exist_ok=True)
        tmp = archive.with_suffix(archive.suffix + '.part')
        request = urllib.request.Request(asset['url'], headers={'User-Agent': 'li2025-reproduction/1'})
        try:
            with urllib.request.urlopen(request, timeout=120) as source, tmp.open('xb') as target:
                for block in iter(lambda: source.read(1024 * 1024), b''):
                    target.write(block)
            os.replace(tmp, archive)
        except Exception:
            tmp.unlink(missing_ok=True)
            raise
    if archive.stat().st_size != asset['bytes'] or digest(archive) != asset['sha256']:
        raise ValueError('Archive size or SHA-256 mismatch')
    n = extract_checked(archive, records)
    verify(records)
    print(f'Extracted and verified {n} files from {archive.name}.')


if __name__ == '__main__':
    main()
