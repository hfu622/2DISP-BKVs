"""Fetch only seven pinned author files; reuse exact cached bytes, never overwrite.

Source licenses are not inferred from neighboring projects. The unlicensed
author wrapper stays in ignored external/ instead of being redistributed.
"""
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

ROOT=Path(__file__).resolve().parents[1]


def prepare():
    lock=json.loads((ROOT/'data/f057_source_lock.json').read_text(encoding='utf-8'))
    for row in lock['files']:
        target=ROOT/row['destination']
        if target.exists():
            if hashlib.sha256(target.read_bytes()).hexdigest()!=row['sha256']:
                raise ValueError(f'Existing source differs; preserved without overwrite: {target}')
            print('REUSE',row['destination'])
            continue
        repo=row['repository'].removeprefix('https://github.com/')
        url=f"https://raw.githubusercontent.com/{repo}/{row['commit']}/{row['path']}"
        with urlopen(Request(url,headers={'User-Agent':'2DISP-BKVs-source-reproduction'}),timeout=45) as response:
            data=response.read(4*1024*1024)
        if hashlib.sha256(data).hexdigest()!=row['sha256']:
            raise ValueError(f'Pinned source digest mismatch: {url}')
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(data)
        print('FETCH',row['destination'])


if __name__=='__main__':prepare()
