"""Run F057's published fixed-order vector baseline, retaining author poses.

No browser GA, extra rotation search, repair, or post-compaction is added.
The author wrapper and SVGnest must be fetched at the recorded commits first.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import time

from dataset import load_xml, sha256
from verify_layout import verify

ROOT = Path(__file__).resolve().parents[1]
AUTHOR = ROOT / 'external/Raster-Convolution-Search-and-CTM-Bench'
AUTHOR_COMMIT = 'c7aaa7e8c31db15f4f7df7673b2cecfc362b21cd'
SVG_COMMIT = '1248dc21efd3f90d1aa52ba5785e27e5217ed2c9'
ALGORITHM = 'SVGNest-F057-fixed-order'


def run(name, tag):
    assert name in ['han', 'poly3b', 'poly4b'], 'Audited fixed-zero configurations only'
    folder = ROOT / 'official_esicup/results' / tag / name
    if folder.exists():
        raise FileExistsError(folder)
    sources = []
    for repo, expected, url, files in [
        (AUTHOR, AUTHOR_COMMIT, 'https://github.com/Hermione-HKX/Raster-Convolution-Search-and-CTM-Bench',
         ['code/svgnest_albano_runner.js', 'results/svgnest_public_xml.csv']),
        (AUTHOR / 'external/SVGnest', SVG_COMMIT, 'https://github.com/Jack000/SVGnest',
         ['util/clipper.js', 'util/geometryutil.js', 'util/placementworker.js'])]:
        if (repo / '.git').exists():
            assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip() == expected
        for file in files:
            lock=json.loads((ROOT/'data/f057_source_lock.json').read_text(encoding='utf-8'))
            locked=next(r for r in lock['files'] if r['repository']==url and r['commit']==expected and r['path']==file)
            assert locked['sha256']==sha256(repo / file),f'Run prepare_f057_sources.py or check source line endings: {file}'
            sources.append(dict(url=url, commit=expected, path=file, sha256=sha256(repo / file)))
    node = shutil.which('node')
    if not node:
        raise FileNotFoundError('Node.js is required')
    xml = ROOT / 'official_esicup/xml' / f'{name}.xml'
    instance, manifest = load_xml(xml)
    assert manifest['coordinate_origin'] == 'up-left'
    pieces, copies = [], []
    for item in instance['items']:
        assert {a % 360 for a in item['orientation']['rotation']['angles']} == {0.0}
        vertices = item['shape']['data']
        minx, maxy = min(p[0] for p in vertices), max(p[1] for p in vertices)
        raw_normalized = [dict(x=x-minx, y=maxy-y) for x, y in vertices]
        width = max(p['x'] for p in raw_normalized)
        height = max(p['y'] for p in raw_normalized)
        area = float(manifest['pieces'][item['id']]['area_exact'])
        for copy in range(item['demand']):
            copies.append(dict(item_id=item['id'], minx=minx, maxy=maxy))
            pieces.append(dict(id=len(copies)-1, name=f"{item['id']}_{copy}", vertices=raw_normalized,
                               rotation=0, allowed_rotations=[0], area=area, width=width, height=height))
    height, length = manifest['fixed_width'], manifest['board_length']
    payload = dict(bin=[dict(x=0,y=0),dict(x=length,y=0),dict(x=length,y=height),dict(x=0,y=height)],
                   parts=sorted(pieces, key=lambda p:p['width'], reverse=True))
    folder.mkdir(parents=True)
    (folder / 'payload.json').write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8')
    command = [node, str(AUTHOR / 'code/svgnest_albano_runner.js')]
    record = dict(instance=name, algorithm=ALGORITHM, status='running', run_count=1, seed=None,
                  source_snapshot_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  source_files={f'scripts/{p}':sha256(ROOT/'scripts'/p) for p in
                                ['run_svgnest_reference.py','dataset.py','verify_layout.py']},
                  author_sources=sources, node_version=subprocess.check_output([node,'--version'],text=True).strip(),
                  xml_path=xml.relative_to(ROOT).as_posix(), xml_sha256=sha256(xml),
                  input_path=(folder/'payload.json').relative_to(ROOT).as_posix(), input_sha256=sha256(folder/'payload.json'),
                  command=command, dataset_id='esicup-official-xml', budget_seconds=600,
                  budget_kind='project_watchdog_for_single_deterministic_constructor', paper_budget_seconds=None,
                  paper_protocol='F057 TableIV source CSV vector baseline; width descending, stable XML ties, first allowed angle, one placement-worker evaluation',
                  paper_statistical_run_count=1, placeholder=False, manual_role='official_angle_policy_algorithm_test',
                  strict_full_paper_protocol_reproduction=False,
                  started_utc=datetime.now(timezone.utc).isoformat())
    target=folder/'run.json'
    def save(): target.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    save()
    start=time.perf_counter()
    result=subprocess.run(command,input=json.dumps(payload),text=True,capture_output=True,cwd=AUTHOR,timeout=600)
    record.update(exit_code=result.returncode,wall_seconds=time.perf_counter()-start,finished_utc=datetime.now(timezone.utc).isoformat())
    (folder/'stdout.log').write_text(result.stdout,encoding='utf-8')
    (folder/'stderr.txt').write_text(result.stderr,encoding='utf-8')
    raw=json.loads(result.stdout)
    (folder/'author_output.json').write_text(json.dumps(raw,indent=2)+'\n',encoding='utf-8')
    assert raw['nfp_errors']==0 and raw['bins_used']==1 and raw['items_placed_first_bin']==len(copies),raw
    placements=[]
    for p in raw['placements']:
        copy=copies[p['id']]
        assert p['rotation']==0
        placements.append(dict(item_id=copy['item_id'], transformation=dict(rotation=0,
                              translation=[p['x']-copy['minx'],height-copy['maxy']-p['y']])))
    output=dict(algorithm=ALGORITHM,instance=instance,seed=None,solve_seconds=record['wall_seconds'],
                solution=dict(strip_width=raw['occupied_length'],density=raw['strip_window_utilization_percent']/100,
                              layout=dict(placed_items=placements)),
                stats=dict(paper_run_status='completed_single_author_reference_pass',completed_variants=1,
                           evaluations=1,nfp_pairs=raw['nfp_pairs']),
                assumptions=['Unmodified F057 author JavaScript wrapper and pinned SVGnest placement worker',
                             'Common y-reflection converts the complete raw layout to canonical XML coordinates; no individual reflection or repair',
                             'Source supplementary CSV supplies per-instance results supporting the aggregate paper TableIV',
                             'Single deterministic baseline; this is not SVGnest browser GA or a new paper-algorithm reproduction'],
                strict_full_paper_protocol_reproduction=False)
    out=folder/'output'/f'final_{name}.json';out.parent.mkdir()
    out.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    check=verify(xml,out)
    (folder/'verification.json').write_text(json.dumps(check,indent=2)+'\n',encoding='utf-8')
    record.update(status='verified_with_tolerance' if check['passed_with_declared_tolerances'] else 'invalid_layout',
                  solver_status='completed_single_author_reference_pass',algorithm_completed=True,
                  utilization_percent=check['utilization_percent_recomputed'],length=check['length'],
                  assumptions=output['assumptions'])
    save()
    print(json.dumps(dict(instance=name,status=record['status'],length=record['length'],percent=record['utilization_percent'],
                          max_overlap=check['max_pair_overlap_area'],max_outside=check['max_outside_area'])))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--instance',choices=['han','poly3b','poly4b'],required=True)
    p.add_argument('--tag',required=True)
    a=p.parse_args();run(a.instance,a.tag)
