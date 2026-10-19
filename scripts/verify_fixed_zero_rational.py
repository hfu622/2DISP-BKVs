"""Exact dyadic-rational geometry audit for simple polygons at exactly0 radians.

No buffering, snapping, polygon repair, contour scaling or coordinate rounding.
XML vertices and the encoded binary64 translations are interpreted exactly.
Used independently when GEOS raises on otherwise valid, near-contact polygons.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from dataset import load_xml


def cross(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def signed2(poly):
    return sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(poly,poly[1:]+poly[:1]))


def area(poly):
    return abs(signed2(poly))/2 if len(poly)>2 else Q(0)


def triangles(points):
    poly=list(points)
    if poly[0]==poly[-1]:poly.pop()
    for i,(a,b) in enumerate(zip(poly,poly[1:]+poly[:1])):
        assert a!=b,'Zero-length polygon edge'
        for j,(c,d) in enumerate(zip(poly,poly[1:]+poly[:1])):
            if j<=i or j==(i+1)%len(poly) or i==(j+1)%len(poly):continue
            boxes_meet=all(max(min(a[k],b[k]),min(c[k],d[k]))<=min(max(a[k],b[k]),max(c[k],d[k])) for k in (0,1))
            assert not (boxes_meet and cross(a,b,c)*cross(a,b,d)<=0 and cross(c,d,a)*cross(c,d,b)<=0),'Non-simple original polygon'
    if signed2(poly)<0:poly.reverse()
    while len(poly)>3:
        removed=False
        for i,b in enumerate(poly):
            a,c=poly[i-1],poly[(i+1)%len(poly)]
            if cross(a,b,c)==0 and all(min(a[k],c[k])<=b[k]<=max(a[k],c[k]) for k in (0,1)):
                poly.pop(i);removed=True;break
        if not removed:break
    result=[]
    while len(poly)>3:
        for i,b in enumerate(poly):
            a,c=poly[i-1],poly[(i+1)%len(poly)]
            if cross(a,b,c)<=0:continue
            rest=[p for j,p in enumerate(poly) if j not in ((i-1)%len(poly),i,(i+1)%len(poly))]
            if any(cross(a,b,p)>=0 and cross(b,c,p)>=0 and cross(c,a,p)>=0 for p in rest):continue
            result.append([a,b,c]);poly.pop(i);break
        else:raise ValueError('Exact ear clipping failed; do not repair the polygon')
    if len(poly)==3 and area(poly)>0:result.append(poly)
    assert sum(map(area,result))==area(points)
    return result


def clipped(subject,clipper):
    out=list(subject)
    for a,b in zip(clipper,clipper[1:]+clipper[:1]):
        source=out;out=[]
        if not source:break
        for p,q in zip(source,source[1:]+source[:1]):
            fp,fq=cross(a,b,p),cross(a,b,q)
            if fp>=0:out.append(p)
            if (fp<0 and fq>0) or (fp>0 and fq<0):
                t=fp/(fp-fq)
                out.append((p[0]+t*(q[0]-p[0]),p[1]+t*(q[1]-p[1])))
    return out


def bounds(poly):
    return min(p[0] for p in poly),min(p[1] for p in poly),max(p[0] for p in poly),max(p[1] for p in poly)


def overlap_boxes(a,b):
    return a[0]<b[2] and b[0]<a[2] and a[1]<b[3] and b[1]<a[3]


def rational(value):
    return dict(numerator=str(value.numerator),denominator=str(value.denominator),value=float(value))


def verify_exact(xml,output):
    original,manifest=load_xml(xml)
    raw=json.loads(output.read_text(encoding='utf-8'));instance=raw.get('instance',raw)
    assert instance['name'].casefold()==original['name'].casefold()
    height=Q(original['strip_height']);actual_height=Q(instance['strip_height'])
    assert actual_height==height,'This certificate requires the original fixed width'
    assert instance.get('min_item_separation',0)==0
    length=Q(raw['solution']['strip_width']);assert length>0
    expected={i['id']:i for i in original['items']}
    assert len(instance['items'])==len(expected)
    assert {i['id'] for i in instance['items']}==set(expected)
    for item in instance['items']:
        source=expected[item['id']]
        assert item['shape']==source['shape'] and item['demand']==source['demand']
        angles=item.get('allowed_orientations',item.get('orientation',{}).get('rotation',{}).get('angles'))
        assert angles is not None
        assert {Q(a)%360 for a in angles}=={Q(a)%360 for a in source['orientation']['rotation']['angles']}
    placed=raw['solution']['layout']['placed_items']
    assert Counter(p['item_id'] for p in placed)==Counter({i['id']:i['demand'] for i in original['items']})
    local={i['id']:[(Q(x),Q(y)) for x,y in i['shape']['data']] for i in original['items']}
    assert all(Q(0) in {Q(a)%360 for a in i['orientation']['rotation']['angles']} for i in original['items'])
    triangulations={k:triangles(v) for k,v in local.items()}
    container=[(Q(0),Q(0)),(length,Q(0)),(length,height),(Q(0),height)]
    polys=[];tris=[];boxes=[];outside=[];excess=[]
    for pose in placed:
        transform=pose['transformation'];assert transform['rotation']==0.,'Only exact0-radian rotations are supported'
        dx,dy=map(Q,transform['translation']);key=pose['item_id']
        shifted=lambda points:[(x+dx,y+dy) for x,y in points]
        poly=shifted(local[key]);ts=[shifted(t) for t in triangulations[key]]
        box=bounds(poly);polys.append(poly);tris.append(ts);boxes.append(box)
        outside.append(sum(area(t)-area(clipped(t,container)) for t in ts))
        excess.append(max(Q(0),-box[0],-box[1],box[2]-length,box[3]-height))
    overlaps=[];max_overlap=Q(0)
    for i,a in enumerate(tris):
        for j,b in enumerate(tris[:i]):
            if not overlap_boxes(boxes[i],boxes[j]):continue
            value=sum(area(clipped(x,y)) for x in a for y in b if overlap_boxes(bounds(x),bounds(y)))
            max_overlap=max(max_overlap,value)
            if value:overlaps.append(dict(i=i,j=j,area=rational(value)))
    total=Q(manifest['total_area_exact']);linear_tol=max(height,1)/10**6;area_tol=max(total,1)/10**8
    max_outside=max(outside);max_excess=max(excess)
    errors=[]
    if max_overlap>area_tol:errors.append('pair overlap exceeds declared tolerance')
    if max_outside>area_tol or max_excess>linear_tol:errors.append('boundary violation exceeds declared tolerance')
    return dict(instance=original['name'],method='Exact fractions; XML-derived vertices and encoded dyadic translations; ear triangulation and exact convex clipping; exactly0 radians only',
        xml_sha256=hashlib.sha256(xml.read_bytes()).hexdigest(),output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        passed_with_declared_tolerances=not errors,errors=errors,piece_count=len(placed),length=float(length),
        utilization_percent_recomputed=float(100*total/(height*length)),width_expansion_ratio=0.,
        max_pair_overlap_area=rational(max_overlap),max_outside_area=rational(max_outside),
        max_boundary_excess=rational(max_excess),pair_area_tolerance=rational(area_tol),linear_tolerance=rational(linear_tol),
        zero_area_violation=max_overlap==0 and max_outside==0,positive_pairs=overlaps,
        coordinates_modified=False,polygon_repair_applied=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--xml',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True);args=parser.parse_args()
    result=verify_exact(args.xml,args.output)
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['passed_with_declared_tolerances'] else 1)
