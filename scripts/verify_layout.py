"""Independent f64 GEOS verification against XML-derived original contours.

Never uses jagua's simplified geometry or its reported density. Edge contact is
allowed. Both zero-tolerance measurements and declared numerical tolerances are
reported; passing with tolerance is not a mathematical feasibility proof.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
from shapely.geometry import Polygon, box
from shapely.affinity import affine_transform
from dataset import load_xml, sha256


def verify(xml, output, relative_linear_tolerance=1e-6, relative_area_tolerance=1e-8,
           allowed_width_expansion_ratio=0.0):
    original, manifest = load_xml(xml)
    raw = json.loads(Path(output).read_text(encoding='utf-8'))
    instance = raw.get('instance', raw)
    solution = raw['solution']
    errors = []
    actual_height = instance['strip_height']
    nominal_height = original['strip_height']
    width_ratio = actual_height / nominal_height - 1
    if instance['name'].casefold() != original['name'].casefold():
        errors.append('instance name changed')
    if not math.isfinite(actual_height) or width_ratio < -1e-12 or width_ratio > allowed_width_expansion_ratio + 1e-12:
        errors.append('fixed strip width changed beyond declared allowance')
    if instance['items'] != json.loads(json.dumps(original['items'])):
        # Legacy jagua outputs use allowed_orientations and min_quality=null.
        expected = {i['id']: i for i in original['items']}
        if len(instance['items']) != len(expected):
            errors.append('item type count changed')
        for item in instance['items']:
            e = expected.get(item['id'])
            if e is None or item['demand'] != e['demand'] or item['shape'] != e['shape']:
                errors.append('input contours / quantities changed')
            else:
                angles = item.get('allowed_orientations', item.get('orientation', {}).get('rotation', {}).get('angles'))
                # XML up-left parsing can reverse angles and retain -360 for 0.
                # Compare legal rotation sets modulo a full turn; list ordering
                # and duplicate full-turn representatives do not change legality.
                canon = lambda aa: sorted(set(round(float(a) % 360, 5) % 360 for a in aa))
                if angles is None or canon(angles) != canon(e['orientation']['rotation']['angles']):
                    errors.append('input allowed rotations changed')
    height = manifest['fixed_width']
    length = solution['strip_width']
    if not math.isfinite(length) or length <= 0:
        raise ValueError('invalid strip length')
    linear_tol = relative_linear_tolerance * max(height, 1.0)
    area_tol = relative_area_tolerance * max(manifest['total_area'], 1.0)
    container = box(0, 0, length, actual_height)
    nominal_container = box(0, 0, length, height)
    expected = {i['id']: i for i in original['items']}
    placed = solution['layout']['placed_items']
    counts = Counter(p['item_id'] for p in placed)
    if counts != Counter({i['id']: i['demand'] for i in original['items']}):
        errors.append('missing / duplicate / unknown items')
    polygons, max_outside_area, max_boundary_excess = [], 0.0, 0.0
    nominal_outside_area, nominal_boundary_excess = 0.0, 0.0
    for p in placed:
        if p['item_id'] not in expected:
            continue
        item = expected[p['item_id']]
        t = p['transformation']
        angle, (x, y) = t['rotation'], t['translation']
        if not all(math.isfinite(v) for v in [angle, x, y]):
            errors.append('nonfinite transformation')
            continue
        deg = math.degrees(angle) % 360
        if min(abs((deg-a+180)%360-180) for a in item['orientation']['rotation']['angles']) > 1e-4:
            errors.append('forbidden rotation')
        c, s = math.cos(angle), math.sin(angle)
        poly = affine_transform(Polygon(item['shape']['data']), [c, -s, s, c, x, y])
        if not poly.is_valid:
            errors.append('invalid original polygon')
        minx, miny, maxx, maxy = poly.bounds
        excess = max(0, -minx, -miny, maxx-length, maxy-actual_height)
        nominal_boundary_excess = max(nominal_boundary_excess, -minx, -miny, maxx-length, maxy-height)
        nominal_outside_area = max(nominal_outside_area, poly.difference(nominal_container).area)
        max_boundary_excess = max(max_boundary_excess, excess)
        max_outside_area = max(max_outside_area, poly.difference(container).area)
        polygons.append(poly)
    intersections, max_overlap_area = [], 0.0
    for i, a in enumerate(polygons):
        for j in range(i):
            b = polygons[j]
            if a.intersects(b):
                area = a.intersection(b).area
                max_overlap_area = max(max_overlap_area, area)
                if area > area_tol:
                    intersections.append({'i': i, 'j': j, 'overlap_area': area})
    if intersections:
        errors.append('positive-area overlap exceeds tolerance')
    if max_boundary_excess > linear_tol or max_outside_area > area_tol:
        errors.append('original contour protrudes beyond strip tolerance')
    return {'instance': manifest['instance'], 'xml_sha256': manifest['xml_sha256'],
        'output_sha256': sha256(output), 'passed_with_declared_tolerances': not errors,
        'zero_area_violation': max_overlap_area == 0 and max_outside_area == 0,
        'errors': sorted(set(errors)), 'piece_count': len(placed), 'length': length,
        'utilization_percent_recomputed': 100*manifest['total_area']/(height*length),
        'utilization_percent_actual_width': 100*manifest['total_area']/(actual_height*length),
        'nominal_fixed_width': height, 'actual_fixed_width': actual_height,
        'width_expansion_ratio': width_ratio,
        'allowed_width_expansion_ratio': allowed_width_expansion_ratio,
        'max_nominal_boundary_excess': nominal_boundary_excess,
        'max_nominal_outside_area': nominal_outside_area,
        'reported_solver_density': solution.get('density'),
        'linear_tolerance': linear_tol, 'pair_area_tolerance': area_tol,
        'max_boundary_excess': max_boundary_excess, 'max_outside_area': max_outside_area,
        'max_pair_overlap_area': max_overlap_area, 'overlaps': intersections,
        'boundary_discrepancies': manifest['boundary_discrepancies'],
        'method': 'original XML vertices, f64 GEOS/Shapely, proper rigid transforms, all pairs'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xml', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--report', required=True, type=Path)
    args = p.parse_args()
    result = verify(args.xml, args.output)
    args.report.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['instance','passed_with_declared_tolerances','errors','utilization_percent_recomputed']}))
    raise SystemExit(0 if result['passed_with_declared_tolerances'] else 1)
