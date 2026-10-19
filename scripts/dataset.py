"""Read original ESICUP XML; never infer rotations from an instance name.

The source coordinate frame is reflected to Cartesian coordinates when the XML
declares up-left. Rotation signs are reflected together with coordinates.
Only one solid polygon per piece and a rectangular board are supported. Reject
other geometry rather than quietly dropping components or inflating the board.
"""
import argparse
from collections import Counter
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_xml(path):
    path = Path(path)
    root = ET.parse(path).getroot()
    for element in root.iter():
        element.tag = element.tag.split('}')[-1]
    origin = root.findtext('coordinatesOrigin')
    if origin not in ('up-left', 'down-left'):
        raise ValueError(f'Unsupported coordinate origin: {origin}')
    sign = -1 if origin == 'up-left' else 1
    polygons = {}
    boundary_discrepancies = []
    for poly in root.findall('./polygons/polygon'):
        segments = sorted(poly.findall('./lines/segment'), key=lambda s: int(s.attrib['n']))
        vertices = [(Decimal(s.attrib['x0']), Decimal(s.attrib['y0'])) for s in segments]
        if len(vertices) != int(poly.attrib['nVertices']):
            raise ValueError('Vertex count mismatch')
        for i, s in enumerate(segments):
            end = (Decimal(s.attrib['x1']), Decimal(s.attrib['y1']))
            nxt = vertices[(i + 1) % len(vertices)]
            gap = max(abs(end[0]-nxt[0]), abs(end[1]-nxt[1]))
            if gap:
                # The original Jakobs2 polygon8 has three inconsistent endpoints.
                # Preserve all declared start vertices; never silently snap/repair.
                if path.stem.lower() != 'jakobs2' or poly.attrib['id'] != 'polygon8' or gap > Decimal('0.001'):
                    raise ValueError(f'Disconnected boundary: {poly.attrib["id"]}, gap={gap}')
                boundary_discrepancies.append({'polygon': poly.attrib['id'], 'segment': s.attrib['n'],
                    'endpoint': list(map(str, end)), 'next_start': list(map(str, nxt)), 'gap': str(gap)})
        polygons[poly.attrib['id']] = vertices

    def component(piece):
        comps = piece.findall('./component')
        if len(comps) != 1 or comps[0].attrib['type'] != '0':
            raise ValueError('Unsupported multiple components / holes')
        c = comps[0]
        x, y = Decimal(c.attrib['xOffset']), Decimal(c.attrib['yOffset'])
        return [(vx + x, (vy + y) * sign) for vx, vy in polygons[c.attrib['idPolygon']]]

    boards = root.findall('./problem/boards/piece')
    if len(boards) != 1 or int(boards[0].attrib['quantity']) != 1:
        raise ValueError('Expected one board')
    board = component(boards[0])
    xs, ys = {p[0] for p in board}, {p[1] for p in board}
    if len(board) != 4 or len(xs) != 2 or len(ys) != 2 or set(board) != {(x, y) for x in xs for y in ys}:
        raise ValueError('Nonrectangular board')
    height = max(ys) - min(ys)
    items, identities, total_area = [], [], Decimal(0)
    for index, piece in enumerate(root.findall('./problem/lot/piece')):
        angles = piece.findall('./orientation/enumeration')
        if not angles or len(piece.find('orientation')) != len(angles):
            raise ValueError('Unsupported continuous / missing orientation')
        raw_angles = [Decimal(e.attrib['angle']) for e in angles]
        # Decimal remainder preserves sign, hence the explicit addition of 360.
        solver_angles = sorted({float((a * sign % 360 + 360) % 360) for a in raw_angles})
        vertices = component(piece)
        area = abs(sum(vertices[i][0] * vertices[(i+1) % len(vertices)][1] -
                       vertices[(i+1) % len(vertices)][0] * vertices[i][1]
                       for i in range(len(vertices)))) / 2
        quantity = int(piece.attrib['quantity'])
        if area <= 0 or quantity <= 0:
            raise ValueError('Degenerate item / invalid demand')
        total_area += area * quantity
        items.append({'id': index, 'demand': quantity,
                      'orientation': {'rotation': {'mode': 'discrete', 'angles': solver_angles}},
                      'shape': {'type': 'simple_polygon', 'data': [[float(x), float(y)] for x, y in vertices]}})
        identities.append({'id': index, 'piece_id': piece.attrib['id'], 'demand': quantity,
                           'xml_angles': [float(a) for a in raw_angles], 'solver_angles': solver_angles,
                           'area_exact': str(area), 'vertices': len(vertices)})
    name = path.stem.lower()
    instance = {'name': name, 'strip_height': float(height), 'min_item_separation': 0.0, 'items': items}
    manifest = {'instance': name, 'xml_sha256': sha256(path), 'coordinate_origin': origin,
                'fixed_width': float(height), 'board_length': float(max(xs)-min(xs)),
                'total_area_exact': str(total_area), 'total_area': float(total_area),
                'piece_count': sum(i['demand'] for i in items), 'type_count': len(items),
                'pieces': identities, 'reflection_allowed': False, 'width_inflation': 0.0,
                'contour_policy': 'declared x0/y0 vertices in segment order',
                'boundary_discrepancies': boundary_discrepancies}
    return instance, manifest


def prepare(nesting, out, manifest_path):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    files = sorted(Path(nesting).glob('*.xml'))
    if not files:
        raise ValueError('No XML instances')
    manifests = []
    for path in files:
        instance, manifest = load_xml(path)
        (out / f'{instance["name"]}.json').write_text(json.dumps(instance, indent=2)+'\n', encoding='utf-8')
        legacy = dict(instance)
        legacy['items'] = [{k:v for k,v in i.items() if k != 'orientation'} | {
            'allowed_orientations': i['orientation']['rotation']['angles']} for i in instance['items']]
        (out / f'{instance["name"]}.legacy.json').write_text(json.dumps(legacy, indent=2)+'\n', encoding='utf-8')
        manifests.append(manifest)
    target = Path(manifest_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(manifests, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'instances': len(manifests), 'piece_counts': {m['instance']: m['piece_count'] for m in manifests}}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nesting', type=Path, required=True)
    parser.add_argument('--out', type=Path, default=Path('inputs'))
    parser.add_argument('--manifest', type=Path, default=Path('data/instances.json'))
    args = parser.parse_args()
    prepare(args.nesting, args.out, args.manifest)
