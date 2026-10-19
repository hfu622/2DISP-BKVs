//! f64 simple-polygon decomposition and segment arrangements.
use anyhow::{Result, ensure};
use geo::{Area, Coord, LineString, Polygon, TriangulateEarcut};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct P {
    pub x: f64,
    pub y: f64,
}
impl P {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    pub fn shift(self, d: P) -> Self {
        Self::new(self.x + d.x, self.y + d.y)
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Segment {
    pub a: P,
    pub b: P,
}
pub fn cross(a: P, b: P, c: P) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}
pub fn signed_area(v: &[P]) -> f64 {
    (0..v.len())
        .map(|i| v[i].x * v[(i + 1) % v.len()].y - v[(i + 1) % v.len()].x * v[i].y)
        .sum::<f64>()
        * 0.5
}
pub fn polygon(v: &[P]) -> Polygon<f64> {
    Polygon::new(
        LineString::new(v.iter().map(|p| Coord { x: p.x, y: p.y }).collect()),
        vec![],
    )
}
pub fn vertices(poly: &Polygon<f64>) -> Vec<P> {
    let mut v: Vec<_> = poly.exterior().0.iter().map(|c| P::new(c.x, c.y)).collect();
    if v.len() > 1 && v.first() == v.last() {
        v.pop();
    }
    v
}
pub fn edges(v: &[P]) -> Vec<Segment> {
    (0..v.len())
        .map(|i| Segment {
            a: v[i],
            b: v[(i + 1) % v.len()],
        })
        .collect()
}
pub fn clean(v: &[P]) -> Result<Vec<P>> {
    ensure!(
        v.iter().all(|p| p.x.is_finite() && p.y.is_finite()),
        "non-finite polygon"
    );
    let mut out = Vec::new();
    for &p in v {
        if out.last() != Some(&p) {
            out.push(p);
        }
    }
    if out.len() > 1 && out.first() == out.last() {
        out.pop();
    }
    // Remove only exactly collinear vertices. Near-collinearity is not simplified.
    loop {
        if out.len() <= 3 {
            break;
        }
        let remove = (0..out.len()).find(|&i| {
            cross(
                out[(i + out.len() - 1) % out.len()],
                out[i],
                out[(i + 1) % out.len()],
            ) == 0.0
        });
        if let Some(i) = remove {
            out.remove(i);
        } else {
            break;
        }
    }
    ensure!(
        out.len() >= 3 && signed_area(&out).abs() > 0.0,
        "degenerate polygon"
    );
    if signed_area(&out) < 0.0 {
        out.reverse();
    }
    Ok(out)
}
pub fn convex_hull(v: &[P]) -> Vec<P> {
    let mut points = v.to_vec();
    points.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    points.dedup();
    if points.len() < 3 {
        return points;
    }
    let mut lower = Vec::new();
    for &p in &points {
        while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], p) <= 0.0 {
            lower.pop();
        }
        lower.push(p);
    }
    let mut upper = Vec::new();
    for &p in points.iter().rev() {
        while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], p) <= 0.0 {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}
pub fn triangulate(v: &[P]) -> Result<Vec<[P; 3]>> {
    let v = clean(v)?;
    let p = polygon(&v);
    // The external contour can touch itself at an articulation vertex after
    // joining pieces. Earcut handles these weakly simple contour degeneracies.
    let result: Vec<_> = p
        .earcut_triangles()
        .into_iter()
        .map(|t| {
            [
                P::new(t.v1().x, t.v1().y),
                P::new(t.v2().x, t.v2().y),
                P::new(t.v3().x, t.v3().y),
            ]
        })
        .filter(|t| signed_area(t).abs() > 0.0)
        .collect();
    let area = result.iter().map(|t| signed_area(t).abs()).sum::<f64>();
    ensure!(
        !result.is_empty() && (area - p.unsigned_area()).abs() <= p.unsigned_area() * 1e-10,
        "triangulation does not preserve contour area"
    );
    Ok(result)
}

/// Strict simple-polygon implementation retained for focused decomposition tests.
pub fn ear_clipping_simple(v: &[P]) -> Result<Vec<[P; 3]>> {
    let v = clean(v)?;
    let mut active: Vec<_> = (0..v.len()).collect();
    let mut result = Vec::with_capacity(v.len() - 2);
    while active.len() > 3 {
        let ear = (0..active.len()).find(|&k| {
            let a = v[active[(k + active.len() - 1) % active.len()]];
            let b = v[active[k]];
            let c = v[active[(k + 1) % active.len()]];
            cross(a, b, c) > 0.0
                && !active.iter().enumerate().any(|(j, &i)| {
                    j != k
                        && j != (k + active.len() - 1) % active.len()
                        && j != (k + 1) % active.len()
                        && cross(a, b, v[i]) >= 0.0
                        && cross(b, c, v[i]) >= 0.0
                        && cross(c, a, v[i]) >= 0.0
                })
        });
        let k =
            ear.ok_or_else(|| anyhow::anyhow!("ear clipping failed: polygon must be simple"))?;
        result.push([
            v[active[(k + active.len() - 1) % active.len()]],
            v[active[k]],
            v[active[(k + 1) % active.len()]],
        ]);
        active.remove(k);
    }
    result.push([v[active[0]], v[active[1]], v[active[2]]]);
    Ok(result)
}
pub fn intersection(a: Segment, b: Segment) -> Option<P> {
    let r = P::new(a.b.x - a.a.x, a.b.y - a.a.y);
    let s = P::new(b.b.x - b.a.x, b.b.y - b.a.y);
    let det = r.x * s.y - r.y * s.x;
    if det == 0.0 {
        return None;
    } // Collinear overlap endpoints are emitted separately.
    let q = P::new(b.a.x - a.a.x, b.a.y - a.a.y);
    let t = (q.x * s.y - q.y * s.x) / det;
    let u = (q.x * r.y - q.y * r.x) / det;
    if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
        Some(P::new(a.a.x + t * r.x, a.a.y + t * r.y))
    } else {
        None
    }
}
pub fn vertical(e: Segment, x: f64) -> Option<P> {
    if e.a.x == e.b.x {
        return None;
    }
    let t = (x - e.a.x) / (e.b.x - e.a.x);
    (0.0..=1.0)
        .contains(&t)
        .then(|| P::new(x, e.a.y + t * (e.b.y - e.a.y)))
}
pub fn horizontal(e: Segment, y: f64) -> Option<P> {
    if e.a.y == e.b.y {
        return None;
    }
    let t = (y - e.a.y) / (e.b.y - e.a.y);
    (0.0..=1.0)
        .contains(&t)
        .then(|| P::new(e.a.x + t * (e.b.x - e.a.x), y))
}
