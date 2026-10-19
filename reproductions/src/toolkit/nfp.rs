//! Concave Minkowski differences via convex decomposition and Boolean union.
//! Reference point: translation of the unshifted input polygon, not a vertex.
use super::polygon::*;
use crate::model::Budget;
use anyhow::Result;
use geo::{Area, BooleanOps, Buffer, Contains, MultiPolygon, Point};
use std::{
    collections::{BTreeMap, VecDeque},
    sync::Arc,
};

pub struct NoFit {
    pub region: MultiPolygon<f64>,
    pub boundary: Vec<Segment>,
    pub contacts: Vec<P>,
    /// Feasible isolated points/segments erased by a regularized union.
    /// A point is represented by a segment with equal endpoints.
    pub zero_strata: Vec<Segment>,
}
impl NoFit {
    /// Retained allocations, including spare vector capacity and hole-ring
    /// descriptors. This excludes allocator bookkeeping and build temporaries.
    /// Saturation makes an unrepresentable estimate fail a bounded-cache check.
    pub fn allocated_bytes(&self) -> usize {
        use std::mem::size_of;
        let mut bytes = size_of::<Self>().saturating_add(2 * size_of::<usize>());
        bytes = bytes.saturating_add(
            self.boundary
                .capacity()
                .saturating_mul(size_of::<Segment>()),
        );
        bytes = bytes.saturating_add(self.contacts.capacity().saturating_mul(size_of::<P>()));
        bytes = bytes.saturating_add(
            self.zero_strata
                .capacity()
                .saturating_mul(size_of::<Segment>()),
        );
        bytes = bytes.saturating_add(
            self.region
                .0
                .capacity()
                .saturating_mul(size_of::<geo::Polygon>()),
        );
        for polygon in &self.region.0 {
            bytes = bytes.saturating_add(
                polygon
                    .exterior()
                    .0
                    .capacity()
                    .saturating_mul(size_of::<geo::Coord>()),
            );
            // Polygon exposes its interior rings as a slice; its hidden Vec
            // capacity is unavailable. Count the descriptors actually retained.
            bytes = bytes.saturating_add(
                polygon
                    .interiors()
                    .len()
                    .saturating_mul(size_of::<geo::LineString>()),
            );
            for ring in polygon.interiors() {
                bytes =
                    bytes.saturating_add(ring.0.capacity().saturating_mul(size_of::<geo::Coord>()));
            }
        }
        bytes
    }
}
/// Open forbidden intervals on a horizontal translation ray. Boundary-only
/// contact and inaccessible feasible holes do not become forbidden intervals.
pub fn horizontal_intervals(nfp: &NoFit, y: f64) -> Vec<(f64, f64)> {
    let mut xs = Vec::new();
    for &edge in &nfp.boundary {
        if let Some(p) = horizontal(edge, y) {
            xs.push(p.x);
        }
        if edge.a.y == y {
            xs.push(edge.a.x);
        }
        if edge.b.y == y {
            xs.push(edge.b.x);
        }
    }
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    xs.windows(2)
        .filter_map(|pair| {
            nfp.region
                .contains(&Point::new((pair[0] + pair[1]) * 0.5, y))
                .then_some((pair[0], pair[1]))
        })
        .collect()
}
/// Vertical counterpart for a reachable downward slide.
pub fn vertical_intervals(nfp: &NoFit, x: f64) -> Vec<(f64, f64)> {
    let mut ys = Vec::new();
    for &edge in &nfp.boundary {
        if let Some(p) = vertical(edge, x) {
            ys.push(p.y);
        }
        if edge.a.x == x {
            ys.push(edge.a.y);
        }
        if edge.b.x == x {
            ys.push(edge.b.y);
        }
    }
    ys.sort_by(f64::total_cmp);
    ys.dedup();
    ys.windows(2)
        .filter_map(|pair| {
            nfp.region
                .contains(&Point::new(x, (pair[0] + pair[1]) * 0.5))
                .then_some((pair[0], pair[1]))
        })
        .collect()
}
pub fn build(a: &[P], b: &[P], budget: &Budget) -> Result<Option<NoFit>> {
    let ta = triangulate(a)?;
    let tb = triangulate(b)?;
    let mut level = Vec::new();
    for x in &ta {
        for y in &tb {
            if budget.expired() {
                return Ok(None);
            }
            let sums: Vec<_> = x
                .iter()
                .flat_map(|p| y.iter().map(move |q| P::new(p.x - q.x, p.y - q.y)))
                .collect();
            level.push(MultiPolygon(vec![polygon(&convex_hull(&sums))]));
        }
    }
    // Balanced union bounds intermediate growth better than a left-fold union.
    while level.len() > 1 {
        let mut next = Vec::with_capacity(level.len().div_ceil(2));
        let mut it = level.into_iter();
        while let Some(x) = it.next() {
            if budget.expired() {
                return Ok(None);
            }
            next.push(if let Some(y) = it.next() {
                x.union(&y)
            } else {
                x
            });
        }
        level = next;
    }
    let region = level.pop().unwrap_or_else(MultiPolygon::empty);
    let boundary = boundaries(&region);
    // Retain vertex-contact events that Boolean union can erase at exact-fit
    // degeneracies. Feasibility is decided by original polygons, not Contains.
    let mut contacts: Vec<P> = a
        .iter()
        .flat_map(|p| b.iter().map(move |q| P::new(p.x - q.x, p.y - q.y)))
        .collect();
    let Some(zero_strata) = erased_feasible_strata(a, b, &region, budget)? else {
        return Ok(None);
    };
    contacts.extend(zero_strata.iter().flat_map(|s| [s.a, s.b]));
    Ok(Some(NoFit {
        region,
        boundary,
        contacts,
        zero_strata,
    }))
}

/// Complete original vertex-edge contact arrangement. Overlap status is
/// constant on each open segment after every support intersection is split.
/// Original-contour predicates retain valid 0D/1D strata inside the regularized
/// NFP; outer/hole boundaries are already retained by the regular union.
fn erased_feasible_strata(
    a: &[P],
    b: &[P],
    region: &MultiPolygon<f64>,
    budget: &Budget,
) -> Result<Option<Vec<Segment>>> {
    let convex = |v: &[P]| {
        (0..v.len()).all(|i| cross(v[i], v[(i + 1) % v.len()], v[(i + 2) % v.len()]) >= 0.)
            || (0..v.len()).all(|i| cross(v[i], v[(i + 1) % v.len()], v[(i + 2) % v.len()]) <= 0.)
    };
    if convex(a) && convex(b) {
        return Ok(Some(Vec::new()));
    }
    let mut traces = Vec::new();
    for &v in a {
        for e in edges(b) {
            traces.push(Segment {
                a: P::new(v.x - e.a.x, v.y - e.a.y),
                b: P::new(v.x - e.b.x, v.y - e.b.y),
            });
        }
    }
    for &v in b {
        for e in edges(a) {
            traces.push(Segment {
                a: P::new(e.a.x - v.x, e.a.y - v.y),
                b: P::new(e.b.x - v.x, e.b.y - v.y),
            });
        }
    }
    anyhow::ensure!(
        traces.len() <= 10_000,
        "raw NFP contact arrangement exceeds explicit limit"
    );
    let fixed = polygon(a);
    let mut out = Vec::new();
    let valid = |p: P| {
        region.contains(&Point::new(p.x, p.y))
            && fixed
                .intersection(&polygon(&b.iter().map(|q| q.shift(p)).collect::<Vec<_>>()))
                .unsigned_area()
                == 0.
    };
    for (index, &s) in traces.iter().enumerate() {
        if budget.expired() {
            return Ok(None);
        }
        let dx = s.b.x - s.a.x;
        let dy = s.b.y - s.a.y;
        let norm = dx * dx + dy * dy;
        if norm == 0. {
            continue;
        }
        let mut cuts = vec![0., 1.];
        let parameter = |p: P| ((p.x - s.a.x) * dx + (p.y - s.a.y) * dy) / norm;
        for (j, &t) in traces.iter().enumerate() {
            if j % 1024 == 0 && budget.expired() {
                return Ok(None);
            }
            if index == j {
                continue;
            }
            if let Some(p) = intersection(s, t) {
                cuts.push(parameter(p).clamp(0., 1.));
            }
            // Collinear overlaps change status only at their endpoints.
            for p in [t.a, t.b] {
                let u = parameter(p);
                if (0.0..=1.0).contains(&u) && cross(s.a, s.b, p) == 0. {
                    cuts.push(u)
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.dedup();
        let at = |u: f64| P::new(s.a.x + u * dx, s.a.y + u * dy);
        for &u in &cuts {
            let p = at(u);
            if valid(p) {
                out.push(Segment { a: p, b: p });
            }
        }
        for range in cuts.windows(2) {
            if valid(at((range[0] + range[1]) * 0.5)) {
                out.push(Segment {
                    a: at(range[0]),
                    b: at(range[1]),
                });
            }
        }
        anyhow::ensure!(
            out.len() <= 1_000_000,
            "feasible NFP strata exceed explicit memory limit"
        );
    }
    for s in &mut out {
        if s.a
            .x
            .total_cmp(&s.b.x)
            .then(s.a.y.total_cmp(&s.b.y))
            .is_gt()
        {
            std::mem::swap(&mut s.a, &mut s.b);
        }
    }
    out.sort_by(|a, b| {
        a.a.x
            .total_cmp(&b.a.x)
            .then(a.a.y.total_cmp(&b.a.y))
            .then(a.b.x.total_cmp(&b.b.x))
            .then(a.b.y.total_cmp(&b.b.y))
    });
    out.dedup_by(|a, b| a.a == b.a && a.b == b.b);
    Ok(Some(out))
}
fn boundaries(region: &MultiPolygon<f64>) -> Vec<Segment> {
    let mut boundary = Vec::new();
    for p in &region.0 {
        boundary.extend(edges(&vertices(p)));
        for ring in p.interiors() {
            let mut v: Vec<_> = ring.0.iter().map(|c| P::new(c.x, c.y)).collect();
            if v.first() == v.last() {
                v.pop();
            }
            boundary.extend(edges(&v));
        }
    }
    boundary
}
type Key = (usize, usize, usize, usize);
pub struct Cache {
    entries: BTreeMap<Key, (Arc<NoFit>, usize)>,
    fifo: VecDeque<Key>,
    bytes: usize,
    limit: usize,
    clearance: f64,
    pub hits: usize,
    pub misses: usize,
}
impl Cache {
    pub fn new(limit: usize) -> Self {
        Self {
            entries: BTreeMap::new(),
            fifo: VecDeque::new(),
            bytes: 0,
            limit,
            clearance: 0.0,
            hits: 0,
            misses: 0,
        }
    }
    /// Positive clearance is an explicit numerical modification, scoped to
    /// this cache. The default preserves the paper's contact geometry.
    pub fn with_clearance(limit: usize, clearance: f64) -> Result<Self> {
        anyhow::ensure!(
            clearance.is_finite() && clearance >= 0.0,
            "Invalid NFP clearance"
        );
        let mut cache = Self::new(limit);
        cache.clearance = clearance;
        Ok(cache)
    }
    pub fn bytes(&self) -> usize {
        self.bytes
    }
    pub fn get(
        &mut self,
        key: Key,
        a: &[P],
        b: &[P],
        budget: &Budget,
    ) -> Result<Option<Arc<NoFit>>> {
        if let Some((v, _)) = self.entries.get(&key) {
            self.hits += 1;
            return Ok(Some(v.clone()));
        }
        self.misses += 1;
        let Some(mut nfp) = build(a, b, budget)? else {
            return Ok(None);
        };
        if self.clearance > 0.0 {
            nfp.region = nfp.region.buffer(self.clearance);
            nfp.boundary = boundaries(&nfp.region);
            nfp.contacts.clear();
            nfp.zero_strata.clear();
            if budget.expired() {
                return Ok(None);
            }
        }
        let size = nfp.allocated_bytes();
        let value = Arc::new(nfp);
        if size <= self.limit {
            while size > self.limit.saturating_sub(self.bytes) {
                if let Some(old) = self.fifo.pop_front() {
                    if let Some((_, s)) = self.entries.remove(&old) {
                        self.bytes -= s;
                    }
                } else {
                    break;
                }
            }
            self.bytes += size;
            self.entries.insert(key, (value.clone(), size));
            self.fifo.push_back(key);
        }
        Ok(Some(value))
    }
}
