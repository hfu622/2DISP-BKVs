//! Continuous placement candidates from complete NFP boundaries and strip IFP.
use super::{nfp::NoFit, polygon::*};
use crate::{
    geometry::{Entity, Geometry},
    model::*,
};
use anyhow::{Result, ensure};
use geo::{Area, BooleanOps, BoundingRect, Polygon};
use jagua_rs::collision_detection::hazards::{Hazard, filter::NoFilter};
use std::sync::Arc;

#[derive(Clone, Copy, Debug)]
pub struct Bounds {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
}
impl Bounds {
    pub fn of(v: &[P]) -> Self {
        Self {
            x0: v.iter().map(|p| p.x).fold(f64::INFINITY, f64::min),
            y0: v.iter().map(|p| p.y).fold(f64::INFINITY, f64::min),
            x1: v.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max),
            y1: v.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max),
        }
    }
    pub fn width(self) -> f64 {
        self.x1 - self.x0
    }
    pub fn height(self) -> f64 {
        self.y1 - self.y0
    }
    pub fn shift(self, p: P) -> Self {
        Self {
            x0: self.x0 + p.x,
            y0: self.y0 + p.y,
            x1: self.x1 + p.x,
            y1: self.y1 + p.y,
        }
    }
    pub fn join(self, b: Self) -> Self {
        Self {
            x0: self.x0.min(b.x0),
            y0: self.y0.min(b.y0),
            x1: self.x1.max(b.x1),
            y1: self.y1.max(b.y1),
        }
    }
    pub fn overlap(self, b: Self) -> f64 {
        (self.x1.min(b.x1) - self.x0.max(b.x0)).max(0.0)
            * (self.y1.min(b.y1) - self.y0.max(b.y0)).max(0.0)
    }
    pub fn center(self) -> P {
        P::new((self.x0 + self.x1) * 0.5, (self.y0 + self.y1) * 0.5)
    }
}
pub fn oriented(g: &Geometry, piece: usize, orientation: usize) -> Vec<P> {
    let o = &g.types[g.pieces[piece]].orientations[orientation];
    let (sin, cos) = o.angle.to_radians().sin_cos();
    g.types[g.pieces[piece]]
        .exact_vertices
        .iter()
        .map(|&(x, y)| P::new(cos * x - sin * y - o.min_x, sin * x + cos * y - o.min_y))
        .collect()
}
pub fn world(g: &Geometry, p: &Placement) -> Polygon<f64> {
    polygon(
        &oriented(g, p.piece, p.orientation)
            .iter()
            .map(|v| v.shift(P::new(p.x, p.y)))
            .collect::<Vec<_>>(),
    )
}
pub fn feasible(poly: &Polygon<f64>, others: &[Polygon<f64>], tolerance: f64) -> bool {
    let a = poly.bounding_rect().unwrap();
    others.iter().all(|b| {
        let r = b.bounding_rect().unwrap();
        a.max().x <= r.min().x
            || r.max().x <= a.min().x
            || a.max().y <= r.min().y
            || r.max().y <= a.min().y
            || poly.intersection(b).unsigned_area() <= tolerance
    })
}
/// A jagua edge hit can be boundary contact. Admit that hit only when the
/// complete original f64 contour is inside the strip and has zero intersection
/// area with every stationary contour. Never perturb an exact-fit placement.
pub fn original_contact_clear(g: &Geometry, poly: &Polygon<f64>, others: &[Polygon<f64>]) -> bool {
    let Some(b) = poly.bounding_rect() else {
        return false;
    };
    // Common-frame reversal can turn exact contacts into ~1e-14 area
    // slivers through cancellation. Use the same original-contour area gate
    // as Geometry::validate; Jagua remains an accelerator, not the f64 oracle.
    let eps = g.height * 1e-12;
    b.min().x >= -eps
        && b.min().y >= -eps
        && b.max().x <= g.upper + eps
        && b.max().y <= g.height + eps
        && feasible(poly, others, g.total_area * 1e-12)
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlacementDecision {
    Clear,
    BoundaryContact,
    InvalidPlacement,
    InvalidPrefix,
    StripBounds,
    OriginalOverlap,
    UnresolvedJaguaCollision,
}
pub fn jagua_clear(g: &Geometry, placed: &[Placement]) -> Result<bool> {
    let mut engine = g.engine()?;
    let mut precise = Vec::new();
    for p in placed {
        if !g.fit(p) {
            return Ok(false);
        }
        let shape = Arc::new(g.world(p));
        let poly = world(g, p);
        if !feasible(&poly, &precise, g.total_area * 1e-12)
            || (engine.detect_poly_collision(&shape, &NoFilter)
                && !original_contact_clear(g, &poly, &precise))
        {
            return Ok(false);
        }
        engine.register_hazard(Hazard::new(Entity::Piece(p.piece), shape));
        precise.push(poly);
    }
    Ok(true)
}

/// Keep stationary hazards registered while testing alternative placements.
/// Jagua detects edge hits; original contours distinguish legal boundary contact.
pub struct FixedPlacementOracle<'a> {
    g: &'a Geometry,
    engine: jagua_rs::collision_detection::CDEngine<Entity>,
    polygons: Vec<Polygon<f64>>,
    valid: bool,
}
impl<'a> FixedPlacementOracle<'a> {
    pub fn new(g: &'a Geometry, placed: &[Placement]) -> Result<Self> {
        let mut oracle = Self {
            g,
            engine: g.engine()?,
            polygons: Vec::new(),
            valid: true,
        };
        for p in placed {
            // Floating-point contacts can depend on hazard registration order.
            // An invalid reordered prefix rejects queries; it is not a fatal
            // solver error (shrinking can simply retain the original layout).
            oracle.valid &= oracle.clear(p);
            oracle.register(p);
        }
        Ok(oracle)
    }
    pub fn decision(&self, p: &Placement) -> PlacementDecision {
        if !self.valid {
            return PlacementDecision::InvalidPrefix;
        }
        if p.piece >= self.g.pieces.len()
            || p.orientation >= self.g.types[self.g.pieces[p.piece]].orientations.len()
        {
            return PlacementDecision::InvalidPlacement;
        }
        if !self.g.fit(p) {
            return PlacementDecision::StripBounds;
        }
        let poly = world(self.g, p);
        if !feasible(&poly, &self.polygons, self.g.total_area * 1e-12) {
            return PlacementDecision::OriginalOverlap;
        }
        if !self
            .engine
            .detect_poly_collision(&self.g.world(p), &NoFilter)
        {
            return PlacementDecision::Clear;
        }
        if original_contact_clear(self.g, &poly, &self.polygons) {
            PlacementDecision::BoundaryContact
        } else if !feasible(&poly, &self.polygons, 0.0) {
            PlacementDecision::OriginalOverlap
        } else {
            PlacementDecision::UnresolvedJaguaCollision
        }
    }
    pub fn clear(&self, p: &Placement) -> bool {
        matches!(
            self.decision(p),
            PlacementDecision::Clear | PlacementDecision::BoundaryContact
        )
    }
    /// Only register a candidate after it has passed `clear`.
    pub fn register(&mut self, p: &Placement) {
        self.engine.register_hazard(Hazard::new(
            Entity::Piece(p.piece),
            Arc::new(self.g.world(p)),
        ));
        self.polygons.push(world(self.g, p));
    }
}
pub fn translated_boundary(nfp: &NoFit, offset: P) -> Vec<Segment> {
    nfp.boundary
        .iter()
        .chain(&nfp.zero_strata)
        .map(|s| Segment {
            a: s.a.shift(offset),
            b: s.b.shift(offset),
        })
        .collect()
}

/// Rectangle IFP including all four walls. TOPOS uses a free left/down origin,
/// so its admissible y interval is supplied by the current partial enclosure.
pub fn candidates(
    groups: &[Vec<Segment>],
    ifp: Bounds,
    extra: &[P],
    xcuts: &[f64],
    ycuts: &[f64],
    budget: &Budget,
) -> Result<Option<Vec<P>>> {
    let mut points = vec![
        P::new(ifp.x0, ifp.y0),
        P::new(ifp.x0, ifp.y1),
        P::new(ifp.x1, ifp.y0),
        P::new(ifp.x1, ifp.y1),
    ];
    points.extend_from_slice(extra);
    let xs: Vec<_> = [ifp.x0, ifp.x1]
        .into_iter()
        .chain(xcuts.iter().copied())
        .collect();
    let ys: Vec<_> = [ifp.y0, ifp.y1]
        .into_iter()
        .chain(ycuts.iter().copied())
        .collect();
    for group in groups {
        for &edge in group {
            if budget.expired() {
                return Ok(None);
            }
            points.extend([edge.a, edge.b]);
            for &x in &xs {
                if let Some(p) = vertical(edge, x) {
                    points.push(p);
                }
            }
            for &y in &ys {
                if let Some(p) = horizontal(edge, y) {
                    points.push(p);
                }
            }
        }
    }
    for (i, a) in groups.iter().enumerate() {
        for b in groups.iter().skip(i + 1) {
            for &s in a {
                if budget.expired() {
                    return Ok(None);
                }
                let x0 = s.a.x.min(s.b.x);
                let x1 = s.a.x.max(s.b.x);
                let y0 = s.a.y.min(s.b.y);
                let y1 = s.a.y.max(s.b.y);
                for &t in b {
                    if x1 < t.a.x.min(t.b.x)
                        || t.a.x.max(t.b.x) < x0
                        || y1 < t.a.y.min(t.b.y)
                        || t.a.y.max(t.b.y) < y0
                    {
                        continue;
                    }
                    if let Some(p) = intersection(s, t) {
                        points.push(p);
                    }
                }
                // A resource limit is explicit, never a silently truncated decoder.
                ensure!(
                    points.len() < 2_000_000,
                    "candidate arrangement exceeds two million points"
                );
            }
        }
    }
    points.retain(|p| {
        p.x >= ifp.x0 - 1e-9 && p.x <= ifp.x1 + 1e-9 && p.y >= ifp.y0 - 1e-9 && p.y <= ifp.y1 + 1e-9
    });
    points.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    points.dedup_by(|a, b| a == b);
    Ok(Some(points))
}
/// Piecewise quadratic enclosure metrics attain extrema at interval endpoints
/// or at a stationary point. Include both, rather than sampling an edge grid.
pub fn extended_points(
    edges: &[Segment],
    xcuts: &[f64],
    ycuts: &[f64],
    metric: impl Fn(P) -> f64,
) -> Vec<P> {
    let mut points = Vec::new();
    for &s in edges {
        let dx = s.b.x - s.a.x;
        let dy = s.b.y - s.a.y;
        let mut cuts = vec![0.0, 1.0];
        for &x in xcuts {
            if dx != 0.0 {
                let t = (x - s.a.x) / dx;
                if t > 0.0 && t < 1.0 {
                    cuts.push(t);
                }
            }
        }
        for &y in ycuts {
            if dy != 0.0 {
                let t = (y - s.a.y) / dy;
                if t > 0.0 && t < 1.0 {
                    cuts.push(t);
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.dedup();
        let at = |t: f64| P::new(s.a.x + dx * t, s.a.y + dy * t);
        for pair in cuts.windows(2) {
            let lo = pair[0];
            let hi = pair[1];
            let f0 = metric(at(lo));
            let f1 = metric(at(hi));
            let fm = metric(at((lo + hi) * 0.5));
            let a = 2.0 * (f1 + f0 - 2.0 * fm);
            let b = f1 - f0 - a;
            points.extend([at(lo), at(hi)]);
            if a != 0.0 {
                let u = -b / (2.0 * a);
                if u > 0.0 && u < 1.0 {
                    points.push(at(lo + (hi - lo) * u));
                }
            }
        }
    }
    points
}
