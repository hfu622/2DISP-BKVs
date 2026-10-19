//! Directional NFP penetration and breakpoints of its piecewise-linear profile.
use super::{nfp::NoFit, polygon::*};
use crate::model::Budget;
use geo::{Contains, Point};

/// Squared Euclidean penetration and its gradient with respect to the moving
/// reference point. At medial-axis ties use the first nearest boundary point.
pub fn squared(nfp: &NoFit, p: P) -> (f64, P) {
    if !nfp.region.contains(&Point::new(p.x, p.y)) {
        return (0.0, P::new(0.0, 0.0));
    }
    let mut distance = f64::INFINITY;
    let mut gradient = P::new(0.0, 0.0);
    for e in nfp.boundary.iter().chain(&nfp.zero_strata) {
        let dx = e.b.x - e.a.x;
        let dy = e.b.y - e.a.y;
        let norm = dx * dx + dy * dy;
        let t = if norm > 0.0 {
            ((p.x - e.a.x) * dx + (p.y - e.a.y) * dy) / norm
        } else {
            0.0
        }
        .clamp(0.0, 1.0);
        let vx = p.x - e.a.x - t * dx;
        let vy = p.y - e.a.y - t * dy;
        let value = vx * vx + vy * vy;
        if value < distance {
            distance = value;
            gradient = P::new(2.0 * vx, 2.0 * vy);
        }
    }
    (distance, gradient)
}

pub fn depth(nfp: &NoFit, p: P) -> f64 {
    if !nfp.region.contains(&Point::new(p.x, p.y)) {
        return 0.0;
    }
    let mut distance = f64::INFINITY;
    for &e in nfp.boundary.iter().chain(&nfp.zero_strata) {
        for q in [e.a, e.b] {
            if p.x == q.x {
                distance = distance.min((q.y - p.y).abs());
            }
            if p.y == q.y {
                distance = distance.min((q.x - p.x).abs());
            }
        }
        if let Some(q) = horizontal(e, p.y) {
            distance = distance.min((q.x - p.x).abs());
        }
        if let Some(q) = vertical(e, p.x) {
            distance = distance.min((q.y - p.y).abs());
        }
    }
    distance
}

/// Superset of every slope-change point of min(horizontal, vertical) depth.
/// Brute-force line intersections replace the paper's trapezoidal partitions.
/// This preserves the continuous minimization candidates but not linear complexity.
pub fn breakpoints(
    nfp: &NoFit,
    orthogonal: f64,
    vertical_axis: bool,
    budget: &Budget,
) -> Option<Vec<f64>> {
    let mut points = Vec::new();
    let mut lines = Vec::<(f64, f64)>::new();
    for e in &nfp.boundary {
        if budget.expired() {
            return None;
        }
        let transform = |p: P| if vertical_axis { P::new(p.y, p.x) } else { p };
        let e = Segment {
            a: transform(e.a),
            b: transform(e.b),
        };
        points.extend([e.a.x, e.b.x]);
        if let Some(p) = horizontal(e, orthogonal) {
            points.push(p.x);
            lines.extend([(1.0, -p.x), (-1.0, p.x)]);
        }
        if (e.b.x - e.a.x).abs() > 1e-14 {
            let slope = (e.b.y - e.a.y) / (e.b.x - e.a.x);
            let intercept = e.a.y - slope * e.a.x - orthogonal;
            lines.extend([(slope, intercept), (-slope, -intercept)]);
            if slope.abs() > 1e-14 {
                points.push(-intercept / slope);
            }
        }
    }
    if lines
        .len()
        .checked_mul(lines.len())
        .is_none_or(|n| n > 4_000_000)
    {
        return None;
    }
    for (i, &(a, b)) in lines.iter().enumerate() {
        if budget.expired() {
            return None;
        }
        for &(c, d) in &lines[i + 1..] {
            if (a - c).abs() > 1e-14 {
                let x = (d - b) / (a - c);
                if x.is_finite() {
                    points.push(x);
                }
            }
        }
    }
    points.retain(|x| x.is_finite());
    points.sort_by(f64::total_cmp);
    points.dedup();
    Some(points)
}
