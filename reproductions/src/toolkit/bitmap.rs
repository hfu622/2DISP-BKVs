//! Conservative raster envelopes, scan-line rectangles, and open rectangle NFPs.
use super::{placement::Bounds, polygon::*};
use crate::model::Budget;
use anyhow::{Result, ensure};

#[derive(Clone, Copy, Debug)]
pub struct Rect {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
}

/// Cover triangle portions in each pixel row. Rounded intervals cover the entire
/// polygon, including cells cut by an edge, rather than testing only cell centers.
pub fn scanline(v: &[P], step: f64, budget: &Budget) -> Result<Option<Vec<Rect>>> {
    ensure!(step.is_finite() && step > 0.0, "Invalid bitmap resolution");
    let triangles = triangulate(v)?;
    let bounds = Bounds::of(v);
    let first = (bounds.y0 / step).floor() as i64;
    let last = (bounds.y1 / step).ceil() as i64;
    ensure!(
        last.checked_sub(first)
            .is_some_and(|rows| (0..=4096).contains(&rows)),
        "Bitmap row resource limit exceeded"
    );
    let mut result = Vec::new();
    for row in first..last {
        if budget.expired() {
            return Ok(None);
        }
        let y0 = row as f64 * step;
        let y1 = y0 + step;
        let mut runs = Vec::new();
        for t in &triangles {
            let mut points = Vec::new();
            for &p in t {
                if p.y >= y0 && p.y <= y1 {
                    points.push(p);
                }
            }
            for i in 0..3 {
                for y in [y0, y1] {
                    if let Some(p) = horizontal(
                        Segment {
                            a: t[i],
                            b: t[(i + 1) % 3],
                        },
                        y,
                    ) {
                        points.push(p);
                    }
                }
            }
            if points.len() < 2 {
                continue;
            }
            let a = points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
            let b = points.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
            if b > a {
                runs.push(((a / step).floor() * step, (b / step).ceil() * step));
            }
        }
        runs.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut merged: Vec<(f64, f64)> = Vec::new();
        for (a, b) in runs {
            if let Some(last) = merged.last_mut()
                && a <= last.1
            {
                last.1 = last.1.max(b);
            } else {
                merged.push((a, b));
            }
        }
        result.extend(merged.into_iter().map(|(x0, x1)| Rect { x0, y0, x1, y1 }));
    }
    Ok(Some(result))
}

/// Equation (4.1): union of open Minkowski-difference rectangles. Keeping the
/// rectangles separate avoids expensive polygon Boolean operations.
pub fn append_nfp(
    out: &mut Vec<Rect>,
    fixed: &[Rect],
    moving: &[Rect],
    offset: P,
    gap: f64,
    budget: &Budget,
) -> Result<bool> {
    let count = fixed
        .len()
        .checked_mul(moving.len())
        .and_then(|n| out.len().checked_add(n));
    ensure!(
        count.is_some_and(|n| n <= 1_000_000),
        "Rectangle NFP resource limit exceeded"
    );
    for a in fixed {
        if budget.expired() {
            return Ok(false);
        }
        for b in moving {
            out.push(Rect {
                x0: a.x0 - b.x1 + offset.x - gap,
                x1: a.x1 - b.x0 + offset.x + gap,
                y0: a.y0 - b.y1 + offset.y - gap,
                y1: a.y1 - b.y0 + offset.y + gap,
            });
        }
    }
    Ok(true)
}

/// Sweep the free strip axis. At every opening/closing event, merge the active
/// forbidden intervals on the fixed axis and return its first uncovered point.
/// Open endpoints allow exact contact. This vector sweep replaces the 2-3 tree.
pub fn find_bl(rectangles: &[Rect], ifp: Bounds, budget: &Budget) -> Option<P> {
    if ifp.x1 < ifp.x0 || ifp.y1 < ifp.y0 {
        return None;
    }
    let mut events = vec![ifp.x0];
    events.extend(
        rectangles
            .iter()
            .flat_map(|r| [r.x0, r.x1])
            .filter(|x| *x >= ifp.x0 && *x <= ifp.x1),
    );
    events.sort_by(f64::total_cmp);
    events.dedup();
    for x in events {
        if budget.expired() {
            return None;
        }
        let mut intervals: Vec<_> = rectangles
            .iter()
            .filter(|r| r.x0 < x && x < r.x1)
            .map(|r| (r.y0, r.y1))
            .collect();
        intervals.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut y = ifp.y0;
        for (a, b) in intervals {
            if b <= y {
                continue;
            }
            if a >= y {
                break;
            }
            y = b;
        }
        if y <= ifp.y1 {
            return Some(P::new(x, y));
        }
    }
    None
}
