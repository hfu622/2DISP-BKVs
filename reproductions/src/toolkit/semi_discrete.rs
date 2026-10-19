//! Vertical resolution lines with continuous y and boundary-side labels.
//! Based on Chehrazad, Roose and Wauters, sections 3 and 4.
use super::{placement::Bounds, polygon::*};
use crate::model::Budget;
use anyhow::{Result, ensure};
use geo::{Contains, Point};
use std::collections::{BTreeMap, VecDeque};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Label {
    Middle,
    Left,
    Right,
}
#[derive(Clone, Copy, Debug)]
pub struct Interval {
    pub bottom: f64,
    pub top: f64,
    pub label: Label,
}
impl Interval {
    pub fn new(a: f64, b: f64, label: Label) -> Self {
        Self {
            bottom: a.min(b),
            top: a.max(b),
            label,
        }
    }
    pub fn shifted(self, y: f64) -> Self {
        Self {
            bottom: self.bottom + y,
            top: self.top + y,
            ..self
        }
    }
}
pub struct Shape {
    pub lines: Vec<Vec<Interval>>,
    pub bounds: Bounds,
    pub origin: P,
}
pub type Strip = BTreeMap<usize, Vec<Interval>>;
fn combine(labels: impl Iterator<Item = Label>) -> Option<Label> {
    let mut left = false;
    let mut right = false;
    for label in labels {
        match label {
            Label::Middle => return Some(Label::Middle),
            Label::Left => left = true,
            Label::Right => right = true,
        }
    }
    match (left, right) {
        (true, true) => Some(Label::Middle),
        (true, false) => Some(Label::Left),
        (false, true) => Some(Label::Right),
        _ => None,
    }
}
/// Overlay opposite side labels as M, preserving uncovered subintervals and
/// zero-length vertex events. Joining ignores neither holes nor side labels.
pub fn normalize(raw: &[Interval]) -> Vec<Interval> {
    let mut points: Vec<_> = raw.iter().flat_map(|i| [i.bottom, i.top]).collect();
    points.sort_by(f64::total_cmp);
    points.dedup();
    let mut result = Vec::<Interval>::new();
    for pair in points.windows(2) {
        let middle = (pair[0] + pair[1]) * 0.5;
        let Some(label) = combine(
            raw.iter()
                .filter(|i| i.bottom < middle && middle < i.top)
                .map(|i| i.label),
        ) else {
            continue;
        };
        if let Some(last) = result.last_mut()
            && last.top == pair[0]
            && last.label == label
        {
            last.top = pair[1];
        } else {
            result.push(Interval::new(pair[0], pair[1], label));
        }
    }
    for y in points {
        if result.iter().any(|i| i.bottom <= y && y <= i.top) {
            continue;
        }
        if let Some(label) = combine(
            raw.iter()
                .filter(|i| i.bottom == y && i.top == y)
                .map(|i| i.label),
        ) {
            result.push(Interval::new(y, y, label));
        }
    }
    result.sort_by(|a, b| a.bottom.total_cmp(&b.bottom).then(a.top.total_cmp(&b.top)));
    result
}
pub fn incompatible(a: Interval, b: Interval) -> bool {
    if matches!(
        (a.label, b.label),
        (Label::Left, Label::Right) | (Label::Right, Label::Left)
    ) {
        return false;
    }
    let lo = a.bottom.max(b.bottom);
    let hi = a.top.min(b.top);
    lo < hi
        || (lo == hi
            && (a.bottom < a.top && a.bottom < lo && lo < a.top
                || b.bottom < b.top && b.bottom < lo && lo < b.top))
}
pub fn checking_order(n: usize) -> Vec<usize> {
    if n == 0 {
        return vec![];
    }
    let mut order = vec![0];
    if n > 1 {
        order.push(n - 1);
    }
    let mut queue = VecDeque::from([(0, n - 1)]);
    while let Some((a, b)) = queue.pop_front() {
        if b > a + 1 {
            let mid = (a + b) / 2;
            order.push(mid);
            queue.push_back((a, mid));
            queue.push_back((mid, b));
        }
    }
    order
}
pub fn discretize(vertices: &[P], step: f64, budget: &Budget) -> Result<Option<Shape>> {
    ensure!(step.is_finite() && step > 0.0, "Invalid resolution");
    let source = clean(vertices)?;
    let original = Bounds::of(&source);
    let origin = P::new(original.x0, original.y0);
    let v: Vec<_> = source
        .iter()
        .map(|p| P::new(p.x - origin.x, p.y - origin.y))
        .collect();
    let bounds = Bounds::of(&v);
    let count = (bounds.width() / step).ceil() as usize + 1;
    ensure!(
        count <= 8192,
        "Semi-discrete oriented shape exceeds 8192 resolution lines"
    );
    let poly = polygon(&v);
    let es = edges(&v);
    let mut lines = vec![Vec::<Interval>::new(); count];
    for (i, line) in lines.iter_mut().enumerate() {
        if budget.expired() {
            return Ok(None);
        }
        let x = i as f64 * step;
        let mut ys = Vec::new();
        for &e in &es {
            if let Some(p) = vertical(e, x) {
                ys.push(p.y);
            }
            if e.a.x == x {
                ys.push(e.a.y);
            }
            if e.b.x == x {
                ys.push(e.b.y);
            }
        }
        ys.sort_by(f64::total_cmp);
        ys.dedup();
        for pair in ys.windows(2) {
            let y = (pair[0] + pair[1]) * 0.5;
            // A CCW upward vertical edge has material to its left;
            // a downward edge has material to its right. Interior sections
            // are M. This is the exact side test, with no epsilon probes
            // that can jump across a narrow feature.
            let boundary=combine(es.iter().filter(|e| e.a.x==x && e.b.x==x
                && e.a.y.min(e.b.y)<y && y<e.a.y.max(e.b.y))
                .map(|e|if e.b.y>e.a.y {Label::Left}else{Label::Right}));
            let label=boundary.or_else(||poly.contains(&Point::new(x,y)).then_some(Label::Middle));
            if let Some(label) = label {
                line.push(Interval::new(pair[0], pair[1], label));
            }
        }
        for j in 0..v.len() {
            let p = v[j];
            let prev = v[(j + v.len() - 1) % v.len()];
            let next = v[(j + 1) % v.len()];
            if p.x == x && cross(prev, p, next) > 0.0 {
                if prev.x >= x && next.x >= x {
                    line.push(Interval::new(p.y, p.y, Label::Right));
                } else if prev.x <= x && next.x <= x {
                    line.push(Interval::new(p.y, p.y, Label::Left));
                }
            }
        }
    }
    // Minimal convex-vertex extensions, including edges wholly between lines.
    for j in 0..v.len() {
        if budget.expired() {
            return Ok(None);
        }
        let p = v[j];
        let prev = v[(j + v.len() - 1) % v.len()];
        let next = v[(j + 1) % v.len()];
        if cross(prev, p, next) <= 0.0 {
            continue;
        }
        let grid = p.x / step;
        if grid == grid.floor() {
            continue;
        }
        let lo = grid.floor() as usize;
        let hi = lo + 1;
        for q in [prev, next] {
            let e = Segment { a: p, b: q };
            if let Some(a) = vertical(e, lo as f64 * step) {
                lines[lo].push(Interval::new(a.y, p.y, Label::Right));
                lines[hi].push(Interval::new(p.y, p.y, Label::Left));
            } else if let Some(b) = vertical(e, hi as f64 * step) {
                lines[hi].push(Interval::new(b.y, p.y, Label::Left));
                lines[lo].push(Interval::new(p.y, p.y, Label::Right));
            } else {
                lines[lo].push(Interval::new(q.y, p.y, Label::Right));
                lines[hi].push(Interval::new(q.y, p.y, Label::Left));
            }
        }
    }
    // Close separated sections bounded by the same two edges. Applying this
    // only zero-length vertex intervals, as in the paper default.
    // Section3.2 / Fig12d: extend ONLY the zero-length interval towards the
    // nearest endpoint on the adjacent line, preserving its boundary side.
    // Never fill the same gap with M on both lines: that forbids legal fits.
    for i in 0..count - 1 {
        if budget.expired() {
            return Ok(None);
        }
        let xl = i as f64 * step;
        let xr = (i + 1) as f64 * step;
        let xm = (xl + xr) * 0.5;
        let mut active: Vec<_> = es
            .iter()
            .filter_map(|&e| vertical(e, xm).map(|p| (p.y, e)))
            .collect();
        active.sort_by(|a, b| a.0.total_cmp(&b.0));
        for pair in active.chunks_exact(2) {
            if let (Some(a), Some(b), Some(c), Some(d)) = (
                vertical(pair[0].1, xl),
                vertical(pair[1].1, xl),
                vertical(pair[0].1, xr),
                vertical(pair[1].1, xr),
            ) {
                let al = a.y.min(b.y);
                let ah = a.y.max(b.y);
                let bl = c.y.min(d.y);
                let bh = c.y.max(d.y);
                let bridge = if ah < bl {
                    Some((ah, bl))
                } else if bh < al {
                    Some((bh, al))
                } else {
                    None
                };
                if bridge.is_some() {
                    if al==ah {
                        lines[i].push(Interval::new(al,al.clamp(bl,bh),Label::Right));
                    }
                    if bl==bh {
                        lines[i+1].push(Interval::new(bl,bl.clamp(al,ah),Label::Left));
                    }
                }
            }
        }
    }
    for line in &mut lines {
        *line = normalize(line);
    }
    Ok(Some(Shape {
        lines,
        bounds,
        origin,
    }))
}

/// At a fixed discrete x, monotonically raise the continuous translation y;
/// restart checks after every update so earlier lines stay collision-free.
pub fn find_y(
    shape: &Shape,
    strip: &Strip,
    column: usize,
    minimum: f64,
    maximum: f64,
    budget: &Budget,
) -> Option<f64> {
    let order = checking_order(shape.lines.len());
    let mut y = minimum;
    loop {
        if budget.expired() || y > maximum {
            return None;
        }
        let mut update = false;
        'checking: for &i in &order {
            let Some(filled) = strip.get(&(column + i)) else {
                continue;
            };
            for &local in &shape.lines[i] {
                for &obstacle in filled {
                    if incompatible(local.shifted(y), obstacle) {
                        let next = obstacle.top - local.bottom;
                        // A zero vertex in an interval needs to reach its end.
                        // Subtraction/addition at the same contact can round
                        // inconsistently: next==y while bottom+y is one ulp
                        // below obstacle.top. It does NOT make this column
                        // infeasible. Move to the next representable ordinate
                        // and restart all checks, never skip later free gaps.
                        y = if next<=y {y.next_up()} else {next};
                        update = true;
                        break 'checking;
                    }
                }
            }
        }
        if !update {
            return Some(y);
        }
    }
}
pub fn place(strip: &mut Strip, shape: &Shape, column: usize, y: f64) {
    for (i, line) in shape.lines.iter().enumerate() {
        if line.is_empty() {
            continue;
        }
        let filled = strip.entry(column + i).or_default();
        filled.extend(line.iter().map(|v| v.shifted(y)));
        *filled = normalize(filled);
    }
}
