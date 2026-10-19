//! Stoyan, Zlotnik & Chugay (2012), polygon/fixed-XML-angle specialization.
//! Source Sections 4.1-4.3 and Appendix; unresolved conventions are exported.
use crate::{
    backends::{
        self,
        cbc::{Model, Row, Sense, Variable},
    },
    geometry::{self, Geometry},
    model::*,
    toolkit::{bitmap, placement::*, polygon::*},
};
use anyhow::{Result, ensure};
use geo::{Area, BooleanOps, BoundingRect};

pub const ASSUMPTIONS: &[&str] = &[
    "Primary source: user-supplied DOI10.1057/jors.2011.41, Sections4.1-4.3, Appendix and Table1",
    "Official XML angles override free continuous rotations by user instruction. This binary requires one legal physical orientation per type; no rotation variable is relaxed",
    "Source Table1 lists time-to-best, not a stopping budget. Source DNS stops each gamma after >0.6*lambda identical objectives in each of three neighbourhoods; gamma=0,1,2",
    "Source Note p386: DNS evaluates approximate starting layouts; renewed centres are stored and locally optimized after DNS",
    "Unreported per-instance constants chosen explicitly: strip approximation h=H/128, lambda=50 (source minimum), omega=3, zeta=3; objective equality and linear epsilon-active H*1e-8; raw area-unit phi epsilon H^2*1e-8; line/objective progress threshold H*1e-12",
    "Orthogonal approximation uses conservative scanline rectangles; fixed-angle group placement minimizes x then y using the union of open rectangle differences, with no extra item clearance",
    "Permutation neighbourhood samples retain their centre as one of50 points, followed by49 symmetric constrained-swap Markov-chain samples (32*n proposals each); centre retention is an explicit convention so the carried best cannot worsen, and finite samples are correlated rather than ideal independent uniform samples",
    "Metric uses minimal enclosing vertex circles and original areas; diameter uses a proven bounding-box upper bound rather than an unavailable exact maximum permutation-pair distance",
    "Recurrent centre/radius equations follow p385, mu0.8/mu1=0.6, three-sigma criterion. Unreported U is the metric change of swapping the first distinct types",
    "Convex decomposition is exact triangulation. Fixed-angle Appendix phi is min over convex-component pairs, max over edge separating functions; raw edge coefficients retained",
    "Local search solves printed feasible-direction Eq6-7 with epsilon-active branches and |direction|<=1. At a tied maximum, one branch is selected, with each alternate tied branch tried when stationary; joint tied-branch combinations are not exhaustively certified",
    "For each direction, exact first-exit bounds for all affine max branches and container edges give a feasible line step; full original-contour QA rejects any infeasible move",
    "Unreported numerical termination convention: a feasible non-worsening line step improving length by at most H*1e-8 ends that local procedure; this is not a proof of nonsmooth local optimality",
    "CBC/CLP substitutes the Meszaros interior-point LP solver; solver path and floating-point results need not match the author implementation",
    "A 600s project guard preserves an intermediate approximation if the natural protocol does not complete; it is not an invented paper budget",
    "Missing approximation/sampling/tie details prevent certification as a strict full-protocol reproduction; all phases and interruption state are recorded",
];
#[derive(Clone, Copy, Debug)]
struct Branch {
    nx: f64,
    ny: f64,
    rhs: f64,
}
impl Branch {
    fn at(self, a: P, b: P) -> f64 {
        self.nx * (b.x - a.x) + self.ny * (b.y - a.y) - self.rhs
    }
}
fn branches(a: &[P], b: &[P]) -> Vec<Branch> {
    let mut out = Vec::new();
    for (own, other, sign) in [(a, b, 1.), (b, a, -1.)] {
        for e in edges(own) {
            let nx = e.b.y - e.a.y;
            let ny = e.a.x - e.b.x;
            let edge = nx * e.a.x + ny * e.a.y;
            let vertex = other
                .iter()
                .map(|p| nx * p.x + ny * p.y)
                .fold(f64::INFINITY, f64::min);
            out.push(Branch {
                nx: sign * nx,
                ny: sign * ny,
                rhs: edge - vertex,
            });
        }
    }
    out
}
#[derive(Clone)]
struct Pair {
    i: usize,
    j: usize,
    b: Vec<Branch>,
}
fn pairs(g: &Geometry, budget: &Budget) -> Result<Vec<Pair>> {
    let mut parts = Vec::new();
    for i in 0..g.pieces.len() {
        let mut ts = triangulate(&oriented(g, i, 0))?;
        for t in &mut ts {
            if signed_area(t) < 0. {
                t.reverse();
            }
        }
        parts.push(ts);
    }
    let mut result = Vec::new();
    for i in 0..parts.len() {
        for j in i + 1..parts.len() {
            for a in &parts[i] {
                for b in &parts[j] {
                    if budget.expired() {
                        anyhow::bail!("phi preparation interrupted");
                    }
                    result.push(Pair {
                        i,
                        j,
                        b: branches(a, b),
                    });
                }
            }
        }
    }
    Ok(result)
}
fn valid(g: &Geometry, l: &Layout) -> bool {
    let eps = g.height.max(1.) * 1e-8;
    let mut prior = Vec::<geo::Polygon<f64>>::new();
    if l.placements.len() != g.pieces.len() {
        return false;
    }
    for p in &l.placements {
        let v = oriented(g, p.piece, p.orientation);
        let b = Bounds::of(&v);
        if p.x + b.x0 < -eps
            || p.y + b.y0 < -eps
            || p.y + b.y1 > g.height + eps
            || p.x + b.x1 > l.length + eps
        {
            return false;
        }
        let shape = world(g, p);
        let b = shape.bounding_rect().unwrap();
        for other in &prior {
            let a = other.bounding_rect().unwrap();
            if a.min().x < b.max().x
                && b.min().x < a.max().x
                && a.min().y < b.max().y
                && b.min().y < a.max().y
                && shape.intersection(other).unsigned_area() > g.total_area * 1e-11
            {
                return false;
            }
        }
        prior.push(shape);
    }
    true
}
// Small benchmark polygons: enumerate support pairs/triples for the true MEC.
fn radius(v: &[P]) -> f64 {
    let mut best = f64::INFINITY;
    let mut check = |c: P, r: f64| {
        if r < best && v.iter().all(|p| (p.x - c.x).hypot(p.y - c.y) <= r + 1e-8) {
            best = r;
        }
    };
    for i in 0..v.len() {
        for j in i + 1..v.len() {
            let a = v[i];
            let b = v[j];
            let c = P::new((a.x + b.x) / 2., (a.y + b.y) / 2.);
            check(c, (a.x - c.x).hypot(a.y - c.y));
            for &c in &v[j + 1..] {
                let d = 2. * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
                if d.abs() < 1e-12 {
                    continue;
                }
                let aa = a.x * a.x + a.y * a.y;
                let bb = b.x * b.x + b.y * b.y;
                let cc = c.x * c.x + c.y * c.y;
                let z = P::new(
                    (aa * (b.y - c.y) + bb * (c.y - a.y) + cc * (a.y - b.y)) / d,
                    (aa * (c.x - b.x) + bb * (a.x - c.x) + cc * (b.x - a.x)) / d,
                );
                check(z, (z.x - a.x).hypot(z.y - a.y));
            }
        }
    }
    best
}
struct GridShape {
    columns: Vec<u128>,
    ox: i32,
    oy: i32,
    h: usize,
}
struct Initializer {
    shapes: Vec<GridShape>,
    step: f64,
}
impl Initializer {
    fn new(g: &Geometry, budget: &Budget) -> Result<Self> {
        let mut shapes = Vec::new();
        let step = g.height / 128.;
        for i in 0..g.pieces.len() {
            let v = oriented(g, i, 0);
            let rs = bitmap::scanline(&v, step, budget)?
                .ok_or_else(|| anyhow::anyhow!("approximation interrupted"))?;
            let ox = (rs.iter().map(|r| r.x0).fold(f64::INFINITY, f64::min) / step).round() as i32;
            let oy = (rs.iter().map(|r| r.y0).fold(f64::INFINITY, f64::min) / step).round() as i32;
            let x1 =
                (rs.iter().map(|r| r.x1).fold(f64::NEG_INFINITY, f64::max) / step).round() as i32;
            let y1 =
                (rs.iter().map(|r| r.y1).fold(f64::NEG_INFINITY, f64::max) / step).round() as i32;
            ensure!(
                y1 - oy <= 128,
                "orthogonal approximation exceeds fixed strip height"
            );
            let mut columns = vec![0u128; (x1 - ox) as usize];
            for r in rs {
                let a = (r.x0 / step).round() as i32 - ox;
                let b = (r.x1 / step).round() as i32 - ox;
                let y = (r.y0 / step).round() as i32 - oy;
                for x in a..b {
                    columns[x as usize] |= 1u128 << y;
                }
            }
            shapes.push(GridShape {
                columns,
                ox,
                oy,
                h: (y1 - oy) as usize,
            });
        }
        Ok(Self { shapes, step })
    }
    fn build(&self, g: &Geometry, order: &[usize], budget: &Budget) -> Result<Option<Layout>> {
        let mut poses = Vec::<Placement>::new();
        let mut board = Vec::<u128>::new();
        for &piece in order {
            if budget.expired() {
                return Ok(None);
            }
            let s = &self.shapes[piece];
            let mut chosen = None;
            // All rectangle edges are multiples of h. Lexicographic integer scanning
            // is exactly the same group-variable minimum as the open-rectangle sweep,
            // without allocating every pair's Minkowski rectangles for every sample.
            'scan: for x in 0..=board.len() {
                if budget.expired() {
                    return Ok(None);
                }
                for y in 0..=128 - s.h {
                    if s.columns
                        .iter()
                        .enumerate()
                        .all(|(i, &c)| board.get(x + i).is_none_or(|&b| b & (c << y) == 0))
                    {
                        chosen = Some((x, y));
                        break 'scan;
                    }
                }
            }
            let Some((x, y)) = chosen else {
                return Ok(None);
            };
            board.resize(board.len().max(x + s.columns.len()), 0);
            for (i, &c) in s.columns.iter().enumerate() {
                board[x + i] |= c << y;
            }
            poses.push(Placement {
                piece,
                orientation: 0,
                x: (x as f64 - s.ox as f64) * self.step,
                y: (y as f64 - s.oy as f64) * self.step,
            });
        }
        poses.sort_by_key(|p| p.piece);
        let l = Layout {
            length: geometry::length(g, &poses),
            placements: poses,
        };
        ensure!(
            valid(g, &l),
            "orthogonal approximation failed original-contour validation"
        );
        Ok(Some(l))
    }
}
#[derive(Clone)]
struct Sample {
    order: Vec<usize>,
    l: Layout,
}
#[derive(Clone)]
struct Hood {
    centre: Sample,
    radius: f64,
    mean: f64,
    sd: f64,
}
fn distance(a: &[usize], b: &[usize], rs: &[(f64, f64)], gamma: i32) -> f64 {
    let n = a.len();
    a.iter()
        .zip(b)
        .enumerate()
        .map(|(k, (&i, &j))| {
            let (r, s) = rs[i];
            let (rr, ss) = rs[j];
            ((n - k) as f64 / n as f64).powi(gamma) * ((r - rr).powi(2) + (s - ss).powi(2))
        })
        .sum::<f64>()
        .sqrt()
}
fn diameter(rs: &[(f64, f64)], gamma: i32) -> f64 {
    let range = |j: usize| {
        let f = |p: &(f64, f64)| if j == 0 { p.0 } else { p.1 };
        rs.iter().map(f).fold(f64::NEG_INFINITY, f64::max)
            - rs.iter().map(f).fold(f64::INFINITY, f64::min)
    };
    let q = range(0).powi(2) + range(1).powi(2);
    (q * (0..rs.len())
        .map(|i| ((rs.len() - i) as f64 / rs.len() as f64).powi(gamma))
        .sum::<f64>())
    .sqrt()
}
fn summarize(v: &[Sample], eps: f64) -> (f64, f64, usize) {
    let mean = v.iter().map(|s| s.l.length).sum::<f64>() / v.len() as f64;
    let sd = (v.iter().map(|s| (s.l.length - mean).powi(2)).sum::<f64>() / v.len() as f64).sqrt();
    let count = v
        .iter()
        .map(|a| {
            v.iter()
                .filter(|b| (a.l.length - b.l.length).abs() <= eps)
                .count()
        })
        .max()
        .unwrap_or(0);
    (mean, sd, count)
}
fn sample_hood(
    init: &Initializer,
    g: &Geometry,
    centre: &Sample,
    r: f64,
    gamma: i32,
    rs: &[(f64, f64)],
    rng: &mut Rng,
    budget: &Budget,
    stats: &mut Stats,
) -> Result<Vec<Sample>> {
    let mut out = vec![centre.clone()];
    let mut chain = centre.order.clone();
    let n = chain.len();
    for _ in 0..49 {
        if budget.expired() {
            break;
        }
        for _ in 0..32 * n {
            let i = rng.index(n);
            let j = rng.index(n);
            chain.swap(i, j);
            if distance(&chain, &centre.order, rs, gamma) > r + 1e-9 {
                chain.swap(i, j);
            }
        }
        let Some(l) = init.build(g, &chain, budget)? else {
            break;
        };
        stats.evaluations += 1;
        out.push(Sample {
            order: chain.clone(),
            l,
        });
    }
    Ok(out)
}
// First exit of max_j(a_j+t*b_j)>=0 on the ray t>=0. The union of
// closed half-lines admits exact contact and can reconnect immediately at a tie.
fn ray_exit(lines: &[(f64, f64)]) -> f64 {
    let mut end = 0.0f64;
    let mut opening = f64::INFINITY;
    for &(a, b) in lines {
        if b.abs() < 1e-14 {
            if a >= -1e-9 {
                return f64::INFINITY;
            }
        } else if b > 0. {
            if a >= -1e-9 {
                return f64::INFINITY;
            }
            opening = opening.min(-a / b);
        } else if a >= -1e-9 {
            end = end.max(a.max(0.) / -b);
        }
    }
    if opening <= end + 1e-9 {
        f64::INFINITY
    } else {
        end
    }
}
fn local(
    g: &Geometry,
    start: &Layout,
    ps: &[Pair],
    budget: &Budget,
    stats: &mut Stats,
) -> Result<Layout> {
    let n = g.pieces.len();
    let dim = 2 * n + 1;
    let mut l = start.clone();
    let eps = g.height.max(1.) * 1e-8;
    let phi_eps=g.height.max(1.).powi(2)*1e-8;
    let step_eps=g.height.max(1.)*1e-12;
    loop {
        if budget.expired() {
            return Ok(l);
        }
        let points = l
            .placements
            .iter()
            .map(|p| P::new(p.x, p.y))
            .collect::<Vec<_>>();
        let mut rows = Vec::<(Vec<(usize, f64)>, f64, f64)>::new();
        for (i, p) in l.placements.iter().enumerate() {
            let b = Bounds::of(&oriented(g, p.piece, p.orientation));
            rows.push((vec![(2 * i, 1.)], p.x + b.x0,eps));
            rows.push((vec![(2 * i + 1, 1.)], p.y + b.y0,eps));
            rows.push((vec![(2 * i + 1, -1.)], g.height - p.y - b.y1,eps));
            rows.push((vec![(2 * i, -1.), (2 * n, 1.)], l.length - p.x - b.x1,eps));
        }
        let mut alternatives = Vec::<(usize, Vec<Vec<(usize, f64)>>)>::new();
        for pair in ps {
            let values = pair
                .b
                .iter()
                .map(|b| b.at(points[pair.i], points[pair.j]))
                .collect::<Vec<_>>();
            let best = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            ensure!(best >= -eps * 100., "phi start infeasible");
            if best > phi_eps {
                continue;
            }
            let mut candidates = Vec::new();
            for (b, &v) in pair.b.iter().zip(&values) {
                if (v - best).abs() <= phi_eps {
                    let mut t = vec![
                        (2 * pair.i, -b.nx),
                        (2 * pair.i + 1, -b.ny),
                        (2 * pair.j, b.nx),
                        (2 * pair.j + 1, b.ny),
                    ];
                    t.retain(|(_, v)| *v != 0.);
                    candidates.push((v,t));
                }
            }
            // A smooth feasible subregion must use an actual maximizing branch.
            // Selecting the first epsilon-near branch could choose a negative
            // separating function and produce spurious zero-length line steps.
            candidates.sort_by(|a,b|b.0.total_cmp(&a.0));
            let candidates=candidates.into_iter().map(|(_,t)|t).collect::<Vec<_>>();
            let index = rows.len();
            rows.push((candidates[0].clone(), best,phi_eps));
            if candidates.len() > 1 {
                alternatives.push((index, candidates));
            }
        }
        let choices = std::iter::once(None).chain(
            alternatives
                .iter()
                .flat_map(|(i, ts)| ts.iter().skip(1).map(move |t| Some((*i, t)))),
        );
        let mut chosen = None;
        for choice in choices {
            if budget.expired() {
                return Ok(l);
            }
            let mut m = Model::default();
            m.variables = vec![
                Variable {
                    lower: -1.,
                    upper: 1.,
                    objective: 0.,
                    binary: false
                };
                dim
            ];
            m.variables.push(Variable {
                lower: 0.,
                upper: 1.,
                objective: -1.,
                binary: false,
            });
            m.rows.push(Row {
                terms: vec![(2 * n, -1.), (dim, -1.)],
                sense: Sense::Greater,
                rhs: 0.,
            });
            for (i, (terms, value,active_eps)) in rows.iter().enumerate() {
                if *value > *active_eps {
                    continue;
                }
                let mut t = if choice.is_some_and(|(j, _)| j == i) {
                    choice.unwrap().1.clone()
                } else {
                    terms.clone()
                };
                t.push((dim, -1.));
                m.rows.push(Row {
                    terms: t,
                    sense: Sense::Greater,
                    rhs: 0.,
                });
            }
            stats.backend_calls += 1;
            let Some(s) = backends::cbc::solve(&m, &backends::Settings::default(), budget)? else {
                stats
                    .backend_events
                    .push("Stoyan feasible-direction LP interrupted/unresolved".into());
                return Ok(l);
            };
            if s.values[dim] > 1e-9 {
                chosen = Some(s.values);
                break;
            }
        }
        let Some(d) = chosen else {
            stats.completed_variants += 1;
            return Ok(l);
        };
        let mut t = f64::INFINITY;
        // Every containment inequality remains active along the complete line segment.
        for (terms, value,_) in rows.iter().take(4 * n) {
            let slope = terms.iter().map(|&(i, c)| c * d[i]).sum::<f64>();
            if slope < -1e-12 {
                t = t.min(value.max(0.) / -slope);
            }
        }
        for p in ps {
            let lines = p
                .b
                .iter()
                .map(|b| {
                    (
                        b.at(points[p.i], points[p.j]),
                        b.nx * (d[2 * p.j] - d[2 * p.i]) + b.ny * (d[2 * p.j + 1] - d[2 * p.i + 1]),
                    )
                })
                .collect::<Vec<_>>();
            t = t.min(ray_exit(&lines));
        }
        if !t.is_finite() || t < step_eps {
            stats
                .backend_events
                .push(format!("Stoyan zero first-exit step={t}; tied branch convention unresolved"));
            return Ok(l);
        }
        let mut next = l.clone();
        for (i, p) in next.placements.iter_mut().enumerate() {
            p.x += t * d[2 * i];
            p.y += t * d[2 * i + 1];
        }
        next.length = geometry::length(g, &next.placements);
        if next.length > l.length + step_eps || !valid(g, &next) {
            stats
                .backend_events
                .push("Stoyan feasible-direction candidate rejected by original-contour QA".into());
            return Ok(l);
        }
        if l.length-next.length<=eps{
            stats.completed_variants+=1;
            stats.backend_events.push(format!("Stoyan declared numerical objective stop: improvement={} tolerance={eps}; not a nonsmooth local-optimality proof",l.length-next.length));
            return Ok(if next.length<l.length{next}else{l});
        }
        stats.accepted_moves += 1;
        l = next;
    }
}
pub fn solve(g: &Geometry, budget: &Budget, seed: u64, stats: &mut Stats) -> Result<Layout> {
    ensure!(
        g.types.iter().all(|t| {
            let a = t.orientations[0].angle;
            t.orientations.iter().all(|o| {
                (o.angle - a)
                    .rem_euclid(360.)
                    .min((a - o.angle).rem_euclid(360.))
                    < 1e-6
            })
        }),
        "Stoyan fixed-angle specialization requires one physical orientation"
    );
    let init = Initializer::new(g, budget)?;
    let mut rng = Rng::new(seed);
    let rs = (0..g.pieces.len())
        .map(|i| (radius(&oriented(g, i, 0)), g.types[g.pieces[i]].area))
        .collect::<Vec<_>>();
    let mut pool = Vec::<Sample>::new();
    for _ in 0..50 {
        if budget.expired() {
            break;
        }
        let mut order = (0..g.pieces.len()).collect::<Vec<_>>();
        rng.shuffle(&mut order);
        let Some(l) = init.build(g, &order, budget)? else {
            break;
        };
        stats.evaluations += 1;
        pool.push(Sample { order, l });
    }
    ensure!(
        !pool.is_empty(),
        "DNS produced no complete initial approximation"
    );
    pool.sort_by(|a, b| a.l.length.total_cmp(&b.l.length));
    let mut best = pool[0].l.clone();
    let mut stored = vec![pool[0].clone()];
    let eps = g.height.max(1.) * 1e-8;
    let mut hoods = Vec::new();
    let beta = diameter(&rs, 0);
    let mut newest = (pool[0].clone(), pool[0].clone());
    for s in pool.iter().take(3) {
        let v = sample_hood(&init, g, s, beta * 0.25, 0, &rs, &mut rng, budget, stats)?;
        if v.is_empty() {
            break;
        }
        let (mean, sd, _) = summarize(&v, eps);
        let chosen = v
            .iter()
            .min_by(|a, b| a.l.length.total_cmp(&b.l.length))
            .unwrap()
            .clone();
        if chosen.l.length < best.length {
            best = chosen.l.clone();
            newest = (chosen.clone(), s.clone());
        }
        stored.push(chosen);
        hoods.push(Hood {
            centre: s.clone(),
            radius: beta,
            mean,
            sd,
        });
    }
    let mut older = newest.clone();
    let mut gammas = 0;
    let mut iterations = 0;
    let mut third_length = best.length;
    for gamma in 0..3 {
        if budget.expired() || hoods.len() != 3 {
            break;
        }
        let beta = diameter(&rs, gamma);
        let mut global_radius = beta;
        if gamma > 0 {
            let order = &hoods[0].centre.order;
            let mut swapped = order.clone();
            let j = (1..order.len())
                .find(|&j| rs[order[j]] != rs[order[0]])
                .unwrap_or(0);
            swapped.swap(0, j);
            let u = distance(order, &swapped, &rs, gamma);
            for h in &mut hoods {
                h.radius = u.max(eps);
            }
        }
        let mut stage_iteration = 0;
        loop {
            if budget.expired() {
                break;
            }
            let previous = hoods.clone();
            let promising = previous
                .iter()
                .min_by(|a, b| (a.mean - 3. * a.sd).total_cmp(&(b.mean - 3. * b.sd)))
                .unwrap()
                .centre
                .clone();
            let incumbent = if newest.0.l.length < older.0.l.length - eps {
                newest.clone()
            } else {
                older.clone()
            };
            let centres = [incumbent.0.clone(), incumbent.1.clone(), promising];
            if stage_iteration > 0 {
                if newest.0.l.length >= older.0.l.length - eps
                    && newest.0.l.length >= third_length - eps
                {
                    global_radius *= 0.6;
                } else if newest.0.l.length < older.0.l.length - eps {
                    global_radius = (global_radius / 0.6).min(beta);
                }
            }
            let mut converged = true;
            let mut next = Vec::new();
            let mut stage_best = None::<(Sample, Sample)>;
            for (i, centre) in centres.into_iter().enumerate() {
                let mut h = previous[i].clone();
                h.centre = centre;
                if stage_iteration > 0 {
                    h.radius = if previous[i].mean - 3. * previous[i].sd >= incumbent.0.l.length {
                        h.radius * 0.8
                    } else {
                        (h.radius / 0.8).min(global_radius)
                    };
                }
                let v = sample_hood(
                    &init, g, &h.centre, h.radius, gamma, &rs, &mut rng, budget, stats,
                )?;
                if v.len() != 50 {
                    converged = false;
                    break;
                }
                let (mean, sd, count) = summarize(&v, eps);
                converged &= count > 30;
                let candidate = v
                    .iter()
                    .min_by(|a, b| a.l.length.total_cmp(&b.l.length))
                    .unwrap()
                    .clone();
                if stage_best
                    .as_ref()
                    .is_none_or(|s| candidate.l.length < s.0.l.length)
                {
                    stage_best = Some((candidate.clone(), h.centre.clone()));
                }
                if candidate.l.length < h.centre.l.length - eps {
                    stored.push(candidate.clone());
                }
                if candidate.l.length < best.length {
                    best = candidate.l.clone();
                }
                h.mean = mean;
                h.sd = sd;
                next.push(h);
            }
            iterations += 1;
            stats.completed_generations = Some(iterations);
            println!(
                "Stoyan DNS gamma={gamma} iteration={iterations} stored={} L={} converged={converged}",
                stored.len(),
                best.length
            );
            if next.len() != 3 {
                break;
            }
            hoods = next;
            third_length = older.0.l.length;
            older = newest;
            newest = stage_best.unwrap();
            stage_iteration += 1;
            if converged {
                gammas += 1;
                break;
            }
        }
    }
    stats.population_initialized = Some(pool.len());
    stats.backend_events.push(format!(
        "DNS phases={gammas}/3, initial samples={}, retained renewed centres={}",
        pool.len(),
        stored.len()
    ));
    let mut optimized = 0;
    if !budget.expired() {
        let ps = pairs(g, budget)?;
        stored.sort_by(|a, b| a.l.length.total_cmp(&b.l.length));
        for s in stored {
            if budget.expired() {
                break;
            }
            let l = local(g, &s.l, &ps, budget, stats)?;
            optimized += 1;
            if l.length < best.length {
                best = l;
            }
            println!("Stoyan local={optimized} bestL={}", best.length);
        }
    }
    stats
        .backend_events
        .push(format!("locally optimized centre starts={optimized}"));
    stats.paper_variant = Some(
        "Stoyan2012 DNS + modified feasible directions; official fixed-angle policy adaptation"
            .into(),
    );
    stats.paper_run_status = Some(
        if budget.expired() {
            "completed_initial_approximations_DNS_project_truncation"
        } else if stats.completed_variants!=optimized {
            "partial_Stoyan2012_unresolved_local_steps"
        } else {
            "completed_Stoyan2012_documented_conventions"
        }
        .into(),
    );
    ensure!(valid(g, &best), "Stoyan final original-contour check");
    Ok(best)
}
#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> Geometry {
        let instance=serde_json::from_value(serde_json::json!({"name":"fixture","strip_height":6,"items":[{"id":0,"demand":2,"shape":{"type":"simple_polygon","data":[[0,0],[2,0],[2,2],[0,2]]},"orientation":{"rotation":{"mode":"discrete","angles":[0]}}}]})).unwrap();
        let input = crate::input::LoadedInput {
            instance,
            item_labels: Default::default(),
            source_path: Default::default(),
            exact_vertices: std::collections::BTreeMap::from([(
                0,
                vec![(0., 0.), (2., 0.), (2., 2.), (0., 2.)],
            )]),
            source_up_left: false,
        };
        Geometry::from_input(&input).unwrap()
    }
    fn mixed_fixture() -> Geometry {
        let shapes=vec![vec![(0.,0.),(2.,0.),(2.,2.),(0.,2.)],vec![(0.,0.),(3.,0.),(0.,1.)],vec![(0.,0.),(3.,0.),(3.,1.),(1.,1.),(1.,3.),(0.,3.)]];
        let items=shapes.iter().enumerate().map(|(id,v)|serde_json::json!({"id":id,"demand":if id<2{2}else{1},"shape":{"type":"simple_polygon","data":v},"orientation":{"rotation":{"mode":"discrete","angles":[0]}}})).collect::<Vec<_>>();
        let input=crate::input::LoadedInput{instance:serde_json::from_value(serde_json::json!({"name":"mixed","strip_height":6,"items":items})).unwrap(),item_labels:Default::default(),source_path:Default::default(),exact_vertices:shapes.into_iter().enumerate().map(|(i,v)|(i as u64,v)).collect(),source_up_left:false};Geometry::from_input(&input).unwrap()
    }
    #[test]
    fn phi_and_ray_cover_contacts_and_disconnected_convex_free_region() {
        let a = vec![
            P::new(0., 0.),
            P::new(2., 0.),
            P::new(2., 1.),
            P::new(0., 1.),
        ];
        let bs = branches(&a, &a);
        for (x, y) in [(0., 0.), (1., 0.), (2., 0.), (3., 0.), (0., 1.), (-2., 0.)] {
            let phi = bs
                .iter()
                .map(|b| b.at(P::new(0., 0.), P::new(x, y)))
                .fold(f64::NEG_INFINITY, f64::max);
            let overlap = polygon(&a)
                .intersection(&polygon(
                    &a.iter()
                        .map(|p| P::new(p.x + x, p.y + y))
                        .collect::<Vec<_>>(),
                ))
                .unsigned_area();
            assert_eq!(phi < -1e-8, overlap > 1e-8);
        }
        assert!((ray_exit(&[(2., -1.), (-4., 1.)]) - 2.).abs() < 1e-9);
        assert!(ray_exit(&[(2., -1.), (-2., 1.)]).is_infinite());
        assert_eq!(ray_exit(&[(0., -1.)]), 0.);
    }
    #[test]
    fn mec_metric_and_dns_stopping_use_geometry_and_multiplicity() {
        let r = radius(&[P::new(0., 0.), P::new(4., 0.), P::new(0., 3.)]);
        assert!((r - 2.5).abs() < 1e-9);
        let rs = vec![(1., 2.), (3., 4.), (1., 2.)];
        assert_eq!(distance(&[0, 1, 2], &[2, 1, 0], &rs, 0), 0.);
        assert!(distance(&[0, 1, 2], &[1, 0, 2], &rs, 2) > 0.);
    }
    #[test]
    fn integer_group_placement_matches_open_rectangle_sweep_and_direction_improves() {
        let g = fixture();
        let budget = Budget::new(8.);
        let init = Initializer::new(&g, &budget).unwrap();
        let l = init.build(&g, &[0, 1], &budget).unwrap().unwrap();
        let a = bitmap::scanline(&oriented(&g, 0, 0), g.height / 128., &budget)
            .unwrap()
            .unwrap();
        let mut forbidden = Vec::new();
        bitmap::append_nfp(
            &mut forbidden,
            &a,
            &a,
            P::new(l.placements[0].x, l.placements[0].y),
            0.,
            &budget,
        )
        .unwrap();
        let q = bitmap::find_bl(
            &forbidden,
            Bounds {
                x0: 0.,
                y0: 0.,
                x1: g.upper,
                y1: g.height - a.iter().map(|r| r.y1).fold(0., f64::max),
            },
            &budget,
        )
        .unwrap();
        assert!((l.placements[1].x - q.x).abs() < 1e-8 && (l.placements[1].y - q.y).abs() < 1e-8);
        let start = Layout {
            length: 6.,
            placements: vec![
                Placement {
                    piece: 0,
                    orientation: 0,
                    x: 0.,
                    y: 0.,
                },
                Placement {
                    piece: 1,
                    orientation: 0,
                    x: 4.,
                    y: 0.,
                },
            ],
        };
        let ps = pairs(&g, &budget).unwrap();
        let mut stats = Stats::default();
        let result = local(&g, &start, &ps, &budget, &mut stats).unwrap();
        assert!(valid(&g, &result));
        assert!(result.length < 5.99);
        assert!(stats.backend_calls > 0 && stats.accepted_moves > 0);
    }
    #[test]
    fn full_dns_reaches_three_weight_stops_then_optimizes_saved_centres() {
        let g = mixed_fixture();
        let mut stats = Stats::default();
        let l = solve(&g, &Budget::new(15.), 2026101506, &mut stats).unwrap();
        assert!(valid(&g, &l));
        assert_eq!(stats.population_initialized, Some(50));
        assert!(stats.completed_generations.unwrap()>=3);
        assert!(stats.backend_calls > 0);
        assert!(
            stats
                .backend_events
                .iter()
                .any(|s| s.contains("DNS phases=3/3"))
        );
        assert_eq!(
            stats.paper_run_status.as_deref(),
            Some("completed_Stoyan2012_documented_conventions")
        );
    }
}
