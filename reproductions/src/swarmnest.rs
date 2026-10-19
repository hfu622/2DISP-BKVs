//! SwarmNest's direct-polygon (D-function) branch, required for Han.
use crate::{
    geometry::{self, Geometry},
    model::*,
    toolkit::{placement::*, polygon::*},
};
use anyhow::Result;

pub const ASSUMPTIONS: &str = "SwarmNest direct-polygon branch for Han, not the raster branch. Paper specifies ten heuristic seeds but omits their complete list and numerical PSO/scan settings. Documented choices: population12=10 heuristic+2 random; ascending/descending area,width,height,perimeter,rectangularity; c1=c2=2; zero initial velocity; omega .9, beta .99, omega_min .4; velocity bound n/2; scan step H/58 (one original Han unit); two compaction passes, step H/256; five explicitly listed directions despite text saying eight; one attempt per each of three local operators per personal/global best per generation. Original XML angles, quantities and contours. Han's published per-instance rotations are not explicit: fixed0 XML adaptation, comparison to77.19% is conditional. Original3600s budget exceeds30min, so this is the user-requested600s placeholder, full experiment deferred. No claim of strict fully specified original implementation equivalence.";

pub fn inverse_ranks(x: &[f64]) -> Vec<usize> {
    let mut ii: Vec<_> = (0..x.len()).collect();
    ii.sort_by(|&a, &b| x[a].total_cmp(&x[b]).then(a.cmp(&b)));
    let mut out = vec![0; x.len()];
    for (r, i) in ii.into_iter().enumerate() {
        out[i] = r
    }
    out
}
struct Particle {
    x: Vec<f64>,
    v: Vec<f64>,
    personal: Vec<f64>,
    length: f64,
}

/// Length of coincident edge intervals, including the fixed strip boundaries.
fn contact(g: &Geometry, p: &Placement, placed: &[Placement]) -> f64 {
    let poly = oriented(g, p.piece, p.orientation)
        .into_iter()
        .map(|q| q.shift(P::new(p.x, p.y)))
        .collect::<Vec<_>>();
    let mut target = vec![
        Segment {
            a: P::new(0., 0.),
            b: P::new(g.upper, 0.),
        },
        Segment {
            a: P::new(0., g.height),
            b: P::new(g.upper, g.height),
        },
        Segment {
            a: P::new(0., 0.),
            b: P::new(0., g.height),
        },
    ];
    for q in placed {
        let v = oriented(g, q.piece, q.orientation)
            .into_iter()
            .map(|a| a.shift(P::new(q.x, q.y)))
            .collect::<Vec<_>>();
        target.extend(edges(&v));
    }
    let mut total = 0.;
    let tol = g.height * 1e-9;
    for a in edges(&poly) {
        let dx = a.b.x - a.a.x;
        let dy = a.b.y - a.a.y;
        let length = dx.hypot(dy);
        if length == 0. {
            continue;
        }
        for b in &target {
            if cross(a.a, a.b, b.a).abs() > tol * length
                || cross(a.a, a.b, b.b).abs() > tol * length
            {
                continue;
            }
            let dot = |p: P| ((p.x - a.a.x) * dx + (p.y - a.a.y) * dy) / length;
            total += (dot(b.a).max(dot(b.b)).min(length) - dot(b.a).min(dot(b.b)).max(0.)).max(0.);
        }
    }
    total
}
fn evaluate(
    g: &Geometry,
    ids: &[usize],
    x: &[f64],
    budget: &Budget,
    stats: &mut Stats,
) -> Result<Option<Layout>> {
    let order: Vec<_> = inverse_ranks(x).into_iter().map(|i| ids[i]).collect();
    let step = g.height / 58.;
    let mut placed = Vec::<Placement>::new();
    for &piece in &order {
        let t = &g.types[g.pieces[piece]];
        let mut start = (0usize, 0usize);
        for p in &placed {
            let q = &g.types[g.pieces[p.piece]];
            if g.pieces[p.piece] == g.pieces[piece]
                || (q.area < t.area
                    && q.orientations[0].width < t.orientations[0].width
                    && q.orientations[0].height < t.orientations[0].height)
            {
                start = start.max(((p.x / step).round() as usize, (p.y / step).round() as usize));
            }
        }
        let oracle = FixedPlacementOracle::new(g, &placed)?;
        let mut choice = None;
        'scan: for col in start.0..=(g.upper / step).ceil() as usize {
            for row in (if col == start.0 { start.1 } else { 0 })..=58 {
                if budget.expired() {
                    return Ok(None);
                }
                let mut score = f64::NEG_INFINITY;
                for o in 0..t.orientations.len() {
                    let p = Placement {
                        piece,
                        orientation: o,
                        x: col as f64 * step,
                        y: row as f64 * step,
                    };
                    stats.collision_queries += 1;
                    if oracle.clear(&p) {
                        let value = contact(g, &p, &placed);
                        if value > score {
                            score = value;
                            choice = Some(p);
                        }
                    }
                }
                if choice.is_some() {
                    break 'scan;
                }
            }
        }
        let Some(p) = choice else { return Ok(None) };
        placed.push(p);
    }
    let mut layout = Layout {
        length: geometry::length(g, &placed),
        placements: placed,
    };
    for _ in 0..2 {
        for i in 0..layout.placements.len() {
            let others: Vec<_> = layout
                .placements
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, p)| p.clone())
                .collect();
            let oracle = FixedPlacementOracle::new(g, &others)?;
            for (dx, dy) in [(-1., 0.), (0., 1.), (0., -1.), (-1., 1.), (-1., -1.)] {
                loop {
                    if budget.expired() {
                        stats.evaluations += 1;
                        return Ok(Some(layout));
                    }
                    let mut p = layout.placements[i].clone();
                    p.x += dx * g.height / 256.;
                    p.y += dy * g.height / 256.;
                    stats.collision_queries += 1;
                    if !oracle.clear(&p) {
                        break;
                    }
                    layout.placements[i] = p;
                    layout.length = geometry::length(g, &layout.placements);
                    stats.accepted_moves += 1;
                }
            }
        }
    }
    stats.evaluations += 1;
    Ok(Some(layout))
}
pub fn solve(g: &Geometry, budget: &Budget, seed: u64, stats: &mut Stats) -> Result<Layout> {
    let mut ids: Vec<_> = (0..g.pieces.len()).collect();
    ids.sort_by(|&a, &b| {
        g.types[g.pieces[a]]
            .area
            .total_cmp(&g.types[g.pieces[b]].area)
            .then(a.cmp(&b))
    });
    let mut rng = Rng::new(seed);
    let mut particles = Vec::new();
    let mut best = g.shelf(&ids);
    let mut global = Vec::new();
    for i in 0..12 {
        if budget.expired() {
            break;
        }
        let mut order: Vec<_> = (0..ids.len()).collect();
        if i < 10 {
            let value = |rank: usize| {
                let t = &g.types[g.pieces[ids[rank]]];
                let o = &t.orientations[0];
                match i / 2 {
                    0 => t.area,
                    1 => o.width,
                    2 => o.height,
                    3 => t.perimeter,
                    _ => t.area / (o.width * o.height),
                }
            };
            order.sort_by(|&a, &b| value(a).total_cmp(&value(b)).then(a.cmp(&b)));
            if i % 2 == 1 {
                order.reverse()
            }
        } else {
            rng.shuffle(&mut order)
        }
        let x: Vec<_> = order.into_iter().map(|r| r as f64 + 1.).collect();
        let Some(layout) = evaluate(g, &ids, &x, budget, stats)? else {
            break;
        };
        if global.is_empty() || layout.length <= best.length {
            best = layout.clone();
            global = x.clone();
        }
        particles.push(Particle {
            v: vec![0.; x.len()],
            personal: x.clone(),
            x,
            length: layout.length,
        });
    }
    stats.population_initialized = Some(particles.len());
    stats.completed_generations = Some(0);
    let mut omega: f64 = 0.9;
    while particles.len() == 12 && !budget.expired() {
        let mut complete = true;
        for p in &mut particles {
            // Figure7 pseudocode updates g after each particle, so the next
            // particle must see that newly found global best in this iteration.
            let target = global.clone();
            for (j, &target_j) in target.iter().enumerate() {
                p.v[j] = (omega * p.v[j]
                    + 2. * rng.unit() * (p.personal[j] - p.x[j])
                    + 2. * rng.unit() * (target_j - p.x[j]))
                    .clamp(-(ids.len() as f64) / 2., ids.len() as f64 / 2.);
                p.x[j] += p.v[j];
            }
            let Some(layout) = evaluate(g, &ids, &p.x, budget, stats)? else {
                complete = false;
                break;
            };
            if layout.length <= p.length {
                p.personal = p.x.clone();
                p.length = layout.length;
            }
            if layout.length <= best.length {
                best = layout;
                global = p.x.clone();
            }
        }
        if !complete {
            break;
        }
        for index in 0..=particles.len() {
            let source = if index == particles.len() {
                global.clone()
            } else {
                particles[index].personal.clone()
            };
            let mut ranks = inverse_ranks(&source);
            for op in 0..3 {
                // Section4.3 says consecutive item, not a circular neighbor.
                let a = rng.index(if op == 1 { ranks.len() - 1 } else { ranks.len() });
                let b = if op == 1 {
                    a + 1
                } else {
                    rng.index(ranks.len())
                };
                // The different-item restriction applies to interchange only.
                // Inserting past an identical endpoint can still move all
                // intervening different items (e.g. A,B,C,A -> B,C,A,A).
                if op != 2 && g.pieces[ids[ranks[a]]] == g.pieces[ids[ranks[b]]] {
                    continue;
                }
                let mut next = ranks.clone();
                if op == 2 {
                    let v = next.remove(a);
                    next.insert(b, v);
                } else {
                    next.swap(a, b);
                }
                let x: Vec<_> = next.iter().map(|&r| r as f64 + 1.).collect();
                let Some(layout) = evaluate(g, &ids, &x, budget, stats)? else {
                    complete = false;
                    break;
                };
                let incumbent = if index == particles.len() {
                    best.length
                } else {
                    particles[index].length
                };
                if layout.length < incumbent {
                    ranks = next;
                    stats.accepted_moves += 1;
                    if index < particles.len() {
                        particles[index].personal = x.clone();
                        particles[index].length = layout.length;
                    }
                    if layout.length < best.length {
                        best = layout;
                        global = x;
                    }
                }
            }
            if !complete {
                break;
            }
        }
        if !complete {
            break;
        }
        *stats.completed_generations.as_mut().unwrap() += 1;
        omega = (omega * 0.99).max(0.4);
    }
    stats.paper_run_status = Some(
        if particles.len() == 12 {
            "swarmnest_direct_polygon_budget_stop"
        } else {
            "swarmnest_initialization_interrupted"
        }
        .into(),
    );
    stats.backend_events.push(ASSUMPTIONS.into());
    Ok(best)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn extraction_is_inverse_ranks_not_sorted_indices() {
        assert_eq!(inverse_ranks(&[3.8, 0.2, 8.0, 2.4]), vec![2, 0, 3, 1]);
    }
    #[test]
    fn published_figure7_solution_extraction_anchor() {
        // Original Figure7 reports the one-based sequence4,2,6,3,1,5.
        assert_eq!(inverse_ranks(&[1.80,-0.99,3.01,-0.72,-1.20,2.15]),vec![3,1,5,2,0,4]);
    }
}
