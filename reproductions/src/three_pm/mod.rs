//! Three dedicated phases: relax-and-fix DBM, refined-grid VND, local facet MIP.
pub(crate) mod models;
use self::models::{Dot, Movement, Pose, Shape};
use crate::{
    backends::{Settings, cbc},
    geometry::{self, Geometry},
    model::*,
    toolkit::{
        placement::{Bounds, oriented},
        polygon::P,
    },
};
use anyhow::{Result, ensure};

struct Prepared {
    shapes: Vec<Shape>,
    offsets: Vec<P>,
    scale: f64,
    height: f64,
    demand: Vec<usize>,
}
fn prepare(g: &Geometry, budget: &Budget) -> Result<Option<Prepared>> {
    ensure!(g.pieces.len() <= 512, "A080 piece cap");
    let scale = 1.0;
    let mut shapes = Vec::new();
    let mut offsets = Vec::new();
    let mut demand = vec![0; g.types.len()];
    for &t in &g.pieces {
        demand[t] += 1;
    }
    for (t, kind) in g.types.iter().enumerate() {
        let Some(piece) = g.pieces.iter().position(|&v| v == t) else {
            continue;
        };
        for (o, a) in kind.orientations.iter().enumerate() {
            if budget.expired() {
                return Ok(None);
            }
            if a.height > g.height {
                continue;
            }
            let v = oriented(g, piece, o);
            let b = Bounds::of(&v);
            let mut vertices: Vec<P> = v
                .iter()
                .map(|p| P::new((p.x - b.x0) * scale, (p.y - b.y0) * scale))
                .collect();
            let reference=crate::config::get().three_pm.reference_vertex_index;
            ensure!(reference<vertices.len(),"reference vertex outside original contour");
            vertices.rotate_left(reference);
            shapes.push(Shape {
                kind: t,
                orientation: o,
                vertices,
                clearance: 0.0,
            });
            offsets.push(P::new(b.x0, b.y0));
        }
    }
    Ok(Some(Prepared {
        shapes,
        offsets,
        scale,
        height: g.height,
        demand,
    }))
}
fn to_dots(g: &Geometry, p: &Prepared, layout: &Layout) -> Option<Vec<Dot>> {
    layout
        .placements
        .iter()
        .map(|q| {
            let s = p
                .shapes
                .iter()
                .position(|s| s.kind == g.pieces[q.piece] && s.orientation == q.orientation)?;
            Some(Dot {
                shape: s,
                x: ((q.x + p.offsets[s].x) * p.scale).max(0.0),
                y: ((q.y + p.offsets[s].y) * p.scale).max(0.0),
            })
        })
        .collect()
}
fn export(g: &Geometry, p: &Prepared, dots: &[Dot]) -> Option<Layout> {
    let mut available: Vec<Vec<usize>> = (0..g.types.len())
        .map(|t| {
            g.pieces
                .iter()
                .enumerate()
                .filter(|(_, s)| **s == t)
                .map(|(i, _)| i)
                .collect()
        })
        .collect();
    let mut placements = Vec::new();
    for d in dots {
        let s = &p.shapes[d.shape];
        let piece = available[s.kind].pop()?;
        let x = d.x / p.scale - p.offsets[d.shape].x;
        let y = d.y / p.scale - p.offsets[d.shape].y;
        if x < -4.0 * g.gap || y < -4.0 * g.gap {
            return None;
        }
        placements.push(Placement {
            piece,
            orientation: s.orientation,
            // Recover a mathematically zero coordinate from negative solver ulps.
            // The original-contour validation gate checks the clamped candidate.
            x: x.max(0.0),
            y: y.max(0.0),
        });
    }
    if available.iter().any(|a| !a.is_empty()) {
        return None;
    }
    placements.sort_by_key(|v| v.piece);
    let layout = Layout {
        length: geometry::length(g, &placements),
        placements,
    };
    g.validate(&layout).ok().map(|_| layout)
}
fn selected(dots: &[Dot], s: &cbc::Solution) -> Vec<Dot> {
    dots.iter()
        .enumerate()
        .filter(|(i, _)| s.values[i + 1] > 0.5)
        .map(|(_, d)| *d)
        .collect()
}
/// Build a feasible seed from retained dots, adding only conflict-free candidates.
/// This prepares a MIP start; CBC still optimizes the full phase model.
fn mip_start(
    p: &Prepared,
    dots: &[Dot],
    model: &cbc::Model,
    demand: &[usize],
    previous: &[Dot],
    movement: Option<&Movement>,
    budget: &Budget,
) -> Option<Vec<f64>> {
    let mut values = vec![0.0; model.variables.len()];
    values[0] = model.variables[0].upper;
    let mut count = vec![0; demand.len()];
    let mut adjacent = vec![Vec::new(); dots.len()];
    for r in &model.rows {
        if budget.expired() {
            return None;
        }
        if matches!(r.sense, cbc::Sense::Less)
            && r.rhs == 1.0
            && r.terms.len() >= 2
            && r.terms
                .iter()
                .all(|(i, a)| *i > 0 && *a == 1.0 && model.variables[*i].binary)
        {
            for (k,&(a,_)) in r.terms.iter().enumerate(){for &(b,_) in &r.terms[k+1..]{
                adjacent[a-1].push(b-1);adjacent[b-1].push(a-1);
            }}
        }
    }
    for (i, d) in dots.iter().enumerate() {
        if model.variables[i + 1].lower == 1.0 || previous.iter().any(|v| same(v, d)) {
            values[i + 1] = 1.0;
            count[p.shapes[d.shape].kind] += 1;
        }
    }
    let mut ids: Vec<_> = (0..dots.len()).collect();
    ids.sort_by(|&a, &b| {
        dots[a]
            .x
            .total_cmp(&dots[b].x)
            .then(dots[a].y.total_cmp(&dots[b].y))
            .then(a.cmp(&b))
    });
    for i in ids {
        if budget.expired() {
            return None;
        }
        let t = p.shapes[dots[i].shape].kind;
        if count[t] >= demand[t]
            || values[i + 1] > 0.5
            || movement.is_some_and(|w| w.window[i])
            || adjacent[i].iter().any(|&j| values[j + 1] > 0.5)
        {
            continue;
        }
        values[i + 1] = 1.0;
        count[t] += 1;
    }
    (count == demand && model.feasible(&values)).then_some(values)
}
#[allow(clippy::too_many_arguments)]
fn solve_dbm(
    p: &Prepared,
    dots: &[Dot],
    upper: f64,
    demand: &[usize],
    fixed: &[usize],
    movement: Option<&Movement>,
    previous: &[Dot],
    cfg: &Settings,
    budget: &Budget,
    stats: &mut Stats,
) -> Result<Option<Vec<Dot>>> {
    let Some(model) = models::dotted(
        &p.shapes, dots, p.height, upper, demand, fixed, movement, budget,
    )?
    else {
        stats
            .backend_events
            .push("A080 DBM interrupted or model resource cap".into());
        return Ok(None);
    };
    stats.evaluations += model.rows.len();
    if budget.expired() {
        return Ok(None);
    }
    stats.backend_calls += 1;
    let start = mip_start(p, dots, &model, demand, previous, movement, budget);
    stats.backend_events.push(format!(
        "A080 DBM feasible MIP start submitted={}",
        start.is_some()
    ));
    let Some(s) = cbc::solve_with_start(&model, cfg, budget, start.as_deref())? else {
        return Ok(None);
    };
    stats.backend_events.push(format!(
        "A080 actual DBM {} solution optimal={} cols={} rows={} exact_epigraph_thresholds={} cpsat_probe_count={} native_dbm_nodes={} native_dbm_probes={}",
        cfg.solver_name(),
        s.optimal,
        model.variables.len(),
        model.rows.len(),
        s.exact_epigraph_thresholds,
        s.cpsat_probe_count,s.native_dbm_nodes,s.native_dbm_probes
    ));
    if !s.optimal { stats.nonoptimal_backend_calls += 1; }
    Ok(Some(selected(dots, &s)))
}
fn same(a: &Dot, b: &Dot) -> bool {
    a.shape == b.shape && a.x == b.x && a.y == b.y
}
/// Generate only types with positive partial demand, preserving global shape IDs.
/// Retained dots keep the previous subproblem solution available after rebuilding.
fn partial_grid(
    p: &Prepared,
    upper: f64,
    demand: &[usize],
    old: &[Dot],
    budget: &Budget,
    stats: &mut Stats,
) -> Result<Option<Vec<Dot>>> {
    let ids: Vec<_> = p
        .shapes
        .iter()
        .enumerate()
        .filter(|(_, s)| demand[s.kind] > 0)
        .map(|(i, _)| i)
        .collect();
    let shapes: Vec<_> = ids.iter().map(|&i| p.shapes[i].clone()).collect();
    let retained: Vec<_> = old
        .iter()
        .map(|d| Dot {
            shape: ids.iter().position(|&i| i == d.shape).unwrap(),
            x: d.x,
            y: d.y,
        })
        .collect();
    let Some(mut dots) = models::grid(&shapes, p.height, upper, crate::config::get().three_pm.initial_spacing, &retained, budget)? else {
        return Ok(None);
    };
    stats
        .backend_events
        .push(format!("3PM original gmin={} dots={}", crate::config::get().three_pm.initial_spacing, dots.len()));
    for d in &mut dots {
        d.shape = ids[d.shape];
    }
    Ok(Some(dots))
}
/// Source's printed largest Theta/qt selection (despite its homogeneity claim).
fn add_batch(full: &[usize], partial: &mut [usize], size: usize) {
    let mut batch = vec![0; full.len()];
    for _ in 0..size {
        let t = (0..full.len())
            .filter(|&t| partial[t] < full[t])
            .max_by(|&a, &b| {
                (batch[a] as f64 / full[a] as f64)
                    .total_cmp(&(batch[b] as f64 / full[b] as f64))
                    .then(b.cmp(&a))
            });
        let Some(t) = t else {
            break;
        };
        partial[t] += 1;
        batch[t] += 1;
    }
}
/// A feasible grid-aligned append bound only limits L; it never supplies the
/// algorithm's result. New rectangles start strictly beyond the previous W.
fn append_upper(p:&Prepared,old:&[usize],new:&[usize],mu:f64)->f64{
    let mut x=mu;let mut y=0.0;let mut column_end=mu;
    let mut pieces=Vec::new();
    for t in 0..new.len(){for _ in old[t]..new[t]{
        let s=p.shapes.iter().filter(|s|s.kind==t).min_by(|a,b|Bounds::of(&a.vertices).width().total_cmp(&Bounds::of(&b.vertices).width())).unwrap();
        pieces.push(s);
    }}
    pieces.sort_by(|a,b|Bounds::of(&b.vertices).height().total_cmp(&Bounds::of(&a.vertices).height()));
    for s in pieces{
        let dx=models::spacing_axis(&s.vertices,true,2.0);let dy=models::spacing_axis(&s.vertices,false,2.0);
        let b=Bounds::of(&s.vertices);
        let reference=s.vertices[0];
        let yy=((y+reference.y)/dy).ceil()*dy-reference.y;
        if yy+b.height()>p.height{ x=column_end;y=0.;}
        let xx=if mu>0.0 && x<=mu{(((mu+reference.x)/dx).floor()+1.)*dx-reference.x}else{((x+reference.x)/dx).ceil()*dx-reference.x};
        let yy=((y+reference.y)/dy).ceil()*dy-reference.y;
        column_end=column_end.max(xx+b.width());y=yy+b.height();
    }
    column_end
}
fn constructive(
    p: &Prepared,
    upper: f64,
    cfg: &Settings,
    budget: &Budget,
    stats: &mut Stats,
) -> Result<Option<Vec<Dot>>> {
    let mut partial = vec![0; p.demand.len()];
    let mut previous = Vec::new();
    let mut mu0 = 0.0;
    let mut mu = 0.0;
    let source_sigma = models::batch_size(p.demand.iter().sum());
    let sigma = source_sigma;
    stats
        .backend_events
        .push(format!("3PM constructive sigma={sigma}"));
    while partial != p.demand && !budget.expired() {
        let before=partial.clone();
        add_batch(&p.demand, &mut partial, sigma);
        // Engineering candidate-domain bound for a growing partial instance.
        // Two grid units per piece reserve coarse-grid rounding room; the fixed
        // overall upper bound remains authoritative in the mathematical model.
        let column_bound: f64 = partial
            .iter()
            .enumerate()
            .map(|(t, &q)| {
                let width = p
                    .shapes
                    .iter()
                    .filter(|s| s.kind == t)
                    .map(|s| Bounds::of(&s.vertices).x1 + s.clearance)
                    .fold(f64::INFINITY, f64::min);
                if q == 0 {
                    0.0
                } else {
                    q as f64 * (width + 2.0)
                }
            })
            .sum();
        let grid_upper = upper.min(column_bound.max(mu + 2.0)).min(append_upper(p,&before,&partial,mu));
        let Some(dots) = partial_grid(p, grid_upper, &partial, &previous, budget, stats)? else {
            stats
                .backend_events
                .push("A080 partial-demand constructive grid cap/budget stop".into());
            return Ok(None);
        };
        let mut fixed = Vec::new();
        let mut alpha = vec![0; p.demand.len()];
        let window: Vec<_> = dots.iter().map(|d| {let rx=d.x+p.shapes[d.shape].vertices[0].x;mu>mu0 && rx>=mu0 && rx<=mu}).collect();
        let prev: Vec<_> = dots
            .iter()
            .map(|d| previous.iter().any(|v| same(d, v)))
            .collect();
        for (i, d) in dots.iter().enumerate() {
            if prev[i] {
                if d.x+p.shapes[d.shape].vertices[0].x < mu0 {
                    fixed.push(i);
                } else if window[i] {
                    alpha[p.shapes[d.shape].kind] += 1;
                }
            }
        }
        let movement = Movement {
            window,
            previous: prev,
            alpha,
        };
        let Some(next) = solve_dbm(
            p,
            &dots,
            upper,
            &partial,
            &fixed,
            Some(&movement),
            &previous,
            cfg,
            budget,
            stats,
        )?
        else {
            return Ok(None);
        };
        previous = next;
        mu0 = mu;
        mu = previous
            .iter()
            .map(|d| d.x + Bounds::of(&p.shapes[d.shape].vertices).x1)
            .fold(0.0, f64::max);
        stats.backend_events.push(format!("3PM partial={} L={mu}; grid upper={grid_upper}",partial.iter().sum::<usize>()));
        stats.completed_generations = Some(stats.completed_generations.unwrap_or(0) + 1);
    }
    Ok((partial == p.demand).then_some(previous))
}
/// Restrict candidate dots to the union of each same-type incumbent neighborhood.
/// The undefined Fig.6 s is interpreted as gmin; side4*gmin is disclosed.
fn neighborhood(p: &Prepared, grid: &[Dot], old: &[Dot], kind: usize, step: f64) -> Vec<Dot> {
    let mut chosen = vec![false; grid.len()];
    for center in old {
        let t = p.shapes[center.shape].kind;
        let cr=p.shapes[center.shape].vertices[0];
        let delta=|d:&Dot| {let r=p.shapes[d.shape].vertices[0];((d.x+r.x-center.x-cr.x).abs(),(d.y+r.y-center.y-cr.y).abs())};
        let mut ids: Vec<_> = grid
            .iter()
            .enumerate()
            .filter(|(_, d)| p.shapes[d.shape].kind == t)
            .map(|(i, _)| i)
            .collect();
        match kind {
            1 => ids.retain(|&i| {
                delta(&grid[i]).0 <= crate::config::get().three_pm.first_square_half_side_factor * step
                    && delta(&grid[i]).1 <= crate::config::get().three_pm.first_square_half_side_factor * step
            }),
            2 | 3 => {
                let distance = |i: usize| {
                    if kind == 2 {
                        delta(&grid[i]).0
                    } else {
                        delta(&grid[i]).1
                    }
                };
                ids.sort_by(|&a, &b| distance(a).total_cmp(&distance(b)).then(a.cmp(&b)));
                // Section3.2/Fig6: md bounds EACH reference-point rectangle.
                // Exclude a complete distance boundary instead of cutting a
                // row/column into a nonrectangular subset.
                let md=crate::config::get().three_pm.maximum_dots_per_reference;
                if ids.len()>md {
                    let boundary=distance(ids[md]);
                    ids.retain(|&i|distance(i)<boundary);
                }
            }
            _ => unreachable!(),
        }
        for i in ids {
            chosen[i] = true;
        }
    }
    let mut dots: Vec<_> = grid
        .iter()
        .enumerate()
        .filter(|(i, _)| chosen[*i])
        .map(|(_, d)| *d)
        .collect();
    for &d in old {
        if !dots.iter().any(|v| same(&d, v)) {
            dots.push(d);
        }
    }
    dots
}
#[cfg(test)]mod neighborhood_tests{
    use super::*;
    #[test]fn md_limits_each_reference_region_not_the_union(){
        let mut shapes=Vec::new();for kind in 0..2{shapes.push(Shape{kind,orientation:0,clearance:0.,vertices:vec![P::new(0.,0.),P::new(1.,0.),P::new(1.,1.),P::new(0.,1.)]});}
        let p=Prepared{shapes,offsets:vec![P::new(0.,0.);2],scale:1.,height:100.,demand:vec![1,1]};
        let old=vec![Dot{shape:0,x:0.,y:0.},Dot{shape:1,x:1.,y:1.}];let mut grid=Vec::new();
        for shape in 0..2{for x in 0..100{for y in 0..60{grid.push(Dot{shape,x:x as f64,y:y as f64});}}}
        for k in [2,3]{let dots=neighborhood(&p,&grid,&old,k,1.);assert!(dots.len()>3000);assert!(old.iter().all(|o|dots.iter().any(|d|same(o,d))));for shape in 0..2{assert!(dots.iter().filter(|d|d.shape==shape).count()<=3000);}}
    }
    #[test]fn md_preserves_a_full_rectangle_distance_boundary(){
        let p=Prepared{shapes:vec![Shape{kind:0,orientation:0,clearance:0.,vertices:vec![P::new(0.,0.),P::new(1.,0.),P::new(0.,1.)]}],offsets:vec![P::new(0.,0.)],scale:1.,height:2000.,demand:vec![1]};
        let old=vec![Dot{shape:0,x:0.,y:0.}];let mut grid=Vec::new();
        for x in 0..2{for y in 0..2000{grid.push(Dot{shape:0,x:x as f64,y:y as f64});}}
        let dots=neighborhood(&p,&grid,&old,2,1.);assert_eq!(dots.len(),2000);assert!(dots.iter().all(|d|d.x==0.));
    }
}
fn improvement(
    g: &Geometry,
    p: &Prepared,
    best: &mut Layout,
    cfg: &Settings,
    budget: &Budget,
    stats: &mut Stats,
) -> Result<bool> {
    let mut step = crate::config::get().three_pm.initial_spacing;
    let mut all_solved = true;
    // Only this invocation's optimal, original-geometry-validated DBMs can
    // certify a repeated domain. No persisted result or interrupted model is
    // reused. A different upper bound is deliberately a different key.
    let mut certificates = VndCertificates::default();
    while step >= crate::config::get().three_pm.minimum_spacing {
        if budget.expired() {
            return Ok(false);
        }
        let Some(mut old) = to_dots(g, p, best) else {
            return Ok(false);
        };
        let upper = ((best.length + 4.0 * g.gap) * p.scale).min(g.upper * p.scale);
        let actual_step = step;
        let Some(grid) = models::grid(&p.shapes, p.height, upper, step, &old, budget)? else {
            return Ok(false);
        };
        stats
            .backend_events
            .push(format!("3PM VND gmin={step} dots={}", grid.len()));
        let mut k = 1;
        while k <= 3 {
            if budget.expired() {
                return Ok(false);
            }
            let dots = neighborhood(p, &grid, &old, k, actual_step);
            stats.backend_events.push(format!(
                "A080 VND neighborhood{k} grid_min={step} dots={}",
                dots.len()
            ));
            if certificates.no_improvement(&dots,upper,best.length,1e-8*g.height) {
                stats.backend_events.push(format!("A080 VND neighborhood{k} reused same-run optimal-domain certificate; no improving neighbour"));
                k+=1;
                continue;
            }
            let nonoptimal_before=stats.nonoptimal_backend_calls;
            let next = solve_dbm(
                p,
                &dots,
                upper,
                &p.demand,
                &[],
                None,
                &old,
                cfg,
                budget,
                stats,
            )?;
            all_solved &= next.is_some() && stats.nonoptimal_backend_calls==nonoptimal_before;
            let validated=next.as_ref().and_then(|d|export(g,p,d));
            certificates.remember(&dots,upper,validated.as_ref().map(|l|l.length),
                next.is_some() && stats.nonoptimal_backend_calls==nonoptimal_before);
            if let Some(layout) = validated && layout.length < best.length - 1e-8 * g.height {
                *best = layout;
                let Some(next_old) = to_dots(g, p, best) else {
                    return Ok(false);
                };
                old = next_old;
                stats.accepted_moves += 1;
                k = 1;
                continue;
            }
            k += 1;
        }
        step *= crate::config::get().three_pm.refinement;
    }
    Ok(all_solved)
}

#[derive(Default)]
struct VndCertificates {
    entries: std::collections::BTreeMap<(u64,Vec<(usize,u64,u64)>),f64>,
}
impl VndCertificates {
    fn key(dots:&[Dot],upper:f64)->(u64,Vec<(usize,u64,u64)>) {
        let bits=|x:f64|if x==0. {0}else{x.to_bits()};
        let mut positions:Vec<_>=dots.iter().map(|d|(d.shape,bits(d.x),bits(d.y))).collect();
        positions.sort_unstable();
        // Keep duplicates: changing multiplicity must not certify a new model.
        (bits(upper),positions)
    }
    fn remember(&mut self,dots:&[Dot],upper:f64,validated_length:Option<f64>,optimal:bool) {
        if optimal && let Some(length)=validated_length && length.is_finite() {
            // Bound memory without changing the search: forgetting only causes
            // a later model to be solved again.
            if self.entries.len()>=16 { self.entries.clear(); }
            self.entries.insert(Self::key(dots,upper),length);
        }
    }
    fn no_improvement(&self,dots:&[Dot],upper:f64,best:f64,tolerance:f64)->bool {
        let (bound,domain)=Self::key(dots,upper);
        self.entries.iter().any(|((old_bound,old_domain),&length)|{
            if length<best-tolerance || f64::from_bits(*old_bound)<f64::from_bits(bound){return false;}
            // With unchanged VND demands/geometry, a subset of a certified
            // domain and a tighter length upper bound cannot improve its lower
            // bound. Compare multisets so duplicate columns remain meaningful.
            let mut i=0;
            for candidate in &domain{
                while i<old_domain.len() && old_domain[i]<*candidate{i+=1;}
                if i==old_domain.len() || old_domain[i]!=*candidate{return false;}
                i+=1;
            }
            true
        })
    }
}

#[cfg(test)]mod certificate_tests {
    use super::*;
    #[test]fn certificates_require_optimality_and_original_geometry_validation(){
        let d=[Dot{shape:0,x:0.,y:1.},Dot{shape:1,x:2.,y:3.}];
        let mut c=VndCertificates::default();
        c.remember(&d,10.,Some(7.),false);assert!(!c.no_improvement(&d,10.,7.,1e-8));
        c.remember(&d,10.,None,true);assert!(!c.no_improvement(&d,10.,7.,1e-8));
        c.remember(&d,10.,Some(7.),true);assert!(c.no_improvement(&d,10.,7.,1e-8));
        assert!(!c.no_improvement(&d,10.,8.,1e-8)); // A certified improving optimum cannot be skipped.
        assert!(c.no_improvement(&d,10.,6.,1e-8));
    }
    #[test]fn certificates_accept_subsets_tighter_bounds_and_column_reordering_but_not_new_candidates(){
        let a=Dot{shape:0,x:0.,y:1.};let b=Dot{shape:1,x:2.,y:3.};
        let mut c=VndCertificates::default();c.remember(&[a,b],10.,Some(7.),true);
        assert!(c.no_improvement(&[b,Dot{x:-0.,..a}],10.,7.,1e-8));
        assert!(c.no_improvement(&[a,b],9.,7.,1e-8));
        assert!(c.no_improvement(&[a],9.,7.,1e-8));
        assert!(!c.no_improvement(&[a,b],11.,7.,1e-8));
        assert!(!c.no_improvement(&[a,b,Dot{x:2.5,..b}],10.,7.,1e-8));
        assert!(!c.no_improvement(&[a,b,b],10.,7.,1e-8));
        assert!(!c.no_improvement(&[a,Dot{shape:0,..b}],10.,7.,1e-8));
    }
}
fn local_lambda(shape:&Shape,n:usize)->f64{
    let b=Bounds::of(&shape.vertices);
    if n<13{b.width().max(b.height())}else if n<=20{0.5*b.width().max(b.height())}else{0.5*b.width().min(b.height())}
}
fn compaction_start(p:&Prepared,poses:&[Pose],m:&cbc::Model)->Option<Vec<f64>>{
    let base=1+2*poses.len();let mut values=vec![0.;m.variables.len()];values[0]=m.variables[0].lower;
    for (i,pose) in poses.iter().enumerate(){
        values[1+2*i]=pose.x;values[2+2*i]=pose.y;
        values[0]=values[0].max(pose.x+Bounds::of(&p.shapes[pose.shape].vertices).x1);
    }
    let mut allowed=vec![true;m.variables.len()];
    for row in &m.rows{
        let binary:Vec<_>=row.terms.iter().filter(|&&(i,_)|i>=base).copied().collect();
        if binary.len()==1 && row.terms.iter().any(|&(i,_)|i<base){
            let (z,a)=binary[0];let lhs=row.terms.iter().filter(|&&(i,_)|i<base).map(|&(i,a)|a*values[i]).sum::<f64>()+a;
            allowed[z]&=match row.sense{cbc::Sense::Less=>lhs<=row.rhs+1e-9,cbc::Sense::Greater=>lhs>=row.rhs-1e-9,cbc::Sense::Equal=>(lhs-row.rhs).abs()<=1e-9};
        }
    }
    for row in &m.rows{
        if row.rhs==1. && row.terms.iter().all(|&(i,a)|i>=base && a==1.){
            let &(z,_)=row.terms.iter().find(|&&(i,_)|allowed[i])?;values[z]=1.;
        }
    }
    m.feasible(&values).then_some(values)
}
#[cfg(test)]mod compaction_seed_tests{
    use super::*;
    #[test]fn feasible_start_uses_one_true_facet_and_rejects_an_overlapping_seed(){
        let shape=Shape{kind:0,orientation:0,clearance:0.,vertices:vec![P::new(0.,0.),P::new(1.,0.),P::new(1.,1.),P::new(0.,1.)]};
        let p=Prepared{shapes:vec![shape],offsets:vec![P::new(0.,0.)],scale:1.,height:1.,demand:vec![2]};
        let poses=[Pose{shape:0,x:0.,y:0.},Pose{shape:0,x:2.,y:0.}];
        let model=models::compact(&p.shapes,&poses,1.,3.,&[3.;2],&Budget::new(2.)).unwrap().unwrap();
        let start=compaction_start(&p,&poses,&model).unwrap();assert!(model.feasible(&start));assert_eq!(start[0],3.);
        let mut cfg=Settings::default();cfg.cbc_preprocess=false;cfg.cbc_presolve=false;
        let s=cbc::solve_with_start(&model,&cfg,&Budget::new(5.),Some(&start)).unwrap().unwrap();
        assert!(s.optimal && (s.objective-2.).abs()<1e-8);
        let overlap=[poses[0],Pose{shape:0,x:0.5,y:0.}];assert!(compaction_start(&p,&overlap,&model).is_none());
    }
}
fn compaction(
    g: &Geometry,
    p: &Prepared,
    best: &mut Layout,
    cfg: &Settings,
    budget: &Budget,
    stats: &mut Stats,
) -> Result<bool> {
    while !budget.expired() {
        let Some(old) = to_dots(g, p, best) else {
            return Ok(false);
        };
        let poses: Vec<_> = old
            .iter()
            .map(|d| Pose {
                shape: d.shape,
                x: d.x,
                y: d.y,
            })
            .collect();
        let lambda: Vec<_> = old
            .iter()
            .map(|d| {
                local_lambda(&p.shapes[d.shape],old.len())
            })
            .collect();
        let upper = (best.length + 4.0 * g.gap) * p.scale;
        let Some(model) = models::compact(&p.shapes, &poses, p.height, upper, &lambda, budget)?
        else {
            stats
                .backend_events
                .push("A080 facet MIP resource/budget stop".into());
            return Ok(false);
        };
        if budget.expired() {
            return Ok(false);
        }
        stats.backend_calls += 1;
        stats.evaluations += model.rows.len();
        stats.backend_events.push(format!(
            "A080 actual facet-disjunctive local MIP cols={} rows={}",
            model.variables.len(),
            model.rows.len()
        ));
        let start=compaction_start(p,&poses,&model);
        stats.backend_events.push(format!("A080 continuous feasible incumbent submitted={}",start.is_some()));
        let Some(s) = cbc::solve_with_start(&model, cfg, budget,start.as_deref())? else {
            return Ok(false);
        };
        if !s.optimal { stats.nonoptimal_backend_calls += 1; }
        stats.backend_events.push(format!("A080 facet {} solution optimal={}",cfg.solver_name(),s.optimal));
        let dots: Vec<_> = old
            .iter()
            .enumerate()
            .map(|(i, d)| Dot {
                shape: d.shape,
                x: s.values[1 + 2 * i],
                y: s.values[2 + 2 * i],
            })
            .collect();
        let Some(layout) = export(g, p, &dots) else {
            stats
                .backend_events
                .push("A080 facet MIP original-geometry rejection; incumbent retained".into());
            return Ok(false);
        };
        if layout.length >= best.length - 1e-8 * g.height {
            return Ok(s.optimal);
        }
        *best = layout;
        stats.accepted_moves += 1;
    }
    Ok(false)
}
pub fn solve(g: &Geometry, budget: &Budget, cfg: &Settings, stats: &mut Stats) -> Result<Layout> {
    stats.effective_algorithm_controls=Some(serde_json::to_value(&crate::config::get().three_pm)?);
    let mut best = g.shelf(&geometry::initial_order(g, OrderRule::Area));
    stats.paper_variant = Some("three_phase_DBM_VND_convex_component_facet_MIP".into());
    let prepared = match prepare(g, budget) {
        Ok(Some(p)) => p,
        other => {
            stats.backend_events.push(format!(
                "A080 preparation stopped: {}",
                match other {
                    Err(e) => e.to_string(),
                    _ => "deadline".into(),
                }
            ));
            stats.paper_run_status = Some("preparation_stop_shelf_retained".into());
            return Ok(best);
        }
    };
    let constructed = constructive(
        &prepared,
        (best.length + 4.0 * g.gap) * prepared.scale,
        cfg,
        budget,
        stats,
    )
    .unwrap_or_else(|e| {
        stats
            .backend_events
            .push(format!("A080 constructive model error: {e}"));
        None
    });
    let construction_ok = constructed
        .as_ref()
        .and_then(|d| export(g, &prepared, d))
        .is_some();
    if let Some(layout) = constructed.as_ref().and_then(|d| export(g, &prepared, d))
        && layout.length < best.length
    {
        best = layout;
        stats.accepted_moves += 1;
    }
    if construction_ok {
        stats.completed_variants += 1;
    } else {
        stats.paper_run_status = Some("construction_incomplete_no_three_phase_result".into());
        return Ok(best);
    }
    let improved = improvement(g, &prepared, &mut best, cfg, budget, stats).unwrap_or_else(|e| {
        stats
            .backend_events
            .push(format!("A080 VND model error: {e}"));
        false
    });
    if improved {
        stats.completed_variants += 1;
    }
    let compacted = compaction(g, &prepared, &mut best, cfg, budget, stats).unwrap_or_else(|e| {
        stats
            .backend_events
            .push(format!("A080 compaction model error: {e}"));
        false
    });
    if compacted {
        stats.completed_variants += 1;
    }
    stats.paper_run_status = Some(
        if construction_ok && improved && compacted && stats.nonoptimal_backend_calls==0 {
            "completed_three_phase_source_described_flow"
        } else {
            "partial_three_phase_flow_validated_incumbent"
        }
        .into(),
    );
    Ok(best)
}

/// Recover a live full-grid DBM incumbent only after byte-identical model replay.
/// This does not run or complete the3PM controller and cannot certify optimality.
pub fn validate_full_grid_checkpoint(
    g:&Geometry,step:f64,upper:f64,original_mps:&[u8],values:&[f64],budget:&Budget,
)->Result<Layout>{
    let p=prepare(g,budget)?.ok_or_else(||anyhow::anyhow!("checkpoint preparation expired"))?;
    let dots=models::grid(&p.shapes,p.height,upper,step,&[],budget)?
        .ok_or_else(||anyhow::anyhow!("checkpoint grid generation expired"))?;
    let model=models::dotted(&p.shapes,&dots,p.height,upper,&p.demand,&[],None,budget)?
        .ok_or_else(||anyhow::anyhow!("checkpoint model generation expired"))?;
    ensure!(model.mps()?.as_bytes()==original_mps,"checkpoint column mapping rejected: MPS differs from live solver model");
    ensure!(model.feasible(values),"checkpoint violates original submitted model");
    let chosen:Vec<_>=dots.iter().enumerate().filter(|(i,_)|values[i+1]>0.5).map(|(_,d)|*d).collect();
    export(g,&p,&chosen).ok_or_else(||anyhow::anyhow!("checkpoint original-contour export failed"))
}

/// Offline single-model diagnostic: preserve every original continuous bound.
/// It neither runs nor completes the3PM controller, and never touches its job.
pub fn audit_compaction_from_bounds(
    g:&Geometry,bounds:&[[f64;2]],cfg:&Settings,budget:&Budget,model_path:&std::path::Path,
)->Result<Option<(Layout,bool,usize,usize)>>{
    let p=prepare(g,budget)?.ok_or_else(||anyhow::anyhow!("compaction audit preparation expired"))?;
    let n=g.pieces.len();ensure!(bounds.len()==2*n+1 && bounds[0][0]==0.,"original continuous-column count/bound differs");
    ensure!(p.demand.iter().all(|&q|q==1) && p.shapes.len()==n,"audit requires unit demand and one orientation");
    let upper=bounds[0][1];let mut poses=Vec::new();let mut lambda=Vec::new();
    for (i,&kind) in g.pieces.iter().enumerate(){
        let shape=p.shapes.iter().position(|s|s.kind==kind).unwrap();let b=Bounds::of(&p.shapes[shape].vertices);
        let lam=local_lambda(&p.shapes[shape],n);
        let recover=|limits:[f64;2],global:f64|->Result<f64>{
            let [lo,hi]=limits;ensure!(lo>=0. && hi>=lo && hi<=global,"original bound outside board: {limits:?}, global={global:.17}, lambda={lam:.17}");
            let center=if lo>0.{lo+lam}else if hi<global{hi-lam}else{
                let low=(global-lam).max(0.);let high=lam.min(global);ensure!(low<=high,"inconsistent clipped window");(low+high)*0.5
            };
            ensure!(center>=0. && center<=global && (center-lam).max(0.)==lo && (center+lam).min(global)==hi,"recovered window is not byte-exact");
            Ok(center)
        };
        poses.push(Pose{shape,x:recover(bounds[1+2*i],upper-b.width())?,y:recover(bounds[2+2*i],p.height-b.height())?});lambda.push(lam);
    }
    // Recover any grid incumbent which induces these identical clipped windows.
    // This only supplies the already-known feasible start to the diagnostic.
    // The production controller has its real incumbent and does not do this.
    let grid=models::grid(&p.shapes,p.height,upper,1.,&[],budget)?.ok_or_else(||anyhow::anyhow!("audit seed grid expired"))?;
    let restricted:Vec<_>=grid.into_iter().filter(|d|{
        let kind=p.shapes[d.shape].kind;let i=g.pieces.iter().position(|&k|k==kind).unwrap();let b=Bounds::of(&p.shapes[d.shape].vertices);let lam=lambda[i];
        [(d.x-lam).max(0.),(d.x+lam).min(upper-b.width())]==bounds[1+2*i]
            && [(d.y-lam).max(0.),(d.y+lam).min(p.height-b.height())]==bounds[2+2*i]
    }).collect();
    let seed_model=models::dotted(&p.shapes,&restricted,p.height,upper,&p.demand,&[],None,budget)?.ok_or_else(||anyhow::anyhow!("audit seed model expired"))?;
    let mut seed_cfg=Settings::default();seed_cfg.finite_dbm=true;
    let seed=cbc::solve(&seed_model,&seed_cfg,budget)?.ok_or_else(||anyhow::anyhow!("no feasible incumbent matching original windows"))?;
    let selected=selected(&restricted,&seed);
    poses=g.pieces.iter().map(|&kind|{let d=selected.iter().find(|d|p.shapes[d.shape].kind==kind).unwrap();Pose{shape:d.shape,x:d.x,y:d.y}}).collect();
    let Some(m)=models::compact(&p.shapes,&poses,p.height,upper,&lambda,budget)?else{return Ok(None);};
    for (v,b) in m.variables.iter().zip(bounds){ensure!([v.lower,v.upper]==*b,"original continuous bounds changed");}
    ensure!(!model_path.exists(),"refuse to replace previous audit model");std::fs::write(model_path,m.mps()?)?;
    let start=compaction_start(&p,&poses,&m).ok_or_else(||anyhow::anyhow!("recovered compaction start rejected"))?;
    let Some(s)=cbc::solve_with_start(&m,cfg,budget,Some(&start))?else{return Ok(None);};
    let dots:Vec<_>=poses.iter().enumerate().map(|(i,p)|Dot{shape:p.shape,x:s.values[1+2*i],y:s.values[2+2*i]}).collect();
    let layout=export(g,&p,&dots).ok_or_else(||anyhow::anyhow!("compaction audit original geometry rejected"))?;
    Ok(Some((layout,s.optimal,m.variables.len(),m.rows.len())))
}
