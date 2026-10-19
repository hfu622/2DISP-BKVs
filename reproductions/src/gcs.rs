//! Elkeran (2013), the published zero-clustering Dighe configuration.
//! The implementation deliberately rejects other instances: pair clustering is
//! not implemented here and must never be silently replaced by a generic search.
use crate::{
    geometry::{self, Geometry},
    model::*,
    toolkit::{
        nfp::{Cache, NoFit},
        penetration,
        placement::*,
        polygon::*,
    },
};
use anyhow::{Result, ensure};
use std::sync::Arc;

pub const ASSUMPTIONS: &str = "Published zero-clustering Dighe branch only, not a general positive-clustering GCS implementation. Resolved author/default and project choices are embedded in reproduction_configuration, with alternate ambiguous controls available through solve_with_config. Original-contour and Jagua gates plus independent export verification; no generic postcompactor. One trial does not complete ten-trial source statistics; source input mapping remains separately conditional.";

struct Engine<'a> {
    g: &'a Geometry,
    cache: Cache,
    cfg: &'a crate::config::Gcs,
}
impl<'a> Engine<'a> {
    fn pair(
        &mut self,
        a: &Placement,
        b: &Placement,
        budget: &Budget,
    ) -> Result<Option<Arc<NoFit>>> {
        self.cache.get(
            (
                self.g.pieces[a.piece],
                a.orientation,
                self.g.pieces[b.piece],
                b.orientation,
            ),
            &oriented(self.g, a.piece, a.orientation),
            &oriented(self.g, b.piece, b.orientation),
            budget,
        )
    }
    fn depth(&mut self, a: &Placement, b: &Placement, budget: &Budget) -> Result<f64> {
        let Some(nfp) = self.pair(a, b, budget)? else {
            return Ok(f64::INFINITY);
        };
        Ok(penetration::squared(&nfp, P::new(b.x - a.x, b.y - a.y))
            .0
            .sqrt())
    }
    fn cost(
        &mut self,
        p: &Placement,
        layout: &[Placement],
        weights: &[Vec<f64>],
        budget: &Budget,
    ) -> Result<f64> {
        let mut value = 0.;
        for q in layout {
            if q.piece != p.piece {
                value += weights[p.piece][q.piece] * self.depth(q, p, budget)?;
            }
        }
        Ok(value)
    }
    fn losses(
        &mut self,
        layout: &[Placement],
        budget: &Budget,
    ) -> Result<Vec<(usize, usize, f64)>> {
        let mut out = Vec::new();
        for i in 0..layout.len() {
            for j in 0..i {
                out.push((
                    layout[i].piece,
                    layout[j].piece,
                    self.depth(&layout[j], &layout[i], budget)?,
                ));
            }
        }
        Ok(out)
    }
    fn inside(&self, p: &Placement, length: f64) -> bool {
        let b = Bounds::of(&oriented(self.g, p.piece, p.orientation));
        let eps = self.g.height * 1e-12;
        p.x + b.x0 >= -eps
            && p.y + b.y0 >= -eps
            && p.x + b.x1 <= length + eps
            && p.y + b.y1 <= self.g.height + eps
    }
    /// Section5 / Eq9: deterministic minimum-penetration placement of a
    /// protruding polygon on the complete NFP/IFR contact arrangement.
    fn relocate(
        &mut self,
        index: usize,
        layout: &[Placement],
        length: f64,
        budget: &Budget,
        stats: &mut Stats,
    ) -> Result<Option<Placement>> {
        let old = &layout[index];
        let piece = old.piece;
        let o = old.orientation;
        let b = Bounds::of(&oriented(self.g, piece, o));
        let ifr = Bounds {
            x0: (-b.x0).max(0.),
            y0: (-b.y0).max(0.),
            x1: length - b.x1,
            y1: self.g.height - b.y1,
        };
        if ifr.x1 < ifr.x0 || ifr.y1 < ifr.y0 {
            return Ok(None);
        }
        let mut groups = Vec::new();
        let mut extra = Vec::new();
        for q in layout {
            if q.piece == piece {
                continue;
            }
            let Some(n) = self.pair(q, old, budget)? else {
                return Ok(None);
            };
            let shift = P::new(q.x, q.y);
            groups.push(translated_boundary(&n, shift));
            extra.extend(
                n.zero_strata
                    .iter()
                    .flat_map(|e| [e.a.shift(shift), e.b.shift(shift)]),
            );
        }
        let Some(points) = candidates(&groups, ifr, &extra, &[], &[], budget)? else {
            return Ok(None);
        };
        let weights = vec![vec![1.; layout.len()]; layout.len()];
        let mut best: Option<(Placement, f64)> = None;
        for point in points {
            if budget.expired() {
                return Ok(None);
            }
            let p = Placement {
                piece,
                orientation: o,
                x: point.x.clamp(ifr.x0, ifr.x1),
                y: point.y.clamp(ifr.y0, ifr.y1),
            };
            let value = self.cost(&p, layout, &weights, budget)?;
            stats.evaluations += 1;
            if best.as_ref().is_none_or(|(_, v)| value < *v) {
                best = Some((p, value));
            }
        }
        Ok(best.map(|(p, _)| p))
    }
    fn initial(
        &mut self,
        order: &[usize],
        budget: &Budget,
        stats: &mut Stats,
    ) -> Result<Option<Layout>> {
        let mut placed = Vec::new();
        for &piece in order {
            let mut choice: Option<Placement> = None;
            for o in 0..self.g.types[self.g.pieces[piece]].orientations.len() {
                if budget.expired() {
                    return Ok(None);
                }
                let moving = oriented(self.g, piece, o);
                let bb = Bounds::of(&moving);
                let ifr = Bounds {
                    x0: (-bb.x0).max(0.),
                    y0: (-bb.y0).max(0.),
                    x1: self.g.upper - bb.x1,
                    y1: self.g.height - bb.y1,
                };
                if ifr.y1 < ifr.y0 {
                    continue;
                }
                let probe = Placement {
                    piece,
                    orientation: o,
                    x: 0.,
                    y: 0.,
                };
                let mut groups = Vec::new();
                let mut extra = Vec::new();
                for q in &placed {
                    let Some(nfp) = self.pair(q, &probe, budget)? else {
                        return Ok(None);
                    };
                    let shift = P::new(q.x, q.y);
                    groups.push(translated_boundary(&nfp, shift));
                    extra.extend(nfp.contacts.iter().map(|p| p.shift(shift)));
                }
                let Some(points) = candidates(&groups, ifr, &extra, &[], &[], budget)? else {
                    return Ok(None);
                };
                let oracle = FixedPlacementOracle::new(self.g, &placed)?;
                for point in points {
                    if budget.expired() {
                        return Ok(None);
                    }
                    let p = Placement {
                        piece,
                        orientation: o,
                        x: point.x.max(ifr.x0).min(ifr.x1),
                        y: point.y.max(ifr.y0).min(ifr.y1),
                    };
                    stats.collision_queries += 1;
                    if oracle.clear(&p) {
                        if choice
                            .as_ref()
                            .is_none_or(|q| p.x.total_cmp(&q.x).then(p.y.total_cmp(&q.y)).is_lt())
                        {
                            choice = Some(p)
                        }
                        break;
                    }
                }
            }
            let Some(p) = choice else { return Ok(None) };
            placed.push(p);
        }
        stats.evaluations += 1;
        Ok(Some(Layout {
            length: geometry::length(self.g, &placed),
            placements: placed,
        }))
    }
    fn cuckoo(
        &mut self,
        piece: usize,
        o: usize,
        layout: &[Placement],
        w: &[Vec<f64>],
        length: f64,
        budget: &Budget,
        rng: &mut Rng,
        stats: &mut Stats,
    ) -> Result<Option<(Placement, f64)>> {
        let bb = Bounds::of(&oriented(self.g, piece, o));
        let xmax = length - bb.x1;
        let ymax = self.g.height - bb.y1;
        let xmin = (-bb.x0).max(0.);
        let ymin = (-bb.y0).max(0.);
        if xmax < xmin || ymax < ymin {
            return Ok(None);
        }
        let random = |rng: &mut Rng| Placement {
            piece,
            orientation: o,
            x: xmin + (xmax - xmin) * rng.unit(),
            y: ymin + (ymax - ymin) * rng.unit(),
        };
        let mut nests = Vec::new();
        for _ in 0..self.cfg.nests {
            if budget.expired() {
                return Ok(None);
            }
            let p = random(rng);
            let f = self.cost(&p, layout, w, budget)?;
            nests.push((p, f));
            stats.evaluations += 1;
        }
        let sigma=levy_sigma(self.cfg.levy_beta);
        for _ in 0..self.cfg.generations {
            if budget.expired() {
                return Ok(None);
            }
            let mut p = nests[rng.index(self.cfg.nests)].0.clone();
            let step = normal(rng) * sigma
                / normal(rng).abs().max(f64::MIN_POSITIVE).powf(1. / self.cfg.levy_beta);
            let direction = rng.unit() * std::f64::consts::TAU;
            let x=p.x+step*direction.cos();let y=p.y+step*direction.sin();
            if self.cfg.boundary=="resample" && (x<xmin||x>xmax||y<ymin||y>ymax) {p=random(rng);}
            else {p.x=x.clamp(xmin,xmax);p.y=y.clamp(ymin,ymax);}
            let f = self.cost(&p, layout, w, budget)?;
            stats.evaluations += 1;
            let j = rng.index(self.cfg.nests);
            if f < nests[j].1 {
                nests[j] = (p, f)
            }
            nests.sort_by(|a, b| a.1.total_cmp(&b.1));
            let discard=(self.cfg.abandonment*self.cfg.nests as f64).ceil() as usize;
            for item in nests.iter_mut().skip(self.cfg.nests-discard) {
                if budget.expired() {
                    return Ok(None);
                }
                let p = random(rng);
                let f = self.cost(&p, layout, w, budget)?;
                *item = (p, f);
                stats.evaluations += 1;
            }
        }
        Ok(nests.into_iter().min_by(|a, b| a.1.total_cmp(&b.1)))
    }
    fn minimize(
        &mut self,
        layout: &mut [Placement],
        length: f64,
        budget: &Budget,
        rng: &mut Rng,
        stats: &mut Stats,
    ) -> Result<bool> {
        let n = layout.len();
        if layout.iter().any(|p| !self.inside(p, length)) {
            return Ok(false);
        }
        let mut w = vec![vec![1.; n]; n];
        let mut stagnant = 0;
        let mut best = f64::INFINITY;
        while stagnant < self.cfg.stagnation && !budget.expired() {
            let mut order: Vec<_> = (0..n).collect();
            rng.shuffle(&mut order);
            for index in order {
                let mut f = self.cost(&layout[index], layout, &w, budget)?;
                let piece = layout[index].piece;
                for o in 0..self.g.types[self.g.pieces[piece]].orientations.len() {
                    let Some((p, cost)) =
                        self.cuckoo(piece, o, layout, &w, length, budget, rng, stats)?
                    else {
                        if budget.expired() {
                            return Ok(false);
                        } else {
                            continue;
                        }
                    };
                    if cost < f {
                        layout[index] = p;
                        f = cost;
                        stats.accepted_moves += 1;
                    }
                }
            }
            let pairs = self.losses(layout, budget)?;
            let value: f64 = pairs.iter().map(|p| p.2).sum();
            if value == 0. && layout.iter().all(|p| self.inside(p, length)) {
                return Ok(true);
            }
            if value < best {
                best = value;
                stagnant = 0;
            }
            let max = pairs.iter().map(|p| p.2).fold(0., f64::max);
            if !max.is_finite() {
                return Ok(false);
            }
            for (i, j, f) in pairs {
                w[i][j] += f / max;
                w[j][i] = w[i][j];
            }
            stagnant += 1;
        }
        Ok(false)
    }
}
fn normal(rng: &mut Rng) -> f64 {
    (-2. * rng.unit().max(f64::MIN_POSITIVE).ln()).sqrt()
        * (rng.unit() * std::f64::consts::TAU).cos()
}
fn levy_sigma(beta:f64)->f64 {
    // Lanczos Gamma on positive arguments >= 0.5; Mantegna scale.
    fn gamma(z:f64)->f64 {
        let c=[676.5203681218851,-1259.1392167224028,771.3234287776531,
            -176.6150291621406,12.507343278686905,-0.13857109526572012,
            9.984369578019572e-6,1.5056327351493116e-7];
        let z=z-1.;let mut x=0.9999999999998099;
        for (i,a) in c.iter().enumerate(){x+=a/(z+i as f64+1.);}
        let t=z+7.5;(2.*std::f64::consts::PI).sqrt()*t.powf(z+0.5)*(-t).exp()*x
    }
    (gamma(1.+beta)*(std::f64::consts::PI*beta/2.).sin()
        /(gamma((1.+beta)/2.)*beta*2f64.powf((beta-1.)/2.))).powf(1./beta)
}
fn next_length(current:f64,best:f64,feasible:bool,cfg:&crate::config::Gcs)->f64 {
    if feasible {current*(1.-cfg.shrink)}
    else {(if cfg.expansion_base=="best" {best}else{current})*(1.+cfg.expand)}
}

fn initial_orders(g: &Geometry) -> Vec<Vec<usize>> {
    (0..5)
        .map(|rule| {
            let mut order: Vec<_> = (0..g.pieces.len()).collect();
            let score = |p: usize| {
                let t = &g.types[g.pieces[p]];
                let o = &t.orientations[0];
                let area = t.area;
                let hull = polygon(&convex_hull(&oriented(g, p, 0)));
                use geo::Area;
                match rule {
                    0 => area,
                    1 => o.width,
                    2 => o.height,
                    3 => hull.unsigned_area() / area,
                    _ => area / (o.width * o.height),
                }
            };
            order.sort_by(|&a, &b| score(b).total_cmp(&score(a)).then(a.cmp(&b)));
            order
        })
        .collect()
}

/// Deterministic preprocessing diagnostic; no stochastic optimization or BKV
/// layout is consumed. Compare with original Table3 before timed experiments.
pub fn initial_densities(g: &Geometry, budget: &Budget) -> Result<Vec<f64>> {
    let mut e = Engine {
        g,
        cache: Cache::new(64 * 1024 * 1024),
        cfg: &crate::config::get().gcs,
    };
    let mut stats = Stats::default();
    let mut rows = Vec::new();
    for order in initial_orders(g) {
        let layout = e
            .initial(&order, budget, &mut stats)?
            .expect("complete initializer diagnostic");
        rows.push(100. * g.total_area / (g.height * layout.length));
    }
    Ok(rows)
}

/// Restore all poses from the common reversed working frame. The original
/// item outlines and per-item rotation choices are unchanged on export.
pub fn reverse_global_y(g: &Geometry, layout: &Layout) -> Layout {
    Layout {
        length: layout.length,
        placements: layout
            .placements
            .iter()
            .map(|p| {
                let mut q = p.clone();
                q.y = g.height - g.orientation(p).height - p.y;
                q
            })
            .collect(),
    }
}

pub fn solve(
    g: &Geometry,
    name: &str,
    budget: &Budget,
    seed: u64,
    stats: &mut Stats,
) -> Result<Layout> {
    solve_with_config(g,name,budget,seed,stats,&crate::config::get().gcs)
}
pub fn solve_with_config(g:&Geometry,name:&str,budget:&Budget,seed:u64,stats:&mut Stats,cfg:&crate::config::Gcs)->Result<Layout>{
    cfg.validate()?;
    stats.effective_algorithm_controls=Some(serde_json::to_value(cfg)?);
    ensure!(
        matches!(name.to_ascii_lowercase().as_str(), "dighe1" | "dighe2"),
        "This implementation supports only the published zero-clustering Dighe cases"
    );
    let mut e = Engine {
        g,
        cache: Cache::new(64 * 1024 * 1024),
        cfg,
    };
    let mut rng = Rng::new(seed);
    let order: Vec<_> = (0..g.pieces.len()).collect();
    let mut best = g.shelf(&order);
    for (rule, order) in initial_orders(g).into_iter().enumerate() {
        if let Some(layout) = e.initial(&order, budget, stats)? {
            stats.backend_events.push(format!(
                "initializer_rule={rule}; utilization_percent={:.12}",
                100. * g.total_area / (g.height * layout.length)
            ));
            if layout.length < best.length {
                best = layout;
            }
        }
    }
    let lower = g.total_area / g.height;
    let mut current = best.placements.clone();
    let mut length = best.length * (1.-cfg.shrink);
    let mut outer = 0;
    let mut domain_extension_recorded=false;
    while !budget.expired() {
        if length>g.upper && !domain_extension_recorded {
            stats.backend_events.push(format!("GCS horizontal validation domain may expand: initial={} attempted={length}; fixed strip width={} unchanged",g.upper,g.height));
            domain_extension_recorded=true;
        }
        let mut all_inside = true;
        for i in 0..current.len() {
            if !e.inside(&current[i], length) {
                if let Some(p) = e.relocate(i, &current, length, budget, stats)? {
                    current[i] = p;
                } else {
                    all_inside = false;
                    break;
                }
            }
        }
        let mut feasible = false;
        if all_inside && e.minimize(&mut current, length, budget, &mut rng, stats)? {
            let trial = Layout {
                length: geometry::length(g, &current),
                placements: current.clone(),
            };
            if trial.length <= length + g.height * 1e-9 && g.validate_search_layout(&trial).is_ok() {
                feasible = true;
                if length < best.length {
                    best = Layout {
                        length,
                        placements: trial.placements,
                    };
                }
            }
        }
        // Follow the source prose's current-container controller. Do not clamp
        // to the area bound: that changes the sequence of attempted lengths.
        let next=next_length(length,best.length,feasible,cfg);
        if outer<128 {stats.backend_events.push(format!("GCS controller step={outer} length={length} feasible={feasible} best={} next={next}; first128 steps retained",best.length));}
        length = next;
        outer += 1;
    }
    stats.nfp_cache_bytes = e.cache.bytes();
    stats.paper_run_status = Some(
        if best.length <= lower + g.height * 1e-9 {
            "gcs_area_bound_reached"
        } else {
            "gcs_budget_stop"
        }
        .into(),
    );
    stats.backend_events.push(format!(
        "{ASSUMPTIONS}; outer_iterations={outer}; cache_hits={}; cache_misses={}",
        e.cache.hits, e.cache.misses
    ));
    Ok(best)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]fn search_length_can_exceed_initial_oracle_domain_without_widening_strip()->Result<()>{
        let g=rectangles();
        let poses=vec![Placement{piece:0,orientation:0,x:0.,y:0.},
            Placement{piece:1,orientation:0,x:2.,y:0.},
            Placement{piece:2,orientation:0,x:g.upper+1.,y:0.}];
        let layout=Layout{length:geometry::length(&g,&poses),placements:poses};
        assert!(g.validate(&layout).is_err());
        g.validate_search_layout(&layout)?;
        let mut outside=layout.clone();outside.placements[2].y=0.1;
        assert!(g.validate_search_layout(&outside).is_err());
        let mut overlap=layout.clone();overlap.placements[1].x=0.5;
        assert!(g.validate_search_layout(&overlap).is_err());
        assert_eq!(g.height,1.);
        Ok(())
    }
    #[test]fn levy_scale_and_conflicting_length_controllers_are_explicit(){
        assert!((levy_sigma(1.5)-0.6965745025576967).abs()<1e-14);
        assert!((levy_sigma(1.)-1.).abs()<1e-14);
        let mut c=crate::config::get().gcs.clone();
        assert!((next_length(90.,100.,false,&c)-90.9).abs()<1e-12);
        c.expansion_base="best".into();
        assert_eq!(next_length(90.,100.,false,&c),101.);
        assert_eq!(next_length(100.,100.,true,&c),96.);
    }
    fn rectangles() -> Geometry {
        use crate::input::LoadedInput;
        use std::collections::BTreeMap;
        let vertices = vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.)];
        let instance=serde_json::from_value(serde_json::json!({"name":"rectangles","strip_height":1.,"min_item_separation":0.,"items":[{"id":0,"demand":3,"orientation":{"rotation":{"mode":"discrete","angles":[0.]}},"shape":{"type":"simple_polygon","data":vertices}}]})).unwrap();
        Geometry::from_input(&LoadedInput {
            instance,
            item_labels: BTreeMap::new(),
            source_path: Default::default(),
            exact_vertices: BTreeMap::from([(0, vertices)]),
            source_up_left: false,
        })
        .unwrap()
    }
    #[test]fn configured_cuckoo_counts_and_boundary_modes_are_executed()->Result<()>{
        let g=rectangles();
        for boundary in ["clip","resample"] {for beta in [1.,1.5] {
            let mut cfg=crate::config::get().gcs.clone();cfg.nests=4;cfg.generations=3;cfg.abandonment=0.5;
            cfg.boundary=boundary.into();cfg.levy_beta=beta;cfg.validate()?;
            let mut e=Engine{g:&g,cache:Cache::new(1024),cfg:&cfg};let mut stats=Stats::default();
            let (p,cost)=e.cuckoo(0,0,&[],&vec![vec![1.;3];3],3.,&Budget::new(1.),&mut Rng::new(17),&mut stats)?.unwrap();
            assert!(p.x>=0. && p.x<=2. && p.y==0.);assert_eq!(cost,0.);
            assert_eq!(stats.evaluations,4+3*(1+2));
        }}Ok(())
    }
    #[test]
    fn protruding_piece_is_relocated_into_isolated_contact_gap() -> Result<()> {
        let g = rectangles();
        let mut e = Engine {
            g: &g,
            cache: Cache::new(1024 * 1024),
            cfg: &crate::config::get().gcs,
        };
        let poses: Vec<_> = [0., 2., 4.]
            .into_iter()
            .enumerate()
            .map(|(piece, x)| Placement {
                piece,
                orientation: 0,
                x,
                y: 0.,
            })
            .collect();
        let moved = e
            .relocate(2, &poses, 3., &Budget::new(2.), &mut Stats::default())?
            .unwrap();
        assert!((moved.x - 1.).abs() < 1e-12);
        assert_eq!(moved.y, 0.);
        let mut impossible = poses.clone();
        assert!(!e.minimize(
            &mut impossible,
            0.9,
            &Budget::new(2.),
            &mut Rng::new(1),
            &mut Stats::default()
        )?);
        assert!(
            e.relocate(2, &poses, 0.9, &Budget::new(2.), &mut Stats::default())?
                .is_none()
        );
        Ok(())
    }
    #[test]
    fn published_initialization_metrics_and_exported_original_contours() -> Result<()> {
        for (name, rule, reported) in [("dighe1", 0, 68.78), ("dighe2", 1, 68.70)] {
            let input = crate::input::load(
                &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join(format!("../inputs/{name}.json")),
            )?;
            let canonical = Geometry::from_input(&input)?;
            let working = Geometry::from_input(&input.opposite_y_frame())?;
            let mut e = Engine {
                g: &working,
                cache: Cache::new(64 * 1024 * 1024),
                cfg: &crate::config::get().gcs,
            };
            let layout = e
                .initial(
                    &initial_orders(&working)[rule],
                    &Budget::new(10.),
                    &mut Stats::default(),
                )?
                .unwrap();
            let measured = 100. * working.total_area / (working.height * layout.length);
            assert!(
                (measured * 100.).round() / 100. == reported,
                "{name}: {measured}"
            );
            let restored = reverse_global_y(&canonical, &layout);
            use geo::{Area, BooleanOps, BoundingRect};
            let polygons: Vec<_> = restored
                .placements
                .iter()
                .map(|p| world(&canonical, p))
                .collect();
            let mut max_overlap = 0_f64;
            let mut max_outside = 0_f64;
            for (i, p) in polygons.iter().enumerate() {
                let b = p.bounding_rect().unwrap();
                max_outside = max_outside
                    .max(-b.min().x)
                    .max(-b.min().y)
                    .max(b.max().y - canonical.height);
                for q in &polygons[..i] {
                    max_overlap = max_overlap.max(p.intersection(q).unsigned_area());
                }
            }
            eprintln!("{name} restored overlap={max_overlap:e}; outside={max_outside:e}");
            canonical.validate(&restored)?;
        }
        Ok(())
    }
    #[test]
    fn euclidean_penalty_is_not_axis_or_area_penalty() -> Result<()> {
        let square = vec![
            P::new(0., 0.),
            P::new(2., 0.),
            P::new(2., 2.),
            P::new(0., 2.),
        ];
        let nfp = crate::toolkit::nfp::build(&square, &square, &Budget::new(2.))?.unwrap();
        assert!((penetration::squared(&nfp, P::new(1.5, 1.5)).0.sqrt() - 0.5).abs() < 1e-8);
        assert_eq!(penetration::squared(&nfp, P::new(2., 0.)).0, 0.);
        Ok(())
    }
}
