//! Sykora thesis Algorithm2: HS2 constructive insertion, n_ini=12.
//! Old pair slice decisions are fixed; coordinates remain free at every insertion.
use crate::{
    backends::{Settings, cbc},
    geometry::{self, Geometry},
    model::*,
    toolkit::{nfp, placement::*, polygon::*},
};
use anyhow::{Result, ensure};
use geo::{Area, BooleanOps};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
struct Cell {
    vertices: Vec<P>,
    planes: Vec<(f64, f64, f64)>,
}
fn clip(v: &[P], a: f64, b: f64, c: f64) -> Vec<P> {
    if v.is_empty() {
        return vec![];
    }
    let mut out = Vec::new();
    let eps = 1e-10;
    for edge in edges(v) {
        let f = a * edge.a.x + b * edge.a.y - c;
        let h = a * edge.b.x + b * edge.b.y - c;
        if f <= eps {
            out.push(edge.a);
        }
        if (f < -eps && h > eps) || (f > eps && h < -eps) {
            let t = f / (f - h);
            out.push(P::new(
                edge.a.x + t * (edge.b.x - edge.a.x),
                edge.a.y + t * (edge.b.y - edge.a.y),
            ));
        }
    }
    out.dedup_by(|a, b| (a.x - b.x).hypot(a.y - b.y) < eps);
    out
}
fn cell(v: Vec<P>) -> Option<Cell> {
    if v.is_empty() {
        return None;
    }
    let mut v = v;
    if signed_area(&v) < 0.0 {
        v.reverse();
    }
    let mut planes = Vec::new();
    if v.len() >= 3 && signed_area(&v).abs() > 1e-12 {
        for e in edges(&v) {
            let a = e.b.y - e.a.y;
            let b = e.a.x - e.b.x;
            let n = a.hypot(b);
            if n > 1e-12 {
                planes.push((a / n, b / n, (a * e.a.x + b * e.a.y) / n));
            }
        }
    } else {
        let b = Bounds::of(&v);
        planes.extend([
            (1., 0., b.x1),
            (-1., 0., -b.x0),
            (0., 1., b.y1),
            (0., -1., -b.y0),
        ]);
        if v.len() > 1 {
            let a = v.last().unwrap().y - v[0].y;
            let b = v[0].x - v.last().unwrap().x;
            let c = a * v[0].x + b * v[0].y;
            if a.hypot(b) > 1e-12 {
                planes.push((a, b, c));
                planes.push((-a, -b, -c));
            }
        }
    }
    Some(Cell {
        vertices: v,
        planes,
    })
}
fn rectangle(x: f64, y: f64) -> Vec<P> {
    vec![P::new(-x, -y), P::new(x, -y), P::new(x, y), P::new(-x, y)]
}

/// Close the largest convex concavity repeatedly (thesis2.2.2). A single closure
/// is used per update; disjoint closures commute. Reject a chain intersecting
/// the original NFP interior instead of accepting a false exterior region.
fn slices(n: &nfp::NoFit, upper: f64, height: f64) -> Result<Vec<Cell>> {
    ensure!(
        n.region.0.len() == 1,
        "CHS audited Poly branch requires connected NFP"
    );
    let poly = &n.region.0[0];
    let mut v: Vec<_> = poly.exterior().0.iter().map(|p| P::new(p.x, p.y)).collect();
    if v.first() == v.last() {
        v.pop();
    }
    if signed_area(&v) < 0. {
        v.reverse();
    }
    let mut cells = Vec::new();
    // Inner feasible polygons are supported; the selected Poly input has none.
    for hole in poly.interiors() {
        let h: Vec<_> = hole
            .0
            .iter()
            .take(hole.0.len() - 1)
            .map(|p| P::new(p.x, p.y))
            .collect();
        for t in triangulate(&h)? {
            if let Some(c) = cell(t.to_vec()) {
                cells.push(c);
            }
        }
    }
    while (polygon(&v).unsigned_area() - polygon(&convex_hull(&v)).unsigned_area()).abs() > 1e-9 {
        let mut best: Option<(usize, usize, Vec<P>)> = None;
        for start in 0..v.len() {
            for count in 3..v.len() {
                let q: Vec<_> = (0..count).map(|k| v[(start + k) % v.len()]).collect();
                if signed_area(&q) >= -1e-10 {
                    continue;
                }
                if (convex_hull(&q).len() != q.len())
                    && q.windows(3).any(|p| cross(p[0], p[1], p[2]) > 1e-10)
                {
                    continue;
                }
                if polygon(&q).intersection(&polygon(&v)).unsigned_area() > 1e-8 {
                    continue;
                }
                if (polygon(&q).unsigned_area() - polygon(&convex_hull(&q)).unsigned_area()).abs()
                    > 1e-8
                {
                    continue;
                }
                if best.as_ref().is_none_or(|(_, c, _)| count > *c) {
                    best = Some((start, count, q));
                }
            }
        }
        let Some((start, count, q)) = best else {
            anyhow::bail!("could not close NFP concavity faithfully");
        };
        cells.push(cell(q).unwrap());
        let len = v.len();
        let remove: Vec<_> = (1..count - 1).map(|k| (start + k) % len).collect();
        v = v
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !remove.contains(i))
            .map(|(_, p)| p)
            .collect();
    }
    let hull = convex_hull(&v);
    let bounds = Bounds::of(&hull);
    // Each nonhorizontal convex-hull edge defines a horizontally bounded outer slice.
    for e in edges(&hull) {
        if (e.a.y - e.b.y).abs() < 1e-12 {
            continue;
        }
        let a = -(e.b.y - e.a.y);
        let b = -(e.a.x - e.b.x);
        let c = a * e.a.x + b * e.a.y;
        let mut q = clip(&rectangle(upper, height), a, b, c);
        q = clip(&q, 0., 1., e.a.y.max(e.b.y));
        q = clip(&q, 0., -1., -e.a.y.min(e.b.y));
        if let Some(c) = cell(q) {
            cells.push(c);
        }
    }
    for (a, b, c) in [(0., 1., bounds.y0), (0., -1., -bounds.y1)] {
        if let Some(c) = cell(clip(&rectangle(upper, height), a, b, c)) {
            cells.push(c);
        }
    }
    for z in &n.zero_strata {
        if let Some(c) = cell(vec![z.a, z.b]) {
            cells.push(c);
        }
    }
    ensure!(!cells.is_empty(), "no HS2 feasible slices");
    Ok(cells)
}
fn addrow(m: &mut cbc::Model, terms: Vec<(usize, f64)>, sense: cbc::Sense, rhs: f64) {
    let mut map = BTreeMap::<usize, f64>::new();
    for (i, a) in terms {
        *map.entry(i).or_default() += a;
    }
    let terms = map
        .into_iter()
        .filter(|(_, a)| *a != 0.)
        .collect::<Vec<_>>();
    if !terms.is_empty() {
        m.rows.push(cbc::Row { terms, sense, rhs });
    }
}
fn var(m: &mut cbc::Model, lower: f64, upper: f64, objective: f64, binary: bool) -> usize {
    let i = m.variables.len();
    m.variables.push(cbc::Variable {
        lower,
        upper,
        objective,
        binary,
    });
    i
}
fn weighted_order(g: &Geometry, rng: &mut Rng) -> Vec<usize> {
    let mut remaining = (0..g.pieces.len()).collect::<Vec<_>>();
    let mut out = Vec::new();
    while !remaining.is_empty() {
        let available = (0..g.types.len())
            .filter(|&t| remaining.iter().any(|&p| g.pieces[p] == t))
            .collect::<Vec<_>>();
        let mut q = rng.unit() * available.iter().map(|&t| g.types[t].area).sum::<f64>();
        let mut chosen = *available.last().unwrap();
        for &t in &available {
            q -= g.types[t].area;
            if q <= 0. {
                chosen = t;
                break;
            }
        }
        let i = remaining
            .iter()
            .position(|&p| g.pieces[p] == chosen)
            .unwrap();
        out.push(remaining.remove(i));
    }
    out
}
type Decisions = Vec<(usize, usize, Vec<usize>)>;
fn build_hs2_model(g:&Geometry,order:&[usize],count:usize,upper:f64,
    cached:&mut BTreeMap<(usize,usize),Vec<Cell>>,fixed:&BTreeMap<(usize,usize),usize>,
    budget:&Budget)->Result<(cbc::Model,Decisions)> {
        let mut m = cbc::Model::default();
        var(&mut m, 0., upper, 1., false);
        for &piece in &order[..count] {
            let o = &g.types[g.pieces[piece]].orientations[0];
            let x = var(&mut m, 0., upper - o.width, 0., false);
            var(&mut m, 0., g.height - o.height, 0., false);
            addrow(&mut m, vec![(0, -1.), (x, 1.)], cbc::Sense::Less, -o.width);
        }
        let mut decisions = Vec::new();
        for i in 0..count {
            for j in i + 1..count {
                if budget.expired() {
                    break;
                }
                let key = (g.pieces[order[i]], g.pieces[order[j]]);
                if !cached.contains_key(&key) {
                    let a = oriented(g, order[i], 0);
                    let b = oriented(g, order[j], 0);
                    let n = nfp::build(&a, &b, budget)?
                        .ok_or_else(|| anyhow::anyhow!("CHS NFP deadline"))?;
                    cached.insert(key, slices(&n, upper, g.height)?);
                }
                let cells = &cached[&key];
                let mut zs = Vec::new();
                for k in 0..cells.len() {
                    let isfixed = fixed.get(&(i, j)).copied();
                    let z = var(
                        &mut m,
                        if isfixed == Some(k) { 1. } else { 0. },
                        if isfixed.is_none() || isfixed == Some(k) {
                            1.
                        } else {
                            0.
                        },
                        0.,
                        true,
                    );
                    zs.push(z);
                }
                addrow(
                    &mut m,
                    zs.iter().map(|&z| (z, 1.)).collect(),
                    cbc::Sense::Equal,
                    1.,
                );
                // Eq2.25: all slice support maxima, with no big-M coefficients.
                for c in cells {
                    for &(a, b, _) in &c.planes {
                        let mut terms = vec![
                            (1 + 2 * i, -a),
                            (2 + 2 * i, -b),
                            (1 + 2 * j, a),
                            (2 + 2 * j, b),
                        ];
                        for (&z, d) in zs.iter().zip(cells) {
                            let h = d
                                .vertices
                                .iter()
                                .map(|p| a * p.x + b * p.y)
                                .fold(f64::NEG_INFINITY, f64::max);
                            terms.push((z, -h));
                        }
                        addrow(&mut m, terms, cbc::Sense::Less, 0.);
                    }
                }
                // HS2 bound lifting in both directions. Cell extrema imply the
                // nonnegative coordinate/upper strip bounds of its other participant.
                let oi = &g.types[g.pieces[order[i]]].orientations[0];
                let oj = &g.types[g.pieces[order[j]]].orientations[0];
                for (piece, other, sign) in [(i, oi, 1.), (j, oj, -1.)] {
                    let mut lowx = vec![(1 + 2 * piece, -1.)];
                    let mut lowy = vec![(2 + 2 * piece, -1.)];
                    let mut highx = vec![(0, -1.), (1 + 2 * piece, 1.)];
                    let mut highy = vec![(2 + 2 * piece, 1.)];
                    let own = if piece == i { oi } else { oj };
                    let opposite = if piece == i { oj } else { oi };
                    let _ = other;
                    for (&z, c) in zs.iter().zip(cells) {
                        let b = Bounds::of(&c.vertices);
                        let (dx0, dx1, dy0, dy1) = if sign > 0. {
                            (b.x0, b.x1, b.y0, b.y1)
                        } else {
                            (-b.x1, -b.x0, -b.y1, -b.y0)
                        };
                        lowx.push((z, (-dx1).max(0.)));
                        lowy.push((z, (-dy1).max(0.)));
                        highx.push((z, (opposite.width + dx0 - own.width).max(0.)));
                        highy.push((z, (opposite.height + dy0 - own.height).max(0.)));
                    }
                    addrow(&mut m, lowx, cbc::Sense::Less, 0.);
                    addrow(&mut m, lowy, cbc::Sense::Less, 0.);
                    addrow(&mut m, highx, cbc::Sense::Less, -own.width);
                    addrow(&mut m, highy, cbc::Sense::Less, g.height - own.height);
                }
                decisions.push((i, j, zs));
            }
        }
    Ok((m,decisions))
}

fn warm_slice(m:&cbc::Model,cells:&[Cell],zs:&[usize],fixed:Option<usize>,dx:f64,dy:f64)->Result<usize> {
    ensure!(cells.len()==zs.len() && zs.iter().all(|&z|z<m.variables.len()),"HS2 slice/variable mapping mismatch");
    if let Some(k)=fixed {
        ensure!(k<zs.len() && m.variables[zs[k]].lower==1. && m.variables[zs[k]].upper==1.,
            "HS2 inherited slice disagrees with original model bounds");
        // An existing pair has no new discrete decision. Inherit its already
        // fixed slice; the full original-model check below still validates the
        // coordinates. A stricter reclassification can reject an incumbent
        // accepted by the original MIP's numerical feasibility tolerance.
        return Ok(k);
    }
    zs.iter().enumerate().find(|(k,z)|m.variables[**z].upper>0.
        && cells[*k].planes.iter().all(|&(a,b,c)|a*dx+b*dy<=c+1e-8*(1.+c.abs())))
        .map(|(k,_)|k).ok_or_else(||anyhow::anyhow!("HS2 new pair has no admissible warm-start slice"))
}

fn build_warm_start(g:&Geometry,order:&[usize],count:usize,last:Option<&Layout>,
    m:&cbc::Model,decisions:&Decisions,cached:&BTreeMap<(usize,usize),Vec<Cell>>,
    fixed:&BTreeMap<(usize,usize),usize>)->Result<Vec<f64>> {
        let mut poses=last.map(|l|l.placements.clone()).unwrap_or_default();
        while poses.len()<count {
            let index=poses.len();
            poses.push(Placement{piece:order[index],orientation:0,
                x:geometry::length(g,&poses),y:0.});
        }
        let mut warm=vec![0.;m.variables.len()];
        warm[0]=geometry::length(g,&poses);
        for (i,p) in poses.iter().enumerate(){warm[1+2*i]=p.x;warm[2+2*i]=p.y;}
        for (i,j,zs) in decisions {
            let dx=poses[*j].x-poses[*i].x;let dy=poses[*j].y-poses[*i].y;
            let cells=&cached[&(g.pieces[order[*i]],g.pieces[order[*j]])];
            let k=warm_slice(m,cells,zs,fixed.get(&(*i,*j)).copied(),dx,dy)
                .map_err(|e|anyhow::anyhow!("HS2 warm start pair {i},{j}: {e}"))?;
            warm[zs[k]]=1.;
        }
        ensure!(m.feasible(&warm),"HS2 warm start violates the original MIP");
    Ok(warm)
}

fn save_failed_warm(g:&Geometry,order:&[usize],count:usize,last:Option<&Layout>,
    m:&cbc::Model,decisions:&Decisions,cached:&BTreeMap<(usize,usize),Vec<Cell>>,
    fixed:&BTreeMap<(usize,usize),usize>,error:&anyhow::Error)->Result<()> {
    if let Some(path)=std::env::var_os("PACKING_CBC_DIAGNOSTIC_DIR") {
        let path=std::path::PathBuf::from(path).join(format!("warm-failure-{count:04}"));
        std::fs::create_dir_all(&path)?;
        std::fs::write(path.join("model.mps"),m.mps()?)?;
        let mut candidate=last.map(|l|l.placements.clone()).unwrap_or_default();
        while candidate.len()<count {
            let piece=order[candidate.len()];
            let x=geometry::length(g,&candidate);
            candidate.push(Placement{piece,orientation:0,x,y:0.});
        }
        let pairs=decisions.iter().map(|(i,j,zs)|{
            let cells=&cached[&(g.pieces[order[*i]],g.pieces[order[*j]])];
            let dx=candidate[*j].x-candidate[*i].x;let dy=candidate[*j].y-candidate[*i].y;
            serde_json::json!({"pair":[i,j],"fixed":fixed.get(&(*i,*j)),"variables":zs,"relative_translation":[dx,dy],
                "cells":cells.iter().map(|c|serde_json::json!({"vertices":c.vertices.iter().map(|p|[p.x,p.y]).collect::<Vec<_>>(),"planes":c.planes})).collect::<Vec<_>>()})
        }).collect::<Vec<_>>();
        std::fs::write(path.join("candidate.json"),serde_json::to_vec_pretty(&serde_json::json!({
            "error":error.to_string(),"count":count,"order":order,
            "previous_placements":last.map(|l|&l.placements),"pairs":pairs,
            "candidate_placements":candidate,
            "candidate_rule":"preserve prior coordinates and append new pieces at current right edge",
            "coordinates_modified_by_diagnostic":false}))?)?;
    }
    Ok(())
}

pub fn solve(g: &Geometry, budget: &Budget, seed: u64, stats: &mut Stats) -> Result<Layout> {
    solve_with_config(g,budget,seed,stats,&crate::config::get().chs2)
}
pub fn solve_with_config(g:&Geometry,budget:&Budget,seed:u64,stats:&mut Stats,controls:&crate::config::Chs)->Result<Layout>{
    ensure!(controls.initial_pieces>0 && controls.mip_seconds.is_finite() && controls.mip_seconds>0.,"invalid CHS2 controls");
    stats.effective_algorithm_controls=Some(serde_json::to_value(controls)?);
    ensure!(
        g.types.iter().all(|t| t.orientations.len() == 1),
        "this CHS configuration is fixed-angle only"
    );
    let mut rng = Rng::new(seed);
    let order = weighted_order(g, &mut rng);
    // CBC's integer preprocessing can overrun its internal elapsed-time limit
    // on the retained HS2 insertion models. Submit the original formulation
    // directly; do not change slices, coordinates, constraints or min-L.
    let cfg = Settings { cbc_preprocess: controls.cbc_preprocess,
        cbc_presolve: controls.cbc_presolve, ..Settings::default() };
    let upper = order
        .iter()
        .map(|&p| g.types[g.pieces[p]].orientations[0].width)
        .sum::<f64>();
    let mut cached = BTreeMap::<(usize, usize), Vec<Cell>>::new();
    let mut fixed = BTreeMap::<(usize, usize), usize>::new();
    let mut last: Option<Layout> = None;
    for count in controls.initial_pieces.min(order.len())..=order.len() {
        if budget.expired() {
            break;
        }
        let (m,decisions)=build_hs2_model(g,&order,count,upper,&mut cached,&fixed,budget)?;
        if budget.expired() {
            break;
        }
        stats.backend_calls += 1;
        let call = Budget::new(budget.remaining().min(controls.mip_seconds));
        // A solver start is an incumbent of this exact HS2 model, not a
        // replacement packing. Keep old coordinates/slices and append each new
        // piece at the current right edge, then encode its admissible slices.
        // CBC still minimizes strip length for the original per-MIP100s cap.
        let warm=match build_warm_start(g,&order,count,last.as_ref(),&m,&decisions,&cached,&fixed) {
            Ok(warm)=>warm,
            Err(error)=>{
                save_failed_warm(g,&order,count,last.as_ref(),&m,&decisions,&cached,&fixed,&error)?;
                return Err(error);
            }
        };
        stats.backend_events.push(format!("HS2 insertion{count}: original-model feasible MIP start L={}; cap={}s",warm[0],controls.mip_seconds));
        let Some(sol) = cbc::solve_with_start(&m, &cfg, &call,Some(&warm))? else {
            stats
                .backend_events
                .push(format!("HS2 insertion{count}: no feasible MIP incumbent"));
            break;
        };
        if !sol.optimal {stats.nonoptimal_backend_calls+=1;}
        for (i, j, zs) in decisions {
            let k = zs.iter().position(|&z| sol.values[z] > 0.5).unwrap();
            fixed.insert((i, j), k);
        }
        let placements = order[..count]
            .iter()
            .enumerate()
            .map(|(i, &piece)| Placement {
                piece,
                orientation: 0,
                x: sol.values[1 + 2 * i].max(0.),
                y: sol.values[2 + 2 * i].max(0.),
            })
            .collect::<Vec<_>>();
        let l = Layout {
            length: geometry::length(g, &placements),
            placements,
        };
        stats.completed_generations = Some(count);
        stats.backend_events.push(format!(
            "HS2 insertion{count}: L={} optimal={}; old slice variables fixed, coordinates free",
            l.length, sol.optimal
        ));
        println!("CHS2 inserted={count} L={} optimal={}",l.length,sol.optimal);
        last = Some(l);
    }
    let complete = last
        .as_ref()
        .is_some_and(|l| l.placements.len() == g.pieces.len());
    stats.paper_run_status = Some(
        if complete {
            "completed_CHS2"
        } else {
            "CHS2_incomplete_no_full_algorithm_result"
        }
        .into(),
    );
    stats.paper_variant =
        Some(format!("HS2 Algorithm2; n_ini{}; per-MIP{}s; one area-weighted random type sequence",controls.initial_pieces,controls.mip_seconds));
    if complete {
        let l = last.unwrap();
        Ok(l)
    } else {
        anyhow::bail!(
            "CHS2 incomplete after {:?} placed items; no shelf relabeling",
            stats.completed_generations
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn retained_line_model()->(cbc::Model,Cell,Vec<f64>) {
        let c=cell(vec![P::new(1.,6.),P::new(2.,5.6)]).unwrap();
        let mut m=cbc::Model::default();
        var(&mut m,-248.,248.,0.,false);var(&mut m,-40.,40.,0.,false);var(&mut m,1.,1.,0.,true);
        for &(a,b,_) in &c.planes {
            let support=c.vertices.iter().map(|p|a*p.x+b*p.y).fold(f64::NEG_INFINITY,f64::max);
            addrow(&mut m,vec![(0,a),(1,b),(2,-support)],cbc::Sense::Less,0.);
        }
        (m,c,vec![0.999999940395347,6.,1.])
    }
    #[test]
    fn inherited_degenerate_contact_is_not_reclassified_with_a_stricter_tolerance(){
        let (m,c,values)=retained_line_model();
        assert!(m.feasible(&values));
        assert!(warm_slice(&m,std::slice::from_ref(&c),&[2],None,values[0],values[1]).is_err());
        assert_eq!(warm_slice(&m,&[c],&[2],Some(0),values[0],values[1]).unwrap(),0);
    }
    #[test]
    fn inheritance_does_not_relax_original_model_or_allow_a_different_fixed_slice(){
        let (mut m,c,mut values)=retained_line_model();values[0]=0.9;
        assert!(!m.feasible(&values));
        assert_eq!(warm_slice(&m,std::slice::from_ref(&c),&[2],Some(0),values[0],values[1]).unwrap(),0);
        assert!(warm_slice(&m,std::slice::from_ref(&c),&[2],Some(1),values[0],values[1]).is_err());
        m.variables[2].lower=0.;
        assert!(warm_slice(&m,&[c],&[2],Some(0),values[0],values[1]).is_err());
    }
    #[test]fn actual_insertion_can_move_old_coordinates_and_fill_a_concave_notch(){
        use std::collections::BTreeMap;
        let a=vec![(0.,0.),(2.,0.),(2.,1.),(0.,1.)];
        let b=vec![(0.,0.),(3.,0.),(3.,1.),(1.,1.),(1.,3.),(0.,3.)];
        let instance=serde_json::from_value(serde_json::json!({"name":"chs-notch","strip_height":3.,"min_item_separation":0.,"items":[
            {"id":0,"demand":1,"orientation":{"rotation":{"mode":"discrete","angles":[0.]}},"shape":{"type":"simple_polygon","data":a}},
            {"id":1,"demand":1,"orientation":{"rotation":{"mode":"discrete","angles":[0.]}},"shape":{"type":"simple_polygon","data":b}}]})).unwrap();
        let g=Geometry::from_input(&crate::input::LoadedInput{instance,item_labels:BTreeMap::new(),source_path:Default::default(),exact_vertices:BTreeMap::from([(0,a),(1,b)]),source_up_left:false}).unwrap();
        let seed=(1..1000).find(|&seed|weighted_order(&g,&mut Rng::new(seed))[0]==0).unwrap();
        let mut stats=Stats::default();let controls=crate::config::Chs{initial_pieces:1,mip_seconds:5.,cbc_preprocess:false,cbc_presolve:false};
        let l=solve_with_config(&g,&Budget::new(15.),seed,&mut stats,&controls).unwrap();
        assert_eq!(stats.backend_calls,2);assert!((l.length-3.).abs()<1e-8);
        let rectangle=l.placements.iter().find(|p|p.piece==0).unwrap();
        assert!(rectangle.x>=1.-1e-8 && rectangle.y>=1.-1e-8);g.validate(&l).unwrap();
        assert!(stats.backend_events.iter().filter(|e|e.contains("optimal=true")).count()==2);
    }
    #[test]
    fn horizontal_slices_cover_exact_concave_free_space() {
        let a = vec![
            P::new(0., 0.),
            P::new(4., 0.),
            P::new(4., 1.),
            P::new(1., 1.),
            P::new(1., 4.),
            P::new(0., 4.),
        ];
        let b = vec![
            P::new(0., 0.),
            P::new(1., 0.),
            P::new(1., 1.),
            P::new(0., 1.),
        ];
        let n = nfp::build(&a, &b, &Budget::new(2.)).unwrap().unwrap();
        let s = slices(&n, 10., 10.).unwrap();
        for ix in -12..=20 {
            for iy in -12..=20 {
                let p = P::new(ix as f64 * 0.25, iy as f64 * 0.25);
                let in_slice = s.iter().any(|c| {
                    c.planes
                        .iter()
                        .all(|&(a, b, d)| a * p.x + b * p.y <= d + 1e-9)
                });
                let shifted = b.iter().map(|q| q.shift(p)).collect::<Vec<_>>();
                let legal = polygon(&a).intersection(&polygon(&shifted)).unsigned_area() < 1e-10;
                assert_eq!(in_slice, legal, "slice mismatch at {p:?}");
            }
        }
    }
    #[test]
    fn official_poly2b_every_slice_stays_outside_the_forbidden_region(){
        let path=std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../official_esicup/xml/poly2b.xml");
        let input=crate::input::load(&path).unwrap();let g=Geometry::from_input(&input).unwrap();
        let budget=Budget::new(60.);
        for i in 0..g.pieces.len(){for j in i+1..g.pieces.len(){
            let a=oriented(&g,i,0);let b=oriented(&g,j,0);
            let n=nfp::build(&a,&b,&budget).unwrap().unwrap();let cells=slices(&n,200.,40.).unwrap();
            for c in cells{let center=P::new(c.vertices.iter().map(|p|p.x).sum::<f64>()/c.vertices.len() as f64,c.vertices.iter().map(|p|p.y).sum::<f64>()/c.vertices.len() as f64);
                for v in c.vertices.iter().copied().chain(std::iter::once(center)){
                    for t in [0.,0.25,0.5,0.75,1.]{let p=P::new(center.x*t+v.x*(1.-t),center.y*t+v.y*(1.-t));
                        let q=b.iter().map(|v|v.shift(p)).collect::<Vec<_>>();
                        let overlap=polygon(&a).intersection(&polygon(&q)).unsigned_area();
                        assert!(overlap<1e-7,"pieces{i},{j} invalid slice point{p:?} overlap={overlap}");
                    }
                }
            }
        }}
    }
}
