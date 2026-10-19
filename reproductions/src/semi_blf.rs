use crate::{
    geometry::{self, Geometry},
    model::*,
    toolkit::{nfp, placement::*, polygon::*, semi_discrete::*},
};
use anyhow::Result;
use std::collections::HashMap;

/// Exact original NFP events on ONE fixed x-column. A rejected discretized
/// candidate must not make us skip all later continuous y positions there.
fn contact_heights(g:&Geometry,placed:&[Placement],piece:usize,o:usize,x:f64,
    minimum:f64,maximum:f64,cache:&mut nfp::Cache,budget:&Budget)->Result<Option<Vec<f64>>>{
    let moving=oriented(g,piece,o);let mut ys=vec![minimum,maximum];
    for q in placed {
        let Some(n)=cache.get((g.pieces[q.piece],q.orientation,g.pieces[piece],o),
            &oriented(g,q.piece,q.orientation),&moving,budget)? else{return Ok(None);};
        let relx=x-q.x;
        for &e in n.boundary.iter().chain(&n.zero_strata){
            if let Some(p)=vertical(e,relx){ys.push(q.y+p.y);}
            if e.a.x==relx {ys.push(q.y+e.a.y);}
            if e.b.x==relx {ys.push(q.y+e.b.y);}
        }
    }
    ys.retain(|y|y.is_finite() && *y>=minimum && *y<=maximum);
    ys.sort_by(f64::total_cmp);ys.dedup();Ok(Some(ys))
}

#[allow(dead_code)]
fn resolution(g: &Geometry) -> f64 {
    let mut minimum = f64::INFINITY;
    let mut smallest = (f64::INFINITY, 0.0, 1);
    for (t, kind) in g.types.iter().enumerate() {
        let piece = g.pieces.iter().position(|&p| p == t).unwrap();
        for o in 0..kind.orientations.len() {
            let v = oriented(g, piece, o);
            let b = Bounds::of(&v);
            for edge in edges(&v) {
                let dx = (edge.a.x - edge.b.x).abs();
                if dx > g.height * 1e-12 {
                    minimum = minimum.min(dx);
                }
            }
            if kind.area < smallest.0 {
                smallest = (kind.area, b.width(), v.len());
            }
        }
    }
    minimum
        .max(smallest.1 / smallest.2 as f64)
        .max(g.height * 1e-8)
}
pub fn construct(
    g: &Geometry,
    order: &[usize],
    step: f64,
    budget: &Budget,
    stats: &mut Stats,
) -> Result<Option<Layout>> {
    let mut shapes = HashMap::<(usize, usize), Shape>::new();
    let mut strip = Strip::new();
    let mut placed = Vec::<Placement>::new();
    let mut previous = HashMap::<(usize, usize), (usize, f64)>::new();
    let mut original_cache=nfp::Cache::new(64*1024*1024);
    for &piece in order {
        let t = g.pieces[piece];
        let mut best: Option<(Placement, usize, f64)> = None;
        for o in 0..g.types[t].orientations.len() {
            if budget.expired() {
                return Ok(None);
            }
            if let std::collections::hash_map::Entry::Vacant(entry) = shapes.entry((t, o)) {
                let moving = oriented(g, piece, o);
                let shape = match discretize(&moving, step, budget) {
                    Ok(Some(s)) => s,
                    Ok(None) => return Ok(None),
                    Err(e) => {
                        stats
                            .backend_events
                            .push(format!("Semi-discrete conversion incomplete: {e}"));
                        return Ok(None);
                    }
                };
                entry.insert(shape);
            }
            let shape = &shapes[&(t, o)];
            let h = shape.bounds.height();
            let minimum = shape.origin.y.max(0.0);
            let maximum =
                (g.height - h).min(g.height - g.types[t].orientations[o].height + shape.origin.y);
            if maximum < minimum {
                continue;
            }
            let (mut column, mut start) = previous.get(&(t, o)).copied().unwrap_or((0, minimum));
            let end = ((g.upper - shape.bounds.width() - g.gap) / step)
                .floor()
                .max(0.0) as usize;
            while column <= end && column < 500000 && !budget.expired() {
                if let Some((b, _, _)) = &best
                    && column as f64 * step + g.gap + shape.bounds.width()
                        > b.x + g.types[t].orientations[b.orientation].width + g.gap
                {
                    break;
                }
                let Some(y) = find_y(shape, &strip, column, start.max(minimum), maximum, budget)
                else {
                    column += 1;
                    start = minimum;
                    continue;
                };
                let mut accepted = None;
                let mut heights=vec![y];let mut cursor=0;
                while cursor<heights.len() {
                    let event=heights[cursor];cursor+=1;
                    let Some(yy)=find_y(shape,&strip,column,event,maximum,budget) else{continue;};
                    if yy > maximum {
                        continue;
                    }
                    let p = Placement {
                        piece,
                        orientation: o,
                        x: column as f64 * step - shape.origin.x,
                        y: yy - shape.origin.y,
                    };
                    stats.collision_queries += 1;
                    let mut trial = placed.clone();
                    trial.push(p.clone());
                    if jagua_clear(g, &trial)?
                        && feasible(
                            &world(g, &p),
                            &placed.iter().map(|p| world(g, p)).collect::<Vec<_>>(),
                            g.total_area * 1e-12,
                        )
                    {
                        accepted = Some((p, yy));
                        break;
                    }
                    if heights.len()==1 {
                        stats.backend_events.push(format!("SemiDiscreteBLF original oracle rejected piece={piece} orientation={o} column={column} y={yy}; checking later same-column contact events"));
                        let Some(events)=contact_heights(g,&placed,piece,o,
                            column as f64*step-shape.origin.x,yy-shape.origin.y,
                            maximum-shape.origin.y,&mut original_cache,budget)? else{return Ok(None);};
                        heights.extend(events.into_iter().map(|v|v+shape.origin.y).filter(|v|*v>yy));
                    }
                }
                if let Some((p, y)) = accepted {
                    stats.evaluations += 1;
                    if best.as_ref().is_none_or(|(b, _, _)| {
                        (p.x + g.types[t].orientations[o].width)
                            .total_cmp(&(b.x + g.types[t].orientations[b.orientation].width))
                            .then(
                                (p.y + g.types[t].orientations[o].height).total_cmp(
                                    &(b.y + g.types[t].orientations[b.orientation].height),
                                ),
                            )
                            .then(
                                g.types[t].orientations[o]
                                    .angle
                                    .total_cmp(&g.types[t].orientations[b.orientation].angle),
                            )
                            .is_lt()
                    }) {
                        best = Some((p, column, y));
                    }
                    break;
                }
                // Advance only after the same-column original contact events
                // have also been exhausted; do not lose an available low gap.
                column += 1;
                start = minimum;
            }
        }
        if budget.expired() {
            return Ok(None);
        }
        let Some((p, column, y)) = best else {
            return Ok(None);
        };
        let key = (t, p.orientation);
        place(&mut strip, &shapes[&key], column, y);
        previous.insert(key, (column, y));
        placed.push(p);
    }
    let layout = Layout {
        length: geometry::length(g, &placed),
        placements: placed,
    };
    g.validate(&layout)?;
    Ok(Some(layout))
}
/// A single original constructive pass, not a metaheuristic or timed search.
pub fn solve(g: &Geometry, step: f64, budget: &Budget, stats: &mut Stats) -> Result<Layout> {
    stats.effective_algorithm_controls=Some(serde_json::json!({"resolution":step,"zero_vertex_extensions":true,"interior_gap_extensions":false}));
    let mut order: Vec<_> = (0..g.pieces.len()).collect();
    order.sort_by(|&a, &b| {
        let area = |i: usize| {
            let o = &g.types[g.pieces[i]].orientations[0];
            o.width * o.height
        };
        area(b).total_cmp(&area(a)).then(a.cmp(&b))
    });
    let result = construct(g, &order, step, budget, stats)?;
    anyhow::ensure!(
        result.is_some(),
        "semi-discrete constructive pass incomplete; no shelf relabeling"
    );
    stats.completed_variants = 1;
    stats.paper_variant = Some(format!(
        "A107 fixed XML angles; R={step}; bbox-area descending; input-order ties"
    ));
    stats.paper_run_status = Some("completed_single_constructive_pass".into());
    Ok(result.unwrap())
}
