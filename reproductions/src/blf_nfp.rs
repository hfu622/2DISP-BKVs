//! Whitwell thesis chapter9 BLF-NFP polygon branch; tabu/area variant.
use crate::{geometry::{self,Geometry},model::*,toolkit::{nfp,placement::*}};
use anyhow::Result;
use std::collections::VecDeque;
/// Weighted 1/2/3/4/N-Opt, retaining every piece exactly once.
pub fn mutate_order(order: &mut Vec<usize>, types: &[usize], operator: usize, rng: &mut Rng) {
    let n = order.len();
    if n < 2 {
        return;
    }
    if operator == 0 {
        let a = rng.index(n);
        let b = rng.index(n);
        let item = order.remove(a);
        order.insert(b, item);
        return;
    }
    if operator == 1 {
        let a = rng.index(n);
        let candidates: Vec<_> = (0..n)
            .filter(|&b| types[order[a]] != types[order[b]])
            .collect();
        if !candidates.is_empty() {
            order.swap(a, candidates[rng.index(candidates.len())]);
        }
        return;
    }
    let size = if operator == 4 {
        2 + rng.index(n - 1)
    } else {
        operator + 1
    }
    .min(n);
    let mut indices: Vec<_> = (0..n).collect();
    rng.shuffle(&mut indices);
    indices.truncate(size);
    let values: Vec<_> = indices.iter().map(|&i| order[i]).collect();
    for (k, &i) in indices.iter().enumerate() {
        order[i] = values[(k + 1) % size];
    }
}
fn operator(rng: &mut Rng) -> usize {
    let v = rng.unit();
    if v < 0.4 {
        0
    } else if v < 0.7 {
        1
    } else if v < 0.85 {
        2
    } else if v < 0.95 {
        3
    } else {
        4
    }
}
fn construct(g:&Geometry,order:&[usize],budget:&Budget,cache:&mut nfp::Cache,stats:&mut Stats)->Result<Option<Layout>>{
    let mut placed=Vec::<Placement>::new();let step=g.height/256.;
    for &piece in order{
        let t=g.pieces[piece];let mut best=None::<Placement>;
        for orientation in 0..g.types[t].orientations.len(){
            let moving=oriented(g,piece,orientation);let o=&g.types[t].orientations[orientation];
            if o.height>g.height{continue;}
            let mut obstacles=Vec::new();
            for p in &placed{let a=oriented(g,p.piece,p.orientation);
                let Some(n)=cache.get((g.pieces[p.piece],p.orientation,t,orientation),&a,&moving,budget)?else{return Ok(None);};
                obstacles.push((n,p.x,p.y));
            }
            let previous=placed.iter().rev().find(|p|g.pieces[p.piece]==t&&p.orientation==orientation);
            let mut x=previous.map_or(0.,|p|p.x);let mut first_y=previous.map_or(0.,|p|p.y);
            let oracle=FixedPlacementOracle::new(g,&placed)?;
            while x+o.width<=g.upper&&!budget.expired(){
                if best.as_ref().is_some_and(|p|x+o.width>p.x+g.orientation(p).width){break;}
                let mut intervals=Vec::new();
                for (n,px,py) in &obstacles{intervals.extend(nfp::vertical_intervals(n,x-px).into_iter().map(|(a,b)|(a+py,b+py)));}
                intervals.sort_by(|a,b|a.0.total_cmp(&b.0));let mut y=first_y;
                for (a,b) in intervals{if a<y&&y<b{y=b;}}
                if y+o.height<=g.height{
                    let p=Placement{piece,orientation,x,y};stats.collision_queries+=1;
                    if oracle.clear(&p){if best.as_ref().is_none_or(|b|(x+o.width).total_cmp(&(b.x+g.orientation(b).width)).then(y.total_cmp(&b.y)).is_lt()){best=Some(p);}break;}
                }
                x+=step;first_y=0.;
            }
        }
        let Some(p)=best else{return Ok(None);};placed.push(p);
    }
    stats.evaluations+=1;Ok(Some(Layout{length:geometry::length(g,&placed),placements:placed}))
}
pub fn solve(g:&Geometry,budget:&Budget,seed:u64,stats:&mut Stats)->Result<Layout>{
    let mut order=geometry::initial_order(g,OrderRule::Area);let mut cache=nfp::Cache::new(32*1024*1024);let mut rng=Rng::new(seed);
    let Some(mut current)=construct(g,&order,budget,&mut cache,stats)?else{anyhow::bail!("BLF-NFP initial decode incomplete");};
    let mut best=current.clone();let mut memory=VecDeque::from([order.clone()]);let mut iterations=0;
    while !budget.expired(){let op=operator(&mut rng);let mut selected=None::<(Vec<usize>,Layout)>;
        for _ in 0..5{
            let mut trial=None;for _ in 0..100{let mut v=order.clone();mutate_order(&mut v,&g.pieces,op,&mut rng);if !memory.contains(&v){trial=Some(v);break;}}
            let Some(v)=trial else{continue;};let Some(l)=construct(g,&v,budget,&mut cache,stats)?else{break;};
            if l.length<best.length{best=l.clone();stats.accepted_moves+=1;println!("BLF-NFP improvement L={} eval={}",best.length,stats.evaluations);}
            if selected.as_ref().is_none_or(|(_,b)|l.length<b.length){selected=Some((v,l));}
        }
        let Some((v,l))=selected else{break;};order=v;current=l;memory.push_back(order.clone());if memory.len()>25{memory.pop_front();}iterations+=1;
    }
    let _=current;stats.completed_generations=Some(iterations);stats.completed_variants=1;stats.nfp_cache_bytes=cache.bytes();
    stats.paper_variant=Some("Whitwell2005 Chapter9 tabu/area; five neighbors; tabu25; XML-angle adaptation".into());
    stats.paper_run_status=Some("completed_BLF_NFP_tabu_area_budget".into());Ok(best)
}
