//! Cherri et al. three-phase model primitives. No generic family decoder here.
//! Geometry is in bounding-box coordinates; callers supply XML-legal angles only.
use crate::{
    backends::cbc::{Model, Row, Sense, Variable},
    model::Budget,
    toolkit::{placement::Bounds, polygon::*},
};
use anyhow::{Result, ensure};
use geo::{Area, BooleanOps, Contains, Distance, Euclidean};

#[derive(Clone, Debug)]
pub struct Shape {
    pub kind: usize,
    pub orientation: usize,
    pub vertices: Vec<P>,
    /// Numerical pair clearance; zero reproduces the mathematical contact model.
    pub clearance: f64,
}
#[derive(Clone, Copy, Debug)]
pub struct Dot {
    pub shape: usize,
    pub x: f64,
    pub y: f64,
}
#[derive(Clone, Copy, Debug)]
pub struct Pose {
    pub shape: usize,
    pub x: f64,
    pub y: f64,
}
/// Eq.9 window, previous selected dots, and per-type movement counts.
pub struct Movement {
    pub window: Vec<bool>,
    pub previous: Vec<bool>,
    pub alpha: Vec<usize>,
}
fn variable(lower: f64, upper: f64, objective: f64, binary: bool) -> Variable {
    Variable {
        lower,
        upper,
        objective,
        binary,
    }
}
fn row(terms: Vec<(usize, f64)>, sense: Sense, rhs: f64) -> Row {
    Row { terms, sense, rhs }
}
fn validate_shapes(shapes: &[Shape]) -> Result<()> {
    ensure!(!shapes.is_empty() && shapes.len() <= 2048, "A080 shape cap");
    let mut total = 0;
    for s in shapes {
        ensure!(
            s.clearance.is_finite() && s.clearance >= 0.0 && s.clearance <= 1e4,
            "A080 invalid clearance"
        );
        total += s.vertices.len();
        ensure!(
            s.vertices.len() >= 3 && s.vertices.len() <= 256 && total <= 65_536,
            "A080 contour cap"
        );
        ensure!(
            s.vertices.iter().all(|p| p.x.is_finite()
                && p.y.is_finite()
                && p.x.abs() <= 1e9
                && p.y.abs() <= 1e9),
            "A080 invalid contour"
        );
        let b = Bounds::of(&s.vertices);
        ensure!(
            b.x0 >= 0.0 && b.y0 >= 0.0 && signed_area(&s.vertices).abs() > 0.0,
            "A080 expects nonnegative local bounding-box geometry"
        );
    }
    Ok(())
}
/// Algorithm4 verbatim remainder/tie convention; differs from largest-size prose.
pub fn batch_size(n: usize) -> usize {
    batch_size_with_rule(n,&crate::config::get().three_pm.batch_rule)
}
fn batch_size_with_rule(n:usize,rule:&str)->usize {
    if rule=="prose_divisor_priority" && n>5 {
        let upper=12.min(n-1);
        if let Some(s)=(5..=upper).rev().find(|&s|n%s==0){return s;}
        return (5..=upper).max_by_key(|&s|(n%s,s)).unwrap();
    }
    if n <= 5 {
        return n;
    }
    let mut sigma = 5;
    for s in (6..=12.min(n)).rev() {
        if n.is_multiple_of(s) {
            return s;
        }
        if n % s >= n % sigma {
            sigma = s;
        }
    }
    sigma
}
fn spacing(v: &[P], x_axis: bool, minimum: f64) -> f64 {
    let mut a: Vec<_> = v.iter().map(|p| if x_axis { p.x } else { p.y }).collect();
    a.sort_by(f64::total_cmp);
    a.dedup();
    a.windows(2)
        .map(|p| p[1] - p[0])
        .filter(|d| *d > 0.0)
        .fold(f64::INFINITY, f64::min)
        .max(minimum)
}
pub fn spacing_axis(v:&[P],x_axis:bool,minimum:f64)->f64{spacing(v,x_axis,minimum)}
/// First original vertex is the reference, as required by Section2. Dot x/y
/// remain internal bounding-box origins; the actual reference lies at vertex0
/// plus that origin. Figure4's two phases apply to reference coordinates.
pub fn grid(
    shapes: &[Shape],
    height: f64,
    length: f64,
    minimum: f64,
    retained: &[Dot],
    budget: &Budget,
) -> Result<Option<Vec<Dot>>> {
    validate_shapes(shapes)?;
    ensure!(
        [height, length, minimum]
            .iter()
            .all(|x| x.is_finite() && *x > 0.0)
            && height <= 1e9
            && length <= 1e9,
        "A080 grid dimensions"
    );
    let mut dots = Vec::new();
    for (i, s) in shapes.iter().enumerate() {
        if budget.expired() {
            return Ok(None);
        }
        let b = Bounds::of(&s.vertices);
        if b.x1 > length || b.y1 > height {
            continue;
        }
        let dx = spacing(&s.vertices, true, minimum) + s.clearance;
        let dy = spacing(&s.vertices, false, minimum) + s.clearance;
        let reference=s.vertices[0];
        let nx = (length / dx).floor() as usize + 1;
        let ny = (height / dy).floor() as usize + 1;
        if nx.checked_mul(ny).is_none_or(|n| n > 100_000) {
            return Ok(None);
        }
        let mut ys = Vec::new();
        for k in 0..ny {
            for y in [k as f64 * dy-reference.y, height-k as f64 * dy-reference.y] {
                if y >= 0.0 && y + b.y1 <= height {
                    ys.push(y);
                }
            }
        }
        ys.sort_by(f64::total_cmp);
        ys.dedup();
        for x in 0..nx {
            if budget.expired() {
                return Ok(None);
            }
            let x=x as f64*dx-reference.x;
            if x<0. || x+b.x1>length {continue;}
            for &y in &ys {
                dots.push(Dot {
                    shape: i,
                    x,
                    y,
                });
                // A full fine grid is later restricted to the md=3000
                // neighbourhood. Its preparation cap need not equal the cap
                // on variables of an actual dotted-board MIP.
                if dots.len() > 300000 {
                    return Ok(None);
                }
            }
        }
    }
    for &d in retained {
        ensure!(
            d.shape < shapes.len() && d.x.is_finite() && d.y.is_finite(),
            "invalid retained dot"
        );
        let b = Bounds::of(&shapes[d.shape].vertices);
        ensure!(
            d.x >= 0.0 && d.y >= 0.0 && d.x + b.x1 <= length && d.y + b.y1 <= height,
            "retained dot outside IFR"
        );
        if !dots
            .iter()
            .any(|p| p.shape == d.shape && p.x == d.x && p.y == d.y)
        {
            dots.push(d);
        }
        if dots.len() > 300000 {
            return Ok(None);
        }
    }
    Ok(Some(dots))
}
/// DBM Eqs.2-6 plus the partial-demand, movement and fixed-dot Eqs.8-10.
/// Column0 is length; column(d+1) is the type/angle/dot binary decision.
#[allow(clippy::too_many_arguments)] // Separate mathematical model inputs, not CLI options.
pub fn dotted(
    shapes: &[Shape],
    dots: &[Dot],
    height: f64,
    upper: f64,
    demand: &[usize],
    fixed: &[usize],
    movement: Option<&Movement>,
    budget: &Budget,
) -> Result<Option<Model>> {
    validate_shapes(shapes)?;
    ensure!(
        height.is_finite()
            && height > 0.0
            && upper.is_finite()
            && upper > 0.0
            && !demand.is_empty()
            && dots.len() <= 24000,
        "A080 dotted dimensions"
    );
    ensure!(
        shapes.iter().all(|s| s.kind < demand.len()),
        "A080 demand type"
    );
    let mut m = Model::default();
    m.variables.push(variable(0.0, upper, 1.0, false));
    let mut polygons = Vec::new();
    let mut boxes = Vec::new();
    for d in dots {
        if budget.expired() {
            return Ok(None);
        }
        ensure!(
            d.shape < shapes.len()
                && d.x.is_finite()
                && d.y.is_finite()
                && d.x >= 0.0
                && d.y >= 0.0,
            "A080 invalid dot"
        );
        let v: Vec<_> = shapes[d.shape]
            .vertices
            .iter()
            .map(|p| p.shift(P::new(d.x, d.y)))
            .collect();
        let b = Bounds::of(&v);
        ensure!(b.x1 <= upper && b.y1 <= height, "A080 dot outside IFR");
        boxes.push(b);
        polygons.push(polygon(&v));
        m.variables.push(variable(0.0, 1.0, 0.0, true));
        m.rows.push(row(
            vec![(0, -1.0), (m.variables.len() - 1, b.x1)],
            Sense::Less,
            0.0,
        ));
    }
    for (t, &q) in demand.iter().enumerate() {
        let terms: Vec<_> = dots
            .iter()
            .enumerate()
            .filter(|(_, d)| shapes[d.shape].kind == t)
            .map(|(i, _)| (i + 1, 1.0))
            .collect();
        if terms.is_empty() {
            if q > 0 {
                return Ok(None);
            }
        } else {
            m.rows.push(row(terms, Sense::Equal, q as f64));
        }
    }
    for &d in fixed {
        ensure!(d < dots.len(), "invalid fixed dot");
        m.variables[d + 1].lower = 1.0;
    }
    if let Some(w) = movement {
        ensure!(
            w.window.len() == dots.len()
                && w.previous.len() == dots.len()
                && w.alpha.len() == demand.len(),
            "A080 movement dimensions"
        );
        // Eq.9 is indexed by type AND rotation; alpha_t remains per type.
        for (s, shape) in shapes.iter().enumerate() {
            let mut old = 0;
            let mut terms = Vec::new();
            for (i, d) in dots.iter().enumerate() {
                if d.shape == s && w.window[i] {
                    if w.previous[i] {
                        old += 1;
                    }
                    terms.push((i + 1, if w.previous[i] { -1.0 } else { 1.0 }));
                }
            }
            if !terms.is_empty() {
                m.rows.push(row(
                    terms,
                    Sense::Less,
                    w.alpha[shape.kind] as f64 - old as f64,
                ));
            }
        }
    }
    let mut cache=crate::toolkit::nfp::Cache::new(32*1024*1024);
    for i in 0..dots.len() {
        for j in i + 1..dots.len() {
            if budget.expired() || m.rows.len() >= 2_000_000 {
                return Ok(None);
            }
            let kind = shapes[dots[i].shape].kind;
            if kind == shapes[dots[j].shape].kind && demand[kind] <= 1 {
                // The type-demand equality already forbids selecting both.
                // Omit this logically redundant conflict without changing DBM.
                continue;
            }
            let gap = shapes[dots[i].shape]
                .clearance
                .max(shapes[dots[j].shape].clearance);
            let near = boxes[i].x0 < boxes[j].x1 + gap
                && boxes[j].x0 < boxes[i].x1 + gap
                && boxes[i].y0 < boxes[j].y1 + gap
                && boxes[j].y0 < boxes[i].y1 + gap;
            if !near{continue;}
            let si=dots[i].shape;let sj=dots[j].shape;
            let Some(nfp)=cache.get((si,0,sj,0),&shapes[si].vertices,&shapes[sj].vertices,budget)?else{return Ok(None);};
            let point=geo::Point::new(dots[j].x-dots[i].x,dots[j].y-dots[i].y);
            let mut collision=nfp.region.contains(&point);
            if collision && !nfp.zero_strata.is_empty(){
                collision=polygons[i].intersection(&polygons[j]).unsigned_area()>0.;
            }
            if collision || (gap > 0.0 && Euclidean.distance(&polygons[i], &polygons[j]) < gap * 0.99){
                m.rows
                    .push(row(vec![(i + 1, 1.0), (j + 1, 1.0)], Sense::Less, 1.0));
            }
        }
    }
    if !compress_conflicts(&mut m,budget)?{return Ok(None);}
    m.validate()?;
    Ok(Some(m))
}

fn pair_conflict(m:&Model,r:&Row)->bool{
    matches!(r.sense,Sense::Less)&&r.rhs==1.&&r.terms.len()==2&&r.terms.iter().all(|&(i,a)|a==1.&&m.variables[i].binary)
}
/// Replace covered pairwise binary conflicts by clique inequalities. Every
/// emitted clique has all pair edges; every original edge is covered. Thus
/// integer feasible layouts are unchanged, while the LP bound is strengthened.
/// This is model presolve, not a change of geometric tests or search phases.
fn compress_conflicts(m:&mut Model,budget:&Budget)->Result<bool>{
    let n=m.variables.len();let words=n.div_ceil(64);
    let pairs=m.rows.iter().filter(|r|pair_conflict(m,r)).count();if pairs<128{return Ok(true);}
    ensure!(n<=24001,"conflict compression column cap");
    let mut graph=vec![vec![0u64;words];n];
    let mut add=|i:usize,j:usize|{graph[i][j/64]|=1u64<<(j%64);graph[j][i/64]|=1u64<<(i%64);};
    for r in &m.rows{
        if pair_conflict(m,r){add(r.terms[0].0,r.terms[1].0);}
        // A unit-demand equality supplies valid redundant conflicts too.
        if matches!(r.sense,Sense::Equal)&&r.rhs==1.&&r.terms.iter().all(|&(i,a)|a==1.&&m.variables[i].binary){
            for (k,&(i,_)) in r.terms.iter().enumerate(){for &(j,_) in &r.terms[k+1..]{add(i,j);}}
        }
    }
    let mut uncovered=graph.clone();let mut cliques=Vec::new();
    for i in 0..n{loop{
        if budget.expired(){return Ok(false);}
        let next=uncovered[i].iter().enumerate().find_map(|(w,&bits)|{let masked=if w<i/64{0}else if w==i/64{let shift=i%64+1;if shift==64{0}else{bits&(u64::MAX<<shift)}}else{bits};(masked!=0).then(||w*64+masked.trailing_zeros() as usize)});
        let Some(j)=next else{break;};let mut clique=vec![i,j];
        let mut candidates=graph[i].iter().zip(&graph[j]).map(|(a,b)|a&b).collect::<Vec<_>>();
        for (w,bits) in candidates.iter_mut().enumerate(){if w<i/64{*bits=0;}else if w==i/64{let shift=i%64+1;*bits=if shift==64{0}else{*bits&(u64::MAX<<shift)};}}
        loop{let next=candidates.iter().enumerate().find(|(_,b)|**b!=0).map(|(w,b)|w*64+b.trailing_zeros() as usize);let Some(k)=next else{break;};clique.push(k);for (bits,neighbours) in candidates.iter_mut().zip(&graph[k]){*bits&=*neighbours;}}
        clique.sort_unstable();for (p,&a) in clique.iter().enumerate(){for &b in &clique[p+1..]{uncovered[a][b/64]&=!(1u64<<(b%64));uncovered[b][a/64]&=!(1u64<<(a%64));}}
        cliques.push(row(clique.into_iter().map(|i|(i,1.)).collect(),Sense::Less,1.));
    }}
    let variables=&m.variables;m.rows.retain(|r|!(matches!(r.sense,Sense::Less)&&r.rhs==1.&&r.terms.len()==2&&r.terms.iter().all(|&(i,a)|a==1.&&variables[i].binary)));
    println!("3PM equivalent integer conflict compression: {pairs} pair rows -> {} clique rows",cliques.len());m.rows.extend(cliques);Ok(true)
}

/// Continuous local MIP, Eqs.13-22. Convex-component NFP complements replace
/// unsafe disjunctions over extended boundary edges of a nonconvex union.
/// Columns: L, x0,y0,x1,y1,... then binary facet selectors.
pub fn compact(
    shapes: &[Shape],
    poses: &[Pose],
    height: f64,
    upper: f64,
    lambda: &[f64],
    budget: &Budget,
) -> Result<Option<Model>> {
    validate_shapes(shapes)?;
    ensure!(
        !poses.is_empty()
            && poses.len() <= 512
            && lambda.len() == poses.len()
            && height.is_finite()
            && height > 0.0
            && upper.is_finite()
            && upper > 0.0,
        "A080 compaction dimensions"
    );
    let mut m = Model::default();
    m.variables.push(variable(0.0, upper, 1.0, false));
    let mut components = Vec::new();
    for (i, p) in poses.iter().enumerate() {
        if budget.expired() {
            return Ok(None);
        }
        ensure!(
            p.shape < shapes.len()
                && p.x.is_finite()
                && p.y.is_finite()
                && lambda[i].is_finite()
                && lambda[i] >= 0.0,
            "A080 compaction pose"
        );
        let v = &shapes[p.shape].vertices;
        let b = Bounds::of(v);
        ensure!(
            p.x >= 0.0 && p.y >= 0.0 && p.x + b.x1 <= upper + 1e-7 && p.y + b.y1 <= height + 1e-7,
            "A080 compaction incumbent outside IFR"
        );
        m.variables.push(variable(
            0.0f64.max(p.x - lambda[i]),
            (upper - b.x1).min(p.x + lambda[i]),
            0.0,
            false,
        ));
        m.variables.push(variable(
            0.0f64.max(p.y - lambda[i]),
            (height - b.y1).min(p.y + lambda[i]),
            0.0,
            false,
        ));
        m.rows
            .push(row(vec![(0, -1.0), (1 + 2 * i, 1.0)], Sense::Less, -b.x1));
        components.push(compaction_components(v)?);
    }
    for i in 0..poses.len() {
        for j in i + 1..poses.len() {
            if shapes[poses[i].shape].clearance==0. && shapes[poses[j].shape].clearance==0.{
                let a=Bounds::of(&shapes[poses[i].shape].vertices);let b=Bounds::of(&shapes[poses[j].shape].vertices);
                let ax=&m.variables[1+2*i];let ay=&m.variables[2+2*i];let bx=&m.variables[1+2*j];let by=&m.variables[2+2*j];
                // Entire swept bounding boxes are separated under these exact
                // local variable bounds. All component nonoverlap clauses are
                // then redundant, including closed boundary contact.
                if ax.upper+a.x1<=bx.lower+b.x0 || bx.upper+b.x1<=ax.lower+a.x0
                    || ay.upper+a.y1<=by.lower+b.y0 || by.upper+b.y1<=ay.lower+a.y0{continue;}
            }
            for a in &components[i] {
                for b in &components[j] {
                    if budget.expired() || m.variables.len() > 95_000 || m.rows.len() > 450_000 {
                        return Ok(None);
                    }
                    let sums: Vec<_> = a
                        .iter()
                        .flat_map(|p| b.iter().map(move |q| P::new(p.x - q.x, p.y - q.y)))
                        .collect();
                    let hull = convex_hull(&sums);
                    let mut selectors = Vec::new();
                    for e in edges(&hull) {
                        let dx = e.b.x - e.a.x;
                        let dy = e.b.y - e.a.y;
                        let norm = dx.hypot(dy);
                        ensure!(norm > 0.0, "degenerate convex NFP facet");
                        let gap = shapes[poses[i].shape]
                            .clearance
                            .max(shapes[poses[j].shape].clearance);
                        // Integer XML vertices give integer Minkowski vertices.
                        // Keep the source supporting plane, without normalization
                        // and without reconstructing/rounding coefficients later.
                        let raw_rhs = -dy * e.a.x + dx * e.a.y;
                        let exact_integer = gap == 0.0
                            && [e.a.x,e.a.y,e.b.x,e.b.y,dx,dy,raw_rhs].iter()
                                .all(|a| a.is_finite() && a.fract()==0. && a.abs()<(1u64<<52) as f64)
                            && (e.a.x.abs()*dy.abs()+e.a.y.abs()*dx.abs())<(1u64<<52) as f64;
                        let (nx,ny,rhs) = if exact_integer {(-dy,dx,raw_rhs)}
                            else {(-dy/norm,dx/norm,raw_rhs/norm-gap)};
                        let mut terms = vec![
                            (1 + 2 * i, -nx),
                            (2 + 2 * i, -ny),
                            (1 + 2 * j, nx),
                            (2 + 2 * j, ny),
                        ];
                        terms.retain(|(_, a)| *a != 0.0);
                        terms.sort_by_key(|(v, _)| *v);
                        let max: f64 = terms
                            .iter()
                            .map(|(v, a)| {
                                a * if *a > 0.0 {
                                    m.variables[*v].upper
                                } else {
                                    m.variables[*v].lower
                                }
                            })
                            .sum();
                        // Facet-specific finite M from these exact local variable bounds.
                        // Integral M preserves an integral encoded constant; a
                        // conservative ceiling changes only the inactive facet.
                        let big_m = if exact_integer {
                            (max-rhs).max(0.0).ceil()+1.0
                        } else {(max-rhs).max(0.0)};
                        let z = m.variables.len();
                        m.variables.push(variable(0.0, 1.0, 0.0, true));
                        selectors.push((z, 1.0));
                        if big_m > 0.0 {
                            terms.push((z, big_m));
                        }
                        m.rows.push(row(terms, Sense::Less, rhs + big_m));
                    }
                    // Exactly one exterior facet suffices as a certificate.
                    // Every old feasible coordinate vector can choose one of
                    // its valid facets; redundant binary certificates are
                    // removed without changing physical feasible placements.
                    m.rows.push(row(selectors, Sense::Equal, 1.0));
                }
            }
        }
    }
    m.validate()?;
    Ok(Some(m))
}

/// Keep a convex original contour intact. Cyclic vertex identity with its hull
/// proves this exact fast path; concave contours keep the existing triangles.
/// No approximate hull replacement, vertex movement or clearance is introduced.
fn compaction_components(v:&[P])->Result<Vec<Vec<P>>>{
    let original=clean(v)?;let hull=convex_hull(&original);
    if hull.len()==original.len() && hull.iter().position(|p|*p==original[0]).is_some_and(|start|
        (0..original.len()).all(|k|original[k]==hull[(start+k)%hull.len()])){
        Ok(vec![original])
    }else{Ok(triangulate(&original)?.into_iter().map(|t|t.to_vec()).collect())}
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backends::{Settings, cbc};
    #[test]fn integer_diagonal_supporting_planes_are_direct_and_preserve_geometry(){
        let triangle=Shape{kind:0,orientation:0,clearance:0.,vertices:vec![P::new(0.,0.),P::new(3.,0.),P::new(0.,2.)]};
        let shapes=[triangle];let poses=[Pose{shape:0,x:0.,y:0.},Pose{shape:0,x:3.,y:0.}];
        let m=compact(&shapes,&poses,4.,8.,&[8.,8.],&Budget::new(2.)).unwrap().unwrap();
        for row in &m.rows {assert_eq!(row.rhs.fract(),0.);for (_,a) in &row.terms {assert_eq!(a.fract(),0.);}}
        let facets=m.rows.iter().filter(|r|r.terms.iter().any(|&(i,_)|i>=5) && r.terms.iter().any(|&(i,_)|i<5)).collect::<Vec<_>>();
        assert!(facets.iter().any(|r|r.terms.iter().any(|&(i,a)|i==1 && a.abs()==2.) && r.terms.iter().any(|&(i,a)|i==2 && a.abs()==3.)));
        for ix in -4..=8 {for iy in -4..=8 {
            let dx=ix as f64*0.5;let dy=iy as f64*0.5;
            let outside=facets.iter().any(|r| {
                let selector=r.terms.iter().find(|&&(i,_)|i>=5).unwrap().1;
                let lhs=r.terms.iter().filter(|&&(i,_)|i<5).map(|&(i,a)|a*if i==3{dx}else if i==4{dy}else{0.}).sum::<f64>();
                lhs<=r.rhs-selector
            });
            let a=polygon(&shapes[0].vertices);let b=polygon(&shapes[0].vertices.iter().map(|p|p.shift(P::new(dx,dy))).collect::<Vec<_>>());
            assert_eq!(outside,a.intersection(&b).unsigned_area()==0.,"dx={dx},dy={dy}");
        }}
        assert_eq!(batch_size_with_rule(15,"printed_algorithm_4"),8);
        assert_eq!(batch_size_with_rule(15,"prose_divisor_priority"),5);
    }
    fn square() -> Shape {
        Shape {
            kind: 0,
            orientation: 0,
            clearance: 0.0,
            vertices: vec![
                P::new(0.0, 0.0),
                P::new(1.0, 0.0),
                P::new(1.0, 1.0),
                P::new(0.0, 1.0),
            ],
        }
    }
    #[test]fn exact_convex_fast_path_preserves_all_official_poly1a_contours(){
        let path=std::path::Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap().join("official_esicup/xml/poly1a.xml");
        let input=crate::input::load(&path).unwrap();let g=crate::geometry::Geometry::from_input(&input).unwrap();
        let p=super::super::prepare(&g,&Budget::new(3.)).unwrap().unwrap();
        let mut old=0;let mut new=0;let mut intact=0;
        for shape in &p.shapes{
            let parts=compaction_components(&shape.vertices).unwrap();old+=triangulate(&shape.vertices).unwrap().len();new+=parts.len();
            if parts.len()==1{intact+=1;}
            let original=geo::MultiPolygon(vec![polygon(&shape.vertices)]);let mut union=geo::MultiPolygon(vec![]);
            for part in parts{union=union.union(&geo::MultiPolygon(vec![polygon(&part)]));}
            assert_eq!(original.difference(&union).unsigned_area(),0.);
            assert_eq!(union.difference(&original).unsigned_area(),0.);
        }
        assert_eq!(intact,10);assert!(new<old);
    }
    #[test]fn only_guaranteed_separated_motion_envelopes_skip_nonoverlap_clauses(){
        let s=[square()];let p=[Pose{shape:0,x:0.,y:0.},Pose{shape:0,x:0.,y:3.}];
        let m=compact(&s,&p,5.,3.,&[0.25;2],&Budget::new(2.)).unwrap().unwrap();
        assert_eq!(m.variables.len(),5);
        let p=[p[0],Pose{shape:0,x:0.,y:1.5}];
        let m=compact(&s,&p,5.,3.,&[1.;2],&Budget::new(2.)).unwrap().unwrap();
        assert!(m.variables.len()>5);
    }
    #[test]
    fn conflict_cliques_preserve_every_binary_feasible_assignment(){
        for salt in [0,3]{let mut original=Model::default();original.variables=vec![variable(0.,1.,0.,true);18];
            for i in 0..18{for j in i+1..18{if (i+j+salt)%13!=0{original.rows.push(row(vec![(i,1.),(j,1.)],Sense::Less,1.));}}}
            assert!(original.rows.len()>128);let mut compressed=original.clone();assert!(compress_conflicts(&mut compressed,&Budget::new(5.)).unwrap());assert!(compressed.rows.len()<original.rows.len());
            for mask in 0u32..1<<18{let values=(0..18).map(|i|((mask>>i)&1) as f64).collect::<Vec<_>>();assert_eq!(original.feasible(&values),compressed.feasible(&values),"mask={mask}");}
        }
    }
    #[test]fn fine_grid_is_prepared_before_the_global_neighbourhood_cap(){
        let mut shapes=vec![square();15];for (i,s) in shapes.iter_mut().enumerate(){s.kind=i;}
        let dots=grid(&shapes,40.,45.,0.5,&[],&Budget::new(3.)).unwrap().unwrap();
        assert!(dots.len()>24000&&dots.len()<300000);
    }
    #[test]
    fn meshes_keep_board_top_phase_and_sigma_ties() {
        let s = Shape {
            kind: 0,
            orientation: 0,
            clearance: 0.0,
            vertices: vec![
                P::new(0.0, 0.0),
                P::new(10.0, 6.0),
                P::new(4.0, 12.0),
                P::new(0.0, 12.0),
            ],
        };
        let d = grid(&[s], 20.0, 22.0, 2.0, &[], &Budget::new(2.0))
            .unwrap()
            .unwrap();
        let ys: Vec<_> = d.iter().filter(|p| p.x == 0.0).map(|p| p.y).collect();
        assert_eq!(ys, [0.0, 2.0, 6.0, 8.0]);
        assert!(d.iter().any(|p| p.x == 4.0));
        assert_eq!(
            (
                batch_size(5),
                batch_size(13),
                batch_size(15),
                batch_size(16)
            ),
            (5, 7, 8, 8)
        );
        assert!(
            grid(&[square()], 2.0, 1e9, 0.001, &[], &Budget::new(1.0))
                .unwrap()
                .is_none()
        );
    }
    #[test]
    fn vertex_reference_shifts_the_mesh_without_changing_the_contour() {
        let shape=Shape{kind:0,orientation:0,clearance:0.,vertices:vec![P::new(1.,0.),P::new(2.,2.),P::new(0.,2.)]};
        let points=grid(&[shape],5.,6.,2.,&[],&Budget::new(1.)).unwrap().unwrap();
        assert!(!points.is_empty());
        assert!(points.iter().all(|p|((p.x+1.)/2.).fract().abs()<1e-12));
        assert!(points.iter().all(|p|p.x>=0. && p.x+2.<=6. && p.y>=0. && p.y+2.<=5.));
        assert!(points.iter().any(|p|p.x==1. && p.y==0.));
        assert!(!points.iter().any(|p|p.x==0.));
    }
    #[test]
    fn actual_dotted_demand_collision_fixed_and_hamming_models() {
        let s = vec![square()];
        let cfg = Settings::default();
        let d = vec![
            Dot {
                shape: 0,
                x: 0.0,
                y: 0.0,
            },
            Dot {
                shape: 0,
                x: 0.5,
                y: 0.0,
            },
            Dot {
                shape: 0,
                x: 1.0,
                y: 0.0,
            },
        ];
        let m = dotted(&s, &d, 1.0, 3.0, &[2], &[0], None, &Budget::new(2.0))
            .unwrap()
            .unwrap();
        let sol = cbc::solve(&m, &cfg, &Budget::new(3.0)).unwrap().unwrap();
        assert_eq!(sol.values, [2.0, 1.0, 0.0, 1.0]);
        let w = Movement {
            window: vec![true; 3],
            previous: vec![true, false, false],
            alpha: vec![0],
        };
        let m = dotted(&s, &d, 1.0, 3.0, &[2], &[], Some(&w), &Budget::new(2.0))
            .unwrap()
            .unwrap();
        assert!(cbc::solve(&m, &cfg, &Budget::new(3.0)).unwrap().is_none());
    }
    #[test]
    fn actual_disjunctive_compaction_moves_and_excludes_concave_arms() {
        let cfg = Settings::default();
        let s = vec![square()];
        let p = [
            Pose {
                shape: 0,
                x: 0.0,
                y: 0.0,
            },
            Pose {
                shape: 0,
                x: 2.0,
                y: 0.0,
            },
        ];
        let m = compact(&s, &p, 1.0, 3.0, &[3.0; 2], &Budget::new(2.0))
            .unwrap()
            .unwrap();
        let sol = cbc::solve(&m, &cfg, &Budget::new(3.0)).unwrap().unwrap();
        assert!((sol.values[0] - 2.0).abs() < 1e-8);
        // A square pinned in the solid arm of an L shape must be infeasible.
        let l = Shape {
            kind: 1,
            orientation: 0,
            clearance: 0.0,
            vertices: vec![
                P::new(0.0, 0.0),
                P::new(3.0, 0.0),
                P::new(3.0, 1.0),
                P::new(1.0, 1.0),
                P::new(1.0, 3.0),
                P::new(0.0, 3.0),
            ],
        };
        let s = vec![square(), l];
        let p = [
            Pose {
                shape: 1,
                x: 0.0,
                y: 0.0,
            },
            Pose {
                shape: 0,
                x: 1.0,
                y: 0.0,
            },
        ];
        let m = compact(&s, &p, 3.0, 4.0, &[0.0; 2], &Budget::new(2.0))
            .unwrap()
            .unwrap();
        assert!(cbc::solve(&m, &cfg, &Budget::new(3.0)).unwrap().is_none());
        let p = [
            p[0],
            Pose {
                shape: 0,
                x: 1.0,
                y: 1.0,
            },
        ];
        let m = compact(&s, &p, 3.0, 4.0, &[0.0; 2], &Budget::new(2.0))
            .unwrap()
            .unwrap();
        assert!(cbc::solve(&m, &cfg, &Budget::new(3.0)).unwrap().is_some());
    }
}
