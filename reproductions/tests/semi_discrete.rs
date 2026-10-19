use packing_bkv_reproductions::{model::Budget, toolkit::polygon::P};
fn pts(a: &[(f64, f64)]) -> Vec<P> {
    a.iter().map(|&(x, y)| P::new(x, y)).collect()
}
#[test]
fn actual_blf_uses_lowest_continuous_y_before_the_next_column(){
    use packing_bkv_reproductions::{geometry::Geometry,input::LoadedInput,model::Stats,semi_blf};
    use std::collections::BTreeMap;
    let vertices=vec![(0.,0.),(1.,0.),(1.,0.75),(0.,0.75)];
    let instance=serde_json::from_value(serde_json::json!({"name":"blf-compactness","strip_height":1.5,"min_item_separation":0.,"items":[{"id":0,"demand":3,"orientation":{"rotation":{"mode":"discrete","angles":[0.]}},"shape":{"type":"simple_polygon","data":vertices}}]})).unwrap();
    let g=Geometry::from_input(&LoadedInput{instance,item_labels:BTreeMap::new(),source_path:Default::default(),exact_vertices:BTreeMap::from([(0,vertices)]),source_up_left:false}).unwrap();
    let l=semi_blf::construct(&g,&[0,1,2],0.5,&Budget::new(3.),&mut Stats::default()).unwrap().unwrap();
    assert_eq!(l.placements.iter().map(|p|(p.x,p.y)).collect::<Vec<_>>(),[(0.,0.),(0.,0.75),(1.,0.)]);
    assert_eq!(l.length,2.);g.validate(&l).unwrap();
}
#[test]
fn official_poly5b_every_insertion_has_no_earlier_valid_semidiscrete_position(){
    use packing_bkv_reproductions::{geometry::Geometry,input,model::{Placement,Stats},semi_blf,
        toolkit::{nfp,placement::*,polygon::*,semi_discrete::*}};
    let input=input::load(&std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../official_esicup/xml/poly5b.xml")).unwrap();
    let g=Geometry::from_input(&input).unwrap();assert!(g.types.iter().all(|t|t.orientations.len()==1));
    let budget=Budget::new(60.);let layout=semi_blf::solve(&g,0.5,&budget,&mut Stats::default()).unwrap();
    g.validate(&layout).unwrap();
    let mut strip=Strip::new();let mut prior=Vec::<Placement>::new();let mut cache=nfp::Cache::new(64*1024*1024);
    for p in &layout.placements {
        let moving=oriented(&g,p.piece,p.orientation);let shape=discretize(&moving,0.5,&budget).unwrap().unwrap();
        let col=((p.x+shape.origin.x)/0.5).round() as usize;
        assert!((col as f64*0.5-p.x-shape.origin.x).abs()<1e-10);
        let maximum=g.height-shape.bounds.height();
        // Start at column ZERO every time: this independently checks the
        // previous-copy shortcut and the outer BLF decoder, not only validity.
        for column in 0..=col {
            let x=column as f64*0.5-shape.origin.x;
            let mut events=vec![0.,maximum];
            for q in &prior {
                let n=cache.get((g.pieces[q.piece],q.orientation,g.pieces[p.piece],p.orientation),
                    &oriented(&g,q.piece,q.orientation),&moving,&budget).unwrap().unwrap();
                for &e in n.boundary.iter().chain(&n.zero_strata) {
                    if let Some(v)=vertical(e,x-q.x){events.push(q.y+v.y+shape.origin.y);}
                    if e.a.x==x-q.x {events.push(q.y+e.a.y+shape.origin.y);}
                    if e.b.x==x-q.x {events.push(q.y+e.b.y+shape.origin.y);}
                }
            }
            events.sort_by(f64::total_cmp);events.dedup();
            for event in events.into_iter().filter(|y|*y>=0. && *y<=maximum) {
                let Some(y)=find_y(&shape,&strip,column,event,maximum,&budget) else{continue;};
                if column==col && y>=p.y+shape.origin.y-1e-10 {continue;}
                let candidate=Placement{piece:p.piece,orientation:p.orientation,x,y:y-shape.origin.y};
                let original_clear=feasible(&world(&g,&candidate),&prior.iter().map(|q|world(&g,q)).collect::<Vec<_>>(),g.total_area*1e-12);
                let mut trial=prior.clone();trial.push(candidate.clone());
                assert!(!original_clear || !jagua_clear(&g,&trial).unwrap(),"missed compact placement of piece {}: chosen ({},{}), earlier ({},{})",p.piece,p.x,p.y,candidate.x,candidate.y);
            }
        }
        place(&mut strip,&shape,col,p.y+shape.origin.y);prior.push(p.clone());
    }
    assert!(!budget.expired(),"compactness audit must exhaust all insertion positions");
}
#[test]
fn vertex_bridge_is_one_sided_and_retains_legal_boundary_contact(){
    use packing_bkv_reproductions::toolkit::semi_discrete::*;
    let triangle=pts(&[(0.,0.),(1.,2.),(1.,3.)]);
    let shape=discretize(&triangle,1.,&Budget::new(2.)).unwrap().unwrap();
    assert!(shape.lines[0].iter().any(|s|s.bottom==0. && s.top==2. && s.label==Label::Right));
    assert!(shape.lines[1].iter().all(|s|s.bottom>=2. && s.label==Label::Left));
    assert!(!incompatible(Interval::new(0.,2.,Label::Right),Interval::new(0.,2.,Label::Left)));
    // A symmetric M bridge would incorrectly occupy y=0..2 on the right line.
    let thin=pts(&[(0.,0.),(1e-10,0.),(1e-10,1.),(0.,1.)]);
    let s=discretize(&thin,0.5,&Budget::new(1.)).unwrap().unwrap();
    assert!(s.lines[0].iter().any(|s|s.label==Label::Right && s.top==1.));
}

#[test]
fn continuous_y_decoder_finds_the_first_hole_after_multiple_collisions(){
    use packing_bkv_reproductions::toolkit::semi_discrete::*;
    let shape=discretize(&pts(&[(0.,0.),(1.,0.),(1.,0.75),(0.,0.75)]),0.5,&Budget::new(1.)).unwrap().unwrap();
    let mut strip=Strip::new();
    // First y update clears one line but enters a different line's obstacle.
    strip.insert(0,vec![Interval::new(0.,0.5,Label::Middle),Interval::new(2.,4.,Label::Middle)]);
    strip.insert(1,vec![Interval::new(0.5,1.25,Label::Middle)]);
    assert_eq!(find_y(&shape,&strip,0,0.,5.,&Budget::new(1.)),Some(1.25));
    assert_eq!(find_y(&shape,&strip,0,1.3,5.,&Budget::new(1.)),Some(4.));
}
#[test]
fn semi_discrete_side_labels_continuous_y_and_vertex_extensions_are_preserved() {
    use packing_bkv_reproductions::toolkit::semi_discrete::*;
    let r = Interval::new(0.0, 2.0, Label::Right);
    let l = Interval::new(1.0, 3.0, Label::Left);
    assert!(!incompatible(r, l));
    let normalized = normalize(&[r, l]);
    assert!(
        normalized
            .iter()
            .any(|i| i.bottom == 1.0 && i.top == 2.0 && i.label == Label::Middle)
    );
    assert!(incompatible(
        Interval::new(1.5, 1.5, Label::Middle),
        normalized[1]
    ));
    let triangle = pts(&[(0.0, 0.0), (2.6, 0.0), (0.4, 2.3)]);
    let shape = discretize(&triangle, 1.0, &Budget::new(1.0))
        .unwrap()
        .unwrap();
    assert!(shape.lines[0].iter().any(|i| i.top >= 2.3));
    assert!(shape.lines[1].iter().any(|i| i.top >= 2.3));
    let rectangle = pts(&[(0.0, 0.0), (2.0, 0.0), (2.0, 0.75), (0.0, 0.75)]);
    let shape = discretize(&rectangle, 1.0, &Budget::new(1.0))
        .unwrap()
        .unwrap();
    let mut strip = Strip::new();
    place(&mut strip, &shape, 0, 0.0);
    let y = find_y(&shape, &strip, 0, 0.0, 4.0, &Budget::new(1.0)).unwrap();
    assert_eq!(y, 0.75); // Continuous coordinate, not rounded to a y-grid.
    assert_eq!(
        find_y(&shape, &strip, 2, 0.0, 4.0, &Budget::new(1.0)),
        Some(0.0)
    );
}
