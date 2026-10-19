use geo::{Area, BooleanOps, Contains, Point};
use packing_bkv_reproductions::{
    model::Budget,
    toolkit::{nfp, placement::*, polygon::*},
};
fn pts(a: &[(f64, f64)]) -> Vec<P> {
    a.iter().map(|&(x, y)| P::new(x, y)).collect()
}
#[test]
fn penetration_preserves_single_pair_isolated_fit_and_slide() {
    use packing_bkv_reproductions::toolkit::penetration;
    let a = pts(&[
        (0., 0.),
        (10., 0.),
        (10., 4.),
        (8., 4.),
        (8., 2.),
        (2., 2.),
        (2., 8.),
        (8., 8.),
        (8., 6.),
        (10., 6.),
        (10., 10.),
        (0., 10.),
    ]);
    for h in [6., 5.] {
        let b = pts(&[(0., 0.), (6., 0.), (6., h), (0., h)]);
        let n = nfp::build(&a, &b, &Budget::new(5.)).unwrap().unwrap();
        let y = if h == 6. { 2. } else { 2.5 };
        let p = P::new(2., y);
        let poly = polygon(&b.iter().map(|q| q.shift(p)).collect::<Vec<_>>());
        assert_eq!(polygon(&a).intersection(&poly).unsigned_area(), 0.);
        assert_eq!(
            penetration::squared(&n, p).0,
            0.,
            "isolated/slide contact h={h}"
        );
        assert!((penetration::squared(&n, P::new(2.1, y)).0.sqrt() - 0.1).abs() < 1e-8);
    }
}
#[test]
fn euclidean_depth_for_diagonal_boundary_differs_from_axis_depth() {
    use packing_bkv_reproductions::toolkit::penetration;
    let a = pts(&[(-1., 0.), (0., -1.), (1., 0.), (0., 1.)]);
    let n = nfp::build(&a, &a, &Budget::new(1.)).unwrap().unwrap();
    let p = P::new(0.5, 0.5);
    assert!((penetration::squared(&n, p).0.sqrt() - 0.5_f64.sqrt()).abs() < 1e-9);
    assert!((penetration::depth(&n, p) - 1.).abs() < 1e-9);
}
#[test]
fn arrangement_preserves_isolated_exact_fit_cavity() {
    let moving = pts(&[(0., 0.), (1., 0.), (1., 1.), (0., 1.)]);
    let fixed = [
        pts(&[(0., 0.), (3., 0.), (3., 1.), (0., 1.)]),
        pts(&[(0., 2.), (3., 2.), (3., 3.), (0., 3.)]),
        pts(&[(0., 1.), (1., 1.), (1., 2.), (0., 2.)]),
        pts(&[(2., 1.), (3., 1.), (3., 2.), (2., 2.)]),
    ];
    let mut groups = Vec::new();
    let mut extra = Vec::new();
    for a in &fixed {
        let n = nfp::build(a, &moving, &Budget::new(1.)).unwrap().unwrap();
        groups.push(n.boundary);
        extra.extend(n.contacts);
    }
    let points = candidates(
        &groups,
        Bounds {
            x0: 0.,
            y0: 0.,
            x1: 2.,
            y1: 2.,
        },
        &extra,
        &[],
        &[],
        &Budget::new(1.),
    )
    .unwrap()
    .unwrap();
    assert!(points.contains(&P::new(1., 1.)));
    let q = polygon(
        &moving
            .iter()
            .map(|p| p.shift(P::new(1., 1.)))
            .collect::<Vec<_>>(),
    );
    assert!(
        fixed
            .iter()
            .all(|p| q.intersection(&polygon(p)).unsigned_area() == 0.)
    );
    let shifted = polygon(
        &moving
            .iter()
            .map(|p| p.shift(P::new(1.001, 1.)))
            .collect::<Vec<_>>(),
    );
    assert!(
        fixed
            .iter()
            .any(|p| shifted.intersection(&polygon(p)).unsigned_area() > 0.)
    );
}
#[test]
fn decomposition_preserves_concave_area_and_reversed_winding() {
    let v = pts(&[
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 1.0),
        (1.0, 1.0),
        (1.0, 4.0),
        (0.0, 4.0),
    ]);
    for p in [v.clone(), v.into_iter().rev().collect()] {
        let triangles = triangulate(&p).unwrap();
        assert!((triangles.iter().map(|t| signed_area(t).abs()).sum::<f64>() - 7.0).abs() < 1e-12);
        for i in 0..triangles.len() {
            for j in i + 1..triangles.len() {
                assert!(
                    polygon(&triangles[i])
                        .intersection(&polygon(&triangles[j]))
                        .unsigned_area()
                        < 1e-12
                );
            }
        }
    }
}
#[test]
fn concave_nfp_matches_direct_overlap_including_an_inaccessible_hole() {
    let a = pts(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 4.0),
        (8.0, 4.0),
        (8.0, 2.0),
        (2.0, 2.0),
        (2.0, 8.0),
        (8.0, 8.0),
        (8.0, 6.0),
        (10.0, 6.0),
        (10.0, 10.0),
        (0.0, 10.0),
    ]);
    let b = pts(&[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0)]);
    let result = nfp::build(&a, &b, &Budget::new(10.0)).unwrap().unwrap();
    assert!(result.region.0.iter().any(|p| !p.interiors().is_empty()));
    assert!(!result.region.contains(&Point::new(3.0, 3.0)));
    for ix in -10..30 {
        for iy in -10..30 {
            let p = P::new(ix as f64 * 0.47 + 0.123, iy as f64 * 0.39 + 0.087);
            let shifted = polygon(&b.iter().map(|q| q.shift(p)).collect::<Vec<_>>());
            let overlap = polygon(&a).intersection(&shifted).unsigned_area() > 1e-10;
            assert_eq!(
                result.region.contains(&Point::new(p.x, p.y)),
                overlap,
                "translation {p:?}"
            );
        }
    }
}
#[test]
fn arrangements_include_edge_intersections_and_extended_stationary_points() {
    let s = Segment {
        a: P::new(0.0, 1.0),
        b: P::new(2.0, 3.0),
    };
    let t = Segment {
        a: P::new(0.0, 3.0),
        b: P::new(2.0, 1.0),
    };
    let points = candidates(
        &[vec![s], vec![t]],
        Bounds {
            x0: 0.0,
            y0: 0.0,
            x1: 4.0,
            y1: 4.0,
        },
        &[],
        &[],
        &[],
        &Budget::new(1.0),
    )
    .unwrap()
    .unwrap();
    assert!(points.contains(&P::new(1.0, 2.0)));
    let points = extended_points(&[s], &[], &[], |p| (p.x - 1.25).powi(2));
    assert!(points.iter().any(|p| (p.x - 1.25).abs() < 1e-12));
}
#[test]
fn nfp_cache_and_deadline_do_not_return_partial_geometry() {
    let a = pts(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
    assert!(nfp::build(&a, &a, &Budget::new(0.0)).unwrap().is_none());
    let mut cache = nfp::Cache::new(1024);
    for i in 0..20 {
        cache
            .get((i, 0, 0, 0), &a, &a, &Budget::new(1.0))
            .unwrap()
            .unwrap();
        assert!(cache.bytes() <= 1024);
    }
}
