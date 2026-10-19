//! Shared original-contour geometry and jagua feasibility oracle.
use crate::{input::LoadedInput, model::*};
use anyhow::{Context, Result, ensure};
use jagua_rs::{
    collision_detection::{
        CDEConfig, CDEngine,
        hazards::filter::NoFilter,
        hazards::{Hazard, HazardEntity},
    },
    geometry::{
        DTransformation, Transformation,
        fail_fast::SPSurrogateConfig,
        geo_enums::GeoPosition,
        geo_traits::Transformable,
        primitives::{Point, Rect, SPolygon},
    },
    io::ext_repr::{ExtRotation, ExtShape},
};
use std::sync::Arc;

#[derive(Clone)]
pub struct Orientation {
    pub angle: f64,
    pub min_x: f64,
    pub min_y: f64,
    pub width: f64,
    pub height: f64,
    pub shape: Arc<SPolygon>,
}
#[derive(Clone)]
pub struct PieceType {
    pub id: u64,
    pub area: f64,
    pub perimeter: f64,
    pub orientations: Vec<Orientation>,
    pub exact_vertices: Vec<(f64, f64)>,
}
#[derive(Clone)]
pub struct Geometry {
    pub height: f64,
    pub types: Vec<PieceType>,
    pub pieces: Vec<usize>,
    pub total_area: f64,
    pub gap: f64,
    pub upper: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Entity {
    Exterior,
    Piece(usize),
}
impl HazardEntity for Entity {
    fn scope(&self) -> GeoPosition {
        if matches!(self, Self::Exterior) {
            GeoPosition::Exterior
        } else {
            GeoPosition::Interior
        }
    }
    fn is_dynamic(&self) -> bool {
        matches!(self, Self::Piece(_))
    }
}

impl Geometry {
    pub fn from_input(input: &LoadedInput) -> Result<Self> {
        let height = input.instance.strip_height as f64;
        let gap = (height * 2e-6).max(2e-6);
        let mut types = Vec::new();
        let mut pieces = Vec::new();
        let mut total_area = 0.0;
        let mut upper = gap;
        for item in &input.instance.items {
            let ExtShape::SimplePolygon(poly) = &item.base.shape else {
                anyhow::bail!("only original simple polygons are supported");
            };
            let points = &poly.0;
            let exact_vertices = input.exact_vertices[&item.base.id].clone();
            let area = (0..exact_vertices.len())
                .map(|i| {
                    let a = exact_vertices[i];
                    let b = exact_vertices[(i + 1) % exact_vertices.len()];
                    a.0 * b.1 - b.0 * a.1
                })
                .sum::<f64>()
                .abs()
                * 0.5;
            let perimeter = (0..points.len())
                .map(|i| {
                    let a = points[i];
                    let b = points[(i + 1) % points.len()];
                    ((a.0 - b.0) as f64).hypot((a.1 - b.1) as f64)
                })
                .sum();
            let shape = SPolygon::new(points.iter().map(|&(x, y)| Point(x, y)).collect())?;
            let ExtRotation::Discrete { angles } = &item.base.orientation.rotation else {
                anyhow::bail!("enumerated rotations are required");
            };
            let mut orientations = Vec::new();
            // Retain the exact XML-derived degree values, including +/-360 degrees.
            for &angle in angles {
                let mut rotated = shape.transform_clone(
                    &DTransformation::new(angle.to_radians(), (0.0, 0.0)).compose(),
                );
                let b = rotated.bbox();
                let width = (b.x_max - b.x_min) as f64;
                let h = (b.y_max - b.y_min) as f64;
                let min_x = b.x_min as f64;
                let min_y = b.y_min as f64;
                rotated.transform(&Transformation::from_translation((-b.x_min, -b.y_min)));
                orientations.push(Orientation {
                    angle: angle as f64,
                    min_x,
                    min_y,
                    width,
                    height: h,
                    shape: Arc::new(rotated),
                });
            }
            ensure!(
                orientations.iter().any(|o| o.height <= height),
                "piece {} cannot fit with its specified rotations",
                item.base.id
            );
            let count =
                usize::try_from(item.demand).context("piece count exceeds addressable memory")?;
            ensure!(
                pieces.len().saturating_add(count) <= 100_000,
                "too many pieces"
            );
            let max_w = orientations.iter().map(|o| o.width).fold(0.0, f64::max);
            // Reserve one grid-cell rounding allowance per piece for raster fallback.
            upper += (max_w + gap + height / 32.0) * count as f64;
            total_area += area * count as f64;
            pieces.extend(std::iter::repeat_n(types.len(), count));
            types.push(PieceType {
                id: item.base.id,
                area,
                perimeter,
                orientations,
                exact_vertices,
            });
        }
        ensure!(
            !pieces.is_empty() && upper.is_finite() && total_area > 0.0,
            "invalid or empty instance"
        );
        let safe_gap = gap.max(upper * 4.0 * f32::EPSILON as f64);
        upper += (safe_gap - gap) * pieces.len() as f64;
        Ok(Self {
            height,
            types,
            pieces,
            total_area,
            gap: safe_gap,
            upper,
        })
    }
    pub fn orientation(&self, p: &Placement) -> &Orientation {
        &self.types[self.pieces[p.piece]].orientations[p.orientation]
    }
    pub fn world(&self, p: &Placement) -> SPolygon {
        self.orientation(p)
            .shape
            .transform_clone(&Transformation::from_translation((p.x as f32, p.y as f32)))
    }
    pub fn engine(&self) -> Result<CDEngine<Entity>> {
        // Manual bounds enforce the original strip. Expanded oracle bounds avoid
        // treating legal strip-edge contact as a polygon intersection.
        let e = self.gap * 8.0;
        let b = Rect::try_new(
            -e as f32,
            -e as f32,
            (self.upper + e) as f32,
            (self.height + e) as f32,
        )?;
        let outer = Arc::new(SPolygon::new(vec![
            Point(b.x_min, b.y_min),
            Point(b.x_max, b.y_min),
            Point(b.x_max, b.y_max),
            Point(b.x_min, b.y_max),
        ])?);
        Ok(CDEngine::new(
            b.inflate_to_square(),
            vec![Hazard::new(Entity::Exterior, outer)],
            CDEConfig {
                quadtree_depth: 5,
                cd_threshold: 96,
                item_surrogate_config: SPSurrogateConfig {
                    n_pole_limits: [(100, 0.0), (20, 0.75), (10, 0.90)],
                    ff_pole_area_ratio: 0.5,
                    n_ff_piers: 0,
                },
            },
        ))
    }
    pub fn fit(&self, p: &Placement) -> bool {
        let o = self.orientation(p);
        p.x.is_finite()
            && p.y.is_finite()
            && p.x >= 0.0
            && p.y >= 0.0
            && p.y + o.height <= self.height + 1e-7
            && p.x + o.width <= self.upper
    }
    pub fn validate(&self, layout: &Layout) -> Result<()> {
        ensure!(
            layout.placements.len() == self.pieces.len(),
            "incomplete placement"
        );
        let mut seen = vec![false; self.pieces.len()];
        let mut engine = self.engine()?;
        let mut precise = Vec::<geo::Polygon<f64>>::new();
        for p in &layout.placements {
            ensure!(
                p.piece < seen.len() && !seen[p.piece],
                "duplicate or invalid piece index"
            );
            seen[p.piece] = true;
            ensure!(
                p.orientation < self.types[self.pieces[p.piece]].orientations.len(),
                "unlisted rotation"
            );
            ensure!(self.fit(p), "piece {} violates strip bounds", p.piece);
            let shape = Arc::new(self.world(p));
            let jagua_collision = engine.detect_poly_collision(&shape, &NoFilter);
            engine.register_hazard(Hazard::new(Entity::Piece(p.piece), shape));
            let o = self.orientation(p);
            let angle = o.angle.to_radians();
            let (sin, cos) = angle.sin_cos();
            let coords = self.types[self.pieces[p.piece]]
                .exact_vertices
                .iter()
                .map(|&(x, y)| geo::Coord {
                    x: cos * x - sin * y + p.x - o.min_x,
                    y: sin * x + cos * y + p.y - o.min_y,
                })
                .collect::<Vec<_>>();
            let tolerance = self.upper.max(self.height) * 2e-6;
            ensure!(
                coords.iter().all(|c| c.x >= -tolerance
                    && c.y >= -tolerance
                    && c.y <= self.height + tolerance
                    && c.x <= layout.length + tolerance),
                "exported contour out of bounds"
            );
            let polygon = geo::Polygon::new(geo::LineString::new(coords), vec![]);
            ensure!(
                !jagua_collision
                    || crate::toolkit::placement::original_contact_clear(self, &polygon, &precise),
                "piece {} has a jagua collision beyond zero-area in-strip contact",
                p.piece
            );
            use geo::{Area, BooleanOps, BoundingRect};
            for other in &precise {
                let a = polygon.bounding_rect().unwrap();
                let b = other.bounding_rect().unwrap();
                if a.min().x < b.max().x
                    && b.min().x < a.max().x
                    && a.min().y < b.max().y
                    && b.min().y < a.max().y
                {
                    ensure!(
                        polygon.intersection(other).unsigned_area() <= self.total_area * 1e-12,
                        "f64 exported contours overlap for piece {}",
                        p.piece
                    );
                }
            }
            precise.push(polygon);
        }
        let actual = length(self, &layout.placements);
        ensure!(
            (layout.length - actual).abs() <= self.gap,
            "reported strip length mismatch"
        );
        ensure!(
            self.total_area / (self.height * layout.length) <= 1.00001,
            "utilization exceeds one"
        );
        Ok(())
    }
    /// Search container length is independent of the initial finite oracle
    /// domain. Enlarge only its horizontal domain; never the fixed strip width,
    /// rotations, polygons, pair-overlap tolerance or placement coordinates.
    pub fn validate_search_layout(&self, layout: &Layout) -> Result<()> {
        if layout.length <= self.upper {
            return self.validate(layout);
        }
        ensure!(layout.length.is_finite() && layout.length>0.
            && layout.length<(f32::MAX as f64)/16.,"invalid search container length");
        let mut expanded=self.clone();
        expanded.upper=layout.length+self.gap;
        expanded.validate(layout)
    }
    pub fn shelf(&self, order: &[usize]) -> Layout {
        let mut x = 0.0;
        let mut placements = Vec::new();
        for &piece in order {
            let os = &self.types[self.pieces[piece]].orientations;
            let orientation = os
                .iter()
                .enumerate()
                .filter(|(_, o)| o.height <= self.height)
                .min_by(|(_, a), (_, b)| a.width.total_cmp(&b.width))
                .unwrap()
                .0;
            let p = Placement {
                piece,
                orientation,
                x,
                y: 0.0,
            };
            x += os[orientation].width + self.gap;
            placements.push(p);
        }
        let length = length(self, &placements);
        Layout { placements, length }
    }
}
pub fn length(g: &Geometry, ps: &[Placement]) -> f64 {
    ps.iter()
        .map(|p| p.x + g.orientation(p).width)
        .fold(0.0, f64::max)
}

pub fn initial_order(g: &Geometry, rule: OrderRule) -> Vec<usize> {
    let mut order: Vec<_> = (0..g.pieces.len()).collect();
    let score = |i: usize| {
        let t = &g.types[g.pieces[i]];
        let o = &t.orientations[0];
        match rule {
            OrderRule::Area => t.area,
            OrderRule::Height => o.height,
            OrderRule::Width => o.width,
            OrderRule::Perimeter => t.perimeter,
            OrderRule::Input | OrderRule::PaperDefined => -(i as f64),
            OrderRule::Compactness => t.area / (o.width * o.height),
        }
    };
    order.sort_by(|&a, &b| score(b).total_cmp(&score(a)).then(a.cmp(&b)));
    order
}

fn candidates(
    g: &Geometry,
    placed: &[Placement],
    o: &Orientation,
    rule: PlacementRule,
    cap: usize,
) -> Vec<(f64, f64)> {
    let step = g.height / 64.0;
    let end = g.height - o.height;
    let mut xs = vec![0.0, length(g, placed) + g.gap];
    let mut ys = vec![0.0, end.max(0.0)];
    for p in placed {
        let q = g.orientation(p);
        xs.extend([p.x, p.x + q.width + g.gap, (p.x - o.width - g.gap).max(0.0)]);
        ys.extend([
            p.y,
            p.y + q.height + g.gap,
            (p.y - o.height - g.gap).max(0.0),
        ]);
    }
    if matches!(rule, PlacementRule::Raster | PlacementRule::SemiDiscrete) {
        for k in 0..=64 {
            ys.push(k as f64 * step);
        }
    }
    if matches!(rule, PlacementRule::Column) {
        xs = vec![0.0, length(g, placed) + g.gap];
    }
    ys.retain(|y| *y >= 0.0 && *y <= end);
    xs.sort_by(f64::total_cmp);
    ys.sort_by(f64::total_cmp);
    xs.dedup_by(|a, b| (*a - *b).abs() < g.gap * 0.25);
    ys.dedup_by(|a, b| (*a - *b).abs() < g.gap * 0.25);
    let mut out = Vec::new();
    // Vertex-contact events approximate, but do not replace, a complete NFP/CFR.
    if matches!(
        rule,
        PlacementRule::Contact
            | PlacementRule::Hull
            | PlacementRule::Gravity
            | PlacementRule::Potential
    ) {
        for p in placed.iter().rev().take(6) {
            let q = g.orientation(p);
            let stride = (q.shape.vertices().len() / 12).max(1);
            let stride_o = (o.shape.vertices().len() / 12).max(1);
            for a in q.shape.vertices().iter().step_by(stride) {
                for b in o.shape.vertices().iter().step_by(stride_o) {
                    let x = p.x + (a.0 - b.0) as f64;
                    let y = p.y + (a.1 - b.1) as f64;
                    for (dx, dy) in [(g.gap, 0.0), (-g.gap, 0.0), (0.0, g.gap), (0.0, -g.gap)] {
                        if x + dx >= 0.0 && y + dy >= 0.0 && y + dy <= end {
                            out.push((x + dx, y + dy));
                        }
                    }
                }
            }
        }
    }
    for x in xs {
        for &y in &ys {
            out.push((x, y));
        }
    }
    if matches!(rule, PlacementRule::Raster) {
        for p in &mut out {
            p.0 = (p.0 / step).ceil() * step;
            p.1 = (p.1 / step).ceil() * step;
        }
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    out.dedup_by(|a, b| (a.0 - b.0).abs() < g.gap * 0.25 && (a.1 - b.1).abs() < g.gap * 0.25);
    out.truncate(cap);
    out
}

pub fn decode(
    g: &Geometry,
    order: &[usize],
    preferred: &[usize],
    rule: PlacementRule,
    budget: &Budget,
    stats: &mut Stats,
) -> Result<Layout> {
    let mut engine = g.engine()?;
    let mut placed = Vec::with_capacity(order.len());
    for &piece in order {
        let os = &g.types[g.pieces[piece]].orientations;
        let mut best: Option<(Placement, f64)> = None;
        let orientations: Vec<_> = if preferred.is_empty() {
            (0..os.len()).collect()
        } else {
            vec![preferred[piece] % os.len()]
        };
        for orientation in orientations {
            let o = &os[orientation];
            if o.height > g.height {
                continue;
            }
            let mut buffer = (*o.shape).clone();
            for (x, y) in candidates(g, &placed, o, rule, 512) {
                if budget.expired() {
                    break;
                }
                let p = Placement {
                    piece,
                    orientation,
                    x,
                    y,
                };
                if !g.fit(&p) {
                    continue;
                }
                use jagua_rs::geometry::geo_traits::TransformableFrom;
                buffer.transform_from(
                    &o.shape,
                    &Transformation::from_translation((x as f32, y as f32)),
                );
                stats.collision_queries += 1;
                if engine.detect_poly_collision(&buffer, &NoFilter) {
                    continue;
                }
                let end = (x + o.width).max(length(g, &placed));
                let score = match rule {
                    PlacementRule::Gravity => end + y * 0.1,
                    PlacementRule::Potential => end + (y + o.height * 0.5) * 0.15,
                    PlacementRule::Hull => end + (x + o.width) * (y + o.height) / g.height * 0.05,
                    _ => end + y * 1e-4,
                };
                if best.as_ref().is_none_or(|(_, s)| score < *s) {
                    best = Some((p, score));
                }
            }
        }
        let p = best.map(|v| v.0).unwrap_or_else(|| {
            let orientation = os
                .iter()
                .enumerate()
                .filter(|(_, o)| o.height <= g.height)
                .min_by(|(_, a), (_, b)| a.width.total_cmp(&b.width))
                .unwrap()
                .0;
            let mut x = length(g, &placed) + g.gap;
            if matches!(rule, PlacementRule::Raster) {
                let step = g.height / 64.0;
                x = (x / step).ceil() * step;
            }
            Placement {
                piece,
                orientation,
                x,
                y: 0.0,
            }
        });
        ensure!(g.fit(&p), "fallback exceeded conservative strip bound");
        engine.register_hazard(Hazard::new(Entity::Piece(piece), Arc::new(g.world(&p))));
        placed.push(p);
    }
    stats.evaluations += 1;
    Ok(Layout {
        length: length(g, &placed),
        placements: placed,
    })
}

pub fn compact(
    g: &Geometry,
    layout: &Layout,
    passes: usize,
    budget: &Budget,
    stats: &mut Stats,
) -> Result<Layout> {
    let mut result = layout.clone();
    for pass in 0..passes {
        for i in 0..result.placements.len() {
            if budget.expired() {
                return Ok(result);
            }
            let current = result.placements[i].clone();
            let mut engine = g.engine()?;
            for (j, p) in result.placements.iter().enumerate() {
                if i != j {
                    engine
                        .register_hazard(Hazard::new(Entity::Piece(p.piece), Arc::new(g.world(p))));
                }
            }
            let horizontal = pass % 2 == 0;
            let max = if horizontal { current.x } else { current.y };
            // Scan before bisection: feasibility along a concave shape's path is not monotone.
            let mut best = current.clone();
            for k in 0..=24 {
                if budget.expired() {
                    break;
                }
                let mut p = current.clone();
                let v = max * k as f64 / 24.0;
                if horizontal {
                    p.x = v;
                } else {
                    p.y = v;
                }
                stats.collision_queries += 1;
                if g.fit(&p) && !engine.detect_poly_collision(&g.world(&p), &NoFilter) {
                    best = p;
                    break;
                }
            }
            result.placements[i] = best;
        }
    }
    result.length = length(g, &result.placements);
    Ok(result)
}
