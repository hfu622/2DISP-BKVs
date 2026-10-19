//! Read original ESICUP Nesting XML without using its precomputed NFP/IFP data.
use anyhow::{Context, Result, bail, ensure};
use jagua_rs::io::ext_repr::{ExtItem, ExtOrientation, ExtRotation, ExtSPolygon, ExtShape};
use jagua_rs::probs::spp::io::ext_repr::{ExtItem as SPItem, ExtSPInstance};
use roxmltree::Node;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone)]
pub struct LoadedInput {
    pub instance: ExtSPInstance,
    pub item_labels: BTreeMap<u64, String>,
    pub source_path: PathBuf,
    pub exact_vertices: BTreeMap<u64, Vec<(f64, f64)>>,
    pub source_up_left: bool,
}

impl LoadedInput {
    /// A common coordinate-frame reversal, never an independent item reflection.
    /// Callers must map every placement back before exporting original contours.
    pub fn opposite_y_frame(&self) -> Self {
        let mut out = self.clone();
        for item in &mut out.instance.items {
            if let ExtShape::SimplePolygon(poly) = &mut item.base.shape {
                for (_, y) in &mut poly.0 {
                    *y = -*y;
                }
            }
            if let ExtRotation::Discrete { angles } = &mut item.base.orientation.rotation {
                for a in angles {
                    *a = -*a;
                }
            }
        }
        for points in out.exact_vertices.values_mut() {
            for (_, y) in points {
                *y = -*y;
            }
        }
        out.source_up_left = !out.source_up_left;
        out
    }
}

/// A bare name such as `fu` resolves against the project's Nesting directory,
/// independently of the current working directory. Explicit XML/JSON paths work too.
pub fn load(requested: &Path) -> Result<LoadedInput> {
    let path = if requested.is_file() {
        requested.to_path_buf()
    } else if requested.components().count() == 1 {
        let mut path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../Nesting")
            .join(requested);
        if path.extension().is_none() {
            path.set_extension("xml");
        }
        path
    } else {
        requested.to_path_buf()
    };
    let source_path = path
        .canonicalize()
        .with_context(|| format!("cannot open {}", path.display()))?;
    let text = fs::read_to_string(&source_path)?;
    let (instance, item_labels) = match path.extension().and_then(|e| e.to_str()) {
        Some(extension) if extension.eq_ignore_ascii_case("xml") => parse_xml(&text)?,
        Some(extension) if extension.eq_ignore_ascii_case("json") => {
            let instance: ExtSPInstance = serde_json::from_str(&text)?;
            let labels = instance
                .items
                .iter()
                .map(|i| (i.base.id, i.base.id.to_string()))
                .collect();
            (instance, labels)
        }
        _ => bail!("input must be a Nesting XML or jagua SPP JSON file"),
    };
    let exact_vertices = if path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("xml"))
    {
        exact_xml_vertices(&text)?
    } else {
        instance
            .items
            .iter()
            .filter_map(|i| {
                if let ExtShape::SimplePolygon(p) = &i.base.shape {
                    Some((
                        i.base.id,
                        p.0.iter().map(|&(x, y)| (x as f64, y as f64)).collect(),
                    ))
                } else {
                    None
                }
            })
            .collect()
    };
    let source_up_left = path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("xml"))
        && roxmltree::Document::parse(&text)?.descendants().any(|n| {
            n.has_tag_name("coordinatesOrigin") && n.text().unwrap_or("").trim() == "up-left"
        });
    Ok(LoadedInput {
        instance,
        item_labels,
        source_path,
        exact_vertices,
        source_up_left,
    })
}

/// Preserve XML f64 coordinates for independent exported-layout verification.
/// Jagua still receives f32 geometry; this copy is never used to widen rotations.
fn exact_xml_vertices(text: &str) -> Result<BTreeMap<u64, Vec<(f64, f64)>>> {
    let doc = roxmltree::Document::parse(text)?;
    let root = doc.root_element();
    let flip = child(root, "coordinatesOrigin")?
        .text()
        .unwrap_or("")
        .trim()
        == "up-left";
    let polygons = child(root, "polygons")?;
    let lot = child(child(root, "problem")?, "lot")?;
    let value = |node: Node<'_, '_>, key: &str| -> Result<f64> {
        let v: f64 = node
            .attribute(key)
            .context("missing coordinate")?
            .trim()
            .parse()?;
        ensure!(v.is_finite(), "nonfinite coordinate");
        Ok(v)
    };
    let mut result = BTreeMap::new();
    for (i, piece) in lot
        .children()
        .filter(|n| n.has_tag_name("piece"))
        .enumerate()
    {
        let component = child(piece, "component")?;
        let id = component
            .attribute("idPolygon")
            .context("missing polygon")?;
        let polygon = polygons
            .children()
            .find(|p| p.has_tag_name("polygon") && p.attribute("id") == Some(id))
            .context("missing polygon definition")?;
        let dx = value(component, "xOffset")?;
        let dy = value(component, "yOffset")?;
        let mut segments = child(polygon, "lines")?
            .children()
            .filter(|n| n.has_tag_name("segment"))
            .collect::<Vec<_>>();
        segments.sort_by_key(|n| {
            n.attribute("n")
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or(usize::MAX)
        });
        let vertices = segments
            .iter()
            .map(|&n| {
                let x = value(n, "x0")? + dx;
                let y = value(n, "y0")? + dy;
                Ok((x, if flip { -y } else { y }))
            })
            .collect::<Result<Vec<_>>>()?;
        result.insert(i as u64, vertices);
    }
    Ok(result)
}

fn child<'a, 'input>(node: Node<'a, 'input>, name: &str) -> Result<Node<'a, 'input>> {
    node.children()
        .find(|n| n.has_tag_name(name))
        .with_context(|| format!("missing <{name}>"))
}

fn number(node: Node<'_, '_>, name: &str) -> Result<f32> {
    let value: f32 = node
        .attribute(name)
        .with_context(|| format!("missing {name}"))?
        .trim()
        .parse()
        .with_context(|| format!("invalid {name}"))?;
    ensure!(value.is_finite(), "{name} must be finite");
    Ok(value)
}

fn close(a: (f32, f32), b: (f32, f32)) -> bool {
    (a.0 - b.0).abs() <= 1e-5 && (a.1 - b.1).abs() <= 1e-5
}

fn contour(piece: Node<'_, '_>, polygons: Node<'_, '_>, flip_y: bool) -> Result<Vec<(f32, f32)>> {
    let components: Vec<_> = piece
        .children()
        .filter(|n| n.has_tag_name("component"))
        .collect();
    ensure!(
        components.len() == 1,
        "only single-component pieces/boards are supported"
    );
    let component = components[0];
    ensure!(
        component.attribute("type") == Some("0"),
        "holes are unsupported"
    );
    let id = component
        .attribute("idPolygon")
        .context("missing polygon reference")?;
    let polygon = polygons
        .children()
        .find(|n| n.has_tag_name("polygon") && n.attribute("id") == Some(id))
        .with_context(|| format!("unknown polygon {id}"))?;
    let lines = child(polygon, "lines")?;
    let mut segments = lines
        .children()
        .filter(|n| n.has_tag_name("segment"))
        .map(|n| {
            let index: usize = n.attribute("n").context("missing segment index")?.parse()?;
            Ok((
                index,
                (number(n, "x0")?, number(n, "y0")?),
                (number(n, "x1")?, number(n, "y1")?),
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    segments.sort_by_key(|s| s.0);
    ensure!(
        segments.len() >= 3,
        "polygon {id} has fewer than three edges"
    );
    if let Some(count) = polygon.attribute("nVertices") {
        ensure!(
            count.parse::<usize>()? == segments.len(),
            "polygon {id} vertex count mismatch"
        );
    }
    for i in 0..segments.len() {
        // jakobs2 has a 0.001-unit discrepancy between an edge endpoint and
        // the next start. Treat x0/y0 as authoritative vertices, allowing only
        // this source-coordinate rounding tolerance when checking continuity.
        let end = segments[i].2;
        let next = segments[(i + 1) % segments.len()].1;
        ensure!(
            (end.0 - next.0).abs() <= 0.00101 && (end.1 - next.1).abs() <= 0.00101,
            "polygon {id} segments do not form a closed contour"
        );
    }
    let dx = number(component, "xOffset")?;
    let dy = number(component, "yOffset")?;
    Ok(segments
        .iter()
        .map(|(_, (x, y), _)| (x + dx, if flip_y { -(y + dy) } else { y + dy }))
        .collect())
}

fn parse_xml(text: &str) -> Result<(ExtSPInstance, BTreeMap<u64, String>)> {
    let doc = roxmltree::Document::parse(text).context("invalid Nesting XML")?;
    let root = doc.root_element();
    ensure!(root.has_tag_name("nesting"), "expected <nesting> root");
    let origin = child(root, "coordinatesOrigin")?;
    let flip_y = match origin.text().unwrap_or("").trim() {
        "up-left" => true,
        "down-left" | "bottom-left" => false,
        other => bail!("unsupported coordinatesOrigin: {other}"),
    };
    let problem = child(root, "problem")?;
    let polygons = child(root, "polygons")?;
    let boards: Vec<_> = child(problem, "boards")?
        .children()
        .filter(|n| n.has_tag_name("piece"))
        .collect();
    ensure!(boards.len() == 1, "SPP requires one board type");
    let board = contour(boards[0], polygons, flip_y)?;
    let xmin = board.iter().map(|p| p.0).fold(f32::INFINITY, f32::min);
    let xmax = board.iter().map(|p| p.0).fold(f32::NEG_INFINITY, f32::max);
    let ymin = board.iter().map(|p| p.1).fold(f32::INFINITY, f32::min);
    let ymax = board.iter().map(|p| p.1).fold(f32::NEG_INFINITY, f32::max);
    ensure!(
        board.len() == 4
            && xmin < xmax
            && ymin < ymax
            && [(xmin, ymin), (xmin, ymax), (xmax, ymin), (xmax, ymax)]
                .iter()
                .all(|corner| board.iter().any(|p| close(*p, *corner))),
        "SPP board must be an axis-aligned rectangle"
    );
    let mut items = Vec::new();
    let mut labels = BTreeMap::new();
    for (index, piece) in child(problem, "lot")?
        .children()
        .filter(|n| n.has_tag_name("piece"))
        .enumerate()
    {
        let label = piece
            .attribute("id")
            .context("missing piece id")?
            .to_owned();
        ensure!(
            !labels.values().any(|s| s == &label),
            "duplicate piece id {label}"
        );
        let demand = piece
            .attribute("quantity")
            .context("missing quantity")?
            .parse()?;
        let orientation = child(piece, "orientation")?;
        let angles = orientation
            .children()
            .filter(|n| n.is_element())
            .map(|n| {
                ensure!(
                    n.has_tag_name("enumeration"),
                    "only enumerated XML rotations supported"
                );
                let angle = number(n, "angle")?;
                // Conjugating the rotation by a y-axis sign change negates the angle.
                Ok(if flip_y { -angle } else { angle })
            })
            .collect::<Result<Vec<_>>>()?;
        ensure!(!angles.is_empty(), "piece {label} has no allowed rotations");
        let id = index as u64;
        items.push(SPItem {
            base: ExtItem {
                id,
                orientation: ExtOrientation {
                    rotation: ExtRotation::Discrete { angles },
                    reflection_axes: Vec::new(),
                },
                shape: ExtShape::SimplePolygon(ExtSPolygon(contour(piece, polygons, flip_y)?)),
                min_quality: None,
            },
            demand,
        });
        labels.insert(id, label);
    }
    Ok((
        ExtSPInstance {
            name: child(root, "name")?
                .text()
                .unwrap_or("unnamed")
                .trim()
                .to_owned(),
            min_item_separation: 0.0,
            strip_height: ymax - ymin,
            items,
        },
        labels,
    ))
}
