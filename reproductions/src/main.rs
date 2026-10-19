use anyhow::{Context, Result, ensure};
use packing_bkv_reproductions::{gcs, geometry::Geometry, input, model::*};
use serde_json::json;
use std::{path::Path, time::Instant};
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    ensure!(
        args.len() == 4,
        "usage: packing-bkv-reproductions INSTANCE.xml SECONDS SEED OUTPUT.json"
    );
    let seconds: f64 = args[1].parse()?;
    let seed: u64 = args[2].parse()?;
    let input = input::load(Path::new(&args[0]))?;
    let started = Instant::now();
    let budget = Budget::new(seconds);
    let g = Geometry::from_input(&input)?;
    let mut stats = Stats::default();
    let working_input = input.opposite_y_frame();
    let working_geometry = Geometry::from_input(&working_input)?;
    let working_layout = gcs::solve(
        &working_geometry,
        &input.instance.name,
        &budget,
        seed,
        &mut stats,
    )?;
    let layout = gcs::reverse_global_y(&g, &working_layout);
    let tight_length = packing_bkv_reproductions::geometry::length(&g, &layout.placements);
    ensure!(
        tight_length <= layout.length + g.height * 1e-9,
        "contours exceed the reported container length"
    );
    let tight_layout = Layout {
        length: tight_length,
        placements: layout.placements.clone(),
    };
    g.validate(&tight_layout)
        .context("final original-contour/Jagua validation")?;
    let placements:Vec<_>=layout.placements.iter().map(|p|{let o=g.orientation(p);json!({"item_id":g.types[g.pieces[p.piece]].id,"transformation":{"rotation":o.angle.to_radians(),"translation":[p.x-o.min_x,p.y-o.min_y]}})}).collect();
    let mut embedded = serde_json::to_value(&input.instance)?;
    for item in embedded["items"].as_array_mut().unwrap() {
        let id = item["id"].as_u64().unwrap();
        item["shape"]["data"] = serde_json::to_value(&input.exact_vertices[&id])?;
    }
    let doc = json!({"algorithm":"GCS zero-clustering Dighe configuration","strict_full_paper_protocol_reproduction":false,"assumptions":gcs::ASSUMPTIONS,"coordinate_frame":"common y-axis reversal during search; every pose restored to canonical original XML contours on export","tight_contour_length":tight_length,"seed":seed,"time_limit_seconds":seconds,"solve_seconds":started.elapsed().as_secs_f64(),"stats":stats,"instance":embedded,"solution":{"strip_width":layout.length,"density":g.total_area/(g.height*layout.length),"layout":{"placed_items":placements}}});
    let output = Path::new(&args[3]);
    if let Some(p) = output.parent() {
        std::fs::create_dir_all(p)?;
    }
    std::fs::write(output, serde_json::to_vec_pretty(&doc)?)?;
    println!(
        "{} utilization={:.8}% length={:.8} elapsed={:.3}s status={:?}",
        input.instance.name,
        100. * g.total_area / (g.height * layout.length),
        layout.length,
        started.elapsed().as_secs_f64(),
        stats.paper_run_status
    );
    Ok(())
}
