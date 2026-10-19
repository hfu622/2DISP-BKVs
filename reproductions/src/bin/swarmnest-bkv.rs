use anyhow::{Result, ensure};
use packing_bkv_reproductions::{geometry::Geometry, input, model::*, swarmnest};
use serde_json::json;
use std::{path::Path, time::Instant};
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    ensure!(
        args.len() == 4,
        "usage: swarmnest-bkv INSTANCE.xml SECONDS SEED OUTPUT.json"
    );
    let seconds: f64 = args[1].parse()?;
    let seed: u64 = args[2].parse()?;
    let input = input::load(Path::new(&args[0]))?;
    ensure!(
        input.instance.name.eq_ignore_ascii_case("han"),
        "only the Han direct-polygon configuration is audited"
    );
    let start = Instant::now();
    let budget = Budget::new(seconds);
    let g = Geometry::from_input(&input)?;
    let mut stats = Stats::default();
    let layout = swarmnest::solve(&g, &budget, seed, &mut stats)?;
    g.validate(&layout)?;
    let mut embedded = serde_json::to_value(&input.instance)?;
    for item in embedded["items"].as_array_mut().unwrap() {
        let id = item["id"].as_u64().unwrap();
        item["shape"]["data"] = serde_json::to_value(&input.exact_vertices[&id])?;
    }
    let placed:Vec<_>=layout.placements.iter().map(|p|{let o=g.orientation(p);json!({"item_id":g.types[g.pieces[p.piece]].id,"transformation":{"rotation":o.angle.to_radians(),"translation":[p.x-o.min_x,p.y-o.min_y]}})}).collect();
    let result = json!({"algorithm":"SwarmNest direct-polygon Han configuration","strict_full_paper_protocol_reproduction":false,"assumptions":swarmnest::ASSUMPTIONS,"seed":seed,"time_limit_seconds":seconds,"solve_seconds":start.elapsed().as_secs_f64(),"stats":stats,"instance":embedded,"solution":{"strip_width":layout.length,"density":g.total_area/(g.height*layout.length),"layout":{"placed_items":placed}}});
    let output = Path::new(&args[3]);
    if let Some(p) = output.parent() {
        std::fs::create_dir_all(p)?;
    }
    std::fs::write(output, serde_json::to_vec_pretty(&result)?)?;
    println!(
        "Han utilization={:.8}% length={:.8} elapsed={:.3}s",
        100. * g.total_area / (g.height * layout.length),
        layout.length,
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
