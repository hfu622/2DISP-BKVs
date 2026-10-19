use anyhow::{Result, ensure};
use packing_bkv_reproductions::{export, geometry::Geometry, input, model::*, semi_blf};
use std::{path::Path, time::Instant};
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().skip(1).collect();
    ensure!(
        a.len() == 4,
        "usage: semi-bkv INSTANCE.xml RESOLUTION WATCHDOG_SECONDS OUTPUT.json"
    );
    let input = input::load(Path::new(&a[0]))?;
    let start = Instant::now();
    let seconds: f64 = a[2].parse()?;
    let budget = Budget::new(seconds);
    let g = Geometry::from_input(&input)?;
    let mut stats = Stats::default();
    let layout = semi_blf::solve(&g, a[1].parse()?, &budget, &mut stats)?;
    let mut value = export::output(&input, &g, &layout)?;
    value["algorithm"] = "SemiDiscreteBLF".into();
    value["effective_resolution"] = a[1].parse::<f64>()?.into();
    value["strict_full_paper_protocol_reproduction"] = false.into();
    value["seed"] = serde_json::Value::Null;
    value["assumptions"] = packing_bkv_reproductions::config::record()["semi_discrete_blf"]["decisions"].clone();
    value["stats"] = serde_json::to_value(stats)?;
    value["solve_seconds"] = start.elapsed().as_secs_f64().into();
    let p = Path::new(&a[3]);
    if let Some(parent) = p.parent() {
        std::fs::create_dir_all(parent)?;
    }
    std::fs::write(p, serde_json::to_vec_pretty(&value)?)?;
    println!(
        "{} L={} density={} elapsed={}",
        input.instance.name,
        layout.length,
        g.total_area / (g.height * layout.length) * 100.,
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
