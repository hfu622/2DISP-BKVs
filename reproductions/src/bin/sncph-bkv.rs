use anyhow::{Result, ensure};
use packing_bkv_reproductions::{sncph, export, geometry::Geometry, input, model::*};
use std::{path::Path, time::Instant};
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().skip(1).collect();
    ensure!(
        a.len() == 4,
        "usage: sncph-bkv INSTANCE.xml WATCHDOG_SECONDS SEED OUTPUT.json"
    );
    let input = input::load(Path::new(&a[0]))?;
    let start = Instant::now();
    let budget = Budget::new(a[1].parse()?);
    let g = Geometry::from_input(&input)?;
    let mut stats = Stats::default();
    let layout = sncph::solve(&g, &budget, a[2].parse()?, &mut stats)?;
    let mut v = export::output_for_independent_verification(&input, &g, &layout)?;
    v["strict_geometry_diagnostic"]=serde_json::json!(g.validate(&layout).err().map(|e|e.to_string()));
    v["algorithm"] = "SNCPH".into();
    v["seed"] = serde_json::Value::Null;
    v["strict_full_paper_protocol_reproduction"] = false.into();
    v["stats"] = serde_json::to_value(stats)?;
    v["solve_seconds"] = start.elapsed().as_secs_f64().into();
    v["requested_seed_unused"] = a[2].parse::<u64>()?.into();
    v["assumptions"] = packing_bkv_reproductions::config::record()["sncph"]["decisions"].clone();
    let p = Path::new(&a[3]);
    if let Some(parent) = p.parent() {
        std::fs::create_dir_all(parent)?;
    }
    std::fs::write(p, serde_json::to_vec_pretty(&v)?)?;
    println!(
        "{} L={} density={} elapsed={}",
        input.instance.name,
        layout.length,
        g.total_area / (g.height * layout.length) * 100.,
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
