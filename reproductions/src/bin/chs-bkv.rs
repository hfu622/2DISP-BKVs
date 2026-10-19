use anyhow::{Result, ensure};
use packing_bkv_reproductions::{chs, export, geometry::Geometry, input, model::*};
use std::{path::Path, time::Instant};
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().skip(1).collect();
    ensure!(
        a.len() == 4,
        "usage: chs-bkv INSTANCE.xml WATCHDOG_SECONDS SEED OUTPUT.json"
    );
    let input = input::load(Path::new(&a[0]))?;
    let start = Instant::now();
    let budget = Budget::new(a[1].parse()?);
    let g = Geometry::from_input(&input)?;
    let mut stats = Stats::default();
    let layout = chs::solve(&g, &budget, a[2].parse()?, &mut stats)?;
    let mut v = export::output_for_independent_verification(&input, &g, &layout)?;
    v["strict_geometry_diagnostic"]=serde_json::json!(g.validate(&layout).err().map(|e|e.to_string()));
    v["algorithm"] = "CHS2".into();
    v["seed"] = serde_json::Value::Null;
    v["strict_full_paper_protocol_reproduction"] = false.into();
    v["compactness_diagnostic"]=serde_json::json!({
        "objective":"minimum strip length under the selected HS2 slices",
        "all_insertion_mips_proven_optimal":stats.nonoptimal_backend_calls==0,
        "nonoptimal_insertion_mips":stats.nonoptimal_backend_calls,
        "global_packing_optimum_claimed":false});
    v["stats"] = serde_json::to_value(stats)?;
    v["solve_seconds"] = start.elapsed().as_secs_f64().into();
    v["seed"] = a[2].parse::<u64>()?.into();
    v["assumptions"] = packing_bkv_reproductions::config::record()["chs2"]["decisions"].clone();
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
