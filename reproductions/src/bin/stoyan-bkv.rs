use anyhow::{Result, ensure};
use packing_bkv_reproductions::{export, geometry::Geometry, input, model::*, stoyan2012};
use std::{path::Path, time::Instant};
fn main() -> Result<()> {
    let a = std::env::args().skip(1).collect::<Vec<_>>();
    ensure!(
        a.len() == 4,
        "usage: stoyan-bkv INSTANCE.xml GUARD SEED OUTPUT.json"
    );
    let input = input::load(Path::new(&a[0]))?;
    let start = Instant::now();
    let budget = Budget::new(a[1].parse()?);
    let g = Geometry::from_input(&input)?;
    let mut stats = Stats::default();
    let seed = a[2].parse()?;
    let l = stoyan2012::solve(&g, &budget, seed, &mut stats)?;
    let mut v = export::output_for_independent_verification(&input, &g, &l)?;
    v["strict_geometry_diagnostic"] =
        serde_json::json!(g.validate(&l).err().map(|e| e.to_string()));
    v["algorithm"] = "Stoyan2012-DNS".into();
    v["seed"] = seed.into();
    v["stats"] = serde_json::to_value(stats)?;
    v["strict_full_paper_protocol_reproduction"] = false.into();
    v["solve_seconds"] = start.elapsed().as_secs_f64().into();
    v["assumptions"] = serde_json::json!(stoyan2012::ASSUMPTIONS);
    let path = Path::new(&a[3]);
    if let Some(p) = path.parent() {
        std::fs::create_dir_all(p)?;
    }
    std::fs::write(path, serde_json::to_vec_pretty(&v)?)?;
    println!(
        "{} L={} density={} elapsed={}",
        input.instance.name,
        l.length,
        g.total_area / (g.height * l.length) * 100.,
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
