use anyhow::{Result, ensure};
use packing_bkv_reproductions::{blf_nfp, export, geometry::Geometry, input, model::*};
use std::{path::Path, time::Instant};
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().skip(1).collect();
    ensure!(
        a.len() == 4,
        "usage: blf_nfp-bkv INSTANCE.xml WATCHDOG_SECONDS SEED OUTPUT.json"
    );
    let input = input::load(Path::new(&a[0]))?;
    let start = Instant::now();
    let budget = Budget::new(a[1].parse()?);
    let g = Geometry::from_input(&input)?;
    let mut stats = Stats::default();
    let layout = blf_nfp::solve(&g, &budget, a[2].parse()?, &mut stats)?;
    let mut v = export::output_for_independent_verification(&input, &g, &layout)?;
    v["strict_geometry_diagnostic"]=serde_json::json!(g.validate(&layout).err().map(|e|e.to_string()));
    v["algorithm"] = "BLF-NFP".into();
    v["seed"] = serde_json::Value::Null;
    v["strict_full_paper_protocol_reproduction"] = false.into();
    v["stats"] = serde_json::to_value(stats)?;
    v["solve_seconds"] = start.elapsed().as_secs_f64().into();
    v["seed"] = a[2].parse::<u64>()?.into();
    v["assumptions"]=serde_json::json!([
        "Whitwell author thesis Chapter9.3 original tabu/area300s variant; report best40.50 is the corresponding Poly3B column",
        "NFP vertical overlap escape with continuous y, same-copy starts, discrete x advancement; no minimum-box preorientation and XML0degrees only",
        "Five neighbors and25-solution tabu history follow thesis9.3, not the200-solution history of the later2006 trigonometric paper",
        "Unreported x increment chosen strip_height/256; unreported1/2/3/4/N-opt weights0.4/0.3/0.15/0.1/0.05; cyclic multi-swap",
        "One run and seed; original40-run multi-variant statistical protocol not claimed; raw original contours exported without repairs"]);
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
