use anyhow::{Result, ensure};
use packing_bkv_reproductions::{
    backends::Settings, export, geometry::Geometry, input, model::*, three_pm,
};
use std::{path::Path, time::Instant};
fn main() -> Result<()> {
    run(false)
}
pub fn run(gurobi: bool) -> Result<()> {
    run_with_backend(gurobi,false)
}
pub fn run_with_backend(gurobi:bool,finite:bool)->Result<()> {
    run_config(gurobi,finite,false)
}
pub fn run_config(gurobi:bool,finite:bool,smt:bool)->Result<()> {
    let a: Vec<_> = std::env::args().skip(1).collect();
    ensure!(
        a.len() == 3,
        "usage: threepm-bkv INSTANCE.xml WATCHDOG_SECONDS OUTPUT.json"
    );
    let input = input::load(Path::new(&a[0]))?;
    let start = Instant::now();
    let budget = Budget::new(a[1].parse()?);
    let g = Geometry::from_input(&input)?;
    let mut stats = Stats::default();
    let mut settings=Settings::default();
    settings.finite_dbm=finite;
    if smt {settings.smt_python=Some(std::env::var_os("BKV_PYTHON_BIN").map(Into::into)
        .unwrap_or_else(||if cfg!(windows){"F:/anaconda3/python.exe".into()}else{"python3".into()}));}
    if finite {settings.cbc_preprocess=false;settings.cbc_presolve=false;}
    if gurobi {
        settings.gurobi_python=Some(std::env::var_os("BKV_PYTHON_BIN").map(Into::into)
            .unwrap_or_else(||if cfg!(windows){"F:/anaconda3/python.exe".into()}else{"python3".into()}));
        settings.gurobi_threads=std::env::var("BKV_GUROBI_THREADS").unwrap_or_else(|_|"1".into()).parse()?;
        ensure!((1..=6).contains(&settings.gurobi_threads),"Gurobi threads must be between1 and6");
    }
    let layout = three_pm::solve(&g, &budget, &settings, &mut stats)?;
    let mut v = export::output(&input, &g, &layout)?;
    v["algorithm"] = "3PM".into();
    v["seed"] = serde_json::Value::Null;
    v["strict_full_paper_protocol_reproduction"] = false.into();
    v["stats"] = serde_json::to_value(stats)?;
    v["mip_solver"] = settings.solver_name().into();
    v["solve_seconds"] = start.elapsed().as_secs_f64().into();
    v["assumptions"] = packing_bkv_reproductions::config::record()["three_pm"]["decisions"].clone();
    v["backend_configuration"]=serde_json::json!({
        "finite_dbm":finite,"rational_smt":smt,"gurobi":gurobi,
        "cbc_preprocess":settings.cbc_preprocess,"cbc_presolve":settings.cbc_presolve,
        "smt_python":settings.smt_python,"cbc":settings.cbc,
        "solver_threads":if gurobi {settings.gurobi_threads}else{1}});
    let p = Path::new(&a[2]);
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
