// Same three-phase controller; exact rational integer-contour CPM solver.
#[allow(dead_code)]mod cli { include!("threepm-bkv.rs"); }
fn main()->anyhow::Result<()> { cli::run_config(false,true,true) }
