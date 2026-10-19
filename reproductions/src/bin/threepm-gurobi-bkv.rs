// Reuse the identical controller/CLI; only the selected MIP backend differs.
#[allow(dead_code)]
mod cli { include!("threepm-bkv.rs"); }
fn main()->anyhow::Result<()> { cli::run(true) }
