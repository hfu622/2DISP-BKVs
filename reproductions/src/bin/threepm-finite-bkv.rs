// Identical three-phase controller, native exact DBMs plus CBC continuous MIPs.
#[allow(dead_code)]mod cli { include!("threepm-bkv.rs"); }
fn main()->anyhow::Result<()> { cli::run_with_backend(false,true) }
