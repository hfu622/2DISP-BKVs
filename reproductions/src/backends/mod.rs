pub mod cbc;
pub mod gurobi;
pub mod finite_dbm;
pub mod smt;
mod finite_search;
use crate::model::Budget;
use anyhow::{Context, Result};
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Child, Command},
    sync::atomic::{AtomicU64, Ordering},
    time::Duration,
};

static SERIAL: AtomicU64 = AtomicU64::new(0);
struct TempDir(PathBuf);
impl TempDir {
    fn new(base: &Path) -> Result<Self> {
        fs::create_dir_all(base)?;
        let p = base.join(format!(
            "packing-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ));
        // create_dir (not create_dir_all) refuses to reuse someone else's directory.
        fs::create_dir(&p)?;
        Ok(Self(p))
    }
}
impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
struct ManagedChild(Child);
impl Drop for ManagedChild {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}
fn bounded(cmd: &mut Command, budget: &Budget) -> Result<bool> {
    Ok(bounded_outcome(cmd, budget)?.success)
}
#[derive(Debug)]
struct ProcessOutcome {
    success: bool,
    timed_out: bool,
    exit_code: Option<i32>,
}
fn bounded_outcome(cmd: &mut Command, budget: &Budget) -> Result<ProcessOutcome> {
    let mut child = ManagedChild(cmd.spawn().context("cannot launch optimizer")?);
    loop {
        if let Some(status) = child.0.try_wait()? {
            return Ok(ProcessOutcome { success: status.success(), timed_out: false, exit_code: status.code() });
        }
        if budget.expired() {
            return Ok(ProcessOutcome { success: false, timed_out: true, exit_code: None });
        }
        std::thread::sleep(Duration::from_millis(5));
    }
}

pub struct Settings {
    pub finite_dbm: bool,
    pub smt_python: Option<PathBuf>,
    pub cbc_preprocess: bool,
    pub cbc_presolve: bool,
    pub cbc: PathBuf,
    pub temp: PathBuf,
    pub gurobi_python: Option<PathBuf>,
    pub gurobi_threads: u32,
    /// Optimizer-internal fixed-endpoint SAT assistance; no external layouts.
    pub gurobi_cpsat_min_columns: usize,
}
impl Settings {
    pub fn solver_name(&self) -> &'static str {
        if self.smt_python.is_some() { "Rust-CSP/Z3/CBC" } else if self.finite_dbm { "Rust-CSP/CBC" } else if self.gurobi_python.is_some() { "Gurobi" } else { "CBC" }
    }
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            finite_dbm: false,
            smt_python: None,
            cbc_preprocess: true,
            cbc_presolve: true,
            cbc: std::env::var_os("BKV_CBC_BIN").map(PathBuf::from).unwrap_or_else(||
                PathBuf::from(if cfg!(windows) {"F:/packing-deps/cbc-2.10.13/cbc.exe"}else{"cbc"})),
            temp: std::env::var_os("BKV_TMP_DIR").map(PathBuf::from).unwrap_or_else(||
                if cfg!(windows) {PathBuf::from("F:/packing-deps/tmp-bkv")}else{std::env::temp_dir().join("2disp-bkv")}),
            gurobi_python: None,
            gurobi_threads: 1,
            gurobi_cpsat_min_columns: 3000,
        }
    }
}

#[cfg(all(test, unix))]
mod process_tests {
    use super::*;
    #[test]
    fn optimizer_timeout_is_distinct_from_nonzero_exit() {
        let failed = bounded_outcome(Command::new("sh").args(["-c", "exit 7"]), &Budget::new(1.)).unwrap();
        assert!(!failed.success && !failed.timed_out);
        assert_eq!(failed.exit_code, Some(7));
        let timeout = bounded_outcome(Command::new("sleep").arg("2"), &Budget::new(0.03)).unwrap();
        assert!(!timeout.success && timeout.timed_out);
        assert_eq!(timeout.exit_code, None);
    }
}
