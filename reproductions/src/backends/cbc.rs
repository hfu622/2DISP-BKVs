//! Bounded sparse LP/MIP execution with original-column, full-precision output.
//! This module supplies a solver protocol, not a paper reproduction.
use super::{Settings, TempDir, bounded_outcome};
use crate::model::Budget;
use anyhow::{Result, ensure};
use std::{
    fmt::Write,
    fs,
    process::{Command, Stdio},
};

pub fn native_columns(bytes: &[u8], rows: usize, columns: usize) -> Result<Vec<f64>> {
    ensure!(bytes.len() >= 16, "Truncated CBC binary header");
    let actual_rows = i32::from_le_bytes(bytes[..4].try_into()?);
    let actual_columns = i32::from_le_bytes(bytes[4..8].try_into()?);
    ensure!(
        actual_rows >= 0
            && actual_columns >= 0
            && actual_rows as usize == rows
            && actual_columns as usize == columns,
        "CBC binary dimensions differ from submitted model"
    );
    let expected = rows
        .checked_add(columns)
        .and_then(|n| n.checked_mul(16))
        .and_then(|n| n.checked_add(16))
        .ok_or_else(|| anyhow::anyhow!("CBC binary size overflow"))?;
    ensure!(
        bytes.len() == expected,
        "CBC binary length differs from dimensions"
    );
    let start = 16 + rows * 16;
    let values: Vec<_> = bytes[start..start + columns * 8]
        .chunks_exact(8)
        .map(|p| f64::from_le_bytes(p.try_into().unwrap()))
        .collect();
    ensure!(
        values.iter().all(|v| v.is_finite()),
        "Non-finite CBC binary columns"
    );
    Ok(values)
}

const CAP: usize = 256 * 1024 * 1024;
#[derive(Clone, Copy, Debug)]
pub enum Sense {
    Less,
    Equal,
    Greater,
}
#[derive(Clone, Debug)]
pub struct Variable {
    pub lower: f64,
    pub upper: f64,
    pub objective: f64,
    pub binary: bool,
}
#[derive(Clone, Debug)]
pub struct Row {
    pub terms: Vec<(usize, f64)>,
    pub sense: Sense,
    pub rhs: f64,
}
#[derive(Clone, Debug, Default)]
pub struct Model {
    pub variables: Vec<Variable>,
    pub rows: Vec<Row>,
}
#[derive(Debug)]
pub struct Solution {
    pub values: Vec<f64>,
    pub objective: f64,
    /// False allows a validated incumbent returned at CBC's internal time limit.
    pub optimal: bool,
    /// Audit optimizer-only DBM strengthening; zero for CBC/non-DBM models.
    pub exact_epigraph_thresholds: usize,
    pub cpsat_probe_count: usize,
    pub native_dbm_nodes: u64,
    pub native_dbm_probes: usize,
}
impl Model {
    pub fn validate(&self) -> Result<()> {
        ensure!(
            !self.variables.is_empty() && self.variables.len() <= 100_000,
            "CBC column cap"
        );
        ensure!(
            !self.rows.is_empty() && self.rows.len() <= 2_000_000,
            "CBC row cap"
        );
        let finite = |x: f64| x.is_finite() && x.abs() <= 1e12;
        for v in &self.variables {
            ensure!(
                finite(v.lower) && finite(v.upper) && finite(v.objective) && v.lower <= v.upper,
                "invalid CBC variable"
            );
            ensure!(
                !v.binary || (v.lower >= 0.0 && v.upper <= 1.0),
                "binary bounds must lie within [0,1]"
            );
        }
        let mut count = 0usize;
        for row in &self.rows {
            ensure!(finite(row.rhs) && !row.terms.is_empty(), "invalid CBC row");
            count = count
                .checked_add(row.terms.len())
                .ok_or_else(|| anyhow::anyhow!("CBC nonzero overflow"))?;
            ensure!(count <= 8_000_000, "CBC nonzero cap");
            let mut previous = None;
            for &(i, a) in &row.terms {
                ensure!(
                    i < self.variables.len() && finite(a) && a != 0.0,
                    "invalid CBC coefficient"
                );
                // Unique sorted indices make the sparse representation unambiguous.
                ensure!(
                    previous.is_none_or(|p| p < i),
                    "CBC row indices must be increasing"
                );
                previous = Some(i);
            }
        }
        Ok(())
    }
    /// Explicit MPS column order prevents LP first-occurrence reordering.
    pub fn mps(&self) -> Result<String> {
        self.validate()?;
        let mut columns = vec![Vec::new(); self.variables.len()];
        for (r, row) in self.rows.iter().enumerate() {
            for &(v, a) in &row.terms {
                columns[v].push((r, a));
            }
        }
        let mut s = String::from("NAME PACKING FREE\nROWS\n N OBJ\n");
        for (r, row) in self.rows.iter().enumerate() {
            let kind = match row.sense {
                Sense::Less => 'L',
                Sense::Equal => 'E',
                Sense::Greater => 'G',
            };
            writeln!(s, " {kind} R{r}")?;
        }
        s.push_str("COLUMNS\n");
        let mut integer = false;
        let mut marker = 0;
        for (v, (variable, terms)) in self.variables.iter().zip(columns).enumerate() {
            if integer != variable.binary {
                writeln!(
                    s,
                    " M{marker} 'MARKER' '{}'",
                    if variable.binary { "INTORG" } else { "INTEND" }
                )?;
                integer = variable.binary;
                marker += 1;
            }
            // Even a zero objective entry declares the original column explicitly.
            writeln!(s, " V{v} OBJ {:.17e}", variable.objective)?;
            for (r, a) in terms {
                writeln!(s, " V{v} R{r} {a:.17e}")?;
            }
            ensure!(s.len() <= CAP, "CBC model byte cap");
        }
        if integer {
            writeln!(s, " M{marker} 'MARKER' 'INTEND'")?;
        }
        s.push_str("RHS\n");
        for (r, row) in self.rows.iter().enumerate() {
            writeln!(s, " RHS1 R{r} {:.17e}", row.rhs)?;
        }
        s.push_str("BOUNDS\n");
        for (v, variable) in self.variables.iter().enumerate() {
            // Integer markers plus explicit bounds encode binary variables.
            // Coin rejects LO/UP records following a BV record for that column.
            writeln!(
                s,
                " LO BND V{v} {:.17e}\n UP BND V{v} {:.17e}",
                variable.lower, variable.upper
            )?;
        }
        s.push_str("ENDATA\n");
        ensure!(s.len() <= CAP, "CBC model byte cap");
        Ok(s)
    }
    /// Check the original model, including integrality, after presolve recovery.
    pub fn feasible(&self, values: &[f64]) -> bool {
        if self.validate().is_err() || values.len() != self.variables.len() {
            return false;
        }
        for (v, &x) in self.variables.iter().zip(values) {
            if !x.is_finite()
                || x < v.lower - 1e-7
                || x > v.upper + 1e-7
                || (v.binary && (x - x.round()).abs() > 1e-7)
            {
                return false;
            }
        }
        for row in &self.rows {
            let lhs: f64 = row.terms.iter().map(|&(i, a)| a * values[i]).sum();
            let tolerance = 1e-7 * (1.0 + row.rhs.abs());
            if !lhs.is_finite()
                || match row.sense {
                    Sense::Less => lhs > row.rhs + tolerance,
                    Sense::Equal => (lhs - row.rhs).abs() > tolerance,
                    Sense::Greater => lhs < row.rhs - tolerance,
                }
            {
                return false;
            }
        }
        true
    }
}

pub fn solve(model: &Model, cfg: &Settings, budget: &Budget) -> Result<Option<Solution>> {
    solve_with_start(model, cfg, budget, None)
}
/// Submission does not certify that CBC used or improved the feasible start.
pub fn solve_with_start(
    model: &Model,
    cfg: &Settings,
    budget: &Budget,
    start: Option<&[f64]>,
) -> Result<Option<Solution>> {
    if cfg.finite_dbm {
        match super::finite_dbm::solve(model,budget,start)? {
            super::finite_dbm::Attempt::Finished(answer)=>return Ok(answer),
            super::finite_dbm::Attempt::NotApplicable=>{},
        }
    }
    if cfg.smt_python.is_some() {
        match super::smt::solve(model,cfg,budget,start)? {
            super::smt::Attempt::Finished(answer)=>return Ok(answer),
            super::smt::Attempt::NotApplicable=>{},
        }
    }
    if cfg.gurobi_python.is_some() {
        return super::gurobi::solve(model, cfg, budget, start);
    }
    if budget.expired() {
        return Ok(None);
    }
    let text = model.mps()?;
    if budget.expired() {
        return Ok(None);
    }
    let tmp = TempDir::new(&cfg.temp)?;
    let input = tmp.0.join("model.mps");
    let output = tmp.0.join("solution.sol");
    let binary = tmp.0.join("solution.bin");
    let diagnostic = std::env::var_os("PACKING_CBC_DIAGNOSTIC_DIR");
    let log = tmp.0.join("cbc.log");
    fs::write(&input, text)?;
    if let Some(values) = start {
        ensure!(
            model.feasible(values),
            "CBC MIP start violates original model"
        );
        let mut text = String::new();
        for (i, x) in values.iter().enumerate() {
            writeln!(text, "{i} V{i} {x:.17e}")?;
        }
        ensure!(text.len() <= CAP, "CBC MIP start byte cap");
        fs::write(tmp.0.join("start.sol"), text)?;
    }
    let search_seconds=(budget.remaining()-if model.rows.len()>50_000 {
        (2.0+model.rows.len() as f64/100_000.0).min(10.0)
    }else{2.0}).max(0.001);
    let mut cmd = Command::new(&cfg.cbc);
    cmd.arg(&input);
    if start.is_some() {
        cmd.args(["-mipstart", "start.sol"]);
    }
    if !cfg.cbc_preprocess {cmd.args(["-preprocess","off"]);}
    if !cfg.cbc_presolve {cmd.args(["-presolve","off"]);}
    // This deployed CBC is a serial build: -threads is not a supported option.
    // CPU affinity and the serial binary, rather than an ignored CLI flag,
    // enforce the single-core protocol.
    cmd.args(["-timeMode","elapsed","-randomSeed","1","-randomCbcSeed","1"])
        .arg("-seconds")
        // Reserve serialization time for the incumbent after the optimizer's
        // search limit. A large model cannot reliably write both full-precision
        // binary columns and its text report in the former fixed0.1s margin.
        .arg(format!("{search_seconds:.6}"))
        .args([
            "-primalTolerance",
            "1e-9",
            "-integerTolerance",
            "1e-9",
            "-solve",
            "-printingOptions",
            "all",
            "-solution",
        ])
        .arg("solution.sol")
        .arg("-saveSolution")
        .arg("solution.bin")
        .arg("-quit")
        .current_dir(&tmp.0)
        .env("OMP_NUM_THREADS", "1")
        .stdout(if diagnostic.is_some() {
            Stdio::from(fs::File::create(&log)?)
        } else {
            Stdio::null()
        })
        .stderr(if diagnostic.is_some() {
            Stdio::from(fs::File::create(tmp.0.join("cbc.stderr.log"))?)
        } else { Stdio::null() });
    let process = bounded_outcome(&mut cmd, budget)?;
    let finished = process.success;
    let answer=(||->Result<(Option<Solution>,&'static str)>{
    if !finished || !output.is_file() || !binary.is_file() {
        return Ok((None,"process_failed_or_serialized_solution_missing"));
    }
    ensure!(
        fs::metadata(&output)?.len() <= CAP as u64 && fs::metadata(&binary)?.len() <= CAP as u64,
        "CBC response byte cap"
    );
    let status = fs::read_to_string(&output)?;
    let header = status.lines().next().unwrap_or("");
    let optimal = header.starts_with("Optimal");
    if !(optimal || header.starts_with("Stopped on time"))
        || header.to_ascii_lowercase().contains("no integer solution")
    {
        return Ok((None,"status_has_no_integer_incumbent"));
    }
    let bytes = fs::read(&binary)?;
    let Ok(values) = native_columns(&bytes, model.rows.len(), model.variables.len()) else {
        return Ok((None,"native_solution_dimensions_or_values_invalid"));
    };
    if !model.feasible(&values) {
        return Ok((None,"incumbent_violates_original_model"));
    }
    let objective: f64 = model
        .variables
        .iter()
        .zip(&values)
        .map(|(v, x)| v.objective * x)
        .sum();
    let reported = f64::from_le_bytes(bytes[8..16].try_into()?);
    if !objective.is_finite()
        || !reported.is_finite()
        || (objective - reported).abs() > 1e-7 * (1.0 + objective.abs())
    {
        return Ok((None,"objective_serialization_mismatch"));
    }
    Ok((Some(Solution {
        values,
        objective,
        optimal,
        exact_epigraph_thresholds: 0,
        cpsat_probe_count: 0,
        native_dbm_nodes: 0,
        native_dbm_probes: 0,
    }),"accepted_original_model_incumbent"))
    })();
    if let Some(directory)=diagnostic {
        use std::sync::atomic::{AtomicU64,Ordering};
        static CALL:AtomicU64=AtomicU64::new(0);
        let call=CALL.fetch_add(1,Ordering::Relaxed)+1;
        let directory=std::path::PathBuf::from(directory).join(format!("call-{call:04}"));
        fs::create_dir_all(&directory)?;
        let accepted=answer.as_ref().is_ok_and(|(s,_)|s.is_some());
        let mut bytes=0;
        for name in ["model.mps","start.sol","solution.sol","solution.bin","cbc.log","cbc.stderr.log"] {
            // Keep every status/log/solution, but full models only for the
            // first call and rejected calls. This bounds persistent disk use.
            if call>1 && accepted && matches!(name,"model.mps"|"start.sol"){continue;}
            let source=tmp.0.join(name);
            if source.is_file(){bytes+=fs::metadata(&source)?.len();
                ensure!(bytes<=CAP as u64,"CBC diagnostic byte cap");fs::copy(source,directory.join(name))?;}
        }
        fs::write(directory.join("status.json"),serde_json::to_vec_pretty(&serde_json::json!({
            "call":call,"process_success":finished,"accepted":accepted,
            "process_timed_out":process.timed_out,"process_exit_code":process.exit_code,
            "reason":match &answer {Ok((_,r))=>r.to_string(),Err(e)=>e.to_string()},
            "rows":model.rows.len(),"columns":model.variables.len(),"mip_start_supplied":start.is_some(),
            "search_seconds":search_seconds,"time_mode":"elapsed","cbc_threads":1,
            "thread_control":"serial CBC build plus worker CPU affinity; no unsupported threads option",
            "preprocess":cfg.cbc_preprocess,"presolve":cfg.cbc_presolve,
            "randomSeed":1,"randomCbcSeed":1,"original_model_validated":accepted}))?)?;
    }
    if let Ok((None,reason))=&answer {eprintln!("CBC rejected: {reason}; process_success={finished}");}
    answer.map(|(solution,_)|solution)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn binary_model() -> Model {
        Model {
            variables: vec![
                Variable {
                    lower: 0.0,
                    upper: 1.0,
                    objective: 1.0,
                    binary: true
                };
                2
            ],
            rows: vec![Row {
                terms: vec![(0, 2.0), (1, 2.0)],
                sense: Sense::Greater,
                rhs: 3.0,
            }],
        }
    }
    #[test]
    fn rejects_fractional_and_invalid_sparse_models() {
        let mut m = binary_model();
        assert!(!m.feasible(&[0.5, 1.0]));
        assert!(!m.feasible(&[0.0, 1.0]));
        assert!(!m.feasible(&[f64::NAN, 1.0]));
        assert!(m.feasible(&[1.0, 1.0]));
        m.rows[0].terms.push((0, 1.0));
        assert!(m.mps().is_err());
        m = binary_model();
        m.variables[0].upper = f64::INFINITY;
        assert!(m.mps().is_err());
    }
    #[test]
    fn actual_cbc_integer_and_full_precision_continuous_recovery() {
        let cfg = Settings::default();
        assert!(
            cfg.cbc.is_file(),
            "CBC is required for this native protocol test"
        );
        let m = binary_model();
        let s = solve(&m, &cfg, &Budget::new(3.0)).unwrap().unwrap();
        assert!(s.optimal);
        assert_eq!(s.values, [1.0, 1.0]);
        assert_eq!(s.objective, 2.0);
        let seeded = solve_with_start(&m, &cfg, &Budget::new(3.0), Some(&[1.0, 1.0]))
            .unwrap()
            .unwrap();
        assert_eq!(seeded.values, [1.0, 1.0]);
        assert!(solve_with_start(&m, &cfg, &Budget::new(3.0), Some(&[0.5, 1.0])).is_err());
        // Interleaved integer markers and zero objective column exercise MPS order.
        let rhs = 0.123456789123456;
        let m = Model {
            variables: vec![
                Variable {
                    lower: 0.0,
                    upper: 1.0,
                    objective: 0.0,
                    binary: true,
                },
                Variable {
                    lower: 0.0,
                    upper: 1.0,
                    objective: 1.0,
                    binary: false,
                },
                Variable {
                    lower: 0.0,
                    upper: 1.0,
                    objective: 0.0,
                    binary: true,
                },
            ],
            rows: vec![
                Row {
                    terms: vec![(0, 1.0)],
                    sense: Sense::Equal,
                    rhs: 1.0,
                },
                Row {
                    terms: vec![(1, 1.0)],
                    sense: Sense::Greater,
                    rhs,
                },
                Row {
                    terms: vec![(2, 1.0)],
                    sense: Sense::Less,
                    rhs: 0.0,
                },
            ],
        };
        let s = solve(&m, &cfg, &Budget::new(3.0)).unwrap().unwrap();
        assert_eq!(s.values[0], 1.0);
        assert_eq!(s.values[2], 0.0);
        assert!((s.values[1] - rhs).abs() < 1e-14);
        assert!(solve(&m, &cfg, &Budget::new(0.0)).unwrap().is_none());
    }
}
