//! Recognized integer-contour continuous CPMs use exact rational LRA.
use super::{Settings,TempDir,bounded,cbc::{Model,Solution}};
use crate::model::Budget;
use anyhow::{Result,ensure};
use serde::Deserialize;
use std::{fs,process::{Command,Stdio}};
const WORKER:&str=include_str!("smt_compaction.py");
#[derive(Deserialize)]
struct Response{applicable:bool,optimal:bool,values:Option<Vec<f64>>,objective:Option<f64>,status:String,lower:Option<String>,upper:Option<String>,runtime:Option<f64>}
pub enum Attempt{NotApplicable,Finished(Option<Solution>)}
pub fn solve(model:&Model,cfg:&Settings,budget:&Budget,start:Option<&[f64]>)->Result<Attempt>{
    if budget.expired(){return Ok(Attempt::Finished(None));}
    if let Some(values)=start{ensure!(model.feasible(values),"SMT start violates original model");}
    let tmp=TempDir::new(&cfg.temp)?;
    let input=tmp.0.join("model.mps");let initial=tmp.0.join("start.json");let output=tmp.0.join("solution.json");
    fs::write(&input,model.mps()?)?;fs::write(&initial,serde_json::to_vec(&start)?)?;
    let error=tmp.0.join("error.txt");
    let mut cmd=Command::new(cfg.smt_python.as_ref().unwrap());
    cmd.arg("-c").arg(WORKER).arg(&input).arg(&initial).arg(&output).arg((budget.remaining()-1.).max(0.001).to_string())
        .current_dir(&tmp.0).stdout(Stdio::null()).stderr(Stdio::from(fs::File::create(&error)?));
    let finished=bounded(&mut cmd,budget)?;
    if !output.is_file(){ensure!(budget.expired(),"SMT worker failed (success={}): {}",finished,fs::read_to_string(error).unwrap_or_default());return Ok(Attempt::Finished(None));}
    ensure!(fs::metadata(&output)?.len()<256*1024*1024,"SMT response byte cap");
    let response:Response=serde_json::from_slice(&fs::read(output)?)?;
    if !response.applicable{return Ok(Attempt::NotApplicable);}
    println!("SMT continuous optimal={} status={} lower={:?} upper={:?} runtime={:?}",response.optimal,response.status,response.lower,response.upper,response.runtime);
    let Some(values)=response.values else{return Ok(Attempt::Finished(None));};
    ensure!(values.len()==model.variables.len()&&model.feasible(&values),"SMT output violates original float model");
    let objective:f64=model.variables.iter().zip(&values).map(|(v,x)|v.objective*x).sum();
    ensure!(response.objective.is_some_and(|v|v.is_finite()&&(v-objective).abs()<1e-8),"SMT objective differs");
    Ok(Attempt::Finished(Some(Solution{values,objective,optimal:response.optimal,exact_epigraph_thresholds:0,cpsat_probe_count:0,native_dbm_nodes:0,native_dbm_probes:0})))
}

#[cfg(test)]
mod tests{
    use super::*;
    use crate::{three_pm::models::{Shape,Pose,compact},toolkit::polygon::P};
    fn cfg()->Settings{let mut c=Settings::default();c.smt_python=Some("F:/anaconda3/python.exe".into());c}
    fn square()->Shape{Shape{kind:0,orientation:0,clearance:0.,vertices:vec![P::new(0.,0.),P::new(1.,0.),P::new(1.,1.),P::new(0.,1.)]}}
    #[test]fn actual_direct_diagonal_plane_has_exact_fractional_compaction_optimum(){
        let shape=Shape{kind:0,orientation:0,clearance:0.,vertices:vec![P::new(0.,0.),P::new(3.,0.),P::new(0.,2.)]};
        let m=compact(&[shape],&[Pose{shape:0,x:0.,y:0.},Pose{shape:0,x:3.,y:0.}],3.,6.,&[6.,6.],&Budget::new(2.)).unwrap().unwrap();
        let Attempt::Finished(Some(s))=solve(&m,&cfg(),&Budget::new(10.),None).unwrap() else{panic!("direct integer model did not use SMT");};
        assert!(s.optimal);assert_eq!(s.objective,4.5);assert!(m.feasible(&s.values));
    }
    #[test]fn actual_rational_compaction_optimum_and_exact_fractional_windows(){
        let shapes=[square()];
        for (poses,lambda,expected) in [
            ([Pose{shape:0,x:0.,y:0.},Pose{shape:0,x:2.,y:0.}],[3.,3.],2.),
            ([Pose{shape:0,x:0.,y:0.},Pose{shape:0,x:1.25,y:0.}],[0.,0.],2.25)]{
            let m=compact(&shapes,&poses,1.,3.,&lambda,&Budget::new(2.)).unwrap().unwrap();
            let s=super::super::cbc::solve(&m,&cfg(),&Budget::new(10.)).unwrap().unwrap();
            assert!(s.optimal);assert_eq!(s.objective,expected);assert!(m.feasible(&s.values));
        }
    }
    #[test]fn actual_smt_concave_arm_infeasibility_and_notch_feasibility(){
        let l=Shape{kind:1,orientation:0,clearance:0.,vertices:vec![P::new(0.,0.),P::new(3.,0.),P::new(3.,1.),P::new(1.,1.),P::new(1.,3.),P::new(0.,3.)]};
        let shapes=[square(),l];
        for (y,feasible) in [(0.,false),(1.,true)]{
            let poses=[Pose{shape:1,x:0.,y:0.},Pose{shape:0,x:1.,y}];
            let m=compact(&shapes,&poses,3.,4.,&[0.,0.],&Budget::new(2.)).unwrap().unwrap();
            let s=super::super::cbc::solve(&m,&cfg(),&Budget::new(10.)).unwrap();
            assert_eq!(s.is_some(),feasible);if let Some(s)=s{assert!(s.optimal);assert_eq!(s.objective,3.);}
        }
    }
    #[test]fn unsupported_fractional_plane_falls_back_and_deadline_is_unknown(){
        let shapes=[square()];let poses=[Pose{shape:0,x:0.,y:0.},Pose{shape:0,x:2.,y:0.}];
        let mut m=compact(&shapes,&poses,1.,3.,&[3.,3.],&Budget::new(2.)).unwrap().unwrap();
        for rhs in [-0.3,-1.+1e-14] {
            m.rows[0].rhs=rhs;
            assert!(matches!(solve(&m,&cfg(),&Budget::new(10.),None).unwrap(),Attempt::NotApplicable));
        }
        assert!(matches!(solve(&m,&cfg(),&Budget::new(0.),None).unwrap(),Attempt::Finished(None)));
    }
}
