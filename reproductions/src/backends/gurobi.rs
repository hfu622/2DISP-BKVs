//! Optional installed Gurobi adapter. Original MPS and starts stay authoritative.
//! Recognized DBMs gain an exactly equivalent extended length epigraph; candidate
//! domains, geometry, integer layouts and their length objectives are unchanged.
use super::{Settings, TempDir, bounded, cbc::{Model, Solution}};
use crate::model::Budget;
use anyhow::{Result, ensure};
use serde::Deserialize;
use std::{fs, process::{Command, Stdio}};

const WORKER: &str = r#"
import sys,json,os,math,time
import gurobipy as gp
deadline=time.monotonic()+float(sys.argv[4])
env=gp.Env(empty=True)
env.setParam('OutputFlag',0)
env.start()
m=gp.read(sys.argv[1],env=env)
m.Params.OutputFlag=0
m.Params.Threads=int(sys.argv[5])
m.Params.Seed=20261015
m.Params.FeasibilityTol=1e-9
m.Params.IntFeasTol=1e-9
m.Params.MIPGap=0
m.Params.TimeLimit=float(sys.argv[4])
start=json.load(open(sys.argv[2]))
if start is not None:
    for i,value in enumerate(start):
        m.getVarByName('V'+str(i)).Start=value
original_n=m.NumVars
columns=[m.getVarByName('V'+str(i)) for i in range(original_n)]
# Exact extended epigraph for DBMs. The original MPS stays authoritative.
# A partition-demand group has exactly q selected dots. For every right-end
# threshold r, its selected tail is <=q*y[r]. Monotonic y and L>=sum(delta*y)
# imply L>=the largest selected end. Conversely any original integer solution
# extends with y[r]=1 iff r<=its largest selected end. Thus no integer layout,
# objective or candidate domain is removed; only fractional relaxations change.
ends={};groups=[];original_rows=[];indices={v.VarName:i for i,v in enumerate(columns)}
is_dbm=(original_n>=3 and columns[0].VType=='C' and columns[0].Obj==1
        and all(v.VType in ['B','I'] and 0<=v.LB<=v.UB<=1 and v.Obj==0 for v in columns[1:]))
if is_dbm:
    for row in m.getConstrs():
        a=m.getRow(row)
        terms=[(indices[a.getVar(i).VarName],a.getCoeff(i)) for i in range(a.size())]
        original_rows.append((terms,row.Sense,row.RHS))
        if (row.Sense=='<' and row.RHS==0 and len(terms)==2
            and terms[0]==(0,-1.) and terms[1][0]>0 and terms[1][1]>0):
            ends[terms[1][0]]=terms[1][1]
        elif (row.Sense=='=' and row.RHS>=1 and row.RHS.is_integer()
              and all(i>0 and c==1 for i,c in terms)):
            groups.append(([i for i,c in terms],row.RHS))
    partition=[i for g,q in groups for i in g]
    is_dbm=(len(ends)==original_n-1 and len(partition)==original_n-1
            and set(partition)==set(ends))
threshold_count=0
if is_dbm:
    rights=sorted(set(ends.values()))
    # A cap merely disables strengthening; it never truncates original dots.
    if len(rights)*len(groups)<=100000:
        y=m.addVars(len(rights),lb=0,ub=1,name='epigraph')
        m.addConstr(columns[0]>=gp.quicksum(
            (r-(rights[k-1] if k else 0))*y[k] for k,r in enumerate(rights)))
        for k,r in enumerate(rights):
            if k:m.addConstr(y[k]<=y[k-1])
            for group,q in groups:
                tail=[columns[i] for i in group if ends[i]>=r]
                if tail:m.addConstr(gp.quicksum(tail)<=q*y[k])
        threshold_count=len(rights)
def save(payload):
    with open(sys.argv[3]+'.tmp','w') as f:json.dump(payload,f,allow_nan=False)
    os.replace(sys.argv[3]+'.tmp',sys.argv[3])
cpsat_probes=[]
# A solver-internal primal/proof portfolio for large all-integer DBMs. Every
# probe uses this invocation's original MPS only, never a different algorithm's
# layout or a persisted result. At a fixed length the max-end rows simply fix
# forbidden dots to0. SAT supplies a checked MIP start; UNSAT at the immediately
# preceding endpoint gives a valid lower bound. Gurobi still certifies the
# final strengthened MIP, and Rust checks every original submitted row again.
if (is_dbm and original_n>=int(sys.argv[6]) and columns[0].UB<=20
    and start is not None and all(q==1 for g,q in groups)
    and all(rhs.is_integer() and all(a.is_integer() for i,a in terms)
            for terms,sense,rhs in original_rows)):
    try:
        sys.path.insert(0,os.environ.get('BKV_ORTOOLS_PATH','F:/packing-deps/ortools-bkv'))
        import ortools
        from ortools.sat.python import cp_model
        current=[float(x>0.5) for x in start]
        current[0]=max(right*current[i] for i,right in ends.items())
        endpoints=sorted(set(ends.values()))
        while time.monotonic()<deadline-1:
            targets=[r for r in endpoints if r<current[0]]
            if not targets:break
            target=targets[-1]
            cp=cp_model.CpModel();z=[None]+[cp.NewBoolVar(v.VarName) for v in columns[1:]]
            for i,v in enumerate(columns[1:],1):
                if v.LB==1:cp.Add(z[i]==1)
                if v.UB==0 or ends[i]>target:cp.Add(z[i]==0)
            for terms,sense,rhs in original_rows:
                if any(i==0 for i,a in terms):
                    assert sense=='<' and rhs==0 and len(terms)==2 and terms[0]==(0,-1.)
                    continue
                if all(a==1 for i,a in terms) and rhs==1 and sense=='<':
                    cp.AddAtMostOne([z[i] for i,a in terms]);continue
                if all(a==1 for i,a in terms) and rhs==1 and sense=='=':
                    cp.AddExactlyOne([z[i] for i,a in terms]);continue
                expression=sum(int(a)*z[i] for i,a in terms)
                if sense=='<':cp.Add(expression<=int(rhs))
                elif sense=='>':cp.Add(expression>=int(rhs))
                else:cp.Add(expression==int(rhs))
            sat=cp_model.CpSolver()
            sat.parameters.num_search_workers=min(3,int(sys.argv[5]))
            sat.parameters.max_time_in_seconds=min(120.,max(0.001,deadline-time.monotonic()-1.))
            sat.parameters.random_seed=20261015
            status=sat.Solve(cp)
            probe={'length':target,'status':sat.StatusName(status),'seconds':sat.WallTime(),
                   'workers':sat.parameters.num_search_workers,'ortools_version':ortools.__version__}
            cpsat_probes.append(probe)
            if status in [cp_model.OPTIMAL,cp_model.FEASIBLE]:
                candidate=[0.]+[float(sat.Value(v)) for v in z[1:]]
                candidate[0]=max(right*candidate[i] for i,right in ends.items())
                assert candidate[0]<=target
                assert all(v.LB-1e-9<=x<=v.UB+1e-9 for v,x in zip(columns,candidate))
                for terms,sense,rhs in original_rows:
                    lhs=sum(a*candidate[i] for i,a in terms)
                    assert abs(lhs-rhs)<=1e-8 if sense=='=' else lhs<=rhs+1e-8 if sense=='<' else lhs>=rhs-1e-8
                current=candidate
                for v,x in zip(columns,current):v.Start=x
                probe['original_model_rows_checked']=True
                save({'optimal':False,'status':'same_invocation_CP_SAT_feasible_MIP_start',
                    'threads':int(m.Params.Threads),'values':current,'objective':current[0],
                    'exact_epigraph_thresholds':threshold_count,'cpsat_probes':cpsat_probes})
            elif status==cp_model.INFEASIBLE:
                # No original integer layout can have L<=target. Its selected
                # maximum end must be an endpoint, hence at least current[0].
                columns[0].LB=max(columns[0].LB,current[0]);break
            else:break
    except ImportError:
        cpsat_probes.append({'status':'optional_CP_SAT_runtime_unavailable_Gurobi_only'})
m.Params.TimeLimit=max(0.001,deadline-time.monotonic())
def incumbent(model,where):
    if where==gp.GRB.Callback.MIPSOL:
        save({'optimal':False,'status':'incumbent_checkpoint',
              'threads':int(model.Params.Threads),
              'values':model.cbGetSolution(columns),
              'objective':model.cbGet(gp.GRB.Callback.MIPSOL_OBJ),
              'exact_epigraph_thresholds':threshold_count,'cpsat_probes':cpsat_probes})
m.optimize(incumbent)
out={'status':m.Status,'optimal':m.Status==gp.GRB.OPTIMAL,'version':gp.gurobi.version(),
     'threads':int(m.Params.Threads),
     'solutions':m.SolCount,'runtime':m.Runtime,'variables':original_n,'rows':m.NumConstrs,
     'exact_epigraph_thresholds':threshold_count,
     'cpsat_probes':cpsat_probes,
     'bound':m.ObjBound if m.IsMIP and math.isfinite(m.ObjBound) else None}
if m.SolCount:
    out.update(values=[v.X for v in columns],objective=m.ObjVal)
save(out)
"#;

#[derive(Deserialize)]
struct Response {
    threads: u32,
    optimal: bool,
    values: Option<Vec<f64>>,
    objective: Option<f64>,
    #[serde(default)]
    exact_epigraph_thresholds: usize,
    bound: Option<f64>,
    runtime: Option<f64>,
    #[serde(default)]
    cpsat_probes: Vec<serde_json::Value>,
}

pub fn solve(model:&Model,cfg:&Settings,budget:&Budget,start:Option<&[f64]>)->Result<Option<Solution>> {
    ensure!((1..=6).contains(&cfg.gurobi_threads),"Gurobi threads must be between1 and6");
    if budget.expired() { return Ok(None); }
    let text=model.mps()?;
    if let Some(values)=start { ensure!(model.feasible(values),"Gurobi start violates original model"); }
    let tmp=TempDir::new(&cfg.temp)?;
    let input=tmp.0.join("model.mps");
    let initial=tmp.0.join("start.json");
    let output=tmp.0.join("solution.json");
    fs::write(&input,text)?;
    fs::write(&initial,serde_json::to_vec(&start)?)?;
    let python=cfg.gurobi_python.as_ref().unwrap();
    let mut cmd=Command::new(python);
    cmd.arg("-c").arg(WORKER).arg(&input).arg(&initial).arg(&output)
        .arg((budget.remaining()-0.5).max(0.001).to_string())
        .arg(cfg.gurobi_threads.to_string())
        .arg(cfg.gurobi_cpsat_min_columns.to_string())
        .current_dir(&tmp.0).env("OMP_NUM_THREADS",cfg.gurobi_threads.to_string())
        .stdout(Stdio::null()).stderr(Stdio::from(fs::File::create(tmp.0.join("gurobi-error.txt"))?));
    // ManagedChild has stopped before bounded returns. An atomic checkpoint
    // can preserve a feasible incumbent if the hard watchdog interrupts final
    // solver shutdown. It is always nonoptimal and independently revalidated.
    let _finished=bounded(&mut cmd,budget)?;
    if !output.is_file() { return Ok(None); }
    ensure!(fs::metadata(&output)?.len()<=256*1024*1024,"Gurobi response byte cap");
    let response:Response=serde_json::from_slice(&fs::read(&output)?)?;
    ensure!(response.threads==cfg.gurobi_threads,"Gurobi worker thread metadata mismatch");
    println!("Gurobi original_cols={} exact_epigraph_thresholds={} optimal={} bound={:?} runtime={:?}",
        model.variables.len(),response.exact_epigraph_thresholds,response.optimal,response.bound,response.runtime);
    if !response.cpsat_probes.is_empty(){println!("Same-invocation optimizer CP-SAT probes: {}",serde_json::to_string(&response.cpsat_probes)?);}
    let Some(values)=response.values else { return Ok(None); };
    if values.len()!=model.variables.len() || !model.feasible(&values) { return Ok(None); }
    let objective:f64=model.variables.iter().zip(&values).map(|(v,x)|v.objective*x).sum();
    let Some(reported)=response.objective else { return Ok(None); };
    if !objective.is_finite() || !reported.is_finite() || (objective-reported).abs()>1e-7*(1.+objective.abs()) { return Ok(None); }
    Ok(Some(Solution{values,objective,optimal:response.optimal,exact_epigraph_thresholds:response.exact_epigraph_thresholds,cpsat_probe_count:response.cpsat_probes.len(),native_dbm_nodes:0,native_dbm_probes:0}))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backends::cbc::{Variable,Row,Sense};
    #[test]
    fn dbm_extended_epigraph_preserves_every_integer_placement_and_optimum() {
        let mut cfg=Settings::default();cfg.gurobi_python=Some("F:/anaconda3/python.exe".into());
        let ends=[2.,3.,4.,1.];
        let mut model=Model{variables:vec![Variable{lower:0.,upper:5.,objective:1.,binary:false}],rows:Vec::new()};
        for (k,&end) in ends.iter().enumerate(){
            model.variables.push(Variable{lower:0.,upper:1.,objective:0.,binary:true});
            model.rows.push(Row{terms:vec![(0,-1.),(k+1,end)],sense:Sense::Less,rhs:0.});
        }
        for group in [[1,2],[3,4]]{model.rows.push(Row{terms:group.into_iter().map(|i|(i,1.)).collect(),sense:Sense::Equal,rhs:1.});}
        model.rows.push(Row{terms:vec![(1,1.),(4,1.)],sense:Sense::Less,rhs:1.});
        let mut expected=f64::INFINITY;
        for a in 1..=2 { for b in 3..=4 {
            let mut values=vec![0.;5];values[a]=1.;values[b]=1.;values[0]=ends[a-1].max(ends[b-1]);
            let feasible=model.feasible(&values);
            let mut fixed=model.clone();
            for i in 1..5{fixed.variables[i].lower=values[i];fixed.variables[i].upper=values[i];}
            let solution=solve(&fixed,&cfg,&Budget::new(10.),None).unwrap();
            if feasible {
                expected=expected.min(values[0]);let s=solution.unwrap();assert!(s.optimal);assert!((s.objective-values[0]).abs()<1e-12);assert!(model.feasible(&s.values));
                // The extended epigraph must allow slack continuous length too.
                fixed.variables[0].lower=4.5;
                let s=solve(&fixed,&cfg,&Budget::new(10.),None).unwrap().unwrap();assert!(s.optimal);assert_eq!(s.objective,4.5);
            }else{assert!(solution.is_none());}
        }}
        let s=solve(&model,&cfg,&Budget::new(10.),None).unwrap().unwrap();assert!(s.optimal);assert_eq!(s.objective,expected);assert_eq!(expected,3.);
        assert_eq!(s.exact_epigraph_thresholds,4,"The actual MPS integer0/1 columns must activate the extended epigraph");
    }
    #[test]
    fn same_invocation_cpsat_start_and_unsat_bound_preserve_original_optimum() {
        let mut cfg=Settings::default();cfg.gurobi_python=Some("F:/anaconda3/python.exe".into());
        cfg.gurobi_cpsat_min_columns=3;
        let model=Model{variables:vec![
            Variable{lower:0.,upper:5.,objective:1.,binary:false},
            Variable{lower:0.,upper:1.,objective:0.,binary:true},
            Variable{lower:0.,upper:1.,objective:0.,binary:true},
            Variable{lower:0.,upper:1.,objective:0.,binary:true},
            Variable{lower:0.,upper:1.,objective:0.,binary:true}],
            rows:vec![
                Row{terms:vec![(0,-1.),(1,2.)],sense:Sense::Less,rhs:0.},
                Row{terms:vec![(0,-1.),(2,3.)],sense:Sense::Less,rhs:0.},
                Row{terms:vec![(0,-1.),(3,4.)],sense:Sense::Less,rhs:0.},
                Row{terms:vec![(0,-1.),(4,1.)],sense:Sense::Less,rhs:0.},
                Row{terms:vec![(1,1.),(2,1.)],sense:Sense::Equal,rhs:1.},
                Row{terms:vec![(3,1.),(4,1.)],sense:Sense::Equal,rhs:1.},
                Row{terms:vec![(1,1.),(4,1.)],sense:Sense::Less,rhs:1.}]};
        let s=solve(&model,&cfg,&Budget::new(15.),Some(&[4.,1.,0.,1.,0.])).unwrap().unwrap();
        assert!(s.optimal);assert_eq!(s.objective,3.);assert!(model.feasible(&s.values));
        assert_eq!(s.cpsat_probe_count,2,"SAT at3, then UNSAT at2 must be recorded from this actual invocation");
    }
    #[test]
    fn native_gurobi_interleaved_binary_and_f64_columns() {
        let mut cfg=Settings::default();cfg.gurobi_python=Some("F:/anaconda3/python.exe".into());
        let value=0.123456789123456;
        let model=Model{variables:vec![
            Variable{lower:0.,upper:1.,objective:0.,binary:true},
            Variable{lower:0.,upper:1.,objective:1.,binary:false},
            Variable{lower:0.,upper:1.,objective:0.,binary:true}],
            rows:vec![Row{terms:vec![(0,1.)],sense:Sense::Equal,rhs:1.},
                      Row{terms:vec![(1,1.)],sense:Sense::Greater,rhs:value},
                      Row{terms:vec![(2,1.)],sense:Sense::Equal,rhs:0.}]};
        for threads in [1,6] {
            cfg.gurobi_threads=threads;
            let s=solve(&model,&cfg,&Budget::new(10.),Some(&[1.,1.,0.])).unwrap().unwrap();
            assert!(s.optimal);assert_eq!(s.values[0],1.);assert_eq!(s.values[2],0.);
            assert!((s.values[1]-value).abs()<1e-14);assert!((s.objective-value).abs()<1e-14);
        }
        assert!(solve(&model,&cfg,&Budget::new(0.),None).unwrap().is_none());
        assert!(solve(&model,&cfg,&Budget::new(1.),Some(&[0.,0.,0.])).is_err());
        for threads in [0,7] { cfg.gurobi_threads=threads;assert!(solve(&model,&cfg,&Budget::new(0.),None).is_err()); }
    }
}
