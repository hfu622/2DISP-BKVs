//! Exact optimizer for unit-demand DBMs; continuous/local MIPs fall back to CBC.
use super::{cbc::{Model,Row,Sense,Solution},finite_search::{Answer,Solver}};
use crate::model::Budget;
use anyhow::{Result,ensure};
use std::time::{Duration,Instant};

pub enum Attempt { NotApplicable, Finished(Option<Solution>) }
struct Problem { groups:Vec<Vec<usize>>,ends:Vec<f64>,allowed:Vec<bool>,conflicts:Vec<u64>,words:usize }
fn length_row(row:&Row)->Option<(usize,f64)>{
    if matches!(row.sense,Sense::Less) && row.rhs==0. && row.terms.len()==2
        && row.terms[0]==(0,-1.) && row.terms[1].0>0 && row.terms[1].1>0.
    {Some(row.terms[1])}else{None}
}
fn recognize(m:&Model,budget:&Budget)->Option<Problem>{
    let n=m.variables.len();if !(3..=24001).contains(&n){return None;}
    let l=&m.variables[0];
    if l.binary || l.objective!=1. || l.lower<0. || m.variables[1..].iter().any(|v|!v.binary || v.objective!=0. || ![0.,1.].contains(&v.lower) || ![0.,1.].contains(&v.upper)){return None;}
    let mut ends=vec![f64::NAN;n];let mut groups=Vec::new();
    for row in &m.rows{
        if let Some((i,right))=length_row(row){if ends[i].is_finite(){return None;}ends[i]=right;}
        else if matches!(row.sense,Sense::Equal) && row.rhs==1. && row.terms.iter().all(|&(i,a)|i>0 && a==1.){
            groups.push(row.terms.iter().map(|&(i,_)|i).collect::<Vec<_>>());
        }
    }
    if groups.is_empty() || groups.len()>512 || ends[1..].iter().any(|r|!r.is_finite()){return None;}
    let mut owner=vec![usize::MAX;n];
    for (g,group) in groups.iter().enumerate(){for &i in group{if owner[i]!=usize::MAX{return None;}owner[i]=g;}}
    if owner[1..].contains(&usize::MAX){return None;}
    let words=n.div_ceil(64);let mut conflicts=vec![0u64;n*words];
    let mut allowed:Vec<_>=m.variables.iter().map(|v|v.upper==1.).collect();allowed[0]=false;
    for group in &groups{
        let fixed:Vec<_>=group.iter().filter(|&&i|m.variables[i].lower==1.).copied().collect();
        if !fixed.is_empty(){for &i in group{allowed[i]&=fixed.len()==1 && i==fixed[0];}}
    }
    for row in &m.rows{
        if budget.expired(){return None;}
        if length_row(row).is_some(){continue;}
        if row.terms.iter().any(|&(i,_)|i==0){return None;}
        let group=owner[row.terms[0].0];
        if row.terms.iter().all(|&(i,_)|owner[i]==group){
            // Exactly one candidate in this group is selected. A type-local
            // movement/fixing row becomes an exact unary allowed-value test.
            if row.rhs.fract()!=0. || row.terms.iter().any(|&(_,a)|a.fract()!=0.){return None;}
            for &i in &groups[group]{
                let a=row.terms.binary_search_by_key(&i,|&(j,_)|j).ok().map_or(0.,|k|row.terms[k].1);
                allowed[i]&=match row.sense{Sense::Less=>a<=row.rhs,Sense::Equal=>a==row.rhs,Sense::Greater=>a>=row.rhs};
            }
        }else if matches!(row.sense,Sense::Less) && row.rhs==1. && row.terms.iter().all(|&(_,a)|a==1.){
            // A binary at-most-one clique is exactly the conjunction of all
            // its pair conflicts. No fractional or geometric relaxation here.
            for (p,&(i,_)) in row.terms.iter().enumerate(){for &(j,_) in &row.terms[p+1..]{
                conflicts[i*words+j/64]|=1u64<<(j%64);conflicts[j*words+i/64]|=1u64<<(i%64);
            }}
        }else{return None;}
    }
    Some(Problem{groups,ends,allowed,conflicts,words})
}
fn solution(m:&Model,p:&Problem,selected:&[usize],optimal:bool,nodes:u64,probes:usize)->Result<Solution>{
    let mut values=vec![0.;m.variables.len()];values[0]=m.variables[0].lower;
    for &i in selected{values[i]=1.;values[0]=values[0].max(p.ends[i]);}
    ensure!(m.feasible(&values),"Native DBM assignment violates original submitted model");
    Ok(Solution{objective:values[0],values,optimal,exact_epigraph_thresholds:0,cpsat_probe_count:0,native_dbm_nodes:nodes,native_dbm_probes:probes})
}
pub fn solve(m:&Model,budget:&Budget,start:Option<&[f64]>)->Result<Attempt>{
    m.validate()?;
    let Some(p)=recognize(m,budget)else{return Ok(Attempt::NotApplicable);};
    let mut best=if let Some(values)=start{
        ensure!(m.feasible(values),"Native DBM received an invalid original-model MIP start");
        let candidate=solution(m,&p,&(1..values.len()).filter(|&i|values[i]>0.5).collect::<Vec<_>>(),false,0,0)?;
        (candidate.objective<=m.variables[0].upper).then_some(candidate)
    }else{None};
    let mut endpoints=vec![m.variables[0].lower];
    endpoints.extend(p.ends[1..].iter().copied().filter(|&r|r>=m.variables[0].lower && r<=m.variables[0].upper));
    endpoints.sort_by(f64::total_cmp);endpoints.dedup();
    let mut nodes=0;let mut probes=0;
    loop{
        if budget.expired(){return Ok(Attempt::Finished(best));}
        let target=if let Some(ref s)=best{
            let Some(&r)=endpoints.iter().rev().find(|&&r|r<s.objective)else{
                let mut s=best.unwrap();s.optimal=true;s.native_dbm_nodes=nodes;s.native_dbm_probes=probes;
                return Ok(Attempt::Finished(Some(s)));
            };r
        }else{m.variables[0].upper};
        let masks=p.groups.iter().map(|group|{let mut mask=vec![0u64;p.words];for &i in group{
            if p.allowed[i] && p.ends[i]<=target{mask[i/64]|=1u64<<(i%64);}
        }mask}).collect();
        let mut solver=Solver::new(m.variables.len(),p.conflicts.clone(),masks,Instant::now()+Duration::from_secs_f64(budget.remaining().max(0.)));
        let answer=solver.solve();nodes+=solver.nodes;probes+=1;
        let state=match &answer{Answer::Sat(_)=>"SAT",Answer::Unsat=>"UNSAT",Answer::Unknown=>"UNKNOWN"};
        println!("Native original DBM target={target} state={state} nodes={} negative_states={}",solver.nodes,solver.negative_states());
        match answer{
            Answer::Sat(selected)=>best=Some(solution(m,&p,&selected,false,nodes,probes)?),
            Answer::Unsat=>{
                if let Some(ref mut s)=best{s.optimal=true;s.native_dbm_nodes=nodes;s.native_dbm_probes=probes;}
                return Ok(Attempt::Finished(best));
            },
            Answer::Unknown=>{
                if let Some(ref mut s)=best{s.native_dbm_nodes=nodes;s.native_dbm_probes=probes;}
                return Ok(Attempt::Finished(best));
            }
        }
    }
}

#[cfg(test)]mod tests{
    use super::*;
    use super::super::{Settings,cbc::{self,Variable}};
    fn variable(lower:f64,upper:f64,objective:f64,binary:bool)->Variable{Variable{lower,upper,objective,binary}}
    fn row(terms:Vec<(usize,f64)>,sense:Sense,rhs:f64)->Row{Row{terms,sense,rhs}}
    fn tiny(case:usize)->Model{
        let mut m=Model{variables:vec![variable(if case%4==0{3.5}else{0.},if case%3==0{5.}else{7.},1.,false)],rows:Vec::new()};
        for i in 1..=6{
            m.variables.push(variable(if case%11==0 && i==1{1.}else{0.},if case%13==0 && i==5{0.}else{1.},0.,true));
            let end=1.+((case*7+i*3)%11) as f64*0.5;
            m.rows.push(row(vec![(0,-1.),(i,end)],Sense::Less,0.));
        }
        for group in [1..=3,4..=6]{m.rows.push(row(group.map(|i|(i,1.)).collect(),Sense::Equal,1.));}
        for i in 1..=3{for j in 4..=6{if (case+i*5+j*7)%4==0{m.rows.push(row(vec![(i,1.),(j,1.)],Sense::Less,1.));}}}
        if case%2==0{m.rows.push(row(vec![(1,1.),(2,1.),(4,1.)],Sense::Less,1.));}
        m.rows.push(row(vec![(1,-1.),(2,1.)],Sense::Less,((case%3) as f64)-1.));
        if case%7==0{m.rows.push(row(vec![(4,1.),(6,-1.)],Sense::Greater,0.));}
        m
    }
    #[test]fn exact_optimizer_matches_exhaustive_original_models_with_local_movement_and_fixed_values(){
        for case in 0..120{
            let m=tiny(case);let mut expected=f64::INFINITY;let mut start=None;
            for a in 1..=3{for b in 4..=6{
                let mut v=vec![0.;7];v[a]=1.;v[b]=1.;v[0]=m.variables[0].lower;
                for i in [a,b]{v[0]=v[0].max(m.rows[i-1].terms[1].1);}
                if m.feasible(&v){expected=expected.min(v[0]);v[0]=m.variables[0].upper;start=Some(v);}
            }}
            let supplied=if case%2==0{start.as_deref()}else{None};
            match solve(&m,&Budget::new(2.),supplied).unwrap(){
                Attempt::NotApplicable=>panic!("recognized model rejected at case{case}"),
                Attempt::Finished(None)=>assert!(!expected.is_finite(),"false infeasibility at case{case}"),
                Attempt::Finished(Some(s))=>{assert!(s.optimal);assert_eq!(s.objective,expected);assert!(m.feasible(&s.values));}
            }
        }
    }
    #[test]fn unsupported_coupled_cardinality_does_not_enter_the_native_gate(){
        let mut m=tiny(1);m.rows.push(row(vec![(1,2.),(4,-1.)],Sense::Greater,0.));
        assert!(matches!(solve(&m,&Budget::new(1.),None).unwrap(),Attempt::NotApplicable));
    }
    #[test]fn native_gate_preserves_continuous_mip_cbc_fallback(){
        let value=0.123456789123456;
        let m=Model{variables:vec![variable(0.,1.,1.,false),variable(0.,1.,0.,false),variable(0.,1.,0.,true)],
            rows:vec![row(vec![(0,1.),(1,-1.)],Sense::Greater,0.),row(vec![(1,1.)],Sense::Greater,value),row(vec![(2,1.)],Sense::Equal,1.)]};
        let mut cfg=Settings::default();cfg.finite_dbm=true;
        let s=cbc::solve(&m,&cfg,&Budget::new(10.)).unwrap().unwrap();
        assert!(s.optimal && m.feasible(&s.values));assert!((s.objective-value).abs()<1e-14);assert_eq!(s.native_dbm_probes,0);
    }
}
