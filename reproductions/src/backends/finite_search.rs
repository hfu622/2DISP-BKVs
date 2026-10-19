//! Exact finite-domain MAC search shared by the DBM backend and offline audit.
use std::{collections::HashSet,time::Instant};
pub(crate) struct Solver { g:usize,w:usize,conflicts:Vec<u64>,masks:Vec<Vec<u64>>,words:Vec<Vec<usize>>,deadline:Instant,pub(crate) nodes:u64,dead:HashSet<Vec<u64>> }
pub(crate) enum Answer { Sat(Vec<usize>), Unsat, Unknown }
impl Solver {
    pub(crate) fn new(n:usize,conflicts:Vec<u64>,masks:Vec<Vec<u64>>,deadline:Instant)->Self{
        let g=masks.len();let w=n.div_ceil(64);
        assert!(n<=24001 && g<=512 && conflicts.len()==n*w);
        let mut union=vec![0u64;w];
        for group in &masks{assert_eq!(group.len(),w);for k in 0..w{assert_eq!(union[k]&group[k],0,"groups must be disjoint");union[k]|=group[k];}}
        let words=masks.iter().map(|v|v.iter().enumerate().filter(|(_,v)|**v!=0).map(|(k,_)|k).collect()).collect();
        Self{g,w,conflicts,masks,words,deadline,nodes:0,dead:HashSet::new()}
    }
    pub(crate) fn solve(&mut self)->Answer{
        let mut dom=vec![0u64;self.w];for mask in &self.masks{for k in 0..self.w{dom[k]|=mask[k];}}
        self.search(dom,&(0..self.g).collect::<Vec<_>>())
    }
    pub(crate) fn negative_states(&self)->usize{self.dead.len()}
    fn candidates(&self,dom:&[u64],i:usize)->Vec<usize>{
        let mut out=Vec::new();
        for &k in &self.words[i]{
            let mut x=dom[k]&self.masks[i][k];while x!=0 {let b=x.trailing_zeros() as usize;out.push(k*64+b);x&=x-1;}
        }out
    }
    fn propagate(&self,dom:&mut [u64],changed:&[usize])->Result<bool,()>{
        let mut queue:Vec<_>=changed.iter().flat_map(|&j|(0..self.g).filter(move |&i|i!=j).map(move |i|(i,j))).collect();
        while let Some((i,j))=queue.pop(){
            if Instant::now()>=self.deadline{return Err(());}
            let mut revised=false;
            for c in self.candidates(dom,i){
                let support=self.words[j].iter().any(|&k|dom[k]&self.masks[j][k]&!self.conflicts[c*self.w+k]!=0);
                if !support {dom[c/64]&=!(1u64<<(c%64));revised=true;}
            }
            if revised {
                if self.words[i].iter().all(|&k|dom[k]&self.masks[i][k]==0){return Ok(false);}
                queue.extend((0..self.g).filter(|&k|k!=i && k!=j).map(|k|(k,i)));
            }
        }Ok(true)
    }
    fn search(&mut self,mut dom:Vec<u64>,changed:&[usize])->Answer{
        self.nodes+=1;
        if Instant::now()>=self.deadline{return Answer::Unknown;}
        if (0..self.g).any(|i|self.words[i].iter().all(|&k|dom[k]&self.masks[i][k]==0)){return Answer::Unsat;}
        match self.propagate(&mut dom,changed){Err(())=>return Answer::Unknown,Ok(false)=>return Answer::Unsat,_=>{}}
        if self.dead.contains(&dom){return Answer::Unsat;}
        let choice=(0..self.g).map(|i|(self.words[i].iter().map(|&k|(dom[k]&self.masks[i][k]).count_ones() as usize).sum::<usize>(),i))
            .filter(|&(n,_)|n>1).min();
        let Some((_,i))=choice else{return Answer::Sat((0..self.g).map(|i|self.candidates(&dom,i)[0]).collect());};
        let mut candidates=self.candidates(&dom,i);
        candidates.sort_by_key(|&c|(0..self.w).map(|k|(dom[k]&!self.masks[i][k]&self.conflicts[c*self.w+k]).count_ones()).sum::<u32>());
        for c in candidates {
            let mut next=dom.clone();for &k in &self.words[i]{next[k]&=!self.masks[i][k];}next[c/64]|=1u64<<(c%64);
            match self.search(next,&[i]){Answer::Unsat=>{},other=>return other}
        }
        if self.dead.len()<20000{self.dead.insert(dom);}Answer::Unsat
    }
}

#[cfg(test)]mod tests{
    use super::*;
    #[test]fn expired_search_is_unknown_even_with_a_simple_domain(){
        let mut solver=Solver::new(3,vec![0;3],vec![vec![2],vec![4]],Instant::now());
        assert!(matches!(solver.solve(),Answer::Unknown));
    }
}
