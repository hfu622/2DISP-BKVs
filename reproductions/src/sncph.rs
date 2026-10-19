//! Baldacci et al. SNH1/2/3, PH1/2/3 and pixel-price SNCPH.
//! Unreported control constants/strip wrapper are explicit project conventions.
use crate::{geometry::{self,Geometry},model::*,toolkit::{bitmap,placement::*}};
use anyhow::{Result,ensure};

pub const ASSUMPTIONS:&[&str]=&[
 "Source SNH1/2/3, PH1/2/3 and SNCPH; actual resolved parameters and unavailable-author decisions are embedded in reproduction_configuration",
 "Original XML contours and rotations; conservative raster, original-contour/Jagua validation; no generic postprocessor",
 "Completing the configured patterns is not a claim to have recovered unpublished author controls or the full paper statistical protocol"];

#[derive(Clone)]struct Raster{kind:usize,orientation:usize,w:usize,h:usize,ox:i32,oy:i32,runs:Vec<(usize,usize,usize)>,pixels:usize}
#[derive(Clone,Copy,Debug)]struct Pose{shape:usize,x:usize,y:usize,piece:usize}
struct Board{w:usize,h:usize,bits:Vec<Vec<u64>>,prefix:Vec<u32>,weights:Vec<f64>,vprefix:Vec<f64>,total:usize,extent:usize,valid_offset:f64}
impl Board{
 fn new(w:usize,h:usize,weights:Vec<f64>)->Self{
  // |alpha+beta-gamma+delta| <= 2+2*w*h+sum(abs(weights)).
  // A derived offset preserves the paper's feasible > infeasible ordering
  // for every validated configuration, rather than assuming 1e12 is enough.
  let valid_offset=3.+2.*(w*h) as f64+weights.iter().map(|v|v.abs()).sum::<f64>();
  let mut s=Self{w,h,bits:vec![vec![0;h.div_ceil(64)];w],prefix:vec![0;(w+1)*(h+1)],weights,vprefix:vec![0.0;w*(h+1)],total:0,extent:0,valid_offset};
  for x in 0..w{let mut a=0.;for y in 0..h{a+=s.weights[x*h+y];s.vprefix[x*(h+1)+y+1]=a;}}s
 }
 fn count_col(&self,x:usize,a:usize,b:usize)->usize{
  if a>=b{return 0;}let mut sum=0;
  for k in a/64..=(b-1)/64{let lo=a.max(k*64)-k*64;let hi=b.min((k+1)*64)-k*64;
   let left=u64::MAX<<lo;let right=if hi==64{u64::MAX}else{(1u64<<hi)-1};sum+=(self.bits[x][k]&left&right).count_ones() as usize;}sum
 }
 fn overlaps(&self,s:&Raster,x:usize,y:usize)->usize{s.runs.iter().map(|&(i,a,b)|self.count_col(x+i,y+a,y+b)).sum()}
 fn fits(&self,s:&Raster,x:usize,y:usize)->bool{x+s.w<=self.w&&y+s.h<=self.h&&s.runs.iter().all(|&(i,a,b)|self.count_col(x+i,y+a,y+b)==0)}
 fn put(&mut self,s:&Raster,x:usize,y:usize,add:bool){
  for &(i,a,b) in &s.runs{for yy in y+a..y+b{let bit=1u64<<(yy%64);if add{self.bits[x+i][yy/64]|=bit;}else{self.bits[x+i][yy/64]&=!bit;}}}
 }
 fn rebuild(&mut self){
  self.prefix.fill(0);self.total=0;self.extent=0;
  for x in 0..self.w{let mut col=0;for y in 0..self.h{let used=(self.bits[x][y/64]>>(y%64))&1;col+=used as u32;self.prefix[(x+1)*(self.h+1)+y+1]=self.prefix[x*(self.h+1)+y+1]+col;}
   if col>0{self.extent=x+1;}self.total+=col as usize;}
 }
 fn rectangle(&self,x0:usize,y0:usize,x1:usize,y1:usize)->usize{
  let at=|x,y|self.prefix[x*(self.h+1)+y] as i64;
  (at(x1,y1)+at(x0,y0)-at(x0,y1)-at(x1,y0)) as usize
 }
 fn free(&self,x0:usize,y0:usize,x1:usize,y1:usize)->usize{(x1-x0)*(y1-y0)-self.rectangle(x0,y0,x1,y1)}
 fn phi(&self,s:&Raster,x:usize,y:usize,top:bool,valid:bool)->f64{
  let alpha=(self.total+s.pixels) as f64/(self.h*self.extent.max(x+s.w)) as f64;
  let beta=self.rectangle(x,y,x+s.w,y+s.h) as f64;
  let gamma=(self.free(0,0,x,self.h)+if top{self.free(x,y+s.h,x+s.w,self.h)}else{self.free(x,0,x+s.w,y)}) as f64;
  let delta=s.runs.iter().map(|&(i,a,b)|self.vprefix[(x+i)*(self.h+1)+y+b]-self.vprefix[(x+i)*(self.h+1)+y+a]).sum::<f64>();
  let q=alpha+beta-gamma+delta;
  if valid{q+self.valid_offset}else{q-self.valid_offset-1e6*self.overlaps(s,x,y) as f64}
 }
}
fn rasters(g:&Geometry,h:usize,budget:&Budget)->Result<Vec<Raster>>{
 let step=g.height/h as f64;let mut out=Vec::new();
 for (kind,t) in g.types.iter().enumerate(){let piece=g.pieces.iter().position(|&p|p==kind).unwrap();
  for orientation in 0..t.orientations.len(){let v=oriented(g,piece,orientation);let b=Bounds::of(&v);
   let ox=(b.x0/step).floor() as i32;let oy=(b.y0/step).floor() as i32;
   let w=(b.x1/step).ceil() as i32-ox;let hh=(b.y1/step).ceil() as i32-oy;if w<=0||hh<=0||hh as usize>h{continue;}
   let Some(rects)=bitmap::scanline(&v,step,budget)? else{anyhow::bail!("raster preparation deadline");};
   let mut cells=vec![vec![false;hh as usize];w as usize];
   for r in rects{let x0=(r.x0/step).round() as i32-ox;let x1=(r.x1/step).round() as i32-ox;let yy=(r.y0/step).round() as i32-oy;
    if yy>=0&&yy<hh{for x in x0.max(0)..x1.min(w){cells[x as usize][yy as usize]=true;}}
   }
   let mut runs=Vec::new();let mut pixels=0;
   for (x,c) in cells.iter().enumerate(){let mut y=0;while y<c.len(){if !c[y]{y+=1;continue;}let lo=y;while y<c.len()&&c[y]{y+=1;}runs.push((x,lo,y));pixels+=y-lo;}}
   out.push(Raster{kind,orientation,w:w as usize,h:hh as usize,ox,oy,runs,pixels});
  }
 }ensure!(!out.is_empty(),"no raster shapes");Ok(out)
}
fn ph3(board:&Board,shapes:&[Raster],initial:Pose,top:bool,budget:&Budget,stats:&mut Stats)->Option<(Pose,f64)>{
 let mut p=initial;let mut q=board.phi(&shapes[p.shape],p.x,p.y,top,board.fits(&shapes[p.shape],p.x,p.y));let mut best=None;
 if board.fits(&shapes[p.shape],p.x,p.y){best=Some((p,q));}
 loop{if budget.expired(){return best;}
  let mut chosen=(p,q);
  let o=shapes.iter().enumerate().filter(|(_,s)|s.kind==shapes[p.shape].kind).map(|(i,_)|i).collect::<Vec<_>>();let oi=o.iter().position(|&i|i==p.shape).unwrap();
  for (dx,dy,dr) in [(1,0,0),(0,1,0),(0,0,1),(-1,0,0),(0,-1,0),(0,0,-1)]{
   let shape=o[(oi as i32+dr).rem_euclid(o.len() as i32) as usize];let x=p.x as i32+dx;let y=p.y as i32+dy;
   if x<0||y<0||x as usize+shapes[shape].w>board.w||y as usize+shapes[shape].h>board.h{continue;}
   let n=Pose{shape,x:x as usize,y:y as usize,..p};let valid=board.fits(&shapes[shape],n.x,n.y);stats.collision_queries+=1;
   let value=board.phi(&shapes[shape],n.x,n.y,top,valid);if value>chosen.1{chosen=(n,value);}
  }
  if chosen.1<=q{break;}p=chosen.0;q=chosen.1;
  if board.fits(&shapes[p.shape],p.x,p.y){best=Some((p,q));}
 }best
}
fn placement(board:&Board,shapes:&[Raster],kind:usize,ph:usize,piece:usize,top:bool,budget:&Budget,stats:&mut Stats)->Option<(Pose,f64)>{
 let mut best=None;
 for (shape,s) in shapes.iter().enumerate().filter(|(_,s)|s.kind==kind){
  if s.h>board.h||s.w>board.w{continue;}
  if ph==3{let initial=Pose{shape,x:0,y:if top{board.h-s.h}else{0},piece};
   if let Some(p)=ph3(board,shapes,initial,top,budget,stats){if best.as_ref().is_none_or(|&(_,q)|p.1>q){best=Some(p);}}continue;}
  for x in 0..=board.w-s.w{
   if budget.expired(){return best;}
   let mut found=None;for offset in 0..=board.h-s.h{let y=if top{board.h-s.h-offset}else{offset};stats.collision_queries+=1;
    if board.fits(s,x,y){let p=Pose{shape,x,y,piece};found=Some((p,board.phi(s,x,y,top,true)));break;}}
   if let Some(p)=found{if best.as_ref().is_none_or(|&(_,q)|p.1>q){best=Some(p);}if ph==1{break;}}
  }
 }best
}
fn complete(board:&mut Board,shapes:&[Raster],remaining:&mut Vec<usize>,poses:&mut Vec<Pose>,g:&Geometry,ph:usize,budget:&Budget,stats:&mut Stats){
 while !remaining.is_empty()&&!budget.expired(){let mut best=None;
  for kind in 0..g.types.len(){let Some(piece)=remaining.iter().copied().find(|&p|g.pieces[p]==kind)else{continue;};
   if let Some(p)=placement(board,shapes,kind,ph,piece,true,budget,stats){if best.as_ref().is_none_or(|&(_,q)|p.1>q){best=Some(p);}}
  }
  let Some((p,_))=best else{break;};board.put(&shapes[p.shape],p.x,p.y,true);board.rebuild();remaining.retain(|&v|v!=p.piece);poses.push(p);stats.evaluations+=1;
 }
}
fn initialize_k(board:&mut Board,shapes:&[Raster],remaining:&mut Vec<usize>,poses:&mut Vec<Pose>,g:&Geometry,budget:&Budget,stats:&mut Stats,k:usize){
 let maxweight=board.weights.iter().copied().fold(f64::NEG_INFINITY,f64::max);let mut positions=vec![Vec::<(Pose,f64)>::new();g.types.len()];
 for (shape,s) in shapes.iter().enumerate(){
  if s.w>board.w||s.h>board.h{continue;}
  let list=&mut positions[s.kind];for x in 0..=board.w-s.w{
   if budget.expired(){return;}
   if list.len()==k && board.valid_offset+1.+maxweight*s.pixels as f64-(x*board.h) as f64<=list.last().unwrap().1{break;}
   for y in 0..=board.h-s.h{let p=Pose{shape,x,y,piece:0};let q=board.phi(s,x,y,true,true);
    if list.len()<k||q>list.last().unwrap().1{list.push((p,q));list.sort_by(|a,b|b.1.total_cmp(&a.1));list.truncate(k);}}
  }
 }
 loop{if budget.expired(){break;}let mut best=None;
  for (kind,list) in positions.iter_mut().enumerate(){let Some(piece)=remaining.iter().copied().find(|&p|g.pieces[p]==kind)else{continue;};
   list.retain(|(p,_)|board.fits(&shapes[p.shape],p.x,p.y));
   for &(mut p,_) in list.iter(){p.piece=piece;let q=board.phi(&shapes[p.shape],p.x,p.y,true,true);if best.as_ref().is_none_or(|&(_,v)|q>v){best=Some((p,q));}}
  }
  let Some((p,_))=best else{break;};board.put(&shapes[p.shape],p.x,p.y,true);board.rebuild();poses.push(p);remaining.retain(|&v|v!=p.piece);stats.evaluations+=1;
 }
}
fn improve(board:&mut Board,shapes:&[Raster],remaining:&mut Vec<usize>,poses:&mut Vec<Pose>,g:&Geometry,ph:usize,budget:&Budget,stats:&mut Stats,max_it:usize){
 let mut stagnant=0;
 while stagnant<=max_it&&!budget.expired(){let before=poses.len();complete(board,shapes,remaining,poses,g,ph,budget,stats);if poses.len()>before{stagnant=0;}
  loop{let mut moved=false;for p in poses.iter_mut(){if budget.expired(){return;}
    let old=*p;board.put(&shapes[p.shape],p.x,p.y,false);board.rebuild();
    if let Some((q,_))=ph3(board,shapes,*p,stagnant%2==0,budget,stats){*p=q;}
    board.put(&shapes[p.shape],p.x,p.y,true);board.rebuild();
    if old.shape!=p.shape||old.x!=p.x||old.y!=p.y{moved=true;stats.accepted_moves+=1;}
   }if !moved||budget.expired(){break;}
  }stagnant+=1;
 }
}
fn layout(g:&Geometry,shapes:&[Raster],poses:&[Pose],h:usize)->Layout{
 let step=g.height/h as f64;
 let placements=poses.iter().map(|p|{let s=&shapes[p.shape];Placement{piece:p.piece,orientation:s.orientation,x:(p.x as f64-s.ox as f64)*step,y:(p.y as f64-s.oy as f64)*step}}).collect::<Vec<_>>();
 Layout{length:geometry::length(g,&placements),placements}
}
// The f32 Jagua bounding box may protrude by rounding while the conservative
// f64 raster and exported original vertices fit exactly. Measure the originals
// before accepting such a pattern, with numerical tolerance only.
fn original_feasible(g:&Geometry,l:&Layout)->bool{
 use geo::{Area,BooleanOps,BoundingRect};
 let eps=g.height.max(1.)*1e-9;let area_eps=g.total_area*1e-12;
 let mut prior=Vec::<geo::Polygon<f64>>::new();
 for p in &l.placements{
  let v=oriented(g,p.piece,p.orientation);let b=Bounds::of(&v);
  if p.x+b.x0 < -eps || p.y+b.y0 < -eps || p.y+b.y1 >g.height+eps || p.x+b.x1>l.length+eps{return false;}
  let shape=world(g,p);let b=shape.bounding_rect().unwrap();
  for other in &prior{let a=other.bounding_rect().unwrap();
   if a.min().x<b.max().x && b.min().x<a.max().x && a.min().y<b.max().y && b.min().y<a.max().y && shape.intersection(other).unsigned_area()>area_eps{return false;}
  }prior.push(shape);
 }true
}
pub fn solve(g:&Geometry,budget:&Budget,_seed:u64,stats:&mut Stats)->Result<Layout>{
 solve_with_config(g,budget,_seed,stats,&crate::config::get().sncph)
}
pub fn solve_with_config(g:&Geometry,budget:&Budget,_seed:u64,stats:&mut Stats,cfg:&crate::config::Sncph)->Result<Layout>{
 cfg.validate()?;
 stats.effective_algorithm_controls=Some(serde_json::to_value(cfg)?);
 let final_h=*cfg.raster_heights.last().unwrap();let final_w=final_h*cfg.panel_height_multiple;
 let mut best=None::<Layout>;let mut complete_patterns=0;let mut generated_patterns=0;
 for variant in 0..12{let ph=variant%3+1;let initial2=variant>=3&&variant<6||variant>=9;let multi=variant>=6;
  let levels=if multi{cfg.raster_heights.clone()}else{vec![final_h]};
  let mut values=vec![1.;final_w*final_h];
  for iteration in 0..cfg.patterns_per_variant{
   if budget.expired(){break;}
   let mut current=Vec::<Placement>::new();let mut final_board=None;
   for &h in &levels{
    if budget.expired(){break;}let w=h*cfg.panel_height_multiple;let shapes=match rasters(g,h,budget){Ok(v)=>v,Err(_) if budget.expired()=>break,Err(e)=>return Err(e)};
    let weights=(0..w).flat_map(|x|{let values=&values;(0..h).map(move|y|values[(x*final_h/h)*final_h+y*final_h/h])}).collect();let mut board=Board::new(w,h,weights);
    let step=g.height/h as f64;let mut poses=Vec::new();let mut remaining=(0..g.pieces.len()).collect::<Vec<_>>();
    for p in &current{let shape=shapes.iter().position(|s|s.kind==g.pieces[p.piece]&&s.orientation==p.orientation).unwrap();let s=&shapes[shape];
     let x=(p.x/step+s.ox as f64).round().max(0.) as usize;let y=(p.y/step+s.oy as f64).round().max(0.) as usize;
     ensure!(board.fits(s,x,y),"coarse-to-fine raster transfer must remain feasible");board.put(s,x,y,true);poses.push(Pose{shape,x,y,piece:p.piece});remaining.retain(|&v|v!=p.piece);}
    board.rebuild();
    if current.is_empty(){if initial2{initialize_k(&mut board,&shapes,&mut remaining,&mut poses,g,budget,stats,cfg.k);}else{complete(&mut board,&shapes,&mut remaining,&mut poses,g,ph,budget,stats);}}
    if initial2||!current.is_empty(){improve(&mut board,&shapes,&mut remaining,&mut poses,g,ph,budget,stats,cfg.max_it);}
    current=layout(g,&shapes,&poses,h).placements;
    if remaining.is_empty(){let l=layout(g,&shapes,&poses,h);
     if g.validate(&l).is_err() && !original_feasible(g,&l){anyhow::bail!("SNCPH conservative raster disagrees with original-contour feasibility");}
     if best.as_ref().is_none_or(|b|l.length<b.length){println!("SNCPH variant={variant} pattern={iteration} resolution={w}x{h} L={}",l.length);best=Some(l);}}
    if h==final_h{final_board=Some(board);}
   }
   if let Some(board)=final_board{generated_patterns+=1;if current.len()==g.pieces.len(){complete_patterns+=1;}
    for x in 0..final_w{for y in 0..final_h{values[x*final_h+y]+=if (board.bits[x][y/64]>>(y%64))&1==1{-cfg.penalty}else{cfg.price};}}
   }else{break;}
  }
  if !budget.expired(){stats.completed_variants+=1;}
 }
 stats.backend_events.push(format!("SNCPH generated_patterns={generated_patterns}; all-demand-complete_patterns={complete_patterns}; configured_patterns_per_variant={}",cfg.patterns_per_variant));
 stats.completed_generations=Some(complete_patterns);stats.paper_variant=Some(format!("SNCPH twelve configurations;{} price-updated patterns each; literal XML angles",cfg.patterns_per_variant));
 stats.paper_run_status=Some(if stats.completed_variants==12{"completed_SNCPH_documented_conventions"}else{"completed_patterns_SNCPH_project_truncation"}.into());
 best.ok_or_else(||anyhow::anyhow!("SNCPH produced no complete valid pattern; no shelf relabeling"))
}
#[cfg(test)]mod tests{use super::*;use geo::{Area,BooleanOps};use crate::toolkit::polygon::*;
#[test]fn nondefault_pattern_count_and_nested_rasters_control_the_full_flow(){
 use std::collections::BTreeMap;
 let vertices=vec![(0.,0.),(1.,0.),(1.,1.),(0.,1.)];
 let instance=serde_json::from_value(serde_json::json!({"name":"sncph-controls","strip_height":1.,"min_item_separation":0.,"items":[{"id":0,"demand":3,"orientation":{"rotation":{"mode":"discrete","angles":[0.]}},"shape":{"type":"simple_polygon","data":vertices}}]})).unwrap();
 let g=Geometry::from_input(&crate::input::LoadedInput{instance,item_labels:BTreeMap::new(),source_path:Default::default(),exact_vertices:BTreeMap::from([(0,vertices)]),source_up_left:false}).unwrap();
 let mut cfg=crate::config::get().sncph.clone();cfg.k=2;cfg.max_it=0;cfg.patterns_per_variant=2;cfg.raster_heights=vec![2,4];cfg.price=0.2;cfg.penalty=0.3;
 let mut stats=Stats::default();let layout=solve_with_config(&g,&Budget::new(3.),123,&mut stats,&cfg).unwrap();
 assert_eq!(stats.completed_variants,12);assert!(stats.completed_generations.is_some_and(|n|n>0 && n<=24));
 assert!(stats.backend_events.iter().any(|e|e.contains("generated_patterns=24;")));
 assert_eq!(layout.length,3.);g.validate(&layout).unwrap();
 let s=Raster{kind:0,orientation:0,w:1,h:1,ox:0,oy:0,runs:vec![(0,0,1)],pixels:1};
 let mut b=Board::new(2,2,vec![1e13;4]);b.put(&s,0,0,true);b.rebuild();
 assert!(b.phi(&s,0,1,true,true)>b.phi(&s,0,0,true,false));
}
#[test]fn ph1_ph2_and_ph3_respect_their_compactness_rules(){
 let shapes=[Raster{kind:0,orientation:0,w:2,h:1,ox:0,oy:0,runs:vec![(0,0,1),(1,0,1)],pixels:2}];
 for mask in 0u64..128 {
  let mut b=Board::new(7,4,vec![1.;28]);
  for x in 0..7 {if mask&(1<<x)!=0 {b.bits[x][0]|=1<<3;}}
  b.rebuild();let budget=Budget::new(1.);let mut stats=Stats::default();
  let per_column=(0..=5).filter_map(|x|(0..=3).rev().find(|&y|b.fits(&shapes[0],x,y)).map(|y|(x,y,b.phi(&shapes[0],x,y,true,true)))).collect::<Vec<_>>();
  let first=placement(&b,&shapes,0,1,0,true,&budget,&mut stats).unwrap();
  assert_eq!((first.0.x,first.0.y),(per_column[0].0,per_column[0].1));
  let best=placement(&b,&shapes,0,2,0,true,&budget,&mut stats).unwrap();
  assert!(per_column.iter().all(|p|p.2<=best.1));
  let (pose,score)=ph3(&b,&shapes,Pose{shape:0,x:5,y:1,piece:0},true,&budget,&mut stats).unwrap();
  for (dx,dy) in [(-1,0),(1,0),(0,-1),(0,1)] {
   let x=pose.x as i32+dx;let y=pose.y as i32+dy;
   if x>=0 && y>=0 && x<=5 && y<=3 {let x=x as usize;let y=y as usize;
    assert!(b.phi(&shapes[0],x,y,true,b.fits(&shapes[0],x,y))<=score);}
  }
 }
}
#[test]fn raster_envelope_and_phi_use_uncovered_pixel_prices(){
 let v=vec![P::new(0.,0.),P::new(2.7,0.),P::new(0.,2.3)];let r=bitmap::scanline(&v,0.5,&Budget::new(1.)).unwrap().unwrap();let mut union=geo::MultiPolygon(vec![]);
 for c in r{union=union.union(&polygon(&[P::new(c.x0,c.y0),P::new(c.x1,c.y0),P::new(c.x1,c.y1),P::new(c.x0,c.y1)]));}
 assert!(polygon(&v).difference(&union).unsigned_area()<1e-9);
 let s=Raster{kind:0,orientation:0,w:1,h:1,ox:0,oy:0,runs:vec![(0,0,1)],pixels:1};let mut b=Board::new(4,4,vec![1.;16]);
 assert!(b.fits(&s,0,0));b.put(&s,0,0,true);b.rebuild();assert!(!b.fits(&s,0,0));assert_eq!(b.rectangle(0,0,1,1),1);assert_eq!(b.free(0,0,2,2),3);
 assert!(b.phi(&s,0,1,true,true)>b.phi(&s,0,0,true,false));
}
}
