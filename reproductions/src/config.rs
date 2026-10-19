//! One compiled, reviewable configuration record; never read a mutable file mid-run.
use anyhow::{Result, ensure};
use serde::{Deserialize,Serialize};
use std::sync::OnceLock;
pub const RECORD: &str = include_str!("../../config/reproduction_decisions.json");
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Sncph {
    pub k: usize, pub max_it: usize, pub patterns_per_variant: usize,
    pub price: f64, pub penalty: f64,
    pub panel_height_multiple: usize, pub raster_heights: Vec<usize>,
}
impl Sncph {pub fn validate(&self)->Result<()> {
    ensure!(self.k>0 && self.k<=10000 && self.max_it<=100 && self.patterns_per_variant>0 && self.patterns_per_variant<=10000,"invalid SNCPH counters");
    ensure!(self.price.is_finite() && self.price>0. && self.price<=100. && self.penalty.is_finite() && self.penalty>0. && self.penalty<=100.,"invalid SNCPH prices");
    ensure!(!self.raster_heights.is_empty() && self.raster_heights.len()<=8 && (1..=16).contains(&self.panel_height_multiple),"invalid SNCPH panel");
    ensure!(self.raster_heights.iter().all(|&h|h>0 && h<=2048) && self.raster_heights.windows(2).all(|p|p[1]>p[0] && p[1]%p[0]==0),"raster levels must nest exactly");
    Ok(())
}}
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Gcs {
    pub nests: usize, pub generations: usize, pub abandonment: f64, pub stagnation: usize,
    pub levy_beta: f64, pub boundary: String, pub expansion_base: String,
    pub shrink: f64, pub expand: f64,
}
impl Gcs {pub fn validate(&self)->Result<()> {
    ensure!(self.nests>0 && self.generations>0 && self.stagnation>0,"invalid GCS counters");
    ensure!(self.levy_beta>0. && self.levy_beta<2. && (0.0..=1.0).contains(&self.abandonment) && self.shrink>0. && self.shrink<1. && self.expand>0. && self.expand.is_finite(),"invalid GCS rates");
    ensure!(["clip","resample"].contains(&self.boundary.as_str()) && ["current","best"].contains(&self.expansion_base.as_str()),"unsupported GCS convention");Ok(())
}}
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ThreePm {
    pub initial_spacing: f64, pub minimum_spacing: f64, pub refinement: f64,
    pub maximum_dots_per_reference: usize, pub first_square_half_side_factor: f64,
    pub reference_vertex_index: usize, pub batch_rule: String,
}
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Chs {pub initial_pieces:usize,pub mip_seconds:f64,
    #[serde(default)] pub cbc_preprocess:bool,
    #[serde(default)] pub cbc_presolve:bool}
#[derive(Clone, Debug, Deserialize)]
pub struct Configuration {pub sncph:Sncph,pub gcs:Gcs,pub three_pm:ThreePm,pub chs2:Chs}
pub fn get()->&'static Configuration {
    static CFG:OnceLock<Configuration>=OnceLock::new();
    CFG.get_or_init(||{let c:Configuration=serde_json::from_str(RECORD).expect("embedded configuration");
        c.sncph.validate().expect("SNCPH configuration");c.gcs.validate().expect("GCS configuration");
        assert!(c.three_pm.initial_spacing>=c.three_pm.minimum_spacing && c.three_pm.minimum_spacing>0. && c.three_pm.refinement>0. && c.three_pm.refinement<1. && c.three_pm.maximum_dots_per_reference>0 && c.three_pm.first_square_half_side_factor>0.);
        assert!(["printed_algorithm_4","prose_divisor_priority"].contains(&c.three_pm.batch_rule.as_str()));
        assert!(c.chs2.initial_pieces>0 && c.chs2.mip_seconds.is_finite() && c.chs2.mip_seconds>0.);c})
}
pub fn record()->serde_json::Value {let _=get();serde_json::from_str(RECORD).unwrap()}
#[cfg(test)]mod tests {use super::*;
    #[test]fn embedded_policy_is_valid_and_invalid_controls_are_rejected(){
        let c=get();assert_eq!(c.sncph.raster_heights,[128,256,512]);
        let mut s=c.sncph.clone();s.k=0;assert!(s.validate().is_err());s.k=3;s.raster_heights=vec![100,256];assert!(s.validate().is_err());
        let mut g=c.gcs.clone();g.levy_beta=2.;assert!(g.validate().is_err());g.levy_beta=1.5;g.expansion_base="unknown".into();assert!(g.validate().is_err());
    }
}
