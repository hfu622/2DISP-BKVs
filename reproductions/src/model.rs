use serde::{Deserialize, Serialize};
use std::time::{Duration, Instant};

#[derive(Clone, Copy, Debug, Serialize)]
pub enum Family {
    Constructive,
    MonteCarlo,
    Local,
    Anneal,
    Tabu,
    Genetic,
    RandomKey,
    Jostle,
    Beam,
    Swarm,
    Cuckoo,
    Whale,
    Shake,
    Guided,
    Coordinate,
    Hybrid,
    Pairwise,
    Parallel,
    Bounds,
}

#[derive(Clone, Copy, Debug, Serialize)]
pub enum PlacementRule {
    NfpBottomLeft,
    ToposNfp,
    RayDrop,
    BottomLeft,
    Contact,
    Gravity,
    Hull,
    Potential,
    Raster,
    SemiDiscrete,
    Column,
}

#[derive(Clone, Copy, Debug, Serialize)]
pub enum OrderRule {
    PaperDefined,
    Area,
    Height,
    Width,
    Perimeter,
    Input,
    Compactness,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq)]
pub enum Backend {
    None,
    Lp,
    Mip,
    Ipopt,
    LocalBfgs,
    CpuParallel,
}

#[derive(Clone, Copy, Debug, Serialize)]
pub struct Spec {
    pub id: &'static str,
    pub title: &'static str,
    pub family: Family,
    pub placement: PlacementRule,
    pub order: OrderRule,
    pub backend: Backend,
    pub neighborhood: u8,
    pub population: usize,
    pub beam_width: usize,
    pub passes: usize,
    pub flow: &'static str,
    pub deviations: &'static str,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Placement {
    pub piece: usize,
    pub orientation: usize,
    // Translation of a rotated, bounding-box-normalized polygon.
    pub x: f64,
    pub y: f64,
}

#[derive(Clone, Debug)]
pub struct Layout {
    pub placements: Vec<Placement>,
    pub length: f64,
}

pub struct Budget {
    start: Instant,
    duration: Duration,
}
impl Budget {
    pub fn new(seconds: f64) -> Self {
        Self {
            start: Instant::now(),
            duration: Duration::from_secs_f64(seconds),
        }
    }
    pub fn expired(&self) -> bool {
        self.start.elapsed() >= self.duration
    }
    pub fn remaining(&self) -> f64 {
        self.duration
            .saturating_sub(self.start.elapsed())
            .as_secs_f64()
    }
    pub fn elapsed(&self) -> f64 {
        self.start.elapsed().as_secs_f64()
    }
}

#[derive(Default, Serialize)]
pub struct Stats {
    pub effective_algorithm_controls: Option<serde_json::Value>,
    pub evaluations: usize,
    pub collision_queries: usize,
    pub accepted_moves: usize,
    pub backend_calls: usize,
    pub nonoptimal_backend_calls: usize,
    pub backend_events: Vec<String>,
    pub paper_run_status: Option<String>,
    pub paper_variant: Option<String>,
    pub completed_variants: usize,
    pub nfp_cache_bytes: usize,
    pub population_initialized: Option<usize>,
    pub completed_generations: Option<usize>,
}

// Small reproducible PRNG, independent of the platform's thread scheduler.
pub struct Rng(u64);
impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed.max(1))
    }
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }
    pub fn index(&mut self, n: usize) -> usize {
        (self.next_u64() % n.max(1) as u64) as usize
    }
    pub fn shuffle<T>(&mut self, a: &mut [T]) {
        for i in (1..a.len()).rev() {
            let j = self.index(i + 1);
            a.swap(i, j);
        }
    }
}
