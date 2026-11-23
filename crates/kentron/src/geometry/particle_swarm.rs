//! A deterministic global search for the minimum of a function over a box: a regular grid
//! sampled first, its best points then moved as a particle swarm with constriction
//! coefficients (Kennedy & Eberhart 1995; Clerc & Kennedy 2002), the random factors from a
//! fixed-seed xorshift generator. The answer is a starting point for a local method.

/// The swarm over a box.
pub(crate) struct ParticleSwarm {
    low: Vec<f64>,
    high: Vec<f64>,
    steps: Vec<f64>,
    particle_count: usize,
}

/// The constriction factor χ and the acceleration coefficients c₁ = c₂ for φ = 4.1.
const CHI: f64 = 0.729_843_788_128_357_3;
const ACCELERATION: f64 = 2.05;

/// Grid samples at most (the steps are widened beyond).
const GRID_LIMIT: f64 = 1e6;

/// Swarm updates before the swarm may stop for having slowed down.
const MIN_UPDATES: usize = 10;

/// Marsaglia's xorshift64* generator, fixed seed.
struct Random(u64);

impl Random {
    /// Uniform in [0, 1).
    fn next(&mut self) -> f64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        (x.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 11) as f64 / (1u64 << 53) as f64
    }
}

struct Particle {
    position: Vec<f64>,
    velocity: Vec<f64>,
    best_position: Vec<f64>,
    best_value: Option<f64>,
}

impl ParticleSwarm {
    pub(crate) fn new(low_border: Vec<f64>, upper_border: Vec<f64>, steps: Vec<f64>, particle_count: usize) -> Self {
        assert!(!low_border.is_empty() && low_border.len() == upper_border.len() && low_border.len() == steps.len(), "a box and a step for each coordinate");
        ParticleSwarm { low: low_border, high: upper_border, steps, particle_count }
    }

    /// The least value of `function` found and where (grid spacing and velocity scale
    /// `steps`, at most `iteration_count − 1` swarm updates).
    pub(crate) fn perform(&mut self, function: &mut dyn FnMut(&[f64]) -> Option<f64>, steps: &[f64], iteration_count: usize) -> (f64, Vec<f64>) {
        let dim = self.steps.len();
        self.steps.copy_from_slice(&steps[..dim]);
        // The box pulled in a little from its borders.
        let (lo, hi): (Vec<f64>, Vec<f64>) = (0..dim)
            .map(|k| {
                let (a, b) = (self.low[k].min(self.high[k]), self.low[k].max(self.high[k]));
                let pull = 1e-6 * (b - a);
                (a + pull, b - pull)
            })
            .unzip();
        let mut step: Vec<f64> = (0..dim).map(|k| self.steps[k].abs().max(1e-6 * (hi[k] - lo[k])).max(f64::MIN_POSITIVE)).collect();
        let counts = |step: &[f64]| -> Vec<usize> { (0..dim).map(|k| ((hi[k] - lo[k]) / step[k]).floor().max(0.0) as usize + 1).collect() };
        let mut count = counts(&step);
        let total = count.iter().map(|&c| c as f64).product::<f64>();
        if total > GRID_LIMIT {
            let widen = (total / GRID_LIMIT).powf(1.0 / dim as f64);
            for s in &mut step {
                *s *= widen;
            }
            count = counts(&step);
        }
        // The grid, keeping the best samples.
        let mut samples: Vec<(f64, Vec<f64>)> = Vec::new();
        let mut index = vec![0usize; dim];
        loop {
            let point: Vec<f64> = (0..dim).map(|k| (lo[k] + step[k] * index[k] as f64).min(hi[k])).collect();
            if let Some(v) = function(&point).filter(|v| !v.is_nan()) {
                samples.push((v, point));
            }
            let mut k = 0;
            while k < dim {
                index[k] += 1;
                if index[k] < count[k] {
                    break;
                }
                index[k] = 0;
                k += 1;
            }
            if k == dim {
                break;
            }
        }
        samples.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut random = Random(0x9e37_79b9_7f4a_7c15);
        let mut particles: Vec<Particle> = (0..self.particle_count.max(1))
            .map(|i| {
                let (position, value) = match samples.get(i) {
                    Some((v, p)) => (p.clone(), Some(*v)),
                    None => ((0..dim).map(|k| lo[k] + (hi[k] - lo[k]) * random.next()).collect(), None),
                };
                let velocity = (0..dim).map(|k| step[k] * (random.next() - 0.5)).collect();
                Particle { best_position: position.clone(), position, velocity, best_value: value }
            })
            .collect();
        let mut best: Option<(f64, Vec<f64>)> = samples.first().cloned();
        for update in 1..iteration_count {
            let mut fastest: f64 = 0.0;
            for particle in &mut particles {
                for k in 0..dim {
                    let (r1, r2) = (random.next(), random.next());
                    let own = particle.best_position[k] - particle.position[k];
                    let global = best.as_ref().map_or(0.0, |b| b.1[k] - particle.position[k]);
                    let v = CHI * (particle.velocity[k] + ACCELERATION * r1 * own + ACCELERATION * r2 * global);
                    let limit = hi[k] - lo[k];
                    let v = v.clamp(-limit, limit);
                    particle.velocity[k] = v;
                    particle.position[k] = (particle.position[k] + v).clamp(lo[k], hi[k]);
                    fastest = fastest.max(v.abs() / step[k]);
                }
                let Some(value) = function(&particle.position).filter(|v| !v.is_nan()) else { continue };
                if particle.best_value.is_none_or(|b| value < b) {
                    particle.best_value = Some(value);
                    particle.best_position.clone_from(&particle.position);
                }
                if best.as_ref().is_none_or(|b| value < b.0) {
                    best = Some((value, particle.position.clone()));
                }
            }
            if update >= MIN_UPDATES && fastest < 1e-3 {
                break;
            }
        }
        best.unwrap_or_else(|| (f64::INFINITY, lo.clone()))
    }
}

/// `b` if `a < b`, else `a` (so `a` when either is NaN or both are zeros).
pub(crate) fn std_max(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}

/// `b` if `b < a`, else `a`.
pub(crate) fn std_min(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_function(x: &[f64]) -> f64 {
        x.iter().enumerate().map(|(i, &v)| (5.0 * v).sin() + 0.1 * v * v + 0.05 * (i + 1) as f64 * v).sum()
    }

    #[test]
    fn one_and_two_dimensions() {
        let mut swarm = ParticleSwarm::new(vec![-2.0], vec![2.0], vec![0.1], 8);
        let (v, p) = swarm.perform(&mut |x| Some(test_function(x)), &[0.1], 100);
        assert!((v - -1.00584).abs() < 1e-3, "{v} at {p:?}");
        assert!((p[0] + 0.3137).abs() < 0.05);
        assert_eq!(v, test_function(&p));
        // Two dimensions: no worse than the grid.
        for particles in [13, 32] {
            let mut grid_best = f64::INFINITY;
            let mut swarm = ParticleSwarm::new(vec![-2.0, -1.0], vec![2.0, 3.0], vec![0.2, 0.2], particles);
            let (v, p) = swarm.perform(&mut |x| {
                let v = test_function(x);
                grid_best = grid_best.min(v);
                Some(v)
            }, &[0.2, 0.2], 100);
            assert!(v <= grid_best && v == test_function(&p));
        }
    }

    #[test]
    fn deeper_basin_and_determinism() {
        let f = |x: &[f64]| Some(-(-(x[0] - 0.2).powi(2) * 50.0).exp() - 2.0 * (-(x[0] - 0.8).powi(2) * 50.0).exp());
        let mut swarm = ParticleSwarm::new(vec![0.0], vec![1.0], vec![0.05], 8);
        let first = swarm.perform(&mut |x| f(x), &[0.05], 100);
        assert!((first.1[0] - 0.8).abs() < 0.05, "{first:?}");
        let second = swarm.perform(&mut |x| f(x), &[0.05], 100);
        assert_eq!(first, second);
    }

    #[test]
    fn bowl_and_holes() {
        let mut swarm = ParticleSwarm::new(vec![-1.0, -1.0], vec![1.0, 1.0], vec![0.1, 0.1], 16);
        let (v, p) = swarm.perform(&mut |x| Some((x[0] - 0.33).powi(2) + (x[1] + 0.21).powi(2)), &[0.1, 0.1], 100);
        assert!((p[0] - 0.33).abs() < 0.1 && (p[1] + 0.21).abs() < 0.1 && v < 0.01);
        // Undefined on the left half: never reported there.
        let mut swarm = ParticleSwarm::new(vec![-1.0], vec![1.0], vec![0.1], 8);
        let (v, p) = swarm.perform(&mut |x| (x[0] > 0.0).then(|| (x[0] - 0.5).powi(2)), &[0.1], 100);
        assert!(p[0] > 0.0 && v == (p[0] - 0.5).powi(2));
    }

    #[test]
    fn std_comparisons() {
        assert_eq!(std_max(1.0, 2.0), 2.0);
        assert!(std_max(f64::NAN, 1.0).is_nan());
        assert_eq!(std_max(1.0, f64::NAN), 1.0);
        assert_eq!(std_max(0.0, -0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(std_max(-0.0, 0.0).to_bits(), (-0.0f64).to_bits());
        assert_eq!(std_min(2.0, 1.0), 1.0);
        assert!(std_min(f64::NAN, 1.0).is_nan());
        assert_eq!(std_min(-0.0, 0.0).to_bits(), (-0.0f64).to_bits());
    }
}
