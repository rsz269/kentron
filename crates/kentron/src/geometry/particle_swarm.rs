//! OCCT's global minimizer `math_PSO` (particle swarm optimization) with its particle pool
//! `math_PSOParticlesPool` and random generator `math_BullardGenerator`, ported.

/// `aBorderDivisor`: the search box is shrunk by this fraction of its size on each side.
const BORDER_DIVISOR: f64 = 1.0e+4;

/// `math_BullardGenerator`: OCCT's small deterministic pseudo-random generator.
#[doc(alias = "math_BullardGenerator")]
pub(crate) struct BullardGenerator {
    state_hi: u32,
    state_lo: u32,
}

impl BullardGenerator {
    /// `math_BullardGenerator(theSeed)`.
    pub(crate) fn new(seed: u32) -> Self {
        BullardGenerator { state_hi: seed, state_lo: seed ^ 0x49616E42 }
    }

    #[doc(alias = "math_BullardGenerator::NextInt")]
    pub(crate) fn next_int(&mut self) -> u32 {
        self.state_hi = (self.state_hi >> 2).wrapping_add(self.state_hi << 2);
        self.state_hi = self.state_hi.wrapping_add(self.state_lo);
        self.state_lo = self.state_lo.wrapping_add(self.state_hi);
        self.state_hi
    }

    /// `NextReal()`: in [0, 1].
    #[doc(alias = "math_BullardGenerator::NextReal")]
    pub(crate) fn next_real(&mut self) -> f64 {
        self.next_int() as f64 / 0xFFFFFFFFu32 as f64
    }
}

/// `std::max(a, b)`: `b` only where `a < b` (unlike `f64::max` for NaNs and signed zeros).
pub(crate) fn std_max(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}

/// `std::min(a, b)`: `b` only where `b < a` (unlike `f64::min` for NaNs and signed zeros).
pub(crate) fn std_min(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

/// `PSO_Particle`.
#[doc(alias = "PSO_Particle")]
struct Particle {
    position: Vec<f64>,
    velocity: Vec<f64>,
    best_position: Vec<f64>,
    distance: f64,
    best_distance: f64,
}

/// `math_PSOParticlesPool`.
#[doc(alias = "math_PSOParticlesPool")]
struct ParticlesPool {
    particles: Vec<Particle>,
}

impl ParticlesPool {
    fn new(particle_count: usize, dimension_count: usize) -> Self {
        let particles = (0..particle_count)
            .map(|_| Particle {
                position: vec![0.0; dimension_count],
                velocity: vec![0.0; dimension_count],
                best_position: vec![0.0; dimension_count],
                distance: f64::MAX,
                best_distance: f64::MAX,
            })
            .collect();
        ParticlesPool { particles }
    }

    /// `GetBestParticle()`: `std::min_element` by distance (the first of equal ones).
    fn best_index(&self) -> usize {
        let mut best = 0;
        for (index, particle) in self.particles.iter().enumerate().skip(1) {
            if particle.distance < self.particles[best].distance {
                best = index;
            }
        }
        best
    }

    /// `GetWorstParticle()`: `std::max_element` by distance (the first of equal ones).
    fn worst_index(&self) -> usize {
        let mut worst = 0;
        for (index, particle) in self.particles.iter().enumerate().skip(1) {
            if self.particles[worst].distance < particle.distance {
                worst = index;
            }
        }
        worst
    }
}

/// `math_PSO`: the minimum of a function of several variables over a box, sought by a
/// particle swarm seeded from the best points of a regular grid.
#[doc(alias = "math_PSO")]
pub(crate) struct ParticleSwarm {
    low_border: Vec<f64>,
    upper_border: Vec<f64>,
    steps: Vec<f64>,
    particle_count: usize,
}

impl ParticleSwarm {
    /// `math_PSO(theFunc, theLowBorder, theUppBorder, theSteps, theNbParticles, theNbIter)`
    /// (the function is passed to [`perform`](Self::perform) instead, and the iteration
    /// count there too, as `Perform` overwrites the constructor's).
    pub(crate) fn new(low_border: Vec<f64>, upper_border: Vec<f64>, steps: Vec<f64>, particle_count: usize) -> Self {
        assert!(low_border.len() == upper_border.len() && low_border.len() == steps.len() && !steps.is_empty());
        ParticleSwarm { low_border, upper_border, steps, particle_count }
    }

    fn dimension(&self) -> usize {
        self.low_border.len()
    }

    /// The box the particles move in: `myLowBorder + (myUppBorder - myLowBorder) /
    /// aBorderDivisor` to `myUppBorder - (...)`, with `math_Vector`'s division by
    /// multiplication with the reciprocal.
    fn shrunk_box(&self) -> (Vec<f64>, Vec<f64>) {
        let reciprocal = 1.0 / BORDER_DIVISOR;
        let min_uv = (0..self.dimension())
            .map(|k| self.low_border[k] + (self.upper_border[k] - self.low_border[k]) * reciprocal)
            .collect();
        let max_uv = (0..self.dimension())
            .map(|k| self.upper_border[k] - (self.upper_border[k] - self.low_border[k]) * reciprocal)
            .collect();
        (min_uv, max_uv)
    }

    /// `math_PSO::Perform(theSteps, theValue, theOutPnt, theNbIter)`: the best value found
    /// and its point. `function` returns `None` where OCCT's `Value` returns false; OCCT then
    /// leaves the value it was writing into unchanged, and so does this.
    #[doc(alias = "math_PSO::Perform")]
    pub(crate) fn perform(&mut self, function: &mut dyn FnMut(&[f64]) -> Option<f64>, steps: &[f64], iteration_count: usize) -> (f64, Vec<f64>) {
        let n = self.dimension();
        let (min_uv, max_uv) = self.shrunk_box();
        self.steps = steps.to_vec();
        let mut pool = ParticlesPool::new(self.particle_count, n);

        // OCCT leaves `aCurrValue` uninitialized; it is always written by the first call,
        // whose point lies inside the box.
        let mut current_value = f64::MAX;
        let mut current_point = min_uv.clone();
        let mut worst = pool.worst_index();
        loop {
            if let Some(value) = function(&current_point) {
                current_value = value;
            }
            if current_value < pool.particles[worst].distance {
                let particle = &mut pool.particles[worst];
                particle.position.copy_from_slice(&current_point);
                particle.best_position.copy_from_slice(&current_point);
                particle.distance = current_value;
                particle.best_distance = current_value;
                worst = pool.worst_index();
            }
            current_point[0] += std_max(self.steps[0], 1.0e-15); // Avoid too small step
            for k in 1..n {
                if current_point[k - 1] > max_uv[k - 1] {
                    current_point[k - 1] = min_uv[k - 1];
                    current_point[k] += self.steps[k];
                } else {
                    break;
                }
            }
            if current_point[n - 1] > max_uv[n - 1] {
                break;
            }
        }
        self.perform_with_given_particles(function, &mut pool, iteration_count)
    }

    /// `math_PSO::performPSOWithGivenParticles`.
    fn perform_with_given_particles(&mut self, function: &mut dyn FnMut(&[f64]) -> Option<f64>, pool: &mut ParticlesPool, iteration_count: usize) -> (f64, Vec<f64>) {
        let n = self.dimension();
        let (min_uv, max_uv) = self.shrunk_box();
        let mut current_point = vec![0.0; n];
        let mut random = BullardGenerator::new(1);
        for particle in pool.particles.iter_mut().take(self.particle_count) {
            for k in 0..n {
                let ksi = random.next_real();
                particle.velocity[k] = self.steps[k] * (ksi - 0.5) * 2.0;
            }
        }

        let best = &pool.particles[pool.best_index()];
        let mut best_global_position = best.position.clone();
        let mut best_global_distance = best.distance;

        let termination_velocity: Vec<f64> = self.steps.iter().map(|step| step * (1.0 / 2048.0)).collect();
        let mut minimal_velocity = vec![0.0; n];

        const RETENTION_WEIGHT: f64 = 0.72900;
        const PERSONAL_WEIGHT: f64 = 1.49445;
        const SOCIAL_WEIGHT: f64 = 1.49445;
        const MIN_STEPS: usize = 16;

        for step in 1..iteration_count {
            minimal_velocity.fill(f64::MAX);
            for particle in pool.particles.iter_mut().take(self.particle_count) {
                let ksi1 = random.next_real();
                let ksi2 = random.next_real();
                for k in 0..n {
                    particle.velocity[k] = particle.velocity[k] * RETENTION_WEIGHT
                        + (particle.best_position[k] - particle.position[k]) * (PERSONAL_WEIGHT * ksi1)
                        + (best_global_position[k] - particle.position[k]) * (SOCIAL_WEIGHT * ksi2);
                    particle.position[k] += particle.velocity[k];
                    particle.position[k] = std_min(std_max(particle.position[k], min_uv[k]), max_uv[k]);
                    current_point[k] = particle.position[k];
                    minimal_velocity[k] = std_min(particle.velocity[k].abs(), minimal_velocity[k]);
                }
                if let Some(value) = function(&current_point) {
                    particle.distance = value;
                }
                if particle.distance < particle.best_distance {
                    particle.best_distance = particle.distance;
                    particle.best_position.copy_from_slice(&particle.position);
                    if particle.distance < best_global_distance {
                        best_global_distance = particle.distance;
                        best_global_position.copy_from_slice(&particle.position);
                    }
                }
            }

            let terminal_velocity_reached = !(0..n).any(|k| minimal_velocity[k] > termination_velocity[k]);
            if terminal_velocity_reached {
                if step > MIN_STEPS {
                    break;
                }
                for particle in pool.particles.iter_mut().take(self.particle_count) {
                    let ksi = random.next_real();
                    for k in 0..n {
                        if particle.position[k] == min_uv[k] || particle.position[k] == max_uv[k] {
                            particle.velocity[k] = if particle.position[k] == min_uv[k] { self.steps[k] * ksi } else { -self.steps[k] * ksi };
                        } else {
                            particle.velocity[k] = self.steps[k] * (ksi - 0.5) * 2.0;
                        }
                    }
                }
            }
        }
        (best_global_distance, best_global_position)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `math_BullardGenerator` against OCCT's (`OCP.math.math_BullardGenerator`).
    #[test]
    fn bullard_generator_matches_occt() {
        let mut random = BullardGenerator::new(1);
        let ints: Vec<u32> = (0..5).map(|_| random.next_int()).collect();
        assert_eq!(ints, [1231121991, 3399545143, 3129986794, 819383309, 408636399]);
        assert_eq!(BullardGenerator::new(7).next_real(), 0.28664293193413015);
    }

    /// Against OCCT 8.1's `math_PSO` on `sum_i sin(5 x_i) + 0.1 x_i^2 + 0.05 i x_i` over
    /// `[-3 + 0.1 i, 2.5 + 0.2 i]` with steps `0.37 i` (1-based `i`): `(dimension, particle
    /// count)` to `(value, point)`.
    #[test]
    fn matches_occt() {
        type Case<'a> = ((usize, usize), (f64, &'a [f64]));
        let cases: [Case; 6] = [
            ((1, 8), (-1.0058056112545743, &[-0.3119593991665256])),
            ((1, 13), (-1.005841302449184, &[-0.3134898526731087])),
            ((1, 32), (-1.0058416144068603, &[-0.313680215903205])),
            ((2, 8), (-2.027412751525781, &[-0.3140727317002643, -0.31546616614938644])),
            ((2, 13), (-1.8244002599376243, &[-0.312427315727073, 0.9309182532981642])),
            ((2, 32), (-2.02698747204527, &[-0.31920277244145717, -0.31740286871291157])),
        ];
        for ((n, particle_count), (value, point)) in cases {
            let index = |k: usize| (k + 1) as f64;
            let low: Vec<f64> = (0..n).map(|k| -3.0 + 0.1 * index(k)).collect();
            let upper: Vec<f64> = (0..n).map(|k| 2.5 + 0.2 * index(k)).collect();
            let steps: Vec<f64> = (0..n).map(|k| 0.37 * index(k)).collect();
            let mut function = |x: &[f64]| {
                let mut f = 0.0;
                for (k, &xk) in x.iter().enumerate() {
                    f += (5.0 * xk).sin() + 0.1 * xk * xk + 0.05 * index(k) * xk;
                }
                Some(f)
            };
            let mut swarm = ParticleSwarm::new(low, upper, steps.clone(), particle_count);
            let (got_value, got_point) = swarm.perform(&mut function, &steps, 100);
            assert_eq!((got_value, got_point.as_slice()), (value, point), "{n} {particle_count}");
        }
    }
}
