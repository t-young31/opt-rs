// Steepest decent optimisation
use crate::coordinates::{Point, Vector3D};
use crate::ff::forcefield::Forcefield;
use crate::molecule::Molecule;
use log::info;

pub struct SteepestDecentOptimiser {
    alpha: f64,                   // Step size to take in the parameter space
    iteration: usize,             // Current iteration this optimiser is on
    max_num_iterations: usize,    // Maximum number of steps to take
    grad_rms_tolerance: f64,      // Convergence tolerance on |∇E|
    energy_history: Vec<f64>,     // Values of the energy along this optimisation
    init_coordinates: Vec<Point>, // Initial coordinates this optimisation starts with
}

impl SteepestDecentOptimiser {
    /// A default steepest decent optimiser
    pub fn default() -> Self {
        SteepestDecentOptimiser {
            alpha: 0.0001,
            iteration: 0,
            max_num_iterations: 500,
            grad_rms_tolerance: 0.1,
            energy_history: Default::default(),
            init_coordinates: Default::default(),
        }
    }

    /// Optimiser with a defined number of maximum steps
    pub fn from_max_iterations(max_num_iterations: usize) -> Self {
        let mut optimiser = SteepestDecentOptimiser::default();
        optimiser.max_num_iterations = max_num_iterations;

        optimiser
    }

    /// Optimise a molecule with a steepest decent step
    pub fn optimise(&mut self, molecule: &mut Molecule, forcefield: &mut dyn Forcefield) {
        self.cache_initial_coordinates(molecule);

        while self.iteration < self.max_num_iterations {
            if self.needs_an_energy_evaluation() {
                self.energy_history
                    .push(forcefield.energy(&molecule.coordinates));
            }

            if self.energy_is_rising() {
                self.alpha /= 2.; // Reduce the step size by half
                                  //                    and re-initialise the coordinates and saved energies
                self.energy_history.clear();
                molecule.coordinates = self.init_coordinates.clone();
            }

            let gradient = forcefield.gradient(&molecule.coordinates);

            if self.converged(gradient) {
                info!("Converged in {} steps", self.iteration);
                break;
            }

            for (i, v) in gradient.iter().enumerate() {
                molecule.coordinates[i].x -= self.alpha * v.x;
                molecule.coordinates[i].y -= self.alpha * v.y;
                molecule.coordinates[i].z -= self.alpha * v.z;
            }

            self.iteration += 1;
            // println!("E = {}", forcefield.energy(&molecule.coordinates));
        }
    }

    /// Has this optimiser converged?
    fn converged(&self, gradient: &Vec<Vector3D>) -> bool {
        SteepestDecentOptimiser::grad_rms(gradient) < self.grad_rms_tolerance
    }

    /// Root mean square on the gradient
    fn grad_rms(gradient: &Vec<Vector3D>) -> f64 {
        let mut sum_squares = 0.;

        for v in gradient.iter() {
            sum_squares += v.length();
        }

        (sum_squares / gradient.len() as f64).sqrt()
    }

    /// Cache the initial coordinates of the molecule so they can be reverted if
    /// required
    fn cache_initial_coordinates(&mut self, molecule: &Molecule) {
        self.init_coordinates = molecule.coordinates.clone();
    }

    /// Should an energy evaluation be performed? Only calculate every few iterations,
    /// as to not be overly slow
    fn needs_an_energy_evaluation(&self) -> bool {
        self.energy_history.len() < 5
        // && self.iteration % 5 == 0
    }

    /// Is the energy rising in the optimisation?
    fn energy_is_rising(&mut self) -> bool {
        let n = self.energy_history.len();

        if n < 2 {
            // Cannot determine if the energy is rising with only one point
            return false;
        }

        self.energy_history[n - 1] > self.energy_history[n - 2]
    }
}
