// Steepest decent optimisation
use log::info;
use crate::coordinates::Vector3D;
use crate::{Forcefield, Molecule};


pub struct SteepestDecentOptimiser{
    alpha:              f64,   // Step size to take in the parameter space
    max_num_steps:      usize, // Maximum number of steps to take
    grad_rms_tolerance: f64,   // Convergence tolerance on |∇E|

    gradient:           Vec<Vector3D> // Current gradient
}

impl SteepestDecentOptimiser {

    /// A default steepest decent optimiser
    pub fn default() -> Self{
        SteepestDecentOptimiser{
            alpha: 0.0001,
            max_num_steps: 100,
            grad_rms_tolerance: 1E-4,
            gradient: Default::default()
        }
    }

    /// Optimise a molecule with a steepest decent step
    pub fn optimise(&mut self,
                    molecule:   &mut Molecule,
                    forcefield: &mut dyn Forcefield){

        self.set_zero_gradient(molecule);

        for i in 0..self.max_num_steps{

            self.gradient = molecule.gradient(forcefield);

            if self.converged(){
                info!("Converged in {} steps", i);
                break;
            }

            // println!("E = {}", molecule.energy(forcefield));

            for (i, v) in self.gradient.iter().enumerate(){
                molecule.coordinates[i].x -= self.alpha * v.x;
                molecule.coordinates[i].y -= self.alpha * v.y;
                molecule.coordinates[i].z -= self.alpha * v.z;
            }

        }
    }

    /// Has this optimiser converged?
    fn converged(&self) -> bool{
        self.grad_rms() < self.grad_rms_tolerance
    }

    /// Root mean square on the gradient
    fn grad_rms(&self) -> f64{

        let mut sum_squares = 0.;

        for v in self.gradient.iter(){
            sum_squares += v.length();
        }

        (sum_squares / self.gradient.len() as f64).sqrt()
    }

    /// Set a zero gradient for which will be updated
    fn set_zero_gradient(&mut self, molecule: &mut Molecule){

        self.gradient.clear();
        for _ in 0..molecule.num_atoms(){
            self.gradient.push(Vector3D::default());
        }
    }

}
