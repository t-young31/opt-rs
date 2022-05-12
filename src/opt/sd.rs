// Steepest decent optimisation
use log::info;
use crate::coordinates::Vector3D;
use crate::{Forcefield, Molecule};


pub struct SteepestDecentOptimiser{
    alpha:              f64,   // Step size to take in the parameter space
    max_num_steps:      usize, // Maximum number of steps to take
    grad_rms_tolerance: f64,   // Convergence tolerance on |∇E|
}

impl SteepestDecentOptimiser {

    /// A default steepest decent optimiser
    pub fn default() -> Self{
        SteepestDecentOptimiser{
            alpha:              0.0001,
            max_num_steps:      1000,
            grad_rms_tolerance: 0.1,
        }
    }

    /// Optimise a molecule with a steepest decent step
    pub fn optimise(&mut self,
                    molecule:   &mut Molecule,
                    forcefield: &mut dyn Forcefield){

        for i in 0..self.max_num_steps{

            let gradient = forcefield.gradient(&molecule.coordinates);

            if self.converged(gradient){
                info!("Converged in {} steps", i);
                break;
            }

            // println!("E = {},  RMS(grad){}", molecule.energy(forcefield), self.grad_rms());

            for (i, v) in gradient.iter().enumerate(){
                molecule.coordinates[i].x -= self.alpha * v.x;
                molecule.coordinates[i].y -= self.alpha * v.y;
                molecule.coordinates[i].z -= self.alpha * v.z;
            }

        }
    }

    /// Has this optimiser converged?
    fn converged(&self, gradient: &Vec<Vector3D>) -> bool{
        SteepestDecentOptimiser::grad_rms(gradient) < self.grad_rms_tolerance
    }

    /// Root mean square on the gradient
    fn grad_rms(gradient: &Vec<Vector3D>) -> f64{

        let mut sum_squares = 0.;

        for v in gradient.iter(){
            sum_squares += v.length();
        }

        (sum_squares / gradient.len() as f64).sqrt()
    }

}
