// Steepest decent optimisation
use log::info;
use crate::coordinates::{Point, Vector3D};
use crate::Molecule;


pub struct SteepestDecentOptimiser{
    alpha:              f64,   // Step size to take in the parameter space
    max_num_steps:      usize, // Maximum number of steps to take
    grad_rms_tolerance: f64,   // Convergence tolerance on |∇E|

    gradient:           Vec<Vector3D> // Current gradient
}

impl SteepestDecentOptimiser {

    /// Optimise a molecule with a steepest decent step
    pub fn optimise(&mut self, molecule: &mut Molecule){

        self.set_zero_gradient(molecule);

        for i in 0..self.max_num_steps{

            if self.converged(){
                info!("Converged in {} steps", i);
                break;
            }

        }

    }

    /// Has this optimiser converged?
    fn converged(&self) -> bool{
        self.grad_rms() < self.grad_rms_tolerance
    }

    /// Root mean square on the gradient
    fn grad_rms(&self) -> f64{
        todo!()
    }

    /// Set a zero gradient for which will be updated
    fn set_zero_gradient(&mut self, molecule: &mut Molecule){

        self.gradient.clear();
        for _ in 0..molecule.num_atoms(){
            self.gradient.push(Vector3D::default());
        }
    }

}
