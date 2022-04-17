use crate::atoms::CartesianCoordinate;
use crate::ff::forcefield::EnergyFunction;
use crate::pairs::distance;


pub struct HarmonicBond{
    pub(crate) i:  usize,
    pub(crate) j:  usize,
    pub(crate) r0: f64,
    pub(crate) k:  f64
}


impl EnergyFunction for HarmonicBond {

    /// Energy: k/2 (r-r_0)^2
    fn energy(&self, coordinates: &Vec<CartesianCoordinate>) -> f64 {
        self.k/2.0 * (distance(self.i, self.j, coordinates) - self.r0).powi(2)
    }

    /// Add the gradient for this term
    fn add_gradient(&self,
                    coordinates:      &Vec<CartesianCoordinate>,
                    current_gradient: &Vec<CartesianCoordinate>) -> &Vec<CartesianCoordinate> {
        todo!()
    }
}
