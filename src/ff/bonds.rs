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
                    x:        &Vec<CartesianCoordinate>,
                    gradient: &mut Vec<CartesianCoordinate>){

        let r = distance(self.i, self.j, x);

        for k in 0..3{
            let val = self.k * (1.0 - self.r0 / r) * (x[self.i][k] - x[self.j][k]);

            gradient[self.i][k] += val;
            gradient[self.j][k] -= val;
        }
    }
}
