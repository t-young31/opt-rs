use crate::coordinates::CartesianCoordinate;
use crate::ff::forcefield::EnergyFunction;
use crate::pairs::distance;

pub struct HarmonicAngleTypeA{
    pub(crate) i:  usize,
    pub(crate) j:  usize,
    pub(crate) k:  usize,
    pub(crate) k_ijk: f64,
    pub(crate) n:  f64
}


impl EnergyFunction for HarmonicAngleTypeA {

    /// Energy: k/2 (r-r_0)^2
    fn energy(&self, coordinates: &Vec<CartesianCoordinate>) -> f64 {

        let theta = (r1.dot(r2) / (distance(r1) * distance(r2))).arccos();

        (self.k_ijk / self.n.powi(2)) * (1.0 - (self.n * theta).cos())
    }

    /// Add the gradient for this term
    fn add_gradient(&self,
                    x:        &Vec<CartesianCoordinate>,
                    gradient: &mut Vec<CartesianCoordinate>){

        todo!()
    }
}
