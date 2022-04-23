use std::collections::HashSet;
use crate::coordinates::Point;
use crate::ff::forcefield::EnergyFunction;
use crate::pairs::distance;


pub struct HarmonicBond{
    pub(crate) i:  usize,
    pub(crate) j:  usize,
    pub(crate) r0: f64,
    pub(crate) k_ij:  f64
}


impl EnergyFunction for HarmonicBond {

    fn involves_idxs(&self, idxs: Vec<usize>) -> bool {
        idxs.len() == 2 && HashSet::from([self.i, self.j]) == HashSet::from_iter(idxs)
    }

    fn force_constant(&self) -> f64 { self.k_ij }

    /// Energy: k/2 (r-r_0)^2
    fn energy(&self, coordinates: &Vec<Point>) -> f64 {
        self.k_ij /2.0 * (distance(self.i, self.j, coordinates) - self.r0).powi(2)
    }

    /// Add the gradient for this term
    fn add_gradient(&self,
                    x:        &Vec<Point>,
                    gradient: &mut Vec<Point>){

        let r = distance(self.i, self.j, x);

        for k in 0..3{
            let val = self.k_ij * (1.0 - self.r0 / r) * (x[self.i][k] - x[self.j][k]);

            gradient[self.i][k] += val;
            gradient[self.j][k] -= val;
        }
    }
}
