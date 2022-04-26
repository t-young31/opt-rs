use std::collections::HashSet;
use crate::coordinates::Point;
use crate::ff::forcefield::EnergyFunction;

pub struct TorsionalDihedral{

    pub(crate) i:  usize,
    pub(crate) j:  usize,
    pub(crate) k:  usize,
    pub(crate) l:  usize,

    pub(crate) phi0:  f64,
    pub(crate) n_phi: f64,
    pub(crate) v_phi: f64
}

impl EnergyFunction for TorsionalDihedral {
    fn involves_idxs(&self, idxs: Vec<usize>) -> bool {
        HashSet::from([self.i, self.j, self.k, self.l]) == HashSet::from_iter(idxs)
    }

    fn force_constant(&self) -> f64 {
        self.v_phi
    }

    fn energy(&self, coordinates: &Vec<Point>) -> f64 {

        todo!();
        let phi = 0.;

        0.5 * self.v_phi * (1. - (self.n_phi * self.phi0).cos()*(self.n_phi * phi).cos())
    }

    fn add_gradient(&self,
                    coordinates:      &Vec<Point>,
                    current_gradient: &mut Vec<Point>) {
        todo!()
    }
}

