use crate::coordinates::{Point, Vector3D};
use crate::ff::forcefield::EnergyFunction;
use crate::pairs::distance;
use std::collections::HashSet;

pub struct HarmonicBond {
    pub(crate) i: usize,
    pub(crate) j: usize,
    pub(crate) r0: f64,
    pub(crate) k_ij: f64,
}

impl EnergyFunction for HarmonicBond {
    fn involves_idxs(&self, idxs: Vec<usize>) -> bool {
        idxs.len() == 2 && HashSet::from([self.i, self.j]) == HashSet::from_iter(idxs)
    }

    fn force_constant(&self) -> f64 {
        self.k_ij
    }

    /// Energy: k/2 (r-r_0)^2
    fn energy(&self, coordinates: &[Point]) -> f64 {
        self.k_ij / 2.0 * (distance(self.i, self.j, coordinates) - self.r0).powi(2)
    }

    /// Add the gradient for this term
    fn add_gradient(&self, coordinates: &[Point], gradient: &mut Vec<Vector3D>) {
        let x_i = coordinates[self.i].x;
        let y_i = coordinates[self.i].y;
        let z_i = coordinates[self.i].z;

        let x_j = coordinates[self.j].x;
        let y_j = coordinates[self.j].y;
        let z_j = coordinates[self.j].z;

        let v0 = x_i - x_j;
        let v1 = (z_i - z_j).powi(2);
        let v2 = (y_i - y_j).powi(2);
        let v3 = v0.powi(2);
        let v4 = v3 + v2 + v1;
        let v5 = v4.sqrt();
        let v6 = self.k_ij * (-self.r0 + v5);

        gradient[self.i].x += v6 * v0 / v5;
        gradient[self.i].y += v6 * (y_i - y_j) / v5;
        gradient[self.i].z += v6 * (z_i - z_j) / v5;
        gradient[self.j].x += v6 * (-x_i + x_j) / v5;
        gradient[self.j].y += v6 * (-y_i + y_j) / v5;
        gradient[self.j].z += v6 * (-z_i + z_j) / v5
    }
}
