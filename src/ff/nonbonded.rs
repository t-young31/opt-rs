use crate::coordinates::{Point, Vector3D};
use crate::ff::forcefield::EnergyFunction;
use crate::pairs::{AtomPair, distance};


pub struct LennardJones12x6{
    pub(crate) i:     usize,
    pub(crate) j:     usize,
    pub(crate) sigma: f64,
    pub(crate) d:     f64
}


impl EnergyFunction for LennardJones12x6 {

    fn involves_idxs(&self, idxs: Vec<usize>) -> bool {
        idxs.len() == 2 && AtomPair{i: self.i, j: self.j} == AtomPair{i: idxs[0], j: idxs[1]}
    }

    fn force_constant(&self) -> f64 { self.d }

    /// Energy: D_ij ((σ/r)^12 - 2(σ/r)^6)
    fn energy(&self, coordinates: &Vec<Point>) -> f64 {
        let r = distance(self.i, self.j, coordinates);

        self.d * ((self.sigma /r).powi(12) - 2. * (self.sigma /r).powi(6))
    }

    /// Add the gradient for this term
    fn add_gradient(&self,
                    coordinates: &Vec<Point>,
                    gradient:    &mut Vec<Vector3D>){

        let x_i = coordinates[self.i].x;
        let y_i = coordinates[self.i].y;
        let z_i = coordinates[self.i].z;

        let x_j = coordinates[self.j].x;
        let y_j = coordinates[self.j].y;
        let z_j = coordinates[self.j].z;

        let v0 = (z_i - z_j).powi(2);
        let v1 = (y_i - y_j).powi(2);
        let v2 = (x_i - x_j).powi(2);
        let v3 = v2 + v1 + v0;
        let v4 = v3.powi(7);
        let v5 = self.sigma.powi(12);
        let v6 = 2. * self.sigma.powi(6);
        let v7 = v3.powi(4);

        gradient[self.i].x += self.d*(v5*(-12.*x_i + 12.*x_j)/v4 - v6*(-6.*x_i + 6.*x_j)/v7);
        gradient[self.i].y += self.d*(v5*(-12.*y_i + 12.*y_j)/v4 - v6*(-6.*y_i + 6.*y_j)/v7);
        gradient[self.i].z += self.d*(v5*(-12.*z_i + 12.*z_j)/v4 - v6*(-6.*z_i + 6.*z_j)/v7);
        gradient[self.j].x += self.d*(v5*(12.*x_i  - 12.*x_j)/v4 - v6*(6.*x_i  - 6.*x_j)/v7);
        gradient[self.j].y += self.d*(v5*(12.*y_i  - 12.*y_j)/v4 - v6*(6.*y_i  - 6.*y_j)/v7);
        gradient[self.j].z += self.d*(v5*(12.*z_i  - 12.*z_j)/v4 - v6*(6.*z_i  - 6.*z_j)/v7);
    }
}
