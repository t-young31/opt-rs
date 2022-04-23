use std::collections::HashSet;
use crate::coordinates::{Point, Vector3D};
use crate::ff::forcefield::EnergyFunction;


trait HarmonicAngle{
    fn i(&self) -> usize;
    fn j(&self) -> usize;
    fn k(&self) -> usize;
}


pub struct HarmonicAngleTypeA{
    pub(crate) i:     usize,
    pub(crate) j:     usize,
    pub(crate) k:     usize,
    pub(crate) k_ijk: f64,
    pub(crate) n:     f64
}


impl EnergyFunction for HarmonicAngleTypeA {
    fn involves_idxs(&self, idxs: Vec<usize>) -> bool { involves_idxs(self, idxs) }

    fn force_constant(&self) -> f64 { self.k_ijk }

    /// Energy: k/n^2 (1-cos(nθ))
    fn energy(&self, coordinates: &Vec<Point>) -> f64 {

        let theta = theta(self, coordinates);

        (self.k_ijk / self.n.powi(2)) * (1.0 - (self.n * theta).cos())
    }

    /// Add the gradient for this term
    fn add_gradient(&self,
                    coordinates: &Vec<Point>,
                    gradient:    &mut Vec<Point>){

        let x_i = coordinates[self.i].x;
        let y_i = coordinates[self.i].y;
        let z_i = coordinates[self.i].z;

        let x_j = coordinates[self.j].x;
        let y_j = coordinates[self.j].y;
        let z_j = coordinates[self.j].z;

        let x_k = coordinates[self.k].x;
        let y_k = coordinates[self.k].y;
        let z_k = coordinates[self.k].z;

        let v0 = (z_i - z_j);
        let v1 = (x_i - x_j);
        let v2 = (y_i - y_j);
        let v3 = (-y_j + y_k);
        let v4 = (-z_j + z_k);
        let v5 = (-x_j + x_k);
        let v6 = (v1*v5 + v2*v3 + v0*v4);
        let v7 = (v1.powi(2) + v2.powi(2) + v0.powi(2));
        let v8 = (v5.powi(2) + v3.powi(2) + v4.powi(2));
        let v9 = (v7*v8);
        let v10 = (v7.sqrt()*v8.sqrt());
        let v11 = (v6/v10);
        let v12 = (self.n*v11.acos());
        let v13 = (-v6.powi(2)/v9 + 1.);

        gradient[self.i].x += -self.k_ijk*((-x_i + x_j)*v6/(v7.powf(1.5)*v8.sqrt()) + v5/v10)*v12.sin()/(self.n*v13.sqrt());
        gradient[self.i].y += -self.k_ijk*((-y_i + y_j)*v6/(v7.powf(1.5)*v8.sqrt()) + v3/v10)*v12.sin()/(self.n*v13.sqrt());
        gradient[self.i].z += -self.k_ijk*((-z_i + z_j)*v6/(v7.powf(1.5)*v8.sqrt()) + v4/v10)*v12.sin()/(self.n*v13.sqrt());
        gradient[self.j].x += -self.k_ijk*(v1*v6/(v7.powf(1.5)*v8.sqrt()) + v5*v6/(v7.sqrt()*v8.powf(1.5)) + (-x_i + 2.*x_j - x_k)/v10)*v12.sin()/(self.n*v13.sqrt());
        gradient[self.j].y += -self.k_ijk*(v2*v6/(v7.powf(1.5)*v8.sqrt()) + v3*v6/(v7.sqrt()*v8.powf(1.5)) + (-y_i + 2.*y_j - y_k)/v10)*v12.sin()/(self.n*v13.sqrt());
        gradient[self.j].z += -self.k_ijk*(v0*v6/(v7.powf(1.5)*v8.sqrt()) + v4*v6/(v7.sqrt()*v8.powf(1.5)) + (-z_i + 2.*z_j - z_k)/v10)*v12.sin()/(self.n*v13.sqrt());
        gradient[self.k].x += -self.k_ijk*(v1/v10 + (x_j - x_k)*v6/(v7.sqrt()*v8.powf(1.5)))*v12.sin()/(self.n*v13.sqrt());
        gradient[self.k].y += -self.k_ijk*(v2/v10 + (y_j - y_k)*v6/(v7.sqrt()*v8.powf(1.5)))*v12.sin()/(self.n*v13.sqrt());
        gradient[self.k].z += -self.k_ijk*(v0/v10 + (z_j - z_k)*v6/(v7.sqrt()*v8.powf(1.5)))*v12.sin()/(self.n*v13.sqrt())
    }
}
impl HarmonicAngle for HarmonicAngleTypeA {
    fn i(&self) -> usize {self.i}
    fn j(&self) -> usize {self.j}
    fn k(&self) -> usize {self.k}
}

pub struct HarmonicAngleTypeB{
    pub(crate) i:     usize,
    pub(crate) j:     usize,
    pub(crate) k:     usize,
    pub(crate) k_ijk: f64,
    pub(crate) c0:    f64,
    pub(crate) c1:    f64,
    pub(crate) c2:    f64
}

impl EnergyFunction for HarmonicAngleTypeB {

    fn involves_idxs(&self, idxs: Vec<usize>) -> bool { involves_idxs(self, idxs) }

    fn force_constant(&self) -> f64 { self.k_ijk }

    /// Energy: k(c0 + c1 cos(θ) + c2 cos(2θ))
    fn energy(&self, coordinates: &Vec<Point>) -> f64 {

       let theta = theta(self, coordinates);

        self.k_ijk * (self.c0 + self.c1 * theta.cos() + self.c2 * (2. * theta).cos())
    }

    /// Add the gradient for this term
    fn add_gradient(&self,
                    coordinates: &Vec<Point>,
                    gradient:    &mut Vec<Point>){

        let x_i = coordinates[self.i].x;
        let y_i = coordinates[self.i].y;
        let z_i = coordinates[self.i].z;

        let x_j = coordinates[self.j].x;
        let y_j = coordinates[self.j].y;
        let z_j = coordinates[self.j].z;

        let x_k = coordinates[self.k].x;
        let y_k = coordinates[self.k].y;
        let z_k = coordinates[self.k].z;

        let v0 = (x_i - x_j);
        let v1 = (y_i - y_j);
        let v2 = (z_i - z_j);
        let v3 = (-z_j + z_k);
        let v4 = (-x_j + x_k);
        let v5 = (-y_j + y_k);
        let v6 = (v0*v4 + v1*v5 + v2*v3);
        let v7 = (v0.powi(2) + v1.powi(2) + v2.powi(2));
        let v8 = (v4.powi(2) + v5.powi(2) + v3.powi(2));
        let v9 = (v7*v8);
        let v10 = (-v6.powi(2)/v9 + 1.);
        let v11 = (v7.sqrt()*v8.sqrt());
        let v12 = (v6/v11);
        let v13 = (2.*v12.acos());

        gradient[self.i].x += self.k_ijk*(self.c1*(-x_i + x_j)*v6/(v7.powf(1.5)*v8.sqrt()) + self.c1*v4/v11 + 2.*self.c2*((-x_i + x_j)*v6/(v7.powf(1.5)*v8.sqrt()) + v4/v11)*v13.sin()/v10.sqrt());
        gradient[self.i].y += self.k_ijk*(self.c1*(-y_i + y_j)*v6/(v7.powf(1.5)*v8.sqrt()) + self.c1*v5/v11 + 2.*self.c2*((-y_i + y_j)*v6/(v7.powf(1.5)*v8.sqrt()) + v5/v11)*v13.sin()/v10.sqrt());
        gradient[self.i].z += self.k_ijk*(self.c1*(-z_i + z_j)*v6/(v7.powf(1.5)*v8.sqrt()) + self.c1*v3/v11 + 2.*self.c2*((-z_i + z_j)*v6/(v7.powf(1.5)*v8.sqrt()) + v3/v11)*v13.sin()/v10.sqrt());
        gradient[self.j].x += self.k_ijk*(self.c1*v0*v6/(v7.powf(1.5)*v8.sqrt()) + self.c1*v4*v6/(v7.sqrt()*v8.powf(1.5)) + self.c1*(-x_i + 2.*x_j - x_k)/v11 + 2.*self.c2*(v0*v6/(v7.powf(1.5)*v8.sqrt()) + v4*v6/(v7.sqrt()*v8.powf(1.5)) + (-x_i + 2.*x_j - x_k)/v11)*v13.sin()/v10.sqrt());
        gradient[self.j].y += self.k_ijk*(self.c1*v1*v6/(v7.powf(1.5)*v8.sqrt()) + self.c1*v5*v6/(v7.sqrt()*v8.powf(1.5)) + self.c1*(-y_i + 2.*y_j - y_k)/v11 + 2.*self.c2*(v1*v6/(v7.powf(1.5)*v8.sqrt()) + v5*v6/(v7.sqrt()*v8.powf(1.5)) + (-y_i + 2.*y_j - y_k)/v11)*v13.sin()/v10.sqrt());
        gradient[self.j].z += self.k_ijk*(self.c1*v2*v6/(v7.powf(1.5)*v8.sqrt()) + self.c1*v3*v6/(v7.sqrt()*v8.powf(1.5)) + self.c1*(-z_i + 2.*z_j - z_k)/v11 + 2.*self.c2*(v2*v6/(v7.powf(1.5)*v8.sqrt()) + v3*v6/(v7.sqrt()*v8.powf(1.5)) + (-z_i + 2.*z_j - z_k)/v11)*v13.sin()/v10.sqrt());
        gradient[self.k].x += self.k_ijk*(self.c1*v0/v11 + self.c1*(x_j - x_k)*v6/(v7.sqrt()*v8.powf(1.5)) + 2.*self.c2*(v0/v11 + (x_j - x_k)*v6/(v7.sqrt()*v8.powf(1.5)))*v13.sin()/v10.sqrt());
        gradient[self.k].y += self.k_ijk*(self.c1*v1/v11 + self.c1*(y_j - y_k)*v6/(v7.sqrt()*v8.powf(1.5)) + 2.*self.c2*(v1/v11 + (y_j - y_k)*v6/(v7.sqrt()*v8.powf(1.5)))*v13.sin()/v10.sqrt());
        gradient[self.k].z += self.k_ijk*(self.c1*v2/v11 + self.c1*(z_j - z_k)*v6/(v7.sqrt()*v8.powf(1.5)) + 2.*self.c2*(v2/v11 + (z_j - z_k)*v6/(v7.sqrt()*v8.powf(1.5)))*v13.sin()/v10.sqrt())

    }
}
impl HarmonicAngle for HarmonicAngleTypeB {
    fn i(&self) -> usize {self.i}
    fn j(&self) -> usize {self.j}
    fn k(&self) -> usize {self.k}
}

fn involves_idxs(angle: &dyn HarmonicAngle, idxs: Vec<usize>) -> bool {
    idxs.len() == 3 && HashSet::from([angle.i(), angle.j(), angle.k()]) == HashSet::from_iter(idxs)
}



/// Value of the angle between three atoms (i, j, k)
#[inline(always)]
fn theta(angle:       &dyn HarmonicAngle,
         coordinates: &Vec<Point>) -> f64{

    let r1: Vector3D = &coordinates[angle.i()] - &coordinates[angle.j()];
    let r2: Vector3D = &coordinates[angle.k()] - &coordinates[angle.j()];

    (r1.dot(&r2) / (r1.length() * r2.length())).acos()
}


/*
   /$$                           /$$
  | $$                          | $$
 /$$$$$$    /$$$$$$   /$$$$$$$ /$$$$$$   /$$$$$$$
|_  $$_/   /$$__  $$ /$$_____/|_  $$_/  /$$_____/
  | $$    | $$$$$$$$|  $$$$$$   | $$   |  $$$$$$
  | $$ /$$| $$_____/ \____  $$  | $$ /$$\____  $$
  |  $$$$/|  $$$$$$$ /$$$$$$$/  |  $$$$//$$$$$$$/
   \___/   \_______/|_______/    \___/ |_______/
 */

#[cfg(test)]
mod tests{

    use super::*;
    use crate::utils::*;

    struct TestHarmonicAngle{i: usize, j: usize, k: usize}

    impl HarmonicAngle for TestHarmonicAngle{
        fn i(&self) -> usize {self.i}
        fn j(&self) -> usize {self.j}
        fn k(&self) -> usize {self.k}
    }

    /// Ensure a bend angle can be calculated correctly. Compared to an Avogadro reference
    #[test]
    fn test_angle_calculation(){

        // Coordinates for a H2O molecule
        let x: Vec<Point> = Vec::from([
            Point{x: -2.17305, y: -1.64172, z: 0.01027},
            Point{x: -1.18373, y: -1.60261, z: -0.00725},
            Point{x: -2.46242, y: -0.77749, z: -0.37700}
        ]);

        let angle = theta(&TestHarmonicAngle{i: 1, j: 0, k: 2}, &x);

        assert!(is_close(angle, 1.8239, 1E-3));
    }

}

