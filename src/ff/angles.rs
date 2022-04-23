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

        gradient[self.i].x += -self.k_ijk*((-x_i + x_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()) + (-x_j + x_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()))*(self.n*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt())).acos()).sin()/(self.n*(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2))) + 1.).sqrt());
        gradient[self.i].y += -self.k_ijk*((-y_i + y_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()) + (-y_j + y_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()))*(self.n*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt())).acos()).sin()/(self.n*(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2))) + 1.).sqrt());
        gradient[self.i].z += -self.k_ijk*((-z_i + z_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()) + (-z_j + z_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()))*(self.n*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt())).acos()).sin()/(self.n*(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2))) + 1.).sqrt());
        gradient[self.j].x += -self.k_ijk*((x_i - x_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()) + (-x_j + x_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).powf(1.5)) + (-x_i + 2.*x_j - x_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()))*(self.n*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt())).acos()).sin()/(self.n*(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2))) + 1.).sqrt());
        gradient[self.j].y += -self.k_ijk*((y_i - y_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()) + (-y_j + y_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).powf(1.5)) + (-y_i + 2.*y_j - y_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()))*(self.n*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt())).acos()).sin()/(self.n*(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2))) + 1.).sqrt());
        gradient[self.j].z += -self.k_ijk*((z_i - z_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()) + (-z_j + z_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).powf(1.5)) + (-z_i + 2.*z_j - z_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()))*(self.n*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt())).acos()).sin()/(self.n*(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2))) + 1.).sqrt());
        gradient[self.k].x += -self.k_ijk*((x_i - x_j)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()) + (x_j - x_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).powf(1.5)))*(self.n*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt())).acos()).sin()/(self.n*(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2))) + 1.).sqrt());
        gradient[self.k].y += -self.k_ijk*((y_i - y_j)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()) + (y_j - y_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).powf(1.5)))*(self.n*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt())).acos()).sin()/(self.n*(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2))) + 1.).sqrt());
        gradient[self.k].z += -self.k_ijk*((z_i - z_j)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt()) + (z_j - z_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).powf(1.5)))*(self.n*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2)).sqrt())).acos()).sin()/(self.n*(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((x_j - x_k).powi(2) + (y_j - y_k).powi(2) + (z_j - z_k).powi(2))) + 1.).sqrt());

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

        gradient[self.i].x += self.k_ijk*(self.c1*(-x_i + x_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + self.c1*(-x_j + x_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + 2.*self.c2*((-x_i + x_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + (-x_j + x_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()))*(2.*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt())).acos()).sin()/(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2))) + 1.).sqrt());
        gradient[self.i].y += self.k_ijk*(self.c1*(-y_i + y_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + self.c1*(-y_j + y_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + 2.*self.c2*((-y_i + y_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + (-y_j + y_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()))*(2.*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt())).acos()).sin()/(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2))) + 1.).sqrt());
        gradient[self.i].z += self.k_ijk*(self.c1*(-z_i + z_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + self.c1*(-z_j + z_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + 2.*self.c2*((-z_i + z_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + (-z_j + z_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()))*(2.*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt())).acos()).sin()/(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2))) + 1.).sqrt());
        gradient[self.j].x += self.k_ijk*(self.c1*(x_i - x_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + self.c1*(-x_j + x_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)) + self.c1*(-x_i + 2.*x_j - x_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + 2.*self.c2*((x_i - x_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + (-x_j + x_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)) + (-x_i + 2.*x_j - x_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()))*(2.*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt())).acos()).sin()/(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2))) + 1.).sqrt());
        gradient[self.j].y += self.k_ijk*(self.c1*(y_i - y_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + self.c1*(-y_j + y_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)) + self.c1*(-y_i + 2.*y_j - y_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + 2.*self.c2*((y_i - y_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + (-y_j + y_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)) + (-y_i + 2.*y_j - y_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()))*(2.*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt())).acos()).sin()/(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2))) + 1.).sqrt());
        gradient[self.j].z += self.k_ijk*(self.c1*(z_i - z_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + self.c1*(-z_j + z_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)) + self.c1*(-z_i + 2.*z_j - z_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + 2.*self.c2*((z_i - z_j)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).powf(1.5)*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + (-z_j + z_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)) + (-z_i + 2.*z_j - z_k)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()))*(2.*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt())).acos()).sin()/(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2))) + 1.).sqrt());
        gradient[self.k].x += self.k_ijk*(self.c1*(x_i - x_j)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + self.c1*(x_j - x_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)) + 2.*self.c2*((x_i - x_j)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + (x_j - x_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)))*(2.*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt())).acos()).sin()/(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2))) + 1.).sqrt());
        gradient[self.k].y += self.k_ijk*(self.c1*(y_i - y_j)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + self.c1*(y_j - y_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)) + 2.*self.c2*((y_i - y_j)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + (y_j - y_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)))*(2.*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt())).acos()).sin()/(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2))) + 1.).sqrt());
        gradient[self.k].z += self.k_ijk*(self.c1*(z_i - z_j)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + self.c1*(z_j - z_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)) + 2.*self.c2*((z_i - z_j)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt()) + (z_j - z_k)*((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).powf(1.5)))*(2.*(((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k))/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2)).sqrt()*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2)).sqrt())).acos()).sin()/(-((x_i - x_j)*(-x_j + x_k) + (y_i - y_j)*(-y_j + y_k) + (z_i - z_j)*(-z_j + z_k)).powi(2)/(((x_i - x_j).powi(2) + (y_i - y_j).powi(2) + (z_i - z_j).powi(2))*((-x_j + x_k).powi(2) + (-y_j + y_k).powi(2) + (-z_j + z_k).powi(2))) + 1.).sqrt());

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

