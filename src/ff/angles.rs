use std::collections::HashSet;
use crate::coordinates::{Point, Vector3D};
use crate::ff::forcefield::EnergyFunction;


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
                    x:        &Vec<Point>,
                    gradient: &mut Vec<Point>){

        todo!()
    }
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
                    x:        &Vec<Point>,
                    gradient: &mut Vec<Point>){

        todo!()
    }
}


trait HarmonicAngle{
    fn i(&self) -> usize;
    fn j(&self) -> usize;
    fn k(&self) -> usize;
}

fn involves_idxs(angle: &dyn HarmonicAngle, idxs: Vec<usize>) -> bool {
    idxs.len() == 3 && HashSet::from([angle.i(), angle.j(), angle.k()]) == HashSet::from_iter(idxs)
}

impl HarmonicAngle for HarmonicAngleTypeA {
    fn i(&self) -> usize {self.i}
    fn j(&self) -> usize {self.j}
    fn k(&self) -> usize {self.k}
}
impl HarmonicAngle for HarmonicAngleTypeB {
    fn i(&self) -> usize {self.i}
    fn j(&self) -> usize {self.j}
    fn k(&self) -> usize {self.k}
}


/// Value of the angle between three atoms (i, j, k)
#[inline(always)]
fn theta(angle:       &dyn HarmonicAngle,
         coordinates: &Vec<Point>) -> f64{

    let r1: Vector3D = &coordinates[angle.i()] - &coordinates[angle.j()];
    let r2: Vector3D = &coordinates[angle.k()] - &coordinates[angle.j()];

    (r1.dot(&r2) / (r1.length() * r2.length())).acos()
}
