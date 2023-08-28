use crate::coordinates::{Point, Vector3D};
use crate::molecule::Molecule;

pub trait Forcefield {
    fn new(molecule: &Molecule) -> Self
    where
        Self: Sized;

    fn energy(&mut self, coordinates: &[Point]) -> f64;

    fn gradient(&mut self, coordinates: &[Point]) -> &Vec<Vector3D>;
}

pub trait EnergyFunction {
    fn involves_idxs(&self, idxs: Vec<usize>) -> bool;

    fn force_constant(&self) -> f64;

    fn energy(&self, coordinates: &[Point]) -> f64;

    fn add_gradient(&self, coordinates: &[Point], gradient: &mut Vec<Vector3D>);
}
