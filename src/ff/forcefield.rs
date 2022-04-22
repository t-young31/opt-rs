use crate::coordinates::Point;
use crate::Molecule;

pub trait Forcefield {

    fn new(molecule: &Molecule) -> Self where Self: Sized;

    fn set_atom_types(&mut self, molecule: &Molecule);

    fn energy(&mut self, coordinates: &Vec<Point>) -> f64;

    fn gradient(&mut self, coordinates: &Vec<Point>) -> &Vec<Point>;
}


pub trait EnergyFunction{

    fn involves_idxs(&self, idxs: Vec<usize>) -> bool;

    fn force_constant(&self) -> f64;

    fn energy(&self, coordinates: &Vec<Point>) -> f64;

    fn add_gradient(&self,
                    coordinates:      &Vec<Point>,
                    current_gradient: &mut Vec<Point>);
}
