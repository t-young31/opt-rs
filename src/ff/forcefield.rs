use crate::atoms::CartesianCoordinate;
use crate::Molecule;

pub trait Forcefield {

    fn new(molecule: &Molecule) -> Self where Self: Sized;

    fn set_atom_types(&mut self, molecule: &Molecule);

    fn energy(&mut self, coordinates: &Vec<CartesianCoordinate>) -> f64;

    fn gradient(&mut self, coordinates: &Vec<CartesianCoordinate>) -> &Vec<CartesianCoordinate>;
}


pub trait EnergyFunction{
    fn energy(&self, coordinates: &Vec<CartesianCoordinate>) -> f64;

    fn add_gradient(&self,
                    coordinates:      &Vec<CartesianCoordinate>,
                    current_gradient: &mut Vec<CartesianCoordinate>);
}
