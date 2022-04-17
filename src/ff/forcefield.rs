use crate::atoms::CartesianCoordinate;
use crate::Molecule;

pub trait Forcefield {

    fn new(molecule: &Molecule) -> Self where Self: Sized;

    fn set_atom_types(&self, molecule: &Molecule);

    fn energy(&self, coordinates: &Vec<CartesianCoordinate>);

    fn gradient(&self, coordinates: &Vec<CartesianCoordinate>);
}
