use crate::{Forcefield, Molecule};
use crate::atoms::CartesianCoordinate;

pub(crate) struct UFF;


impl Forcefield for UFF {
    fn new(molecule: &Molecule) -> Self{
        todo!();
        UFF{}
    }

    fn set_atom_types(&self, molecule: &Molecule) {
        todo!()
    }

    fn energy(&self, coordinates: &Vec<CartesianCoordinate>) {
        todo!()
    }

    fn gradient(&self, coordinates: &Vec<CartesianCoordinate>) {
        todo!()
    }
    // TODO
}


