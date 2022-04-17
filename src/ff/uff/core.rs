use crate::{Forcefield, Molecule};
use crate::atoms::CartesianCoordinate;

pub(crate) struct UFF{
    atom_types: Vec<UFFAtomType>
}


impl Forcefield for UFF {
    pub(crate) fn new(molecule: &Molecule) -> Self{
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


