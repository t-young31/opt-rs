use crate::{Forcefield, Molecule};
use crate::atoms::CartesianCoordinate;
use crate::ff::uff::atom_typing::UFFAtomType;

#[derive(Default, Debug)]
pub(crate) struct UFF{
    atom_types: Vec<UFFAtomType>
}


impl Forcefield for UFF {
    fn new(molecule: &Molecule) -> Self{
        todo!();
        UFF::default()
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


