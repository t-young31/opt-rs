use std::cmp::Ordering;
use crate::{Forcefield, Molecule};
use crate::atoms::CartesianCoordinate;
use crate::ff::uff::atom_typing::UFFAtomType;
use crate::ff::uff::atom_types::ATOM_TYPES;


#[derive(Default, Debug)]
pub(crate) struct UFF{
    atom_types: Vec<UFFAtomType>
}


impl Forcefield for UFF {

    /// Create a new, bespoke, forcefield for a molecule by setting atom types which
    fn new(molecule: &Molecule) -> Self{

        let mut ff = UFF::default();
        ff.set_atom_types(molecule);

        ff
    }

    /// Generate atom types for each atom in a molecule, defined by the connectivity and charge on
    /// each atom. As each atom may have multiple matches this function should set the best match
    /// to each atom
    fn set_atom_types(&mut self, molecule: &Molecule) {

        let mut match_qualities: [usize; ATOM_TYPES.len()] = [0; ATOM_TYPES.len()];

        for atom in molecule.atoms().iter(){

            for (i, atom_type) in ATOM_TYPES.iter().enumerate(){
                match_qualities[i] = atom_type.match_quality(atom, molecule);
            }

            let best_match = match_qualities
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                .map(|(index, _)| index)
                .expect("Failed to find a best match");

            self.atom_types.push(ATOM_TYPES[best_match].clone());
        }
    }

    fn energy(&self, coordinates: &Vec<CartesianCoordinate>) {
        todo!()
    }

    fn gradient(&self, coordinates: &Vec<CartesianCoordinate>) {
        todo!()
    }
    // TODO
}


