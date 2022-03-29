use std::cmp::Ordering::Equal;
use log::info;
use crate::atoms::{Atom, AtomicNumber, CartesianCoordinate};
use crate::Forcefield;

use crate::io::xyz::XYZFile;

#[derive(Default)]
pub struct Molecule{

    coordinates:      Vec<CartesianCoordinate>,
    atomic_numbers:   Vec<AtomicNumber>,
    connectivity:     Connectivity,
    non_bonded_pairs: Vec<NBPair>
}

impl Molecule{

    /// Create a Molecule from a .xyz filename
    /// # Arguments
    ///
    /// * `filename` - Name of the .xyz file defining the strucure
    ///
    /// # Examples
    /// use mors::Molecule;
    /// let mol = Molecule::from_xyz_file("name");
    ///
    pub fn from_xyz_file(filename: &str) -> Self{

        let xyz_file = XYZFile::new(filename).unwrap();

        let mut molecule = Molecule{
            coordinates:      xyz_file.coordinates,
            atomic_numbers:   xyz_file.atomic_numbers,
            connectivity:     Default::default(),
            non_bonded_pairs: Default::default()
        };

        molecule.add_bonds();
        // molecule.add_angles();
        // molecule.add_dihedrals();
        // molecule.add_non_bonded_pairs();

        molecule
    }

    pub fn set_forcefield(&mut self, ff: Forcefield) -> (){
        // TODO
    }

    pub fn optimise(&mut self){
        // TODO
    }

    pub fn write_xyz_file(&self){
        // TODO
    }

    /// Get a copy of set of atoms associated with this molecule
    fn atoms(&self) -> Vec<Atom>{

        let mut atoms: Vec<Atom> = Default::default();

        for (i, atomic_number) in self.atomic_numbers.iter().enumerate(){

            let coord = self.coordinates.get(i)
                                                           .expect("N_atoms != N_coords");

            atoms.push(Atom{idx: i,
                                  atomic_number: atomic_number.clone(),
                                  coordinate:    coord.clone()});
        }
        atoms
    }

    /// Number of atoms in this molecule
    fn num_atoms(&self) -> usize{
        return self.atomic_numbers.len()
    }

    /// Add bonds between atoms based on their interatomic distance. For example H-H is consdiered
    /// bonded if the r(HH) < 0.8 Å, or so. Note that this is not particularly well defined, but
    /// essential to determining the energy with a 'classic' (i.e. non-density based) force-field
    fn add_bonds(&mut self) {

        self.connectivity.bonds.clear();

        for atom_i in self.atoms().iter(){

            let mut neighbours: Vec<(Atom, f64)> = Default::default();

            for atom_j in self.atoms(){

                if atom_i.could_be_bonded_to(&atom_j){
                    let distance = atom_i.distance_to(&atom_j);
                    neighbours.push((atom_j, distance));
                }
            }

            neighbours.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Equal));

            for k in 0..atom_i.maximal_valence(){

                let possible_neighbour = neighbours.get(k);

                match possible_neighbour {
                    Some(atom) => self.connectivity.bonds.push(
                        Bond{i: atom_i.idx, j: atom.0.idx}),

                    None => break // Exhausted all the present neighbours
                }
            }

        }
    }

}

#[derive(Default)]
struct Connectivity{
    /*
    Connectivity of a molecule defined in terms of 'bonds'. Includes information about the pairs,
    triples and quadruples which define the bonds, angle and dihedral components required to
    calculate the total energy of a molecule with a force-field.
    */
    pub bonds:     Vec<Bond>,
    pub angles:    Vec<Angle>,
    pub dihedrals: Vec<Dihedral>
}

#[derive(Default)]
struct NBPair{
    i: usize,
    j: usize
}

#[derive(Default)]
struct Bond{
    i: usize,
    j: usize
}

#[derive(Default)]
struct Angle{
    i: usize,
    j: usize,
    k: usize
}

#[derive(Default)]
struct Dihedral{
    i: usize,
    j: usize,
    k: usize,
    l: usize
}
