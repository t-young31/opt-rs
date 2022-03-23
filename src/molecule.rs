use crate::atoms::{AtomicNumber, CartesianCoordinate};
use crate::Forcefield;

use crate::io::xyz::XYZFile;

#[derive(Default)]
pub struct Molecule{

    coordinates: Vec<CartesianCoordinate>,
    atomic_numbers: Vec<AtomicNumber>,
    connectivity: Connectivity,
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

        molecule
    }

    pub fn set_forcefield(&mut self, ff: Forcefield) -> (){
        // TODO
    }

    pub fn optimise(&mut self) -> (){
        // TODO
    }

    pub fn write_xyz_file(&self) -> (){
        // TODO
    }

}



#[derive(Default)]
struct Connectivity {
    /*
    Connectivity of a molecule defined in terms of 'bonds'. Includes information about the pairs,
    triples and quadruples which define the bonds, angle and dihedral components required to
    calculate the total energy of a molecule with a force-field.
    */
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub dihedrals: Vec<Dihedral>
}


#[derive(Default)]
struct NBPair{
    i: u32,
    j: u32
}

#[derive(Default)]
struct Bond{
    i: u32,
    j: u32
}

#[derive(Default)]
struct Angle{
    i: u32,
    j: u32,
    k: u32
}

#[derive(Default)]
struct Dihedral{
    i: u32,
    j: u32,
    k: u32,
    l: u32
}
