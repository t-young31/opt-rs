use std::cmp::Ordering::Equal;
use std::collections::HashSet;
use log::info;
use crate::atoms::{Atom, AtomicNumber, CartesianCoordinate};
use crate::Forcefield;

use crate::io::xyz::XYZFile;

pub struct Molecule{

    coordinates:      Vec<CartesianCoordinate>,
    atomic_numbers:   Vec<AtomicNumber>,
    connectivity:     Connectivity,
    non_bonded_pairs: Vec<NBPair>
}

impl Molecule{

    /// A blank molecule containing no atoms
    pub fn blank() -> Self{ Molecule::from_atomic_symbols(&[]) }

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

        Molecule::from_atomic_nums_and_coords(xyz_file.atomic_numbers, xyz_file.coordinates)
    }

    /// Construct a molecule from atomic symbols, where all atomic positions are at the origin
    /// and an exception is thrown if the atomic number cannot be created
    pub fn from_atomic_symbols(symbols: &[&str]) -> Self{

        let mut atomic_numbers: Vec<AtomicNumber> = Default::default();
        let mut coords: Vec<CartesianCoordinate> = Default::default();

        for symbol in symbols{

            atomic_numbers.push(AtomicNumber::from_string(symbol).unwrap());
            coords.push(CartesianCoordinate::default());
        }

        Molecule::from_atomic_nums_and_coords(atomic_numbers, coords)
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

    /// Construct a molecule from a set of atomic numbers and coordinates of each atom
    fn from_atomic_nums_and_coords(atomic_numbers: Vec<AtomicNumber>,
                                   coordinates:    Vec<CartesianCoordinate>) -> Self{

        let mut molecule = Molecule{coordinates,
            atomic_numbers,
            connectivity:     Default::default(),
            non_bonded_pairs: Default::default()};

        molecule.add_bonds();
        molecule.add_angles();
        // self.add_dihedrals();
        // self.add_non_bonded_pairs();
        molecule
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

    /// Does this molecule have any associated bonds?
    fn has_bonds(&self) -> bool{ return !self.connectivity.bonds.is_empty() }

    /// Does this molecule contain a particular bond?
    fn contains_bond(&self, bond: &Bond) -> bool{
        self.connectivity.bonds.iter().any(|b| b==bond)
    }

    fn contains_bond_between(&self, i: usize, j: usize) -> bool{
        self.contains_bond(&Bond{i, j})
    }

    fn has_dihedrals(&self) -> bool{ !self.connectivity.dihedrals.is_empty()}

    /// Add bonds between atoms based on their interatomic distance. For example H-H is consdiered
    /// bonded if the r(HH) < 0.8 Å, or so. Note that this is not particularly well defined, but
    /// essential to determining the energy with a 'classic' (i.e. non-density based) force-field
    fn add_bonds(&mut self) {

        self.connectivity.bonds.clear();

        for central_atom in self.atoms().iter(){

            let neighbours = Neighbours::from_atom_and_molecule(central_atom, self);

            for k in 0..central_atom.maximal_valence(){

                match neighbours.get(k) {
                    Some(neighbour) => self.add_bond(central_atom, &neighbour.atom),
                    None => break    // Exhausted all the present neighbours
                }
            }
        }
    }

    /// Add a bond between two atoms, provided it is not already present
    fn add_bond(&mut self, atom_i: &Atom, atom_j: &Atom){

        let new_bond = Bond{i: atom_i.idx, j: atom_j.idx};

        for bond in self.connectivity.bonds.iter(){
            if bond == &new_bond{ return; }
        }

        self.connectivity.bonds.insert(new_bond);
    }

    /// Add angles between triples of mutually bonded atoms. For example,
    ///
    ///     j     k
    ///     /----
    ///    /
    ///   i
    ///
    fn add_angles(&mut self){

        if self.num_atoms() < 3 || !self.has_bonds(){
            return;   // No MM angles in molecules with < 3 atoms or that doesn't have any bonds
        }

        for j in 0..self.num_atoms(){
            for i in j+1..self.num_atoms(){

                if !self.contains_bond_between(i, j){ continue; }

                for k in i+1..self.num_atoms(){

                    if !self.contains_bond_between(j, k){ continue; }

                    self.connectivity.angles.insert(Angle{i, j, k});
                }
            }
        }
    }

    /// Add all the proper dihedrals between quadruples pf bonded atoms
    // TODO: Improper dihedrals
    fn add_dihedrals(&mut self){

        if self.num_atoms() < 4{
            return;   // No dihedrals in molecules with < 4 atoms
        }

        if !self.has_bonds(){
            panic!("Cannot add dihedrals without any bonds present")
        }

        let all_neighbours: Vec<Vec<usize>> = self.atoms()
                .iter()
                .map(|a| a.bonded_neighbour_idxs(&self.connectivity.bonds))
                .collect();

        for bond in self.connectivity.bonds.iter(){
            // Add bonded quadruples

            for neighbour_i in all_neighbours.get(bond.i).unwrap(){
                if neighbour_i == &bond.j{ continue; }

                for neighbour_j in all_neighbours.get(bond.j).unwrap(){
                    if neighbour_j == &bond.i{ continue; }

                    self.connectivity.dihedrals.insert(
                        Dihedral{i: neighbour_i.clone(),
                            j: bond.i,
                            k: bond.j,
                            l: neighbour_j.clone()}
                    );
                } // neighbour_j
            }  // neighbour_i
        } // bond

    }
}

#[derive(Default, Debug)]
struct Neighbour{
    atom:     Atom,
    distance: f64
}

#[derive(Default, Debug)]
struct Neighbours{
    values: Vec<Neighbour>
}

impl Neighbours {

    /// Given an atom contained in a molecule, generate the bonded neighbours of that atom
    /// sorted by the distance to the central atom
    pub fn from_atom_and_molecule(central_atom: &Atom, molecule: &Molecule) -> Self{

        let mut neighbours: Neighbours = Default::default();

        neighbours.add_all_while_guessing_bonds(central_atom, molecule);
        neighbours.sort_by_distance();
        neighbours
    }

    /// Get a neighbour, defined by its index
    pub fn get(&self, index: usize) -> Option<&Neighbour>{ self.values.get(index) }

    /// Add neighbours to a central atom while guessing the bonds
    fn add_all_while_guessing_bonds(&mut self, central_atom: &Atom, molecule: &Molecule){

        for atom in molecule.atoms(){

            if central_atom.could_be_bonded_to(&atom){
                let distance = central_atom.distance_to(&atom);
                self.values.push(Neighbour{atom, distance});
            }
        }
    }

    /// Sort these neighbours by the distance to the central atom
    fn sort_by_distance(&mut self){
        self.values.sort_by(|a, b|
                             a.distance.partial_cmp(&b.distance).unwrap_or(Equal));
    }
}

#[derive(Default)]
struct Connectivity{
    /*
    Connectivity of a molecule defined in terms of 'bonds'. Includes information about the pairs,
    triples and quadruples which define the bonds, angle and dihedral components required to
    calculate the total energy of a molecule with a force-field.
    */
    pub bonds:     HashSet<Bond>,
    pub angles:    HashSet<Angle>,
    pub dihedrals: HashSet<Dihedral>
}

impl Connectivity {

    /// Clear all the connectivity (bonds, angles and dihedrals)
    pub fn clear(&mut self){
        self.bonds.clear();
        self.angles.clear();
        self.dihedrals.clear();
    }

}

#[derive(Default)]
struct NBPair{
    i: usize,
    j: usize
}

#[derive(Default, Debug, Hash)]
pub struct Bond{
    i: usize,
    j: usize
}

impl Bond {

    /// Does this bond contain a particular atom?
    pub fn contains(&self, atom: &Atom) -> bool{
        atom.idx == self.i || atom.idx == self.j
    }

    /// Given the index of an atom that may, or may not be present in this bond, return the other
    /// atom index in the bond, if present.
    pub fn other(&self, idx: usize) -> Option<usize>{

        if      idx == self.i { Some(self.j) }
        else if idx == self.j { Some(self.i) }
        else{                       None     }
    }
}

impl PartialEq for Bond {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i && self.j == other.j || self.i == other.j && self.j == other.i
    }
}

impl Eq for Bond {}

#[derive(Default, Debug, Hash)]
pub struct Angle{
    i: usize,
    j: usize,
    k: usize
}

impl PartialEq for Angle {
    fn eq(&self, other: &Self) -> bool {
        self.j == other.j
            && (self.i == other.i && self.k == other.k || self.i == other.k && self.k == other.i)
    }
}

impl Eq for Angle {}

#[derive(Default, Hash, Debug)]
pub struct Dihedral{
    i: usize,
    j: usize,
    k: usize,
    l: usize
}

impl PartialEq for Dihedral{
    fn eq(&self, other: &Self) -> bool {
        (self.i == other.i && self.j == other.j && self.k == other.k && self.l == other.l)
            || (self.i == other.l && self.j == other.k && self.k == other.j && self.l == other.i)
    }
}

impl Eq for Dihedral {}


/*
   /$$                           /$$
  | $$                          | $$
 /$$$$$$    /$$$$$$   /$$$$$$$ /$$$$$$   /$$$$$$$
|_  $$_/   /$$__  $$ /$$_____/|_  $$_/  /$$_____/
  | $$    | $$$$$$$$|  $$$$$$   | $$   |  $$$$$$
  | $$ /$$| $$_____/ \____  $$  | $$ /$$\____  $$
  |  $$$$/|  $$$$$$$ /$$$$$$$/  |  $$$$//$$$$$$$/
   \___/   \_______/|_______/    \___/ |_______/
 */

#[cfg(test)]
mod tests{

    use super::*;
    use crate::utils::*;

    /// Given a valid simple xyz file, when a Molecule is created, then the connectivity is correct
    #[test]
    fn test_molecule_from_xyz_file(){

        print_methane_xyz_file("tmp_methane_1.xyz");
        let mol = Molecule::from_xyz_file("tmp_methane_1.xyz");
        assert_eq!(mol.num_atoms(), 5);

        // Carbon is tetrahedral
        assert_eq!(mol.connectivity.bonds.len(), 4);

        let expected_bond = Bond{i: 0, j: 1};
        assert!(mol.contains_bond(&expected_bond));

        let not_expected_bond = Bond{i: 0, j: 0};
        assert!(!mol.contains_bond(&not_expected_bond));

        assert_eq!(mol.connectivity.angles.len(), 6);
        assert!(mol.connectivity.angles.iter().any(|b| b==&Angle{i: 1, j: 0, k: 2}));

        remove_file_or_panic("tmp_methane_1.xyz");
    }

    /// Bonds are equivilant irrespective of the ordering
    #[test]
    fn test_bond_equality(){
        assert_eq!(&Bond{i: 0, j: 1}, &Bond{i: 0, j: 1});
        assert_eq!(&Bond{i: 0, j: 1}, &Bond{i: 1, j: 0});
        assert_ne!(&Bond{i: 0, j: 1}, &Bond{i: 0, j: 2});
    }

    /// Given a molecule with no atoms, when bonds are added, then no exception is thrown
    #[test]
    fn test_add_bonds_no_atoms(){

        let mut mol = Molecule::blank();
        assert_eq!(mol.num_atoms(), 0);

        mol.add_bonds();

        assert!(!mol.has_bonds());
    }

    /// Given a molecule with <3 atoms, when angles are added, then no exception is thrown
    #[test]
    fn test_add_angles_no_atoms(){

        let mut mol = Molecule::blank();
        mol.add_angles();

        assert_eq!(mol.connectivity.angles.len(), 0);
    }
    #[test]
    fn test_add_angles_two_atoms(){

        let mut mol = Molecule::from_atomic_symbols(&["C", "H"]);
        mol.add_angles();

        assert_eq!(mol.connectivity.angles.len(), 0);
    }

    /// Given a molecule with 3 atoms but no bonds, no angles are added
    #[test]
    fn test_add_angles_no_bonds(){

        print_water_xyz_file("tmp_water_0.xyz");

        let mut mol = Molecule::from_xyz_file("tmp_water_0.xyz");
        mol.connectivity.clear();
        assert!(!mol.has_bonds());

        mol.add_angles();
        println!("{:?}", mol.connectivity.angles);
        assert!(mol.connectivity.angles.is_empty());

        remove_file_or_panic("tmp_water_0.xyz");
    }

    /// The indices in a triple (angle) are ordered, so equality between two angles requires either
    /// an identical order of atom indices or completely reversed
    #[test]
    fn test_angle_equality(){

        assert_eq!(&Angle { i: 0, j: 1, k: 2 }, &Angle { i: 0, j: 1, k: 2 });
        assert_eq!(&Angle { i: 0, j: 1, k: 2 }, &Angle { i: 2, j: 1, k: 0 });
        assert_ne!(&Angle { i: 0, j: 1, k: 2 }, &Angle { i: 0, j: 3, k: 2 });
    }

    /// A bond of two atoms should contain both atoms and be able to return the other atom index
    #[test]
    fn test_bond_contains_and_other(){

        let atom_i = Atom{idx: 0,
                          atomic_number: AtomicNumber::from_string("C").unwrap(),
                          coordinate: CartesianCoordinate::default()};

        let atom_j = Atom{idx: 3,
                          atomic_number: AtomicNumber::from_string("H").unwrap(),
                          coordinate: CartesianCoordinate::default()};

        let bond = Bond{i: atom_i.idx, j: atom_j.idx};

        assert!(bond.contains(&atom_i));
        assert!(bond.contains(&atom_j));

        let mut different_atom = Atom::default();
        different_atom.idx = 1;
        assert!(!bond.contains(&different_atom));

        assert_eq!(bond.other(0).unwrap(), 3);
        assert_eq!(bond.other(3).unwrap(), 0);
        assert!(bond.other(1).is_none());
    }

    /// Dihedrals should not exist for molecules with <4 atoms
    #[test]
    fn test_no_dihedrals_with_fewer_than_four_atoms(){

        assert!(!Molecule::blank().has_dihedrals());
        assert!(!Molecule::from_atomic_symbols(&["H"]).has_dihedrals());

        print_dihydrogen_xyz_file("tmp_h2_0.xyz");
        assert!(!Molecule::from_xyz_file("tmp_h2_0.xyz").has_dihedrals());
        remove_file_or_panic("tmp_h2_0.xyz");

        print_water_xyz_file("tmp_water_1.xyz");
        assert!(!Molecule::from_xyz_file("tmp_water_1.xyz").has_dihedrals());
        remove_file_or_panic("tmp_water_1.xyz");
    }

    /// Dihedrals are equivilant only with forward and reverse ordering
    #[test]
    fn test_dihedral_equality(){

        assert_eq!(&Dihedral{i: 0, j: 1, k: 2, l: 3}, &Dihedral{i: 0, j: 1, k: 2, l: 3});
        assert_eq!(&Dihedral{i: 0, j: 1, k: 2, l: 3}, &Dihedral{i: 3, j: 2, k: 1, l: 0});
        assert_ne!(&Dihedral{i: 0, j: 1, k: 2, l: 3}, &Dihedral{i: 5, j: 1, k: 2, l: 3});
    }
}
