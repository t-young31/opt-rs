use crate::atoms::{Atom, AtomicNumber};
use crate::connectivity::angles::Angle;
use crate::connectivity::bonds::{Bond, BondOrder};
use crate::connectivity::dihedrals::{ImproperDihedral, ProperDihedral};
use crate::coordinates::{angle_value, Point, Vector3D};
use crate::ff::forcefield::Forcefield;
use crate::ff::rb::core::RB;
use crate::io::xyz::XYZFile;
use crate::opt::sd::SteepestDecentOptimiser;
use crate::pairs::{distance, NBPair};
use crate::utils::is_close;
use log::info;
use rand::Rng;
use std::cmp::Ordering::Equal;
use std::collections::HashSet;
use std::f64::consts::PI;

pub struct Molecule {
    pub charge: i32,
    pub coordinates: Vec<Point>,
    atomic_numbers: Vec<AtomicNumber>,
    pub(crate) connectivity: Connectivity,
    pub(crate) non_bonded_pairs: Vec<NBPair>,
}

impl Molecule {
    /// A blank molecule containing no atoms
    #[cfg(test)]
    pub fn blank() -> Self {
        Molecule::from_atomic_symbols(&[])
    }

    /// Create a Molecule from a .xyz filename
    /// # Arguments
    ///
    /// * `filename` - Name of the .xyz file defining the strucure
    ///
    /// # Examples
    /// use optrs::Molecule;
    /// let mol = Molecule::from_xyz_file("name");
    ///
    pub fn from_xyz_file(filename: &str) -> Self {
        let xyz_file = XYZFile::read(filename).unwrap();

        Molecule::from_atomic_nums_and_coords(xyz_file.atomic_numbers, xyz_file.coordinates)
    }

    /// Construct a molecule from atomic symbols, where all atomic positions are at the origin
    /// and an exception is thrown if the atomic number cannot be created
    pub fn from_atomic_symbols(symbols: &[&str]) -> Self {
        let mut atomic_numbers: Vec<AtomicNumber> = Default::default();
        let mut coords: Vec<Point> = Default::default();

        for symbol in symbols {
            atomic_numbers.push(AtomicNumber::from_string(symbol).unwrap());
            coords.push(Point::default());
        }

        Molecule::from_atomic_nums_and_coords(atomic_numbers, coords)
    }

    /// Construct a molecule from a set of atomic numbers and coordinates of each atom
    fn from_atomic_nums_and_coords(
        atomic_numbers: Vec<AtomicNumber>,
        coordinates: Vec<Point>,
    ) -> Self {
        let mut molecule = Molecule {
            coordinates,
            charge: 0,
            atomic_numbers,
            connectivity: Default::default(),
            non_bonded_pairs: Default::default(),
        };

        molecule.add_bonds();
        molecule.add_angles();
        molecule.add_dihedrals();
        molecule.add_non_bonded_pairs();
        molecule
    }

    /// Energy of this molecule using a FF
    pub fn energy(&self, forcefield: &mut dyn Forcefield) -> f64 {
        forcefield.energy(&self.coordinates)
    }

    /// Gradient of the energy with respect to the positions
    pub fn gradient(&mut self, forcefield: &mut dyn Forcefield) -> Vec<Vector3D> {
        forcefield.gradient(&self.coordinates).clone()
    }

    /// Numerical gradient evaluated using central differences
    pub fn numerical_gradient(&mut self, forcefield: &mut dyn Forcefield) -> Vec<Vector3D> {
        let h: f64 = 1E-6;
        let mut grad: Vec<Vector3D> = Default::default();

        for i in 0..self.num_atoms() {
            let mut vec = Vec::default();
            for (k, _) in ['x', 'y', 'z'].iter().enumerate() {
                let mut coords = self.coordinates.clone();

                coords[i][k] -= h; //  -> -h
                let e_minus = forcefield.energy(&coords);

                coords[i][k] += 2.0 * h; //  -> +h
                let e_plus = forcefield.energy(&coords);

                vec.push((e_plus - e_minus) / (2.0 * h));
            }
            grad.push(Vector3D::from_vector(&vec));
        }

        grad
    }

    /// Optimise the positions of the atoms given a forcefield
    pub fn optimise(&mut self, forcefield: &mut dyn Forcefield) {
        let mut optimiser = SteepestDecentOptimiser::default();
        optimiser.optimise(self, forcefield);
    }

    /// Write a .xyz file for this molecule
    pub fn write_xyz_file(&self, filename: &str) {
        XYZFile::write(filename, self);
    }

    /// Get a copy of set of atoms associated with this molecule
    pub(crate) fn atoms(&self) -> Vec<Atom> {
        let mut atoms: Vec<Atom> = Default::default();

        for (i, atomic_number) in self.atomic_numbers.iter().enumerate() {
            let coord = self.coordinates.get(i).expect("N_atoms != N_coords");

            let neighbours: Vec<usize> = self
                .connectivity
                .bonds
                .iter()
                .filter(|b| b.contains_index(i))
                .map(|b| b.other(i).unwrap())
                .collect();

            atoms.push(Atom {
                idx: i,
                atomic_number: atomic_number.clone(),
                coordinate: coord.clone(),
                bonded_neighbours: neighbours,
                formal_charge: 0.0,
            });
        }

        self.set_formal_charges(&mut atoms);

        atoms
    }

    /// Set the formal charges for a set of atoms given this connectivity
    fn set_formal_charges(&self, atoms: &mut Vec<Atom>) {
        // Assign + to least electronegative atom and - to most
        // distribute evenly between equiv atoms. Check for hypervalency
        // TODO: make this better

        let mut remaining_charge: i32 = self.charge;

        for atom in atoms.iter_mut() {
            let mut n = atom.num_valance_electrons() as f64;
            for bond in self.connectivity.bonds.iter().filter(|b| b.contains(atom)) {
                n -= bond.order.value();
            }

            while n > 2.0 {
                n -= 2.0; // Subtract any possible non-bonded pairs
            }

            if n.abs() > 1E-6 {
                // Has some unassigned electrons
                if remaining_charge > 0 {
                    n += 1.0;
                    remaining_charge -= 1;
                }
                if remaining_charge < 0 {
                    n -= 1.0;
                    remaining_charge += 1;
                }
            }

            atom.formal_charge = n;
        }
    }

    /// Number of atoms in this molecule
    pub(crate) fn num_atoms(&self) -> usize {
        self.atomic_numbers.len()
    }

    /// Bonds within this molecule
    pub(crate) fn bonds(&self) -> &HashSet<Bond> {
        &self.connectivity.bonds
    }

    /// Angles within this molecule
    pub(crate) fn angles(&self) -> &HashSet<Angle> {
        &self.connectivity.angles
    }

    /// Dihedrals within this molecule
    pub(crate) fn proper_dihedrals(&self) -> &HashSet<ProperDihedral> {
        &self.connectivity.proper_dihedrals
    }
    pub(crate) fn improper_dihedrals(&self) -> &HashSet<ImproperDihedral> {
        &self.connectivity.improper_dihedrals
    }

    /// Does this molecule have any associated bonds?
    fn has_bonds(&self) -> bool {
        !self.connectivity.bonds.is_empty()
    }

    /// Does this molecule contain a particular bond?
    fn contains_bond(&self, bond: &Bond) -> bool {
        self.connectivity.bonds.iter().any(|b| b == bond)
    }

    fn contains_bond_between(&self, i: usize, j: usize) -> bool {
        i != j && self.contains_bond(&Bond::from_atom_indices(i, j))
    }

    fn has_dihedrals(&self) -> bool {
        self.connectivity.improper_dihedrals.len() + self.connectivity.proper_dihedrals.len() > 0
    }

    /// Add bonds between atoms based on their interatomic distance. For example H-H is consdiered
    /// bonded if the r(HH) < 0.8 Å, or so. Note that this is not particularly well defined, but
    /// essential to determining the energy with a 'classic' (i.e. non-density based) force-field
    pub(crate) fn add_bonds(&mut self) {
        self.connectivity.bonds.clear();

        for atom_i in self.atoms().iter() {
            for n in Neighbours::from_atom_and_molecule(atom_i, self).iter() {
                self.add_bond(atom_i, &n.atom);
            }
        }

        self.set_bond_orders();
    }

    /// Set bond orders for all bonds in this molecule
    fn set_bond_orders(&mut self) {
        let atoms = self.atoms();

        // Iterate though a vector of bonds as the hashset doesn't support iter_mut
        let mut bonds: Vec<_> = self.connectivity.bonds.clone().into_iter().collect();

        for bond in bonds.iter_mut() {
            bond.set_possible_bond_order(&atoms);
        }

        // Refine to remove any hypervalency
        for atom in atoms.iter() {
            if atom.is_hypervalent(&bonds) && atom.period() == 2 {
                atom.reduce_two_double_bonds_to_aromatic(&mut bonds);
                atom.reduce_triple_bond_to_single(&mut bonds);
            }
        }

        self.connectivity.set_bonds_from_vector_of_bonds(bonds);
    }

    /// Add a bond between two atoms, provided it is not already present
    fn add_bond(&mut self, atom_i: &Atom, atom_j: &Atom) {
        let new_bond = Bond::from_atoms(atom_i, atom_j);
        let mut n_bonds = NBonds { i: 0, j: 0 };

        for bond in self.connectivity.bonds.iter() {
            if bond == &new_bond {
                return;
            }

            if bond.contains(atom_i) {
                n_bonds.i += 1
            }
            if bond.contains(atom_j) {
                n_bonds.j += 1
            }
        }

        if n_bonds.i >= atom_i.maximal_valence() || n_bonds.j >= atom_j.maximal_valence() {
            info!(
                "Not adding bond {}-{}. One exceeded the maximal valance",
                atom_i.idx, atom_j.idx
            );
            return;
        }

        self.connectivity.bonds.insert(new_bond);
    }

    /// Add angles between triples of mutually bonded atoms. For example,
    ///
    ///```text
    ///     j     k
    ///     /----
    ///    /
    ///   i
    /// ```
    pub(crate) fn add_angles(&mut self) {
        if self.num_atoms() < 3 || !self.has_bonds() {
            return; // No MM angles in molecules with < 3 atoms or that doesn't have any bonds
        }

        for j in 0..self.num_atoms() {
            for i in 0..self.num_atoms() {
                if !self.contains_bond_between(i, j) {
                    continue;
                }

                for k in 0..self.num_atoms() {
                    if i == k || !self.contains_bond_between(j, k) {
                        continue;
                    }

                    self.connectivity.angles.insert(Angle { i, j, k });
                }
            }
        }
    }

    /// Add all the dihedrals between quadruples pf bonded atoms
    pub(crate) fn add_dihedrals(&mut self) {
        if self.num_atoms() < 4 {
            return; // No dihedrals in molecules with < 4 atoms
        }

        let all_neighbours: Vec<Vec<usize>> = self
            .atoms()
            .iter()
            .map(|a| a.bonded_neighbours.clone())
            .collect();

        self.add_proper_dihedrals(&all_neighbours);
        self.add_improper_dihedrals(&all_neighbours);
    }

    /// Add proper dihedrals between sequentially bonded dihedrals: i -- j -- k -- l
    fn add_proper_dihedrals(&mut self, all_neighbours: &[Vec<usize>]) {
        for bond in self.connectivity.bonds.iter() {
            for neighbour_i in all_neighbours[bond.pair.i].iter() {
                if neighbour_i == &bond.pair.j {
                    continue;
                }

                for neighbour_j in all_neighbours[bond.pair.j].iter() {
                    if neighbour_j == &bond.pair.i {
                        continue;
                    }

                    let dihedral = ProperDihedral {
                        i: *neighbour_i,
                        j: bond.pair.i,
                        k: bond.pair.j,
                        l: *neighbour_j,
                    };

                    self.connectivity.proper_dihedrals.insert(dihedral);
                }
            }
        }
    }

    /// Add improper dihedrals between a central atom and three others e.g.
    ///           j
    ///          /
    ///     i --c-- k
    ///
    fn add_improper_dihedrals(&mut self, all_neighbours: &[Vec<usize>]) {
        for c in 0..self.num_atoms() {
            if all_neighbours[c].len() != 3 {
                continue; // Impropers are defined for a central atom bonded to exactly 3 others
            }

            let neighbours = &all_neighbours[c];
            let dihedral = ImproperDihedral {
                c,
                i: neighbours[0],
                j: neighbours[1],
                k: neighbours[2],
            };

            self.connectivity.improper_dihedrals.insert(dihedral);
        }
    }

    /// Add all pairs (i, j) of atoms which are non-bonded thus are subject to Lennard Jones
    /// interactions, in a classical MM forcefield
    pub(crate) fn add_non_bonded_pairs(&mut self) {
        self.non_bonded_pairs = Vec::with_capacity(self.num_atoms().pow(2));

        for atom_i in self.atoms().iter() {
            for atom_j in self.atoms().iter() {
                if self.contains_bond_between(atom_i.idx, atom_j.idx) {
                    continue; // Skip bonded pairs
                }

                if atom_i.idx <= atom_j.idx {
                    continue; // Only add unique pairs
                }

                self.non_bonded_pairs
                    .push(NBPair::from_atoms(atom_i, atom_j));
            }
        }
    }

    /// Is an angle close to linear?
    pub(crate) fn angle_is_close_to_linear(&self, i: usize, j: usize, k: usize) -> bool {
        is_close(angle_value(i, j, k, &self.coordinates), PI, 1E-1)
    }

    /// Evaluate the minimum atom-atom distance (Å) in this molecule
    fn minimum_pairwise_distance(&self) -> f64 {
        let mut min_distance = f64::MAX;

        for i in 0..self.num_atoms() {
            for j in i + 1..self.num_atoms() {
                let distance = distance(i, j, &self.coordinates);
                if distance < min_distance {
                    min_distance = distance;
                }
            }
        }

        min_distance
    }

    /// Randomise the coordinates of this molecule within a cubic box with side length chosen to
    /// generate an atomic density of ~0.1 atom/Å^3
    pub fn randomise_coordinates(&mut self) {
        let mut rng = rand::thread_rng();
        let box_side_length = ((self.num_atoms() as f64) / 0.1).powf(1. / 3.);
        let max_num_iterations: usize = 1000;

        for _ in 0..max_num_iterations {
            for coord in self.coordinates.iter_mut() {
                for (k, _) in ['x', 'y', 'z'].iter().enumerate() {
                    coord[k] = rng.gen_range(0.0..box_side_length);
                }
            }

            if self.minimum_pairwise_distance() > 0.5 {
                // Don't want any very small distances
                break;
            }
        }
    }

    /// Construct a zero gradient vector over all atoms
    pub fn blank_gradient(&self) -> Vec<Vector3D> {
        let mut gradient: Vec<Vector3D> = Default::default();
        for _ in 0..self.num_atoms() {
            gradient.push(Vector3D::default());
        }

        gradient
    }

    /// Build a possible 3d structure of a molecule structure using iterative bond
    /// addition and minimising with a repulsion+bonded forcefield
    pub fn build_3d(&mut self) {
        self.randomise_coordinates();

        let all_bonds = self.connectivity.bonds.clone();
        self.connectivity.bonds.clear();

        for bond in all_bonds {
            let mut optimiser = SteepestDecentOptimiser::from_max_iterations(20);
            optimiser.optimise(self, &mut RB::new(self));

            self.connectivity.bonds.insert(bond);
        }

        let mut optimiser = SteepestDecentOptimiser::default();
        optimiser.optimise(self, &mut RB::new(self));
    }
}

#[derive(Default, Debug)]
struct Neighbour {
    atom: Atom,
    distance: f64,
}

#[derive(Default, Debug)]
struct Neighbours {
    values: Vec<Neighbour>,
}

impl Neighbours {
    /// Given an atom contained in a molecule, generate the bonded neighbours of that atom
    /// sorted by the distance to the central atom
    pub fn from_atom_and_molecule(central_atom: &Atom, molecule: &Molecule) -> Self {
        let mut neighbours: Neighbours = Default::default();

        neighbours.add_all_while_guessing_bonds(central_atom, molecule);
        neighbours.sort_by_distance();
        neighbours
    }

    /// Add neighbours to a central atom while guessing the bonds
    fn add_all_while_guessing_bonds(&mut self, central_atom: &Atom, molecule: &Molecule) {
        for atom in molecule.atoms() {
            if central_atom.could_be_bonded_to(&atom) {
                let distance = central_atom.distance_to(&atom);
                self.values.push(Neighbour { atom, distance });
            }
        }
    }

    /// Sort these neighbours by the distance to the central atom
    fn sort_by_distance(&mut self) {
        self.values
            .sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Equal));
    }

    fn iter(&self) -> impl Iterator<Item = &Neighbour> {
        self.values.iter()
    }
}

#[derive(Default)]
pub(crate) struct Connectivity {
    /*
    Connectivity of a molecule defined in terms of 'bonds'. Includes information about the pairs,
    triples and quadruples which define the bonds, angle and dihedral components required to
    calculate the total energy of a molecule with a force-field.
    */
    pub bonds: HashSet<Bond>,
    pub angles: HashSet<Angle>,
    pub proper_dihedrals: HashSet<ProperDihedral>,
    pub improper_dihedrals: HashSet<ImproperDihedral>,
}

impl Connectivity {
    /// Clear all the connectivity (bonds, angles and dihedrals)
    pub fn clear(&mut self) {
        self.bonds.clear();
        self.angles.clear();
        self.proper_dihedrals.clear();
        self.improper_dihedrals.clear();
    }

    /// Given a vector of bonds set the connectivity
    fn set_bonds_from_vector_of_bonds(&mut self, bonds: Vec<Bond>) {
        self.bonds.clear();

        for bond in bonds {
            self.bonds.insert(bond);
        }
    }
}

/// Number of bonds present for two atoms
struct NBonds {
    i: usize,
    j: usize,
}

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
mod tests {
    use super::*;
    use crate::connectivity::bonds::BondOrder;
    use crate::ff::uff::core::UFF;
    use crate::pairs::{distance, AtomPair};
    use crate::utils::*;

    fn proper_dihedral(i: usize, j: usize, k: usize, l: usize) -> ProperDihedral {
        ProperDihedral { i, j, k, l }
    }

    /// Given a valid simple xyz file, when a Molecule is created, then the connectivity is correct
    #[test]
    fn test_molecule_from_xyz_file() {
        print_methane_xyz_file("tmp_methane_1.xyz");
        let mol = Molecule::from_xyz_file("tmp_methane_1.xyz");
        assert_eq!(mol.num_atoms(), 5);

        // Carbon is tetrahedral
        assert_eq!(mol.connectivity.bonds.len(), 4);

        let expected_bond = Bond::from_atom_indices(0, 1);
        assert!(mol.contains_bond(&expected_bond));

        let not_expected_bond = Bond::from_atom_indices(0, 99);
        assert!(!mol.contains_bond(&not_expected_bond));

        assert_eq!(mol.connectivity.angles.len(), 6);
        assert!(mol
            .connectivity
            .angles
            .iter()
            .any(|b| b == &Angle { i: 1, j: 0, k: 2 }));

        assert!(is_close(distance(0, 1, &mol.coordinates), 1.1, 1E-1));

        remove_file_or_panic("tmp_methane_1.xyz");
    }

    /// Bonds are equivilant irrespective of the ordering
    #[test]
    fn test_bond_equality() {
        assert_eq!(
            &Bond::from_atom_indices(0, 1),
            &Bond::from_atom_indices(0, 1)
        );
        assert_eq!(
            &Bond::from_atom_indices(0, 1),
            &Bond::from_atom_indices(1, 0)
        );
        assert_ne!(
            &Bond::from_atom_indices(0, 1),
            &Bond::from_atom_indices(0, 2)
        );
    }

    /// Given a molecule with no atoms, when bonds are added, then no exception is thrown
    #[test]
    fn test_add_bonds_no_atoms() {
        let mut mol = Molecule::blank();
        assert_eq!(mol.num_atoms(), 0);

        mol.add_bonds();

        assert!(!mol.has_bonds());
    }

    /// Test that a bond cannot be constructed between two identical atoms
    /// if one exists then there's the chance to get infinities in the FF energy
    #[test]
    #[should_panic]
    fn test_cannot_create_a_bond_between_identical_atoms() {
        let _ = Bond::from_atom_indices(0, 0);
    }

    /// Given hydrogen trimer with close distance, then only a single bond exists
    #[test]
    fn test_hypervalent_hydrogen() {
        let mut h2 = Molecule::from_atomic_symbols(&["H", "H"]);
        h2.coordinates[1] = Point {
            x: 0.77,
            y: 0.,
            z: 0.,
        };
        h2.add_bonds();
        assert_eq!(h2.connectivity.bonds.len(), 1);

        let mut h3 = Molecule::from_atomic_symbols(&["H", "H", "H"]);
        h3.coordinates[1] = Point {
            x: -0.77,
            y: 0.,
            z: 0.,
        };
        h3.coordinates[2] = Point {
            x: 0.77,
            y: 0.,
            z: 0.,
        };
        h3.add_bonds();

        assert_eq!(h3.connectivity.bonds.len(), 1);
    }

    /// Given a molecule with <3 atoms, when angles are added, then no exception is thrown
    #[test]
    fn test_add_angles_no_atoms() {
        let mut mol = Molecule::blank();
        mol.add_angles();

        assert_eq!(mol.connectivity.angles.len(), 0);
    }
    #[test]
    fn test_add_angles_two_atoms() {
        let mut mol = Molecule::from_atomic_symbols(&["C", "H"]);
        mol.add_angles();

        assert_eq!(mol.connectivity.angles.len(), 0);
    }

    /// Given a molecule with 3 atoms but no bonds, no angles are added
    #[test]
    fn test_add_angles_no_bonds() {
        print_water_xyz_file("tmp_water_0.xyz");

        let mut mol = Molecule::from_xyz_file("tmp_water_0.xyz");
        mol.connectivity.clear();
        assert!(!mol.has_bonds());

        mol.add_angles();
        assert!(mol.connectivity.angles.is_empty());

        remove_file_or_panic("tmp_water_0.xyz");
    }

    /// The indices in a triple (angle) are ordered, so equality between two angles requires either
    /// an identical order of atom indices or completely reversed
    #[test]
    fn test_angle_equality() {
        let angle = Angle { i: 0, j: 1, k: 2 };
        assert_eq!(&angle, &Angle { i: 0, j: 1, k: 2 });
        assert_eq!(angle.j, 1 as usize);
        assert_eq!(angle.k, 2 as usize);

        assert_eq!(&angle, &Angle { i: 2, j: 1, k: 0 });
        assert_ne!(&angle, &Angle { i: 0, j: 3, k: 2 });
    }

    /// A bond of two atoms should contain both atoms and be able to return the other atom index
    #[test]
    fn test_bond_contains_and_other() {
        let atom_i = Atom::from_idx_and_atomic_symbol(0, "C");
        let atom_j = Atom::from_idx_and_atomic_symbol(3, "H");

        let bond = Bond::from_atoms(&atom_i, &atom_j);

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
    fn test_no_dihedrals_with_fewer_than_four_atoms() {
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
    fn test_dihedral_equality() {
        assert_eq!(&proper_dihedral(0, 1, 2, 3), &proper_dihedral(0, 1, 2, 3));
        assert_eq!(&proper_dihedral(0, 1, 2, 3), &proper_dihedral(3, 2, 1, 0));
        assert_ne!(&proper_dihedral(0, 1, 2, 3), &proper_dihedral(5, 1, 2, 3));
    }

    /// Given a collection of atoms that are not bonded, then no connectivity dihedrals exist
    #[test]
    fn test_no_dihedral_no_bonds() {
        let mut mol = Molecule::from_atomic_symbols(&["H", "H", "H", "H"]);
        mol.coordinates[1] = Point {
            x: 999.,
            y: 0.,
            z: 0.,
        };
        mol.coordinates[2] = Point {
            x: 0.,
            y: 999.,
            z: 0.,
        };
        mol.coordinates[3] = Point {
            x: 0.,
            y: 0.,
            z: 999.,
        };

        assert!(!mol.has_bonds());
        assert_eq!(mol.connectivity.proper_dihedrals.len(), 0);
    }

    /// Given a non-bonded H2 pair or atoms then there should be a single NB pair present
    #[test]
    fn test_nb_pair_for_non_bonded_h2() {
        let mut h2 = Molecule::from_atomic_symbols(&["H", "H"]);
        h2.coordinates[1].x = 9.99;
        h2.add_non_bonded_pairs();

        assert_eq!(h2.non_bonded_pairs.len(), 1);
    }

    /// Given a bonded H2 molecule then there should be no NB pairs present
    #[test]
    fn test_nb_pair_for_bonded_h2() {
        let mut h2 = Molecule::from_atomic_symbols(&["H", "H"]);
        h2.coordinates[1].x = 0.77;

        h2.add_bonds();
        h2.add_non_bonded_pairs();

        assert_eq!(h2.non_bonded_pairs.len(), 0);
    }

    /// Given two non-bonded pairs of atoms with different ordering then they should be identical
    #[test]
    fn test_non_bonded_pair_equality() {
        let pair_a = NBPair {
            pair: AtomPair { i: 0, j: 1 },
        };
        let pair_b = NBPair {
            pair: AtomPair { i: 1, j: 0 },
        };

        assert_eq!(pair_a, pair_b);
    }

    /// Test that some bond orders are correct
    #[test]
    fn test_bond_order_water() {
        let filename = "test_tbow.xyz";
        print_water_xyz_file(filename);
        let mol = Molecule::from_xyz_file(filename);

        for bond in mol.bonds().iter() {
            assert_eq!(bond.order, BondOrder::Single);
        }

        remove_file_or_panic(filename);
    }
    #[test]
    fn test_bond_order_ethene() {
        let filename = "test_tboe.xyz";
        print_ethene_xyz_file(filename);
        let mol = Molecule::from_xyz_file(filename);

        for bond in mol.bonds().iter() {
            if bond.contains_index(0) && bond.contains_index(1) {
                assert_eq!(bond.order, BondOrder::Double);
            } else {
                assert_eq!(bond.order, BondOrder::Single);
            }
        }

        remove_file_or_panic(filename);
    }
    #[test]
    fn test_bond_order_h2cnh2() {
        let filename = "test_tbonh2.xyz";
        print_ethene_xyz_file(filename);
        let mut mol = Molecule::from_xyz_file(filename);
        mol.atomic_numbers[1] = AtomicNumber::from_string("N").unwrap();
        mol.connectivity.clear();
        mol.add_bonds();

        let expected_bond = Bond {
            pair: AtomPair { i: 0, j: 1 },
            order: BondOrder::Double,
        };
        assert!(mol.bonds().contains(&expected_bond)); // C=N

        remove_file_or_panic(filename);
    }
    #[test]
    fn test_bond_order_protonated_formaldehyde() {
        let filename = "test_tbopf.xyz";
        print_h2coh_xyz_file(filename);
        let mol = Molecule::from_xyz_file(filename);

        let expected_bond = Bond {
            pair: AtomPair { i: 0, j: 1 },
            order: BondOrder::Double,
        };
        assert!(mol.bonds().contains(&expected_bond));

        let o_atom = &mol.atoms()[1];
        assert_eq!(o_atom.atomic_symbol(), "O");
        assert!(is_very_close(o_atom.formal_charge, 1.0));

        remove_file_or_panic(filename);
    }

    /// Test that a square planar Pd complex has a d8 atom
    #[test]
    fn test_square_planar_complex_has_d8_atom() {
        let filename = "pd_complex_tspchda.xyz";
        print_pdcl2nh3h_xyz_file(filename);
        let mol = Molecule::from_xyz_file(filename);

        assert_eq!(mol.atomic_numbers[0].to_atomic_symbol(), "Pd");
        assert!(mol.atoms()[0].is_d8());

        remove_file_or_panic(filename);
    }

    /// Inserting two identical angles into a molecule's connectivity should not be possible
    #[test]
    fn test_no_identical_angles_present() {
        let mut mol = Molecule::blank();
        mol.connectivity.angles.insert(Angle { i: 0, j: 1, k: 2 });
        mol.connectivity.angles.insert(Angle { i: 2, j: 1, k: 0 });

        assert_eq!(mol.angles().len(), 1);
    }

    /// and similarly for two identical dihedrals
    #[test]
    fn test_no_identical_dihedrals_present() {
        let mut mol = Molecule::blank();
        mol.connectivity
            .proper_dihedrals
            .insert(proper_dihedral(0, 1, 2, 3));
        mol.connectivity
            .proper_dihedrals
            .insert(proper_dihedral(3, 2, 1, 0));

        assert_eq!(mol.proper_dihedrals().len(), 1);
    }

    /// Ensure that a trigonal atom has an improper dihedral forcing the current angle
    #[test]
    fn test_ethene_has_two_improper_dihedrals() {
        let filename = "ethene_m_htid.xyz";
        print_ethene_xyz_file(filename);

        let mol = Molecule::from_xyz_file(filename);
        remove_file_or_panic(filename);

        assert_eq!(mol.improper_dihedrals().len(), 2);
    }

    /// Check that benzene is correctly atom typed with aromatic bonds
    #[test]
    fn test_benzene_has_aromatic_bonds() {
        let filename = "benzene_tbhab.xyz";
        print_benzene_xyz_file(filename);

        let mol = Molecule::from_xyz_file(filename);
        let atoms = mol.atoms();

        fn is_carbon_carbon(b: &Bond, atoms: &Vec<Atom>) -> bool {
            atoms[b.pair.i].atomic_symbol() == "C" && atoms[b.pair.j].atomic_symbol() == "C"
        }

        for bond in mol.bonds().iter() {
            if is_carbon_carbon(bond, &atoms) {
                assert_eq!(bond.order, BondOrder::Aromatic);
            }
        }

        remove_file_or_panic(filename);
    }

    /// Test that triple bonds are defined
    #[test]
    fn test_a_triple_bond_is_present_in_c2h2() {
        let filename = "ch2h_tatbinc.xyz";
        print_c2h2_xyz_file(filename);

        let c2h2 = Molecule::from_xyz_file(filename);
        let cc_bond = c2h2
            .bonds()
            .iter()
            .filter(|b| b.contains_index(0) && b.contains_index(1))
            .next()
            .unwrap();

        assert_eq!(cc_bond.order, BondOrder::Triple);

        remove_file_or_panic(filename)
    }

    #[test]
    fn test_reducing_bond_orders_does_nothing_in_ethene() {
        let filename = "ethene_trbodnie.xyz";
        print_ethene_xyz_file(filename);

        let mol = Molecule::from_xyz_file(filename);
        let mut bonds: Vec<Bond> = mol.bonds().clone().into_iter().collect();

        for atom in mol.atoms().iter() {
            atom.reduce_two_double_bonds_to_aromatic(&mut bonds);
        }

        for bond in bonds.iter() {
            assert!(bond.order == BondOrder::Double || bond.order == BondOrder::Single);
        }

        remove_file_or_panic(filename);
    }

    /// Ensure that molecules can be optimised
    #[test]
    fn test_h2_optimisation() {
        let filename = "th2o.xyz";
        print_dihydrogen_xyz_file(filename);

        let mut h2 = Molecule::from_xyz_file(filename);
        remove_file_or_panic(filename);

        let mut ff = UFF::new(&h2);
        h2.optimise(&mut ff);

        let r = (&h2.coordinates[0] - &h2.coordinates[1]).length();
        assert!(is_close(r, 0.7, 1E-1));

        assert!(is_close(h2.energy(&mut ff), 0.0, 1E-4));

        // Also ensure the gradient is close to zero
        for v in h2.gradient(&mut ff).iter() {
            assert!(v.length() < 0.1);
        }
    }
    #[test]
    fn test_methane_optimisation() {
        let filename = "tmo.xyz";
        print_methane_xyz_file(filename);

        let mut ch4 = Molecule::from_xyz_file(filename);
        remove_file_or_panic(filename);

        let mut ff = UFF::new(&ch4);
        let init_energy = ch4.energy(&mut ff);

        ch4.optimise(&mut ff);

        assert!(init_energy > ch4.energy(&mut ff));
    }

    /// Ensure the atom-atom distance function is correct for a simple water molecule
    #[test]
    fn test_min_atom_atom_distance() {
        let mut mol = Molecule::from_atomic_symbols(&["H", "O", "H"]);

        mol.coordinates[0].x = -1.0;
        mol.coordinates[2].x = 1.1;

        assert!(is_very_close(mol.minimum_pairwise_distance(), 1.0));
    }

    /// Test that randomising coordinates in a molecule both generates coordinates with no
    /// short pairwise distance
    #[test]
    fn test_randomise_coordinates() {
        let mut mol = Molecule::from_atomic_symbols(&["H", "O", "H"]);

        let mut r = 0.;

        for _ in 0..100 {
            mol.randomise_coordinates();
            let new_r = mol.minimum_pairwise_distance();
            assert!(new_r > 0.5);
            assert!(!is_very_close(new_r, r));

            r = new_r;
        }
    }

    /// Test that bond orders can be initialised from their numerical value
    #[test]
    fn test_setting_bond_orders_from_values() {
        assert_eq!(BondOrder::from_value(&1.0), BondOrder::Single);
        assert_eq!(BondOrder::from_value(&4.0), BondOrder::Quadruple);
    }

    /// Cannot set a bond order from a value that is not close to something known
    #[should_panic]
    #[test]
    fn test_setting_bond_orders_from_invalid_values() {
        assert_eq!(BondOrder::from_value(&-1.0), BondOrder::Single);
    }

    /// Test that building a molecule results in reasonable molecules
    #[test]
    fn test_methane_structure_can_be_built() {
        let mut methane = Molecule::from_atomic_symbols(&["C", "H", "H", "H", "H"]);

        for i in 1..5 {
            methane
                .connectivity
                .bonds
                .insert(Bond::from_atom_indices(0, i));
        }

        methane.build_3d();

        let mut ff = UFF::new(&methane);
        assert!(methane.energy(&mut ff) < 5.);
    }
}
