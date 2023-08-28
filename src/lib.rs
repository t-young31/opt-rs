extern crate core;

use crate::connectivity::bonds::{Bond, BondOrder};
use crate::ff::forcefield::Forcefield;
use crate::ff::uff::core::UFF;
use crate::molecule::Molecule;
use crate::utils::IsVeryClose;
use pyo3::prelude::*;

mod atoms;
mod cli;
mod connectivity;
mod coordinates;
mod ff;
mod io;
mod molecule;
mod opt;
mod pairs;
mod utils;

use crate::cli::{run, CommandLineArguments};
use clap::Parser;

fn main() {
    run(CommandLineArguments::parse())
}

/// A Python module implemented in Rust.
#[pymodule]
fn optrs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyMoleculeWrapper>()?;
    Ok(())
}

/// Wrapper class over the rust-implemented Molecule struct
#[pyclass(name = "Molecule")]
struct PyMoleculeWrapper {
    molecule: Molecule,
}

/// Methods for wrapper
#[pymethods]
impl PyMoleculeWrapper {
    /// Generate a wrapped molecule from a xyz file that exists
    #[staticmethod]
    fn from_xyz_file(filename: &str) -> Self {
        PyMoleculeWrapper {
            molecule: Molecule::from_xyz_file(filename),
        }
    }

    /// Generate a wrapped molecule from a set of atomic symbols and set some random coordinates
    #[staticmethod]
    fn from_atomic_symbols(symbols: Vec<&str>) -> Self {
        let mol = Molecule::from_atomic_symbols(&symbols);
        PyMoleculeWrapper { molecule: mol }
    }

    /// Given a vector of bond orders for all pairwise interactions. This function is expecting a
    /// flat array in row major layout defining the bond order between atoms i,j
    fn set_bond_orders(&mut self, bond_orders: Vec<f64>) {
        if bond_orders.len() != self.molecule.num_atoms().pow(2) {
            panic!("Cannot set the bond orders. Must have a flat array with (N_atoms)^2 items in");
        }

        self.molecule.connectivity.clear();
        let mut i = 0;

        for (k, bond_order) in bond_orders.iter().enumerate() {
            let j = k % self.molecule.num_atoms();

            if j <= i || bond_order.is_very_close(&0.) {
                continue; // Only want unique pairs to be added
            }

            let mut bond = Bond::from_atom_indices(i, j);
            bond.order = BondOrder::from_value(bond_order);

            self.molecule.connectivity.bonds.insert(bond);

            if j == 0 {
                i += 1;
            }
        }

        self.molecule.add_angles();
        self.molecule.add_dihedrals();
        self.molecule.add_non_bonded_pairs();
    }

    /// Build the 3d structure of a molecule using iterative addition of bonds into the system
    /// with minimisation using a repulsion+bonded (RB) forcefield after each addition
    fn build_3d(&mut self) {
        if self.molecule.num_atoms() > 1 && self.molecule.bonds().is_empty() {
            panic!(
                "Cannot build a 3d structure without any bonds. \
                   Consider calling set_bond_orders()"
            );
        }

        self.molecule.build_3d()
    }

    /// Optimise the structure using the best FF possible
    pub fn optimise(&mut self) {
        self.molecule.optimise(&mut UFF::new(&self.molecule));
    }

    /// Write a .xyz file of this structure with a defined filename
    pub fn write_xyz_file(&self, filename: &str) {
        self.molecule.write_xyz_file(filename);
    }

    /// Set the coordinates using a flat list of coordinates
    pub fn set_coordinates(&mut self, coordinates: Vec<f64>) {
        if coordinates.len() != 3 * self.molecule.num_atoms() {
            panic!("Cannot set the coordinates. Must have a flat list with length 3xN_atoms");
        }

        for (i, coord) in self.molecule.coordinates.iter_mut().enumerate() {
            for (k, _) in ['x', 'y', 'z'].iter().enumerate() {
                coord[k] = coordinates[3 * i + k];
            }
        }
    }

    pub fn generate_connectivty(&mut self) {
        self.molecule.connectivity.clear();

        self.molecule.add_bonds();
        self.molecule.add_angles();
        self.molecule.add_dihedrals();
        self.molecule.add_non_bonded_pairs();
    }
}
