extern crate core;

use pyo3::prelude::*;
use crate::connectivity::bonds::{Bond, BondOrder};
use crate::ff::forcefield::Forcefield;
use crate::ff::uff::core::UFF;
use crate::molecule::Molecule;
use crate::utils::IsVeryClose;

mod molecule;
mod atoms;
mod connectivity;
mod ff;
mod io;
mod opt;
mod utils;
mod pairs;
mod coordinates;


/// A Python module implemented in Rust.
#[pymodule]
fn mors(_py: Python, m: &PyModule) -> PyResult<()> {

    m.add_class::<PyMoleculeWrapper>()?;
    Ok(())
}

/// Wrapper class over the rust-implemented Molecule struct
#[pyclass (name="Molecule")]
struct PyMoleculeWrapper {
    molecule: Molecule,
}

/// Methods for wrapper
#[pymethods]
impl PyMoleculeWrapper {

    /// Generate a wrapped molecule from a xyz file that exists
    #[staticmethod]
    fn from_xyz_file(filename: &str) -> Self {
        PyMoleculeWrapper { molecule: Molecule::from_xyz_file(filename) }
    }

    /// Generate a wrapped molecule from a set of atomic symbols and set some random coordinates
    #[staticmethod]
    fn from_atomic_symbols(symbols: Vec<&str>) -> Self{

        let mut mol = Molecule::from_atomic_symbols(&symbols);
        mol.randomise_coordinates();

        PyMoleculeWrapper { molecule:  mol}
    }

    /// Given a vector of bond orders for all pairwise interactions. This function is expecting a
    /// flat array in row major layout defining the bond order between atoms i,j
    fn set_bond_orders(&mut self, bond_orders: Vec<f64>){

        if bond_orders.len() != self.molecule.num_atoms().pow(2) {
            panic!("Cannot set the bond orders. Must have a flat array with (N_atoms)^2 items in");
        }

        let mut i = 0;

        for (k, bond_order) in bond_orders.iter().enumerate(){
            let j = k % self.molecule.num_atoms();

            if j <= i || bond_order.is_very_close(&0.){
                continue; // Only want unique pairs to be added
            }

            let mut bond = Bond::from_atom_indices(i, j);
            bond.order = BondOrder::from_value(bond_order);

            self.molecule.connectivity.bonds.insert(bond);

            if j == 0{
                i += 1;
            }
        }

    }

    pub fn optimise(&mut self) {
        self.molecule.optimise(&mut UFF::new(&self.molecule));
        self.molecule.write_xyz_file("opt.xyz");
    }
}
