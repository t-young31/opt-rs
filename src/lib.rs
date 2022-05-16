use pyo3::prelude::*;
use crate::ff::forcefield::Forcefield;
use crate::ff::uff::core::UFF;
use crate::molecule::Molecule;

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

    /// Given a vector of bond orders for all pairwise interactions
    fn set_bond_orders(&self, bond_orders: Vec<f64>){
        // todo!()
    }

    pub fn optimise(&mut self) {
        self.molecule.optimise(&mut UFF::new(&self.molecule));
        self.molecule.write_xyz_file("opt.xyz");
    }
}
