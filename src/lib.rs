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

/// Optimise an xyz file
#[pyfunction]
fn optimise(filename: &str){
    let mut mol = Molecule::from_xyz_file(filename);

    mol.optimise(&mut UFF::new(&mol));
    mol.write_xyz_file("opt.xyz")
}

/// A Python module implemented in Rust.
#[pymodule]
fn mors(_py: Python, m: &PyModule) -> PyResult<()> {
    // m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(optimise, m)?)?;
    Ok(())
}
