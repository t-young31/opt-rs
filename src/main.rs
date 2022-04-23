/*
Notes and warnings:

    - Units are:
        -> Å (distance)
        -> rad (angle)
        -> kcal mol-1 (energy)
 */
extern crate core;

mod molecule;
mod atoms;
mod connectivity;
mod ff;
mod io;
mod opt;
mod utils;
mod pairs;
mod coordinates;

use clap::Parser;

use crate::molecule::Molecule;
use crate::ff::forcefield::Forcefield;
use crate::ff::uff::core::UFF;


#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct CommandLineArguments {
    #[clap(short, long)]
    xyz_filename: String,

    #[clap(short, long)]
    forcefield_name: String,
}


fn main() {

    let args = CommandLineArguments::parse();

    let mut mol = Molecule::from_xyz_file(&args.xyz_filename);

    match args.forcefield_name.as_str() {
        "UFF" => mol.optimise(&UFF::new(&mol)),
        _ => panic!("Cannot set a forcefield. Unknown type")
    }

    mol.write_xyz_file("opt.xyz");


    // 0. Setup
    // 	1. Read in .xyz
    // 	2. Atom type
    // 	3. Set lists of bonded pairs, triples, dihedral-s
    // 	4. Add pairs for non-bonded (with exclusions)
    // 	5. Calculate force constants
    // 1. Optimise
    //  1. Calculate gradient
    //  2. SD step (updating positions)
}
