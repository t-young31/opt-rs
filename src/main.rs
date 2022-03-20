mod molecule;
mod atoms;
mod atom_typing;
mod ff;
mod io;
mod opt;

use crate::molecule::Molecule;
use crate::ff::forcefield::Forcefield;
use clap::Parser;


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
    mol.set_forcefield(Forcefield::from_string(&args.forcefield_name));
    mol.optimise();
    mol.write_xyz_file();

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
