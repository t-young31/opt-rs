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
use crate::ff::rb::core::RB;
use crate::ff::uff::core::UFF;


#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct CommandLineArguments {
    #[clap(index = 1)]
    xyz_filename: String,

    #[clap(short, long, default_value="UFF")]
    forcefield: String,
}


fn run(args: CommandLineArguments) {

    let mut mol = Molecule::from_xyz_file(&args.xyz_filename);

    match args.forcefield.as_str() {
        "UFF" => mol.optimise(&mut UFF::new(&mol)),
        "RB" => mol.optimise(&mut RB::new(&mol)),
        _ => panic!("Cannot set a forcefield. Unknown type")
    }

    mol.write_xyz_file("opt.xyz");
}


fn main(){ run(CommandLineArguments::parse()) }



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

    #[test]
    fn test_simple_cli(){

        let filename = "h2_tscli.xyz";
        print_dihydrogen_xyz_file(filename);

        let args = CommandLineArguments{xyz_filename: filename.to_string(),
                                        forcefield: "UFF".to_string()};
        run(args);

        remove_file_or_panic(filename);
        remove_file_or_panic("opt.xyz")
    }

}
