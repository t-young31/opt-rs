/*
Notes and warnings:

    - Units are:
        -> Å (distance)
        -> rad (angle)
        -> kcal mol-1 (energy)
 */
extern crate core;

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

use clap::Parser;

use crate::cli::{run, CommandLineArguments};

fn main() {
    run(CommandLineArguments::parse())
}
