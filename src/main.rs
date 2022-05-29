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
mod cli;

use clap::Parser;

use crate::cli::{run, CommandLineArguments};

fn main() { run(CommandLineArguments::parse()) }
