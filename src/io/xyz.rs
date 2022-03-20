// xyz file parsing. See https://en.wikipedia.org/wiki/XYZ_file_format for a description
use std::fs::File;
use std::path::Path;
use std::io::{self, BufRead, Write};
use crate::atoms::{CartesianCoordinate, AtomicNumber};


#[derive(Default)]
pub struct XYZFile{

    pub filename:       String,
    file_lines:     Vec<String>,

    pub coordinates:    Vec<CartesianCoordinate>,
    pub atomic_numbers: Vec<AtomicNumber>
}


impl XYZFile {

    /// Create an xyz file and extract coordinates from an xyz file
    pub fn new(filename: &str) -> Self{

        if !filename.ends_with(".xyz") {
            panic!("Cannot read {}. File must end with .xyz", filename)
        }

        let mut file: XYZFile = Default::default();

        if let Ok(lines) = read_lines(filename) {
            for line in lines {
                if let Ok(line) = line {
                    file.add_atom_on_line(line);
                }
            }
        }

        file
    }

    /// Add the atomic number and coordinate for an atom on a line
    fn add_atom_on_line(&mut self, line: String){
        let mut items = line.split_whitespace();

        let atomic_symbol = items.next().expect("Failed to parse atomic symbol");
        self.atomic_numbers.push(AtomicNumber::from_string(atomic_symbol));

        self.coordinates.push(
            CartesianCoordinate{
                x: option_str_to_float(items.next()),
                y: option_str_to_float(items.next()),
                z: option_str_to_float(items.next())
            }
        );
    }
}


/// Convert an optional containing a string to a float
fn option_str_to_float(value: Option<&str>) -> f64{
    value.expect("Failed to find item").parse::<f64>().expect("Failed to convert to float")
}


/// Read a set of file lines into an iterator
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
