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
    pub fn new(filename: &str) -> Result<Self, &'static str>{

        if !filename.ends_with(".xyz") {
            panic!("Cannot read {}. File must end with .xyz", filename)
        }

        let mut xyz_file: XYZFile = Default::default();

        for (i, line) in read_lines(filename).enumerate() {
            if let Ok(line) = line {

                // First two lines are number of atoms and the title
                if i <= 1 || line.is_empty() { continue; }

                xyz_file.add_atom_on_line(line);
            }
        }

        if xyz_file.coordinates.len() == 0{
            return Err("Failed to find any coordinates")
        }

        Ok(xyz_file)
    }

    /// Add the atomic number and coordinate for an atom on a line
    fn add_atom_on_line(&mut self, line: String) -> Result<(), &str>{
        let mut items = line.split_whitespace();

        let atomic_symbol = items.next().expect("Failed to parse atomic symbol");
        self.atomic_numbers.push(AtomicNumber::from_string(atomic_symbol));

        let coord = CartesianCoordinate::from_option_strings(items.next(), items.next(), items.next());

        match coord {
            Ok(c) => self.coordinates.push(c),
            Err(..) => return Err("Failed to create a coordinate"),
        }

        Ok(())
    }
}


/// Read a set of file lines into an iterator
fn read_lines<P>(filename: P) -> io::Lines<io::BufReader<File>>
    where P: AsRef<Path>, {
    let file = File::open(filename).expect("Failed to open file");
    io::BufReader::new(file).lines()
}


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
    fn test_ok_file_read(){
        std::fs::write(format!("tmp1.xyz"),
                       "5\n\n\
                        C     0.00000   0.00000   0.00000\n\
                        H    -0.65860  -0.85220  -0.30120\n\
                        H    -0.45940   0.97110  -0.28590\n\
                        H     0.08440  -0.02940   1.10060\n\
                        H     1.02910  -0.10990  -0.41250\n")
            .expect("Failed to write tmp1.xyz!");

        let file = XYZFile::new("tmp1.xyz").unwrap();

        assert_eq!(file.coordinates.len(), 5);
        assert_eq!(file.atomic_numbers.len(), 5);

        let z_carbon = AtomicNumber::from_string("C");
        let z_hydrogen = AtomicNumber::from_string("H");

        assert!(file.atomic_numbers.contains(&z_carbon));
        assert!(file.atomic_numbers.contains(&z_hydrogen));

        std::fs::remove_file("tmp1.xyz");
    }

    #[test]
    #[should_panic]
    fn test_no_coords_file_read(){
        std::fs::write(format!("tmp2.xyz"),
                       "5\n\n")
            .expect("Failed to write tmp2.xyz!");

        let file = XYZFile::new("tmp2.xyz");
        std::fs::remove_file("tmp2.xyz");
        file.unwrap();
    }


    #[test]
    #[should_panic]
    fn test_wrong_format_file_read(){
        std::fs::write(format!("tmp3.xyz"),
                       "1\n\n\
                        H     0.00000   0.00000\n")
            .expect("Failed to write tmp3.xyz!");

        let file = XYZFile::new("tmp3.xyz");
        std::fs::remove_file("tmp3.xyz");
        file.unwrap();
    }
}
