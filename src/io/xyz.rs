// xyz file parsing. See https://en.wikipedia.org/wiki/XYZ_file_format for a description
use std::fs::File;
use std::path::Path;
use std::io::{self, BufRead, Write};
use log::{info, warn};
use crate::atoms::{CartesianCoordinate, AtomicNumber};


#[derive(Default)]
pub struct XYZFile{

    pub filename:       String,
    file_lines:         Vec<String>,

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

                match xyz_file.add_atom_on_line(line){
                    Ok(..) => info!("Added line"),
                    Err(..) => warn!("Failed to add line")
                }
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

        let atomic_number = AtomicNumber::from_option_string(items.next());

        match atomic_number {
            Ok(a) => self.atomic_numbers.push(a),
            Err(..) => return Err("Failed to parse the atomic symbol"),
        }

        let coord = CartesianCoordinate::from_option_strings(
            items.next(), items.next(), items.next()
        );

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
        print_methane_xyz_file("tmp_methane_2.xyz");

        let file = XYZFile::new("tmp_methane_2.xyz").unwrap();

        assert_eq!(file.coordinates.len(), 5);
        assert_eq!(file.atomic_numbers.len(), 5);

        let z_carbon = AtomicNumber::from_string("C").unwrap();
        let z_hydrogen = AtomicNumber::from_string("H").unwrap();

        assert!(file.atomic_numbers.contains(&z_carbon));
        assert!(file.atomic_numbers.contains(&z_hydrogen));

        std::fs::remove_file("tmp_methane_2.xyz").expect("Failed to remove file");
    }

    #[test]
    #[should_panic]
    fn test_no_coords_file_read(){
        std::fs::write(format!("tmp2.xyz"),
                       "5\n\n")
            .expect("Failed to write tmp2.xyz!");

        let file = XYZFile::new("tmp2.xyz");
        std::fs::remove_file("tmp2.xyz").expect("Failed to remove file");
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
        std::fs::remove_file("tmp3.xyz").expect("Failed to remove file");
        file.unwrap();
    }
}
