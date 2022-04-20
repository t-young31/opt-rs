use std::collections::HashSet;
use std::ops::Index;
use log::{warn};
use crate::connectivity::bonds::Bond;
use crate::coordinates::CartesianCoordinate;


#[derive(Default, Debug, Clone)]
pub struct Atom{
    pub(crate) idx:               usize,
    pub(crate) atomic_number:     AtomicNumber,
    pub(crate) coordinate:        CartesianCoordinate,
    pub(crate) bonded_neighbours: Vec<usize>
}

impl Atom {

    /// Construct an atom from only an atomic index and an atomic symbol
    pub fn from_idx_and_atomic_symbol(idx: usize, atomic_symbol: &str) -> Self{

        Atom{idx,
             atomic_number: AtomicNumber::from_string(atomic_symbol).unwrap(),
             coordinate: CartesianCoordinate::default(),
             bonded_neighbours: Default::default()}
    }

    /// Atomic symbol of this atom
    pub fn atomic_symbol(&self) -> &str{
        self.atomic_number.to_atomic_symbol()
    }

    /// Evaluate the distance to another atom, slowly
    pub fn distance_to(&self, atom: &Atom) -> f64{
        let c0 = &self.coordinate;
        let c1 = &atom.coordinate;

        ( (c0.x - c1.x).powi(2)
        + (c0.y - c1.y).powi(2)
        + (c0.z - c1.z).powi(2)).sqrt()
    }

    /// Could this atom be bonded to another atom, based on the distance
    /// between them and their respective covalent radii?
    pub fn could_be_bonded_to(&self, atom: &Atom) -> bool{

        // Relative tolerance on whether a bond could be present
        let tolerance = 1.3f64;

        let r = self.distance_to(atom);
        let is_identical_atom = r < 1E-8;

        !is_identical_atom && r < tolerance*(self.covalent_radius() + atom.covalent_radius())
    }

    /// Maximum coordination number i.e. valance of this atom
    pub fn maximal_valence(&self) -> usize{ self.atomic_number.maximal_valence() }

    /// Covalent radius of this atom
    fn covalent_radius(&self) -> f64{ self.atomic_number.covalent_radius() }

    /// Number of bonded neighbours
    pub fn num_bonded_neighbours(&self, bonds: &HashSet<Bond>) -> usize{

        bonds.iter().filter(|b| b.contains(self)).count()
    }

    pub fn group(&self) -> usize{self.atomic_number.group()}
}

impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        self.idx == other.idx
    }
}

impl Eq for Atom {}

#[derive(Default, Clone, Debug)]
pub struct AtomicNumber{
    value: usize
}

impl AtomicNumber {

    /// Zero indexed index of this atomic number in the periodic table
    fn index(&self) -> usize{
        return self.value - 1
    }

    /// Create an atomic number, must be present in the elements
    pub fn from_integer(value: usize) -> Result<Self, &'static str>{

        if value <= 0 || value > ELEMENTS.len() {
            return Err("Cannot create an atomic number. Not present in the periodic table");
        }

        Ok(AtomicNumber{value})
    }

    /// Construct an atomic number from an optional string
    pub fn from_option_string(value: Option<&str>) -> Result<Self, &'static str>{

        if value.is_none(){
            return Err("Cannot create atomic number. Optional string was none")
        }

        AtomicNumber::from_string(value.unwrap())
    }

    /// Convert a string (slice) into an atomic number. Panics if not present
    pub fn from_string(value: &str) -> Result<Self, &'static str>{

        if !ELEMENTS.contains(&value){
            return Err("Cannot create an atomic number from a symbol not in the periodic table");
        }

        Ok(AtomicNumber{value: ELEMENTS.iter().position(|&x| x == value).unwrap() + 1})
    }

    /// Convert an atomic number to an atomic symbol
    pub fn to_atomic_symbol(&self) -> &str{
        ELEMENTS.index(self.index())
    }

    /// Covalent radius, used to determine if two atoms are bonded
    pub fn covalent_radius(&self) -> f64{

        let radius = COVALENT_RADII_PICOMETERS.get(self.index());

        if radius.is_none(){
            warn!("Covalent radius for atom {} is not defined. Guessing at 2 Å", self.index());
            return 2.0;
        }

        return radius.unwrap().clone() * PICOMETERS_TO_ANGSTROMS;
    }

    /// The maximum number of bonds that this atom can form under 'reasonable' circumstances
    pub fn maximal_valence(&self) -> usize{

        let valance = MAXIMAL_VALENCIES.get(self.index());

        if valance.is_none(){
            warn!("Did not find a value maximal valence value for atomic number {}", self.value);
            return 6;
        }

        return valance.unwrap().clone();
    }

    /// GMP (generalized Mulliken-Pauling) electronegativity from https://doi.org/10.1021/j100161a070
    pub fn gmp_electronegativity(&self) -> f64{

        match GMP_ELECTRONEGATIVITIES.get(self.index()) {
            Some(value) => value.clone(),
            None => {
                warn!("Using a default value of the electronegativity");
                5.0
            }
        }
    }

    /// Group in the periodic table of this atom
    pub fn group(&self) -> usize{

        let index = self.index();
        let mut x = index;

        if index == 0{ return 1 }

        x += 16;  // Period 1

        if index > 3 { x += 10 }
        if index > 11{ x += 10 }
        if index > 55{ x -= 14 }
        if index > 87{ x -= 14 }

        x % 18 + 1
    }
}

impl PartialEq for AtomicNumber {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}


static ELEMENTS: [&'static str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
    "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
];

static PICOMETERS_TO_ANGSTROMS: f64 = 0.01;

static COVALENT_RADII_PICOMETERS: [f64; 86] = [
    31.,                                                                                                   28.,
    128., 96.,                                                               84.,  76.,  71.,  66.,  57.,  58.,
    166., 141.,                                                             121., 111., 107., 105., 102., 106.,
    102., 203., 176., 170., 160., 153., 139., 161., 152., 150., 124., 132., 122., 122., 120., 119., 120., 116.,
    220., 195., 190., 175., 164., 154., 147., 146., 142., 139., 145., 144., 142., 139., 139., 138., 139., 140.,
    244., 215.,
207., 204., 203., 201., 199., 198., 198., 196., 194., 192., 192., 189., 190., 187.,
                175., 187., 170., 162., 151., 144., 141., 136., 136., 132., 145., 146., 148., 140., 150., 150.
];

static MAXIMAL_VALENCIES: [usize; 38] = [
    1, 0, 1, 2, 3, 4, 5, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7,
    0, 1, 2, 3, 4, 5, 6, 7, 7, 5, 4, 4, 6, 3, 4, 5, 6, 7, 2, 1, 2];


// Values tabulated from https://doi.org/10.1021/acs.jctc.8b00669
static GMP_ELECTRONEGATIVITIES: [f64; 51] = [
    4.53, 9.66, 3.01, 4.88, 5.11, 5.34, 6.90, 8.74, 10.87, 11.04, 2.84, 3.95, 4.06, 4.17,
    5.46, 6.93, 8.56, 9.47, 2.42, 3.23, 3.40, 3.47, 3.65, 3.42, 3.33, 3.76, 4.11, 4.47,
    4.20, 5.11, 3.64, 4.05, 5.19, 6.43, 7.79, 8.51, 2.33, 3.02, 3.83, 3.40, 3.55, 3.47,
    3.29, 3.58, 3.98, 4.32, 4.44, 5.03, 3.51, 3.99, 4.90
];


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

    /// Given an atomic number from a number or string, then they should be identical
    #[test]
    fn test_atom_init(){
        assert_eq!(AtomicNumber::from_integer(1), AtomicNumber::from_string("H"));
    }

    /// Test that the group is correct for some random elements
    #[test]
    fn test_atom_group(){

        fn group(symbol: &str) -> usize{
            AtomicNumber::from_string(symbol).unwrap().group()
        }

        assert_eq!(group("H"), 1);
        assert_eq!(group("He"), 18);
        assert_eq!(group("Be"), 2);
        assert_eq!(group("C"), 14);
        assert_eq!(group("P"), 15);
        assert_eq!(group("Mo"), 6);
        assert_eq!(group("Po"), 16);
        assert_eq!(group("Ds"), 10);
        assert_eq!(group("Cd"), 12);
    }
}
