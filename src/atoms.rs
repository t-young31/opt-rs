use crate::connectivity::bonds::{Bond, BondOrder};
use crate::coordinates::Point;
use crate::utils::is_very_close;
use log::{info, warn};
use std::ops::Index;

#[derive(Default, Debug, Clone)]
pub struct Atom {
    pub(crate) idx: usize,
    pub(crate) atomic_number: AtomicNumber,
    pub(crate) coordinate: Point,
    pub(crate) bonded_neighbours: Vec<usize>,
    pub(crate) formal_charge: f64,
}

impl Atom {
    /// Minimal constructor for an atom with an atomic index of 0
    pub fn from_atomic_symbol(atomic_symbol: &str) -> Self {
        Atom::from_idx_and_atomic_symbol(0, atomic_symbol)
    }

    /// Construct an atom from only an atomic index and an atomic symbol
    pub fn from_idx_and_atomic_symbol(idx: usize, atomic_symbol: &str) -> Self {
        Atom {
            idx,
            atomic_number: AtomicNumber::from_string(atomic_symbol).unwrap(),
            coordinate: Point::default(),
            bonded_neighbours: Default::default(),
            formal_charge: 0.0,
        }
    }

    /// Atomic symbol of this atom
    pub fn atomic_symbol(&self) -> &str {
        self.atomic_number.to_atomic_symbol()
    }

    /// Evaluate the distance to another atom, slowly
    pub fn distance_to(&self, atom: &Atom) -> f64 {
        let c0 = &self.coordinate;
        let c1 = &atom.coordinate;

        ((c0.x - c1.x).powi(2) + (c0.y - c1.y).powi(2) + (c0.z - c1.z).powi(2)).sqrt()
    }

    /// Could this atom be bonded to another atom, based on the distance
    /// between them and their respective covalent radii?
    pub fn could_be_bonded_to(&self, atom: &Atom) -> bool {
        // Relative tolerance on whether a bond could be present
        let tolerance = 1.3f64;

        let r = self.distance_to(atom);
        let is_identical_atom = r < 1E-8;

        !is_identical_atom && r < tolerance * (self.covalent_radius() + atom.covalent_radius())
    }

    /// Maximum coordination number i.e. valance of this atom
    pub fn maximal_valence(&self) -> usize {
        self.atomic_number.maximal_valence()
    }

    /// Covalent radius of this atom
    pub fn covalent_radius(&self) -> f64 {
        self.atomic_number.covalent_radius()
    }

    /// Group that this atom is in within the periodic table
    pub fn group(&self) -> usize {
        self.atomic_number.group()
    }

    /// Period that this atom is in within the periodic table
    pub fn period(&self) -> usize {
        self.atomic_number.period()
    }

    /// Is this atom a d8 transition metal, defined by
    pub fn is_d8(&self) -> bool {
        self.is_a_transition_metal() && is_very_close(self.dn(), 8.0)
    }

    /// Number of d electrons
    pub fn dn(&self) -> f64 {
        match self.is_a_transition_metal() {
            true => self.group() as f64 - self.formal_charge,
            false => 0.0,
        }
    }

    /// Is this atom a transition metal?
    pub fn is_a_transition_metal(&self) -> bool {
        let a = self.atomic_number.value;

        a < 31 && a > 20 || a < 48 && a > 38 || a < 81 && a > 71 || a < 113 && a > 103
    }

    /// This this atom a metal?
    pub fn is_metal(&self) -> bool {
        METALLIC_ELEMENTS.contains(&self.atomic_symbol())
    }

    /// Is this atom a from the main group block in the periodic table?
    pub fn is_main_group(&self) -> bool {
        return MAIN_GROUP_ELEMENTS.contains(&self.atomic_number.to_atomic_symbol());
    }

    /// Number of valance electrons that could participate in bonding
    pub fn num_valance_electrons(&self) -> i32 {
        match self.group() {
            13 => 3,
            14 => 4,
            15 => 5,
            16 => 6,
            17 => 7,
            _ => self.group() as i32,
        }

        // TODO: Metals
    }

    /// Number of unpaired electrons that may be present on this atom
    pub fn num_possible_unpaired_electrons(&self) -> i32 {
        if self.group() < 13 {
            warn!("Unknown number of unpaired electrons. Returning 0");
            return 0;
        }

        let mut num_lone_pair_electrons: i32 = 0;
        match self.group() {
            15 => {
                num_lone_pair_electrons = 2;
            }
            16 => {
                num_lone_pair_electrons = 4;
            }
            17 => {
                num_lone_pair_electrons = 6;
            }
            _ => {}
        }

        self.num_valance_electrons() - self.bonded_neighbours.len() as i32 - num_lone_pair_electrons
    }

    /// Can this atom form multiple bonds with others?
    pub fn can_form_multiple_bonds(&self) -> bool {
        static GROUPS: [usize; 3] = [14, 15, 16];

        GROUPS.contains(&self.group())
    }

    /// Is this atom hypervalent given the surrounding bonds. Only useful for main group elements
    pub fn is_hypervalent(&self, bonds: &[Bond]) -> bool {
        let mut n = self.group() as f64 - 10.;

        for bond in bonds.iter() {
            if bond.contains(self) {
                n += bond.order.value();
            }
        }

        self.is_main_group() && n > 8.0
    }

    /// Attempt to reduce any double bonds to aromatic
    pub fn reduce_two_double_bonds_to_aromatic(&self, bonds: &mut [Bond]) {
        let n_double_bonds: usize = bonds
            .iter()
            .filter(|b| b.contains(self) && b.order == BondOrder::Double)
            .count();

        if n_double_bonds < 2 {
            return; // Cannot reduce double bonds to aromatic
        }

        for bond in bonds.iter_mut() {
            if bond.contains(self) && bond.order == BondOrder::Double {
                bond.order = BondOrder::Aromatic;
            }
        }
    }

    /// If this atom is hypervalent and includes a triple bond reduce the bond order to 1,
    /// i.e. a regular single bond
    pub fn reduce_triple_bond_to_single(&self, bonds: &mut [Bond]) {
        for bond in bonds.iter_mut() {
            if bond.contains(self)
                && bond.order == BondOrder::Triple
                && self.bonded_neighbours.len() > 2
            {
                bond.order = BondOrder::Single;
            }
        }
    }
}

impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        self.idx == other.idx
    }
}

impl Eq for Atom {}

#[derive(Default, Clone, Debug)]
pub struct AtomicNumber {
    value: usize,
}

impl AtomicNumber {
    /// Zero indexed index of this atomic number in the periodic table
    fn index(&self) -> usize {
        self.value - 1
    }

    /// Create an atomic number, must be present in the elements
    pub fn from_integer(value: usize) -> Result<Self, &'static str> {
        if value == 0 || value > ELEMENTS.len() {
            return Err("Cannot create an atomic number. Not present in the periodic table");
        }

        Ok(AtomicNumber { value })
    }

    /// Construct an atomic number from an optional string
    pub fn from_option_string(value: Option<&str>) -> Result<Self, &'static str> {
        if value.is_none() {
            return Err("Cannot create atomic number. Optional string was none");
        }

        AtomicNumber::from_string(value.unwrap())
    }

    /// Convert a string (slice) into an atomic number. Panics if not present
    pub fn from_string(value: &str) -> Result<Self, &'static str> {
        if !ELEMENTS.contains(&value) {
            return Err("Cannot create an atomic number from a symbol not in the periodic table");
        }

        Ok(AtomicNumber {
            value: ELEMENTS.iter().position(|&x| x == value).unwrap() + 1,
        })
    }

    /// Convert an atomic number to an atomic symbol
    pub fn to_atomic_symbol(&self) -> &str {
        ELEMENTS.index(self.index())
    }

    /// Covalent radius, used to determine if two atoms are bonded
    pub fn covalent_radius(&self) -> f64 {
        let radius = COVALENT_RADII_PICOMETERS.get(self.index());

        if radius.is_none() {
            warn!(
                "Covalent radius for atom {} is not defined. Guessing at 2 Å",
                self.index()
            );
            return 2.0;
        }

        *radius.unwrap() * PICOMETERS_TO_ANGSTROMS
    }

    /// The maximum number of bonds that this atom can form under 'reasonable' circumstances
    pub fn maximal_valence(&self) -> usize {
        let valance = MAXIMAL_VALENCIES.get(self.index());

        if valance.is_none() {
            warn!(
                "Did not find a value maximal valence value for atomic number {}",
                self.value
            );
            return 6;
        }

        *valance.unwrap()
    }

    /// GMP (generalized Mulliken-Pauling) electronegativity from https://doi.org/10.1021/j100161a070
    pub fn gmp_electronegativity(&self) -> f64 {
        match GMP_ELECTRONEGATIVITIES.get(self.index()) {
            Some(value) => *value,
            None => {
                warn!("Using a default value of the electronegativity");
                5.0
            }
        }
    }

    /// Group in the periodic table of this atom
    pub fn group(&self) -> usize {
        let index = self.index();
        let mut x = index;

        if index == 0 {
            return 1;
        }

        if (56..=70).contains(&index) || (88..=102).contains(&index) {
            info!("Lanthanides and actinides do not have a defined group. Returning 0");
            return 0;
        }

        x += 16; // Period 1

        if index > 3 {
            x += 10
        }
        if index > 11 {
            x += 10
        }
        if index > 55 {
            x -= 14
        }
        if index > 87 {
            x -= 14
        }

        x % 18 + 1
    }

    /// Period of this atomic number
    pub fn period(&self) -> usize {
        if self.value > 0 && self.value < 3 {
            return 1;
        }
        if self.value > 2 && self.value < 11 {
            return 2;
        }
        if self.value > 10 && self.value < 19 {
            return 3;
        }
        if self.value > 18 && self.value < 37 {
            return 4;
        }
        if self.value > 36 && self.value < 55 {
            return 5;
        }
        if self.value > 54 && self.value < 87 {
            return 6;
        }
        if self.value > 86 && self.value < 119 {
            return 7;
        }

        panic!("Unknown period - retuning 0");
    }

    /// Is this atomic number from the main group?
    pub fn is_main_group(&self) -> bool {
        return MAIN_GROUP_ELEMENTS.contains(&self.to_atomic_symbol());
    }
}

impl PartialEq for AtomicNumber {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

static ELEMENTS: [&str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];

static METALLIC_ELEMENTS: [&str; 92] = [
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf",
    "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
];

static MAIN_GROUP_ELEMENTS: [&str; 36] = [
    "B", "C", "N", "O", "F", "Ne", "Al", "Si", "P", "S", "Cl", "Ar", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "In", "Sn", "Sb", "Te", "I", "Xe", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Nh", "Fl", "Mc",
    "Lv", "Ts", "Og",
];

static PICOMETERS_TO_ANGSTROMS: f64 = 0.01;

static COVALENT_RADII_PICOMETERS: [f64; 86] = [
    31., 28., 128., 96., 84., 76., 71., 66., 57., 58., 166., 141., 121., 111., 107., 105., 102.,
    106., 102., 203., 176., 170., 160., 153., 139., 161., 152., 150., 124., 132., 122., 122., 120.,
    119., 120., 116., 220., 195., 190., 175., 164., 154., 147., 146., 142., 139., 145., 144., 142.,
    139., 139., 138., 139., 140., 244., 215., 207., 204., 203., 201., 199., 198., 198., 196., 194.,
    192., 192., 189., 190., 187., 175., 187., 170., 162., 151., 144., 141., 136., 136., 132., 145.,
    146., 148., 140., 150., 150.,
];

static MAXIMAL_VALENCIES: [usize; 38] = [
    1, 0, 1, 2, 3, 4, 5, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 7, 5, 4, 4, 6, 3, 4,
    5, 6, 7, 2, 1, 2,
];

// Values tabulated from https://doi.org/10.1021/acs.jctc.8b00669
static GMP_ELECTRONEGATIVITIES: [f64; 51] = [
    4.53, 9.66, 3.01, 4.88, 5.11, 5.34, 6.90, 8.74, 10.87, 11.04, 2.84, 3.95, 4.06, 4.17, 5.46,
    6.93, 8.56, 9.47, 2.42, 3.23, 3.40, 3.47, 3.65, 3.42, 3.33, 3.76, 4.11, 4.47, 4.20, 5.11, 3.64,
    4.05, 5.19, 6.43, 7.79, 8.51, 2.33, 3.02, 3.83, 3.40, 3.55, 3.47, 3.29, 3.58, 3.98, 4.32, 4.44,
    5.03, 3.51, 3.99, 4.90,
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
mod tests {
    use super::*;

    /// Given an atomic number from a number or string, then they should be identical
    #[test]
    fn test_atom_init() {
        assert_eq!(
            AtomicNumber::from_integer(1),
            AtomicNumber::from_string("H")
        );
    }

    /// Test that the group is correct for some random elements
    #[test]
    fn test_atom_group() {
        fn group(symbol: &str) -> usize {
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
        assert_eq!(group("Yb"), 0);
    }

    /// Test whether atoms are metals
    #[test]
    fn test_metallic_atoms() {
        assert!(Atom::from_atomic_symbol("Na").is_metal());
        assert!(!Atom::from_atomic_symbol("C").is_metal());
        assert!(Atom::from_atomic_symbol("Pt").is_metal());
        assert!(Atom::from_atomic_symbol("Pt").is_a_transition_metal());
    }

    /// Test main group atoms
    #[test]
    fn test_main_group_atoms() {
        assert!(Atom::from_atomic_symbol("P").is_main_group());
        assert!(AtomicNumber::from_string("P").unwrap().is_main_group());

        assert!(Atom::from_atomic_symbol("C").is_main_group());
        assert!(!Atom::from_atomic_symbol("Pt").is_main_group());
    }

    /// Test that the period is correct for some random elements
    #[test]
    fn test_atom_period() {
        fn period(symbol: &str) -> usize {
            AtomicNumber::from_string(symbol).unwrap().period()
        }

        assert_eq!(period("H"), 1);
        assert_eq!(period("He"), 1);
        assert_eq!(period("Be"), 2);
        assert_eq!(period("C"), 2);
        assert_eq!(period("P"), 3);
        assert_eq!(period("Mo"), 5);
        assert_eq!(period("Po"), 6);
        assert_eq!(period("Ds"), 7);
        assert_eq!(period("Cd"), 5);
        assert_eq!(period("Yb"), 6);
    }

    #[test]
    #[should_panic]
    fn test_impossible_period() {
        let _ = AtomicNumber { value: 200 }.period();
    }

    #[test]
    fn test_num_possible_unpaired_electrons() {
        assert_eq!(
            Atom::from_atomic_symbol("Cl").num_possible_unpaired_electrons(),
            1
        );
    }

    /// Test that a carbon with 5 bonds is hypervalent
    #[test]
    fn test_hypervalency_carbon() {
        let atom = Atom::from_atomic_symbol("C");
        let mut bonds: Vec<Bond> = Default::default();
        for i in 1..=5 {
            bonds.push(Bond::from_atom_indices(0, i))
        }

        assert!(atom.is_hypervalent(&bonds));

        // 4 bonds where one is a double bond is also hypervalent
        bonds.remove(4);
        bonds[3].order = BondOrder::Double;

        assert!(atom.is_hypervalent(&bonds));
    }
}
