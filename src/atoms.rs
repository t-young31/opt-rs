use std::f32::consts::E;
use std::ops::Index;
use std::ptr::eq;
use std::str::FromStr;
use log::{warn};

#[derive(Default, Debug)]
pub struct Atom{
    pub(crate) idx:           usize,
    pub(crate) atomic_number: AtomicNumber,
    pub(crate) coordinate:    CartesianCoordinate
}

impl Atom {

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
        let tolerance = 1.2f64;

        let r = self.distance_to(atom);
        let is_identical_atom = r < 1E-8;

        !is_identical_atom && r < tolerance*(self.covalent_radius() + atom.covalent_radius())
    }

    pub fn maximal_valence(&self) -> usize{ self.atomic_number.maximal_valence() }

    fn covalent_radius(&self) -> f64{ self.atomic_number.covalent_radius() }
}

#[derive(Default, Clone, Debug)]
pub struct CartesianCoordinate{
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl CartesianCoordinate {

    /// Create a cartesian coordinate from a set of optional strings
    pub fn from_option_strings(x: Option<&str>,
                               y: Option<&str>,
                               z: Option<&str>
                               ) -> Result<Self, &'static str>{

        for k in [x, y, z]{
            if k.is_none() {
                return Err("An optional was None")
            }
            if f64::from_str(k.unwrap()).is_err(){
                return Err("Failed to parse float")
            };
        }

        let coord = CartesianCoordinate{
            x: x.unwrap().parse::<f64>().unwrap(),
            y: y.unwrap().parse::<f64>().unwrap(),
            z: z.unwrap().parse::<f64>().unwrap()
        };

        Ok(coord)
    }
}

#[derive(Default, Clone, Debug)]
pub struct AtomicNumber{
    value: usize
}

impl AtomicNumber {

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

        Ok(AtomicNumber{value: ELEMENTS.iter().position(|&x| x == value).unwrap()})
    }

    /// Convert an atomic number to an atomic symbol
    pub fn to_atomic_symbol(&self) -> &str{
        ELEMENTS.index(self.value)
    }

    /// Covalent radius, used to determine if two atoms are bonded
    pub fn covalent_radius(&self) -> f64{

        let radius = COVALENT_RADII_PICOMETERS.get(self.value);

        if radius.is_none(){
            warn!("Covalent radius for atom {} is not defined. Guessing at 2 Å", self.value);
            return 2.0;
        }

        return radius.unwrap().clone() * PICOMETERS_TO_ANGSTROMS;
    }

    /// The maximum number of bonds that this atom can form under 'reasonable' circumstances
    pub fn maximal_valence(&self) -> usize{

        let valance = MAXIMAL_VALENCIES.get(self.value);

        if valance.is_none(){
            warn!("Did not find a value maximal valence value for atom {}", self.value);
            return 6;
        }

        return valance.unwrap().clone();
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
