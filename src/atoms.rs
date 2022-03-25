use std::f32::consts::E;
use std::ops::Index;
use std::str::FromStr;

#[derive(Default)]
pub struct Atom{
    pub(crate) atomic_number: AtomicNumber,
    pub(crate) coordinate: CartesianCoordinate
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

        self.distance_to(atom) < (self.covalent_radius() + atom.covalent_radius())
    }

    fn covalent_radius(&self) -> f64{ self.atomic_number.covalent_radius() }
}


#[derive(Default, Clone)]
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

#[derive(Default, Clone)]
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

        // TODO:
        return 1.0;
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
