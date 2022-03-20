use std::ops::Index;

#[derive(Default)]
pub struct AtomicNumber{
    value: usize
}

impl AtomicNumber {

    /// Create an atomic number, must be present in the elements
    pub fn from_integer(value: usize) -> Self{
        if value <= 0 || value > ELEMENTS.len() {
            panic!("Cannot create an atomic number {}. Not present in the periodic table", value);
        }

        let mut atomic_number: AtomicNumber = Default::default();
        atomic_number.value = value;
        atomic_number
    }

    /// Convert an atomic number to an atomic symbol
    pub fn to_atomic_symbol(&self) -> &str{
        ELEMENTS.index(self.value)
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
