use crate::atoms::{Atom, AtomicNumber};
use crate::ff::uff::atom_typing::Hybridisation::SP;
use crate::Molecule;

/// See Table 1 in J. Am. Chem. Soc. 1992, 114, 25, 10024–10035
/// https://doi.org/10.1021/ja00051a040


#[derive(Default, Debug, Clone, PartialEq)]
pub(crate) struct UFFAtomType{

    pub name:            &'static str,  // Standard name of the type
    pub atomic_symbol:   &'static str,  // Atomic symbol string
    pub bridging:        bool,          // Is this bridging?
    pub aromatic:        bool,          // Is this aromatic?
    pub valency:         usize,         // Number of bonded neighbours
    pub oxidation_state: usize,         // Formal charge
    pub environment:     CoordinationEnvironment,

    pub r:         f64,                 // Bonded distance (Å)
    pub theta:     f64,                 // Angle (radians)
    pub x:         f64,                 // Non-bonded distance (Å)
    pub d:         f64,                 // Non-bonded energy (kcal mol-1)
    pub zeta:      f64,                 // Non-bonded scale
    pub z_eff:     f64,                 // Effective charge (e)
    pub v_phi: f64                  // Torsional potential for an sp3-bonded pair
}


#[derive(Debug, Hash, PartialEq, Clone)]
pub enum CoordinationEnvironment{
    None,
    Linear,
    Bent,
    TrigonalPlanar,
    TrigonalPyramidal,
    SquarePlanar,
    Tetrahedral,
    TrigonalBipyramidal,
    Octahedral,
    Unknown
}

impl Default for CoordinationEnvironment {
    fn default() -> Self { CoordinationEnvironment::None }
}

/// sp 'hybridisation' of a particular element
#[derive(Debug, Hash, PartialEq, Clone)]
pub enum Hybridisation{
    SP3,
    SP2,
    SP,
    None
}


impl UFFAtomType {

    /// How well does an atom within a molecule match this atom type. Larger value <=> better match
    pub(crate) fn match_quality(&self,
                                atom:     &Atom,
                                molecule: &Molecule) -> i64{

        let mut value: i64 = 0;

        value += self.match_quality_atomic_symbol(atom);
        value -= (atom.bonded_neighbours.len() as i64 - self.valency as i64).abs();

        // TODO: Match on more things

        value
    }

    /// How well does the atomic symbol match?
    fn match_quality_atomic_symbol(&self, atom: &Atom) -> i64{

        if atom.atomic_symbol() == self.atomic_symbol{ 2 }
        else { 0 }
    }

    /// GMP (generalized Mulliken-Pauling) electronegativity
    pub fn gmp_electronegativity(&self) -> f64{

        AtomicNumber::from_string(self.atomic_symbol).unwrap().gmp_electronegativity()
    }

    /// Set the coordination environment of this atom type given an atom with bonded neighbours
    pub fn set_coordination_environment(&mut self, atom: &Atom){

        let valency = atom.bonded_neighbours.len();

        match valency {
            0 => {self.environment = CoordinationEnvironment::None;},

            1 => {self.environment = CoordinationEnvironment::Linear;}

            2 => {
                match atom.group() {
                    15 | 16 => {self.environment = CoordinationEnvironment::Bent;},
                    _  => {self.environment = CoordinationEnvironment::Linear;}
                }
            }

            3 => {
                match atom.group() {
                    15 => {self.environment = CoordinationEnvironment::TrigonalPyramidal;},
                    _ => {self.environment = CoordinationEnvironment::TrigonalPlanar;}
                }
            }

            4 => {
                if atom.is_d8() || atom.atomic_symbol() == "Xe"{
                    self.environment = CoordinationEnvironment::SquarePlanar;
                }
                else {self.environment = CoordinationEnvironment::Tetrahedral;}
            }
            5 => {self.environment = CoordinationEnvironment::TrigonalBipyramidal;}

            6 => {self.environment = CoordinationEnvironment::Octahedral;}

            _ => {self.environment = CoordinationEnvironment::Unknown;}
        }
    }

    /// What type of bend is this atom? Defined by its coordination environment
    pub fn bend_type(&self) -> char{

        let type_a_environments = [
            CoordinationEnvironment::Linear,
            CoordinationEnvironment::TrigonalPlanar,
            CoordinationEnvironment::SquarePlanar,
            CoordinationEnvironment::Octahedral
        ];

        match type_a_environments.contains(&self.environment) {
            true  => 'A',
            false => 'B'
        }
    }

    /// n defining the equilibrium bend angle π/n
    pub fn bend_n(&self) -> f64{

        match self.environment {
            CoordinationEnvironment::Linear =>         1.,
            CoordinationEnvironment::TrigonalPlanar => 3.,
            CoordinationEnvironment::SquarePlanar =>   4.,
            CoordinationEnvironment::Octahedral =>     4.,
            _ =>                                       0.
        }
    }

    /// Is this atom type from the main group block of the periodic table?
    pub fn is_main_group(&self) -> bool{
        return Atom::from_atomic_symbol(self.atomic_symbol).is_main_group()
    }

    /// Hybridisation of this atom type
    pub fn hybridisation(&self) -> Hybridisation{

        let atom = Atom::from_atomic_symbol(self.atomic_symbol);

        match atom.group(){
            14  => {
                match self.valency {
                    4 => Hybridisation::SP3,
                    3 => Hybridisation::SP2,
                    2 => Hybridisation::SP,
                    _ => Hybridisation::None
                }
            }
            15 => {
                match self.valency {
                    3 | 4 => Hybridisation::SP3,
                    2 => Hybridisation::SP2,
                    1 => Hybridisation::SP,
                    _ => Hybridisation::None
                }
            }
            16 => {
                match self.valency {
                    2  => Hybridisation::SP3,
                    1 | 3 => Hybridisation::SP2,
                    _ => Hybridisation::None
                }
            }
            _ => Hybridisation::None
        }
    }
}
