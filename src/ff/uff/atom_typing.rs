use crate::atoms::{Atom, AtomicNumber};
use crate::connectivity::bonds::BondOrder;
use crate::coordinates::angle_value;
use crate::molecule::Molecule;
use std::f64::consts::PI;

/// See Table 1 in J. Am. Chem. Soc. 1992, 114, 25, 10024–10035
/// https://doi.org/10.1021/ja00051a040

#[derive(Default, Debug, Clone, PartialEq)]
pub(crate) struct UFFAtomType {
    pub name: &'static str,          // Standard name of the type
    pub atomic_symbol: &'static str, // Atomic symbol string
    pub bridging: bool,              // Is this bridging?
    pub aromatic: bool,              // Is this aromatic?
    pub valency: usize,              // Number of bonded neighbours
    pub oxidation_state: usize,      // Formal charge
    pub environment: CoordinationEnvironment,

    pub r: f64,     // Bonded distance (Å)
    pub theta: f64, // Angle (radians)
    pub x: f64,     // Non-bonded distance (Å)
    pub d: f64,     // Non-bonded energy (kcal mol-1)
    pub zeta: f64,  // Non-bonded scale
    pub z_eff: f64, // Effective charge (e)
    pub v_phi: f64, // Torsional potential for an sp3-bonded pair
}

#[derive(Debug, Hash, PartialEq, Eq, Clone, Default)]
pub enum CoordinationEnvironment {
    #[default]
    None,
    Linear,
    Bent,
    TrigonalPlanar,
    TrigonalPyramidal,
    SquarePlanar,
    Tetrahedral,
    TrigonalBipyramidal,
    Octahedral,
    Unknown,
}

/// sp 'hybridisation' of a particular element
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub enum Hybridisation {
    SP3,
    SP2,
    SP,
    None,
}

impl UFFAtomType {
    /// How well does an atom within a molecule match this atom type. Larger value <=> better match
    pub(crate) fn match_quality(&self, atom: &Atom, molecule: &Molecule) -> f64 {
        let mut value: f64 = 0.;

        if atom.atomic_symbol() == self.atomic_symbol {
            value += 10.;
        }

        if self.name == "O_3_z" {
            return 0.; // Noting is considered a zeolite
        }

        let num_neighbours = atom.bonded_neighbours.len();
        value -= (num_neighbours as f64 - self.valency as f64).abs();

        // Match on current angle
        if num_neighbours > 1 {
            let angle = angle_value(
                atom.bonded_neighbours[0],
                atom.idx,
                atom.bonded_neighbours[1],
                &molecule.coordinates,
            );

            // Subtract at most one for the angle not matching
            value -= (angle - self.theta).abs() / PI;
        }

        let num_aromatic_bonds = molecule
            .bonds()
            .iter()
            .filter(|b| b.contains(atom) && b.order == BondOrder::Aromatic)
            .count();

        if num_aromatic_bonds == 2 && self.aromatic {
            value += 5.;
        };

        // TODO: Match on more things?

        value
    }

    /// GMP (generalized Mulliken-Pauling) electronegativity
    pub fn gmp_electronegativity(&self) -> f64 {
        AtomicNumber::from_string(self.atomic_symbol)
            .unwrap()
            .gmp_electronegativity()
    }

    /// Set the coordination environment of this atom type given an atom with bonded neighbours
    pub fn set_coordination_environment(&mut self, atom: &Atom) {
        let valency = atom.bonded_neighbours.len();

        match valency {
            0 => {
                self.environment = CoordinationEnvironment::None;
            }

            1 => {
                self.environment = CoordinationEnvironment::Linear;
            }

            2 => match atom.group() {
                15 | 16 => {
                    self.environment = CoordinationEnvironment::Bent;
                }
                _ => {
                    self.environment = CoordinationEnvironment::Linear;
                }
            },

            3 => match atom.group() {
                15 => {
                    self.environment = CoordinationEnvironment::TrigonalPyramidal;
                }
                _ => {
                    self.environment = CoordinationEnvironment::TrigonalPlanar;
                }
            },

            4 => {
                if atom.is_d8() || atom.atomic_symbol() == "Xe" {
                    self.environment = CoordinationEnvironment::SquarePlanar;
                } else {
                    self.environment = CoordinationEnvironment::Tetrahedral;
                }
            }
            5 => {
                self.environment = CoordinationEnvironment::TrigonalBipyramidal;
            }

            6 => {
                self.environment = CoordinationEnvironment::Octahedral;
            }

            _ => {
                self.environment = CoordinationEnvironment::Unknown;
            }
        }
    }

    /// What type of bend is this atom? Defined by its coordination environment
    pub fn bend_type(&self) -> char {
        let type_a_environments = [
            CoordinationEnvironment::Linear,
            CoordinationEnvironment::TrigonalPlanar,
            CoordinationEnvironment::SquarePlanar,
            CoordinationEnvironment::Octahedral,
        ];

        match type_a_environments.contains(&self.environment) {
            true => 'A',
            false => 'B',
        }
    }

    /// n defining the equilibrium bend angle π/n
    pub fn bend_n(&self) -> f64 {
        match self.environment {
            CoordinationEnvironment::Linear => 4., // This differs from the UFF paper
            CoordinationEnvironment::TrigonalPlanar => 3.,
            CoordinationEnvironment::SquarePlanar => 4.,
            CoordinationEnvironment::Octahedral => 4.,
            _ => 0.,
        }
    }

    /// Is this atom type from the main group block of the periodic table?
    pub fn is_main_group(&self) -> bool {
        Atom::from_atomic_symbol(self.atomic_symbol).is_main_group()
    }

    /// Hybridisation of this atom type
    pub fn hybridisation(&self) -> Hybridisation {
        let atom = Atom::from_atomic_symbol(self.atomic_symbol);

        match atom.group() {
            14 => match self.valency {
                4 => Hybridisation::SP3,
                3 => Hybridisation::SP2,
                2 => Hybridisation::SP,
                _ => Hybridisation::None,
            },
            15 => match self.valency {
                3 | 4 => Hybridisation::SP3,
                2 => Hybridisation::SP2,
                1 => Hybridisation::SP,
                _ => Hybridisation::None,
            },
            16 => match self.valency {
                2 => Hybridisation::SP3,
                1 | 3 => Hybridisation::SP2,
                _ => Hybridisation::None,
            },
            _ => Hybridisation::None,
        }
    }

    /// Torsional U potential defined for sp2-sp2 central bonds within a dihedral
    pub fn u_phi(&self) -> f64 {
        let period = AtomicNumber::from_string(self.atomic_symbol)
            .unwrap()
            .period();

        match period {
            2 => 2.,
            3 => 1.25,
            4 => 0.7,
            5 => 0.2,
            6 => 0.1,
            _ => 0.,
        }
    }
}
