use crate::atoms::Atom;
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

    pub r:     f64,                     // Bonded distance (Å)
    pub theta: f64,                     // Angle (radians)
    pub x:     f64,                     // Non-bonded distance (Å)
    pub d:     f64,                     // Non-bonded energy (kcal mol-1)
    pub zeta:  f64,                     // Non-bonded scale
    pub z_eff: f64                      // Effective charge (e)
}


impl UFFAtomType {

    /// How well does an atom within a molecule match this atom type. Larger value <=> better match
    pub(crate) fn match_quality(&self,
                                atom:     &Atom,
                                molecule: &Molecule) -> i64{

        let mut value: i64 = 0;

        value += self.match_quality_atomic_symbol(atom);
        value -= (atom.num_bonded_neighbours(molecule.bonds()) as i64 - self.valency as i64).abs();

        // TODO: Match on more things

        value
    }

    /// How well does the atomic symbol match?
    fn match_quality_atomic_symbol(&self, atom: &Atom) -> i64{

        if atom.atomic_symbol() == self.atomic_symbol{ 2 }
        else { 0 }
    }
}


