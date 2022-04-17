/// See Table 1 in J. Am. Chem. Soc. 1992, 114, 25, 10024–10035
/// https://doi.org/10.1021/ja00051a040


#[derive(Default, Debug)]
pub(crate) struct UFFAtomType{

    pub name:            &'static str,  // Standard name of the type
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

