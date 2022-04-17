/// See Table 1 in J. Am. Chem. Soc. 1992, 114, 25, 10024–10035
/// https://doi.org/10.1021/ja00051a040


struct UFFAtomType{
    name:            str,     // Standard name of the type
    bridging:        bool,    // Is this bridging?
    aromatic:        bool,    // Is this aromatic?
    valency:         usize,   // Number of bonded neighbours
    oxidation_state: usize,   // Formal charge

    r:     f64,     // Bonded distance (Å)
    theta: f64,     // Angle (radians)
    x:     f64,     // Non-bonded distance (Å)
    d:     f64,     // Non-bonded energy (kcal mol-1)
    zeta:  f64,     // Non-bonded scale
    z_eff: f64      // Effective charge (e)
}

