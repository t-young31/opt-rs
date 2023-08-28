pub struct InversionCentre {
    pub(crate) name: &'static str,
    pub(crate) k: f64,
    pub(crate) c0: f64,
    pub(crate) c1: f64,
    pub(crate) c2: f64,
}

pub static INVERSION_CENTERS: [InversionCentre; 4] = [
    // InversionCentre{name: "N_3",  k: 0.00000, c0: 0.00000, c1: 0.00000, c2: 0.00000},
    InversionCentre {
        name: "P_3",
        k: 1.08449,
        c0: 12.08079,
        c1: -3.72982,
        c2: -11.93509,
    },
    InversionCentre {
        name: "As_3",
        k: 1.15963,
        c0: 10.57534,
        c1: -2.12537,
        c2: -10.52168,
    },
    InversionCentre {
        name: "Sb_3",
        k: 1.10944,
        c0: 10.67772,
        c1: -1.49926,
        c2: -10.65134,
    },
    InversionCentre {
        name: "Bi_3",
        k: 1.00937,
        c0: 11.10497,
        c1: -0.41220,
        c2: -11.10306,
    },
];
