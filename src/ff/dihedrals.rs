use std::collections::HashSet;
use crate::coordinates::Point;
use crate::ff::forcefield::EnergyFunction;

/// Dihedral angle
///
///     E_i  --- E_j
///                 \  φ
///                  \
///                   E_k ---- E_l
///
pub struct TorsionalDihedral{

    pub(crate) i:  usize,
    pub(crate) j:  usize,
    pub(crate) k:  usize,
    pub(crate) l:  usize,

    pub(crate) phi0:  f64,
    pub(crate) n_phi: f64,
    pub(crate) v_phi: f64
}

impl EnergyFunction for TorsionalDihedral {
    fn involves_idxs(&self, idxs: Vec<usize>) -> bool {
        HashSet::from([self.i, self.j, self.k, self.l]) == HashSet::from_iter(idxs)
    }

    fn force_constant(&self) -> f64 {
        self.v_phi
    }

    /// E = V/2 (1 - cos(n φ_0)cos(n φ))
    fn energy(&self, coordinates: &Vec<Point>) -> f64 {

        let phi = phi(self, coordinates);
        0.5 * self.v_phi * (1. - (self.n_phi * self.phi0).cos()*(self.n_phi * phi).cos())
    }

    fn add_gradient(&self,
                    coordinates:      &Vec<Point>,
                    current_gradient: &mut Vec<Point>) {
        todo!()
    }
}


/// Value of the dihedral
fn phi(dihedral: &TorsionalDihedral, coordinates: &Vec<Point>) -> f64{

    let r_ij = &coordinates[dihedral.i] - &coordinates[dihedral.j];
    let r_lk = &coordinates[dihedral.l] - &coordinates[dihedral.k];
    let mut r_kj = &coordinates[dihedral.k] - &coordinates[dihedral.j];

    let mut v0 = r_ij.cross(&r_kj);
    v0.divide_by(v0.length());

    let mut v1 = (-r_kj.clone()).cross(&r_lk);

    v1.divide_by(v1.length());
    r_kj.divide_by(r_kj.length());

    -(v0.cross(&r_kj).dot(&v1)).atan2(v0.dot(&v1))
}


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
mod tests{

    use super::*;
    use crate::utils::*;

    /// Ensure a torsion can be calculated correctly. Compared to an Avogadro reference
    #[test]
    fn test_torsional_dihedral_calculation(){

        // Coordinates for a H2O2 molecule
        let x: Vec<Point> = Vec::from([
            Point{x: -3.80272, y: 0.11331, z: -0.29090},
            Point{x: -3.75673, y: 1.06948, z: -0.03303},
            Point{x: -2.46930, y: 1.33955, z:  0.03004},
            Point{x: -2.25240, y: 1.31540, z:  0.99712}
        ]);
        let dihedral = TorsionalDihedral{i: 0, j: 1, k: 2, l: 3, n_phi: 0., phi0: 0., v_phi: 0.};

        let expected_phi = -1.7593;  // -100.8 degrees

        assert!(is_close(phi(&dihedral, &x), expected_phi, 1E-3));
    }

}
