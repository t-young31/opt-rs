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

/*
w -> i
x -> j
y -> k
z -> l

vec_xw = self.atoms[w].coord - self.atoms[x].coord
vec_yz = self.atoms[z].coord - self.atoms[y].coord
vec_xy = self.atoms[y].coord - self.atoms[x].coord

vec1, vec2 = np.cross(vec_xw, vec_xy), np.cross(-vec_xy, vec_yz)

for vec in (vec1, vec2, vec_xy):
    vec /= np.linalg.norm(vec)


value = -np.arctan2(np.dot(np.cross(vec1, vec_xy), vec2),
                    np.dot(vec1, vec2))

 */


    fn energy(&self, coordinates: &Vec<Point>) -> f64 {

        todo!();
        let phi = 0.;

        let r_ij = &coordinates[self.i] - &coordinates[self.j];
        let r_lk = &coordinates[self.l] - &coordinates[self.k];
        let r_kj = &coordinates[self.k] - &coordinates[self.j];

        let v0 = r_ij.cross(&r_kj);
        let v1 = -r_kj.cross(&r_lk);

        0.5 * self.v_phi * (1. - (self.n_phi * self.phi0).cos()*(self.n_phi * phi).cos())
    }

    fn add_gradient(&self,
                    coordinates:      &Vec<Point>,
                    current_gradient: &mut Vec<Point>) {
        todo!()
    }
}





