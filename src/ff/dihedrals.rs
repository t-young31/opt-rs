use crate::connectivity::dihedrals::ImproperDihedral;
use crate::coordinates::{Point, Vector3D};
use crate::ff::forcefield::EnergyFunction;
use std::collections::HashSet;

/// Dihedral angle
///
/// ```text
///       i  ---  j
///                 \  φ
///                  \
///                    k ---- l
/// ```
pub struct TorsionalDihedral {
    pub(crate) i: usize,
    pub(crate) j: usize,
    pub(crate) k: usize,
    pub(crate) l: usize,

    pub(crate) phi0: f64,
    pub(crate) n_phi: f64,
    pub(crate) v_phi: f64,
}

impl TorsionalDihedral {
    #[inline(always)]
    fn half_v_phi(&self) -> f64 {
        self.v_phi / 2.0
    }
}

impl EnergyFunction for TorsionalDihedral {
    fn involves_idxs(&self, idxs: Vec<usize>) -> bool {
        HashSet::from([self.i, self.j, self.k, self.l]) == HashSet::from_iter(idxs)
    }

    fn force_constant(&self) -> f64 {
        self.v_phi
    }

    /// E = V/2 (1 - cos(n φ_0)cos(n φ))
    fn energy(&self, coordinates: &[Point]) -> f64 {
        let phi = phi(self, coordinates);
        self.half_v_phi() * (1. - (self.n_phi * self.phi0).cos() * (self.n_phi * phi).cos())
    }

    fn add_gradient(&self, coordinates: &[Point], gradient: &mut Vec<Vector3D>) {
        let x_i = coordinates[self.i].x;
        let y_i = coordinates[self.i].y;
        let z_i = coordinates[self.i].z;

        let x_j = coordinates[self.j].x;
        let y_j = coordinates[self.j].y;
        let z_j = coordinates[self.j].z;

        let x_k = coordinates[self.k].x;
        let y_k = coordinates[self.k].y;
        let z_k = coordinates[self.k].z;

        let x_l = coordinates[self.l].x;
        let y_l = coordinates[self.l].y;
        let z_l = coordinates[self.l].z;

        let v0 = -x_j + x_k;
        let v1 = y_j - y_k;
        let v2 = x_i - x_j;
        let v3 = -z_j + z_k;
        let v4 = -x_k + x_l;
        let v5 = z_i - z_j;
        let v6 = -v2 * v3 + v0 * v5;
        let v7 = v0.powi(2);
        let v8 = -y_j + y_k;
        let v9 = v8.powi(2);
        let v10 = v3.powi(2);
        let v11 = -y_k + y_l;
        let v12 = x_j - x_k;
        let v13 = v12 * v11 - v4 * v1;
        let v14 = y_i - y_j;
        let v15 = (v7 + v9 + v10).sqrt();
        let v16 = v6.powi(2);
        let v17 = v14 * v3 - v8 * v5;
        let v18 = v17.powi(2);
        let v19 = (v2 * v8 - v0 * v14).powi(2);
        let v20 = (v19 + v16 + v18).sqrt();
        let v21 = v13.powi(2);
        let v22 = v15 * v20;
        let v23 = -z_k + z_l;
        let v24 = z_j - z_k;
        let v25 = -v0 * v6 / v22 + v8 * v17 / v22;
        let v26 = v1 * v23 - v11 * v24;
        let v27 = v26.powi(2);
        let v28 = (-v12 * v23 + v4 * v24).powi(2);
        let v29 = -v12 * v23 + v4 * v24;
        let v30 = v2 * v8 - v0 * v14;
        let v31 = v0 * v30 / v22 - v3 * v17 / v22;
        let v32 = (v21 + v28 + v27).sqrt();
        let v33 = -v8 * v30 / v22 + v3 * v6 / v22;
        let v34 = v20 * v32;
        let v35 = -v13 * v25 / v32 - v29 * v31 / v32 - v26 * v33 / v32;
        gradient[self.i].x += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * (v8 * v13 / v34
                    + v24 * v29 / v34
                    + v30
                        * v13
                        * (-(-2. * y_j + 2. * y_k) * v30 / 2. - (2. * z_j - 2. * z_k) * v6 / 2.)
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v6 * v29
                        * (-(-2. * y_j + 2. * y_k) * v30 / 2. - (2. * z_j - 2. * z_k) * v6 / 2.)
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v17
                        * (-(-2. * y_j + 2. * y_k) * v30 / 2. - (2. * z_j - 2. * z_k) * v6 / 2.)
                        * v26
                        / ((v19 + v16 + v18).powf(1.5) * v32))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v13
                    * (-v0 * v24 / v22
                        - v0 * v6
                            * (-(-2. * y_j + 2. * y_k) * v30 / 2.
                                - (2. * z_j - 2. * z_k) * v6 / 2.)
                            / (v15 * (v19 + v16 + v18).powf(1.5))
                        + v8 * v17
                            * (-(-2. * y_j + 2. * y_k) * v30 / 2.
                                - (2. * z_j - 2. * z_k) * v6 / 2.)
                            / (v15 * (v19 + v16 + v18).powf(1.5)))
                    / v32
                    + v29
                        * (v0 * v8 / v22
                            + v0 * v30
                                * (-(-2. * y_j + 2. * y_k) * v30 / 2.
                                    - (2. * z_j - 2. * z_k) * v6 / 2.)
                                / (v15 * (v19 + v16 + v18).powf(1.5))
                            - v3 * v17
                                * (-(-2. * y_j + 2. * y_k) * v30 / 2.
                                    - (2. * z_j - 2. * z_k) * v6 / 2.)
                                / (v15 * (v19 + v16 + v18).powf(1.5)))
                        / v32
                    + v26
                        * (-v9 / v22
                            - v8 * v30
                                * (-(-2. * y_j + 2. * y_k) * v30 / 2.
                                    - (2. * z_j - 2. * z_k) * v6 / 2.)
                                / (v15 * (v19 + v16 + v18).powf(1.5))
                            + v3 * v24 / v22
                            + v3 * v6
                                * (-(-2. * y_j + 2. * y_k) * v30 / 2.
                                    - (2. * z_j - 2. * z_k) * v6 / 2.)
                                / (v15 * (v19 + v16 + v18).powf(1.5)))
                        / v32)
                    * (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
        gradient[self.i].y += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * (v12 * v13 / v34
                    + v3 * v26 / v34
                    + v30
                        * v13
                        * (-(2. * x_j - 2. * x_k) * v30 / 2. - (-2. * z_j + 2. * z_k) * v17 / 2.)
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v6 * v29
                        * (-(2. * x_j - 2. * x_k) * v30 / 2. - (-2. * z_j + 2. * z_k) * v17 / 2.)
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + (-(2. * x_j - 2. * x_k) * v30 / 2. - (-2. * z_j + 2. * z_k) * v17 / 2.)
                        * v17
                        * v26
                        / ((v19 + v16 + v18).powf(1.5) * v32))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v13
                    * (-v0
                        * v6
                        * (-(2. * x_j - 2. * x_k) * v30 / 2.
                            - (-2. * z_j + 2. * z_k) * v17 / 2.)
                        / (v15 * (v19 + v16 + v18).powf(1.5))
                        + v8 * v3 / v22
                        + v8 * (-(2. * x_j - 2. * x_k) * v30 / 2.
                            - (-2. * z_j + 2. * z_k) * v17 / 2.)
                            * v17
                            / (v15 * (v19 + v16 + v18).powf(1.5)))
                    / v32
                    + v29
                        * (v0 * v12 / v22
                            + v0 * v30
                                * (-(2. * x_j - 2. * x_k) * v30 / 2.
                                    - (-2. * z_j + 2. * z_k) * v17 / 2.)
                                / (v15 * (v19 + v16 + v18).powf(1.5))
                            - v10 / v22
                            - v3 * (-(2. * x_j - 2. * x_k) * v30 / 2.
                                - (-2. * z_j + 2. * z_k) * v17 / 2.)
                                * v17
                                / (v15 * (v19 + v16 + v18).powf(1.5)))
                        / v32
                    + v26
                        * (-v12 * v8 / v22
                            - v8 * v30
                                * (-(2. * x_j - 2. * x_k) * v30 / 2.
                                    - (-2. * z_j + 2. * z_k) * v17 / 2.)
                                / (v15 * (v19 + v16 + v18).powf(1.5))
                            + v3 * v6
                                * (-(2. * x_j - 2. * x_k) * v30 / 2.
                                    - (-2. * z_j + 2. * z_k) * v17 / 2.)
                                / (v15 * (v19 + v16 + v18).powf(1.5)))
                        / v32)
                    * (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
        gradient[self.i].z += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * (v0 * v29 / v34
                    + v1 * v26 / v34
                    + v30
                        * (-(-2. * x_j + 2. * x_k) * v6 / 2. - (2. * y_j - 2. * y_k) * v17 / 2.)
                        * v13
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v6 * (-(-2. * x_j + 2. * x_k) * v6 / 2. - (2. * y_j - 2. * y_k) * v17 / 2.)
                        * v29
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + (-(-2. * x_j + 2. * x_k) * v6 / 2. - (2. * y_j - 2. * y_k) * v17 / 2.)
                        * v17
                        * v26
                        / ((v19 + v16 + v18).powf(1.5) * v32))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v13
                    * (-v7 / v22
                        - v0 * v6
                            * (-(-2. * x_j + 2. * x_k) * v6 / 2.
                                - (2. * y_j - 2. * y_k) * v17 / 2.)
                            / (v15 * (v19 + v16 + v18).powf(1.5))
                        + v8 * v1 / v22
                        + v8 * (-(-2. * x_j + 2. * x_k) * v6 / 2.
                            - (2. * y_j - 2. * y_k) * v17 / 2.)
                            * v17
                            / (v15 * (v19 + v16 + v18).powf(1.5)))
                    / v32
                    + v29
                        * (v0
                            * v30
                            * (-(-2. * x_j + 2. * x_k) * v6 / 2.
                                - (2. * y_j - 2. * y_k) * v17 / 2.)
                            / (v15 * (v19 + v16 + v18).powf(1.5))
                            - v1 * v3 / v22
                            - v3 * (-(-2. * x_j + 2. * x_k) * v6 / 2.
                                - (2. * y_j - 2. * y_k) * v17 / 2.)
                                * v17
                                / (v15 * (v19 + v16 + v18).powf(1.5)))
                        / v32
                    + v26
                        * (v0 * v3 / v22
                            - v8 * v30
                                * (-(-2. * x_j + 2. * x_k) * v6 / 2.
                                    - (2. * y_j - 2. * y_k) * v17 / 2.)
                                / (v15 * (v19 + v16 + v18).powf(1.5))
                            + v3 * v6
                                * (-(-2. * x_j + 2. * x_k) * v6 / 2.
                                    - (2. * y_j - 2. * y_k) * v17 / 2.)
                                / (v15 * (v19 + v16 + v18).powf(1.5)))
                        / v32)
                    * (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
        gradient[self.j].x += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * ((y_i - y_k) * v13 / v34
                    + v11 * v30 / v34
                    + (-z_i + z_k) * v29 / v34
                    + (z_k - z_l) * v6 / v34
                    + v30
                        * v13
                        * (-(2. * y_i - 2. * y_k) * v30 / 2. - (-2. * z_i + 2. * z_k) * v6 / 2.)
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v30
                        * v13
                        * (-(-2. * y_k + 2. * y_l) * v13 / 2. - (2. * z_k - 2. * z_l) * v29 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + v6 * v29
                        * (-(2. * y_i - 2. * y_k) * v30 / 2. - (-2. * z_i + 2. * z_k) * v6 / 2.)
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v6 * v29
                        * (-(-2. * y_k + 2. * y_l) * v13 / 2. - (2. * z_k - 2. * z_l) * v29 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + v17
                        * (-(2. * y_i - 2. * y_k) * v30 / 2. - (-2. * z_i + 2. * z_k) * v6 / 2.)
                        * v26
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v17
                        * v26
                        * (-(-2. * y_k + 2. * y_l) * v13 / 2. - (2. * z_k - 2. * z_l) * v29 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5)))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    * (v11 * v25 / v32
                        + (z_k - z_l) * v31 / v32
                        + v13
                            * (-(-2. * y_k + 2. * y_l) * v13 / 2.
                                - (2. * z_k - 2. * z_l) * v29 / 2.)
                            * v25
                            / (v21 + v28 + v27).powf(1.5)
                        + v13
                            * (-v7 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v0 * v8 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v0 * (-z_i + z_k) / v22
                                - v0 * v6
                                    * (-(2. * y_i - 2. * y_k) * v30 / 2.
                                        - (-2. * z_i + 2. * z_k) * v6 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v8 * v17
                                    * (-(2. * y_i - 2. * y_k) * v30 / 2.
                                        - (-2. * z_i + 2. * z_k) * v6 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v6 / v22)
                            / v32
                        + v29
                            * (-(-2. * y_k + 2. * y_l) * v13 / 2.
                                - (2. * z_k - 2. * z_l) * v29 / 2.)
                            * v31
                            / (v21 + v28 + v27).powf(1.5)
                        + v29
                            * (v7 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v0 * (y_i - y_k) / v22
                                - v0 * v3 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v0 * v30
                                    * (-(2. * y_i - 2. * y_k) * v30 / 2.
                                        - (-2. * z_i + 2. * z_k) * v6 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - v3 * v17
                                    * (-(2. * y_i - 2. * y_k) * v30 / 2.
                                        - (-2. * z_i + 2. * z_k) * v6 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - v30 / v22)
                            / v32
                        + v26
                            * (-(-2. * y_k + 2. * y_l) * v13 / 2.
                                - (2. * z_k - 2. * z_l) * v29 / 2.)
                            * v33
                            / (v21 + v28 + v27).powf(1.5)
                        + v26
                            * (-v0 * v8 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v0 * v3 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - (y_i - y_k) * v8 / v22
                                - v8 * v30
                                    * (-(2. * y_i - 2. * y_k) * v30 / 2.
                                        - (-2. * z_i + 2. * z_k) * v6 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + (-z_i + z_k) * v3 / v22
                                + v3 * v6
                                    * (-(2. * y_i - 2. * y_k) * v30 / 2.
                                        - (-2. * z_i + 2. * z_k) * v6 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5)))
                            / v32)
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
        gradient[self.j].y += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * ((-x_i + x_k) * v13 / v34
                    + (x_k - x_l) * v30 / v34
                    + (z_i - z_k) * v26 / v34
                    + v23 * v17 / v34
                    + (-(-2. * x_i + 2. * x_k) * v30 / 2. - (2. * z_i - 2. * z_k) * v17 / 2.)
                        * v30
                        * v13
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + (-(-2. * x_i + 2. * x_k) * v30 / 2. - (2. * z_i - 2. * z_k) * v17 / 2.)
                        * v6
                        * v29
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + (-(-2. * x_i + 2. * x_k) * v30 / 2. - (2. * z_i - 2. * z_k) * v17 / 2.)
                        * v17
                        * v26
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v30
                        * v13
                        * (-(2. * x_k - 2. * x_l) * v13 / 2. - (-2. * z_k + 2. * z_l) * v26 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + v6 * v29
                        * (-(2. * x_k - 2. * x_l) * v13 / 2. - (-2. * z_k + 2. * z_l) * v26 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + (-(2. * x_k - 2. * x_l) * v13 / 2. - (-2. * z_k + 2. * z_l) * v26 / 2.)
                        * v17
                        * v26
                        / (v20 * (v21 + v28 + v27).powf(1.5)))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    * ((x_k - x_l) * v25 / v32
                        + v23 * v33 / v32
                        + v13
                            * (-(2. * x_k - 2. * x_l) * v13 / 2.
                                - (-2. * z_k + 2. * z_l) * v26 / 2.)
                            * v25
                            / (v21 + v28 + v27).powf(1.5)
                        + v13
                            * (-v0 * v8 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v0 * (-(-2. * x_i + 2. * x_k) * v30 / 2.
                                    - (2. * z_i - 2. * z_k) * v17 / 2.)
                                    * v6
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v9 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v8 * (z_i - z_k) / v22
                                + v8 * (-(-2. * x_i + 2. * x_k) * v30 / 2.
                                    - (2. * z_i - 2. * z_k) * v17 / 2.)
                                    * v17
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - v17 / v22)
                            / v32
                        + v29
                            * (-(2. * x_k - 2. * x_l) * v13 / 2.
                                - (-2. * z_k + 2. * z_l) * v26 / 2.)
                            * v31
                            / (v21 + v28 + v27).powf(1.5)
                        + v29
                            * ((-x_i + x_k) * v0 / v22
                                + v0 * v8 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v0 * (-(-2. * x_i + 2. * x_k) * v30 / 2.
                                    - (2. * z_i - 2. * z_k) * v17 / 2.)
                                    * v30
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - v8 * v3 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - (z_i - z_k) * v3 / v22
                                - v3 * (-(-2. * x_i + 2. * x_k) * v30 / 2.
                                    - (2. * z_i - 2. * z_k) * v17 / 2.)
                                    * v17
                                    / (v15 * (v19 + v16 + v18).powf(1.5)))
                            / v32
                        + (-(2. * x_k - 2. * x_l) * v13 / 2. - (-2. * z_k + 2. * z_l) * v26 / 2.)
                            * v26
                            * v33
                            / (v21 + v28 + v27).powf(1.5)
                        + v26
                            * (-(-x_i + x_k) * v8 / v22
                                - v9 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v8 * v3 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v8 * (-(-2. * x_i + 2. * x_k) * v30 / 2.
                                    - (2. * z_i - 2. * z_k) * v17 / 2.)
                                    * v30
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v3 * (-(-2. * x_i + 2. * x_k) * v30 / 2.
                                    - (2. * z_i - 2. * z_k) * v17 / 2.)
                                    * v6
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v30 / v22)
                            / v32)
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
        gradient[self.j].z += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * ((x_i - x_k) * v29 / v34
                    + v4 * v6 / v34
                    + (-y_i + y_k) * v26 / v34
                    + (y_k - y_l) * v17 / v34
                    + v30
                        * (-(2. * x_i - 2. * x_k) * v6 / 2. - (-2. * y_i + 2. * y_k) * v17 / 2.)
                        * v13
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v30
                        * v13
                        * (-(-2. * x_k + 2. * x_l) * v29 / 2. - (2. * y_k - 2. * y_l) * v26 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + v6 * (-(2. * x_i - 2. * x_k) * v6 / 2. - (-2. * y_i + 2. * y_k) * v17 / 2.)
                        * v29
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v6 * v29
                        * (-(-2. * x_k + 2. * x_l) * v29 / 2. - (2. * y_k - 2. * y_l) * v26 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + (-(2. * x_i - 2. * x_k) * v6 / 2. - (-2. * y_i + 2. * y_k) * v17 / 2.)
                        * v17
                        * v26
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + (-(-2. * x_k + 2. * x_l) * v29 / 2. - (2. * y_k - 2. * y_l) * v26 / 2.)
                        * v17
                        * v26
                        / (v20 * (v21 + v28 + v27).powf(1.5)))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    * (v4 * v31 / v32
                        + (y_k - y_l) * v33 / v32
                        + v13
                            * (-(-2. * x_k + 2. * x_l) * v29 / 2.
                                - (2. * y_k - 2. * y_l) * v26 / 2.)
                            * v25
                            / (v21 + v28 + v27).powf(1.5)
                        + v13
                            * (-(x_i - x_k) * v0 / v22
                                - v0 * v3 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v0 * v6
                                    * (-(2. * x_i - 2. * x_k) * v6 / 2.
                                        - (-2. * y_i + 2. * y_k) * v17 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + (-y_i + y_k) * v8 / v22
                                + v8 * v3 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v8 * (-(2. * x_i - 2. * x_k) * v6 / 2.
                                    - (-2. * y_i + 2. * y_k) * v17 / 2.)
                                    * v17
                                    / (v15 * (v19 + v16 + v18).powf(1.5)))
                            / v32
                        + v29
                            * (-(-2. * x_k + 2. * x_l) * v29 / 2.
                                - (2. * y_k - 2. * y_l) * v26 / 2.)
                            * v31
                            / (v21 + v28 + v27).powf(1.5)
                        + v29
                            * (v0 * v3 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v0 * v30
                                    * (-(2. * x_i - 2. * x_k) * v6 / 2.
                                        - (-2. * y_i + 2. * y_k) * v17 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - (-y_i + y_k) * v3 / v22
                                - v10 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v3 * (-(2. * x_i - 2. * x_k) * v6 / 2.
                                    - (-2. * y_i + 2. * y_k) * v17 / 2.)
                                    * v17
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v17 / v22)
                            / v32
                        + (-(-2. * x_k + 2. * x_l) * v29 / 2. - (2. * y_k - 2. * y_l) * v26 / 2.)
                            * v26
                            * v33
                            / (v21 + v28 + v27).powf(1.5)
                        + v26
                            * ((x_i - x_k) * v3 / v22
                                - v8 * v3 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v8 * v30
                                    * (-(2. * x_i - 2. * x_k) * v6 / 2.
                                        - (-2. * y_i + 2. * y_k) * v17 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v10 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v3 * v6
                                    * (-(2. * x_i - 2. * x_k) * v6 / 2.
                                        - (-2. * y_i + 2. * y_k) * v17 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - v6 / v22)
                            / v32)
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
        gradient[self.k].x += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * ((-y_i + y_j) * v13 / v34
                    + (y_j - y_l) * v30 / v34
                    + v5 * v29 / v34
                    + (-z_j + z_l) * v6 / v34
                    + v30
                        * v13
                        * (-(-2. * y_i + 2. * y_j) * v30 / 2. - (2. * z_i - 2. * z_j) * v6 / 2.)
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v30
                        * v13
                        * (-(2. * y_j - 2. * y_l) * v13 / 2. - (-2. * z_j + 2. * z_l) * v29 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + v6 * v29
                        * (-(-2. * y_i + 2. * y_j) * v30 / 2. - (2. * z_i - 2. * z_j) * v6 / 2.)
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v6 * v29
                        * (-(2. * y_j - 2. * y_l) * v13 / 2. - (-2. * z_j + 2. * z_l) * v29 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + (-(-2. * y_i + 2. * y_j) * v30 / 2. - (2. * z_i - 2. * z_j) * v6 / 2.)
                        * v17
                        * v26
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v17
                        * v26
                        * (-(2. * y_j - 2. * y_l) * v13 / 2. - (-2. * z_j + 2. * z_l) * v29 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5)))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    * ((y_j - y_l) * v25 / v32
                        + (-z_j + z_l) * v31 / v32
                        + v13
                            * (-(2. * y_j - 2. * y_l) * v13 / 2.
                                - (-2. * z_j + 2. * z_l) * v29 / 2.)
                            * v25
                            / (v21 + v28 + v27).powf(1.5)
                        + v13
                            * (-v0 * v12 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v0 * v5 / v22
                                - v0 * v6
                                    * (-(-2. * y_i + 2. * y_j) * v30 / 2.
                                        - (2. * z_i - 2. * z_j) * v6 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v12 * v8 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v8 * (-(-2. * y_i + 2. * y_j) * v30 / 2.
                                    - (2. * z_i - 2. * z_j) * v6 / 2.)
                                    * v17
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - v6 / v22)
                            / v32
                        + v29
                            * (-(2. * y_j - 2. * y_l) * v13 / 2.
                                - (-2. * z_j + 2. * z_l) * v29 / 2.)
                            * v31
                            / (v21 + v28 + v27).powf(1.5)
                        + v29
                            * (v0 * v12 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v0 * (-y_i + y_j) / v22
                                + v0 * v30
                                    * (-(-2. * y_i + 2. * y_j) * v30 / 2.
                                        - (2. * z_i - 2. * z_j) * v6 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - v12 * v3 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v3 * (-(-2. * y_i + 2. * y_j) * v30 / 2.
                                    - (2. * z_i - 2. * z_j) * v6 / 2.)
                                    * v17
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v30 / v22)
                            / v32
                        + v26
                            * (-(2. * y_j - 2. * y_l) * v13 / 2.
                                - (-2. * z_j + 2. * z_l) * v29 / 2.)
                            * v33
                            / (v21 + v28 + v27).powf(1.5)
                        + v26
                            * (-v12 * v8 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v12 * v3 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - (-y_i + y_j) * v8 / v22
                                - v8 * v30
                                    * (-(-2. * y_i + 2. * y_j) * v30 / 2.
                                        - (2. * z_i - 2. * z_j) * v6 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v5 * v3 / v22
                                + v3 * v6
                                    * (-(-2. * y_i + 2. * y_j) * v30 / 2.
                                        - (2. * z_i - 2. * z_j) * v6 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5)))
                            / v32)
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
        gradient[self.k].y += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * (v2 * v13 / v34
                    + (-x_j + x_l) * v30 / v34
                    + (-z_i + z_j) * v26 / v34
                    + (z_j - z_l) * v17 / v34
                    + v30
                        * (-(2. * x_i - 2. * x_j) * v30 / 2. - (-2. * z_i + 2. * z_j) * v17 / 2.)
                        * v13
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v30
                        * (-(-2. * x_j + 2. * x_l) * v13 / 2. - (2. * z_j - 2. * z_l) * v26 / 2.)
                        * v13
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + v6 * (-(2. * x_i - 2. * x_j) * v30 / 2.
                        - (-2. * z_i + 2. * z_j) * v17 / 2.)
                        * v29
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v6 * (-(-2. * x_j + 2. * x_l) * v13 / 2.
                        - (2. * z_j - 2. * z_l) * v26 / 2.)
                        * v29
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + (-(2. * x_i - 2. * x_j) * v30 / 2. - (-2. * z_i + 2. * z_j) * v17 / 2.)
                        * v17
                        * v26
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + (-(-2. * x_j + 2. * x_l) * v13 / 2. - (2. * z_j - 2. * z_l) * v26 / 2.)
                        * v17
                        * v26
                        / (v20 * (v21 + v28 + v27).powf(1.5)))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    * ((-x_j + x_l) * v25 / v32
                        + (z_j - z_l) * v33 / v32
                        + (-(-2. * x_j + 2. * x_l) * v13 / 2. - (2. * z_j - 2. * z_l) * v26 / 2.)
                            * v13
                            * v25
                            / (v21 + v28 + v27).powf(1.5)
                        + (-(-2. * x_j + 2. * x_l) * v13 / 2. - (2. * z_j - 2. * z_l) * v26 / 2.)
                            * v29
                            * v31
                            / (v21 + v28 + v27).powf(1.5)
                        + (-(-2. * x_j + 2. * x_l) * v13 / 2. - (2. * z_j - 2. * z_l) * v26 / 2.)
                            * v26
                            * v33
                            / (v21 + v28 + v27).powf(1.5)
                        + v13
                            * (-v0 * v1 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v0 * v6
                                    * (-(2. * x_i - 2. * x_j) * v30 / 2.
                                        - (-2. * z_i + 2. * z_j) * v17 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v8 * v1 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v8 * (-z_i + z_j) / v22
                                + v8 * (-(2. * x_i - 2. * x_j) * v30 / 2.
                                    - (-2. * z_i + 2. * z_j) * v17 / 2.)
                                    * v17
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v17 / v22)
                            / v32
                        + v29
                            * (v2 * v0 / v22
                                + v0 * v1 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v0 * v30
                                    * (-(2. * x_i - 2. * x_j) * v30 / 2.
                                        - (-2. * z_i + 2. * z_j) * v17 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - v1 * v3 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - (-z_i + z_j) * v3 / v22
                                - v3 * (-(2. * x_i - 2. * x_j) * v30 / 2.
                                    - (-2. * z_i + 2. * z_j) * v17 / 2.)
                                    * v17
                                    / (v15 * (v19 + v16 + v18).powf(1.5)))
                            / v32
                        + v26
                            * (-v2 * v8 / v22
                                - v8 * v1 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v8 * v30
                                    * (-(2. * x_i - 2. * x_j) * v30 / 2.
                                        - (-2. * z_i + 2. * z_j) * v17 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v1 * v3 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v3 * v6
                                    * (-(2. * x_i - 2. * x_j) * v30 / 2.
                                        - (-2. * z_i + 2. * z_j) * v17 / 2.)
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - v30 / v22)
                            / v32)
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
        gradient[self.k].z += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * ((-x_i + x_j) * v29 / v34
                    + (x_j - x_l) * v6 / v34
                    + v14 * v26 / v34
                    + (-y_j + y_l) * v17 / v34
                    + (-(-2. * x_i + 2. * x_j) * v6 / 2. - (2. * y_i - 2. * y_j) * v17 / 2.)
                        * v30
                        * v13
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + (-(-2. * x_i + 2. * x_j) * v6 / 2. - (2. * y_i - 2. * y_j) * v17 / 2.)
                        * v6
                        * v29
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + (-(-2. * x_i + 2. * x_j) * v6 / 2. - (2. * y_i - 2. * y_j) * v17 / 2.)
                        * v17
                        * v26
                        / ((v19 + v16 + v18).powf(1.5) * v32)
                    + v30
                        * v13
                        * (-(2. * x_j - 2. * x_l) * v29 / 2. - (-2. * y_j + 2. * y_l) * v26 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + v6 * v29
                        * (-(2. * x_j - 2. * x_l) * v29 / 2. - (-2. * y_j + 2. * y_l) * v26 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + (-(2. * x_j - 2. * x_l) * v29 / 2. - (-2. * y_j + 2. * y_l) * v26 / 2.)
                        * v17
                        * v26
                        / (v20 * (v21 + v28 + v27).powf(1.5)))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    * ((x_j - x_l) * v31 / v32
                        + (-y_j + y_l) * v33 / v32
                        + v13
                            * (-(2. * x_j - 2. * x_l) * v29 / 2.
                                - (-2. * y_j + 2. * y_l) * v26 / 2.)
                            * v25
                            / (v21 + v28 + v27).powf(1.5)
                        + v13
                            * (-(-x_i + x_j) * v0 / v22
                                - v0 * v24 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v0 * (-(-2. * x_i + 2. * x_j) * v6 / 2.
                                    - (2. * y_i - 2. * y_j) * v17 / 2.)
                                    * v6
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v14 * v8 / v22
                                + v8 * v24 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v8 * (-(-2. * x_i + 2. * x_j) * v6 / 2.
                                    - (2. * y_i - 2. * y_j) * v17 / 2.)
                                    * v17
                                    / (v15 * (v19 + v16 + v18).powf(1.5)))
                            / v32
                        + v29
                            * (-(2. * x_j - 2. * x_l) * v29 / 2.
                                - (-2. * y_j + 2. * y_l) * v26 / 2.)
                            * v31
                            / (v21 + v28 + v27).powf(1.5)
                        + v29
                            * (v0 * v24 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v0 * (-(-2. * x_i + 2. * x_j) * v6 / 2.
                                    - (2. * y_i - 2. * y_j) * v17 / 2.)
                                    * v30
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - v14 * v3 / v22
                                - v3 * v24 * v17 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v3 * (-(-2. * x_i + 2. * x_j) * v6 / 2.
                                    - (2. * y_i - 2. * y_j) * v17 / 2.)
                                    * v17
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                - v17 / v22)
                            / v32
                        + (-(2. * x_j - 2. * x_l) * v29 / 2. - (-2. * y_j + 2. * y_l) * v26 / 2.)
                            * v26
                            * v33
                            / (v21 + v28 + v27).powf(1.5)
                        + v26
                            * ((-x_i + x_j) * v3 / v22
                                - v8 * v24 * v30 / ((v7 + v9 + v10).powf(1.5) * v20)
                                - v8 * (-(-2. * x_i + 2. * x_j) * v6 / 2.
                                    - (2. * y_i - 2. * y_j) * v17 / 2.)
                                    * v30
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v3 * v24 * v6 / ((v7 + v9 + v10).powf(1.5) * v20)
                                + v3 * (-(-2. * x_i + 2. * x_j) * v6 / 2.
                                    - (2. * y_i - 2. * y_j) * v17 / 2.)
                                    * v6
                                    / (v15 * (v19 + v16 + v18).powf(1.5))
                                + v6 / v22)
                            / v32)
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
        gradient[self.l].x += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * (v8 * v30 / v34
                    + v24 * v6 / v34
                    + v30
                        * v13
                        * (-(-2. * y_j + 2. * y_k) * v13 / 2. - (2. * z_j - 2. * z_k) * v29 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + v6 * v29
                        * (-(-2. * y_j + 2. * y_k) * v13 / 2. - (2. * z_j - 2. * z_k) * v29 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + v17
                        * (-(-2. * y_j + 2. * y_k) * v13 / 2. - (2. * z_j - 2. * z_k) * v29 / 2.)
                        * v26
                        / (v20 * (v21 + v28 + v27).powf(1.5)))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    * (v8 * v25 / v32
                        + v24 * v31 / v32
                        + v13
                            * (-(-2. * y_j + 2. * y_k) * v13 / 2.
                                - (2. * z_j - 2. * z_k) * v29 / 2.)
                            * v25
                            / (v21 + v28 + v27).powf(1.5)
                        + v29
                            * (-(-2. * y_j + 2. * y_k) * v13 / 2.
                                - (2. * z_j - 2. * z_k) * v29 / 2.)
                            * v31
                            / (v21 + v28 + v27).powf(1.5)
                        + (-(-2. * y_j + 2. * y_k) * v13 / 2. - (2. * z_j - 2. * z_k) * v29 / 2.)
                            * v26
                            * v33
                            / (v21 + v28 + v27).powf(1.5))
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
        gradient[self.l].y += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * (v12 * v30 / v34
                    + v3 * v17 / v34
                    + v30
                        * v13
                        * (-(2. * x_j - 2. * x_k) * v13 / 2. - (-2. * z_j + 2. * z_k) * v26 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + v6 * v29
                        * (-(2. * x_j - 2. * x_k) * v13 / 2. - (-2. * z_j + 2. * z_k) * v26 / 2.)
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + (-(2. * x_j - 2. * x_k) * v13 / 2. - (-2. * z_j + 2. * z_k) * v26 / 2.)
                        * v17
                        * v26
                        / (v20 * (v21 + v28 + v27).powf(1.5)))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    * (v12 * v25 / v32
                        + v3 * v33 / v32
                        + v13
                            * (-(2. * x_j - 2. * x_k) * v13 / 2.
                                - (-2. * z_j + 2. * z_k) * v26 / 2.)
                            * v25
                            / (v21 + v28 + v27).powf(1.5)
                        + v29
                            * (-(2. * x_j - 2. * x_k) * v13 / 2.
                                - (-2. * z_j + 2. * z_k) * v26 / 2.)
                            * v31
                            / (v21 + v28 + v27).powf(1.5)
                        + (-(2. * x_j - 2. * x_k) * v13 / 2. - (-2. * z_j + 2. * z_k) * v26 / 2.)
                            * v26
                            * v33
                            / (v21 + v28 + v27).powf(1.5))
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
        gradient[self.l].z += 0.5
            * self.n_phi
            * self.v_phi
            * (v35
                * (v0 * v6 / v34
                    + v1 * v17 / v34
                    + v30
                        * (-(-2. * x_j + 2. * x_k) * v29 / 2. - (2. * y_j - 2. * y_k) * v26 / 2.)
                        * v13
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + v6 * (-(-2. * x_j + 2. * x_k) * v29 / 2.
                        - (2. * y_j - 2. * y_k) * v26 / 2.)
                        * v29
                        / (v20 * (v21 + v28 + v27).powf(1.5))
                    + (-(-2. * x_j + 2. * x_k) * v29 / 2. - (2. * y_j - 2. * y_k) * v26 / 2.)
                        * v17
                        * v26
                        / (v20 * (v21 + v28 + v27).powf(1.5)))
                / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                    + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2))
                + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34)
                    * (v0 * v31 / v32
                        + v1 * v33 / v32
                        + (-(-2. * x_j + 2. * x_k) * v29 / 2. - (2. * y_j - 2. * y_k) * v26 / 2.)
                            * v13
                            * v25
                            / (v21 + v28 + v27).powf(1.5)
                        + (-(-2. * x_j + 2. * x_k) * v29 / 2. - (2. * y_j - 2. * y_k) * v26 / 2.)
                            * v29
                            * v31
                            / (v21 + v28 + v27).powf(1.5)
                        + (-(-2. * x_j + 2. * x_k) * v29 / 2. - (2. * y_j - 2. * y_k) * v26 / 2.)
                            * v26
                            * v33
                            / (v21 + v28 + v27).powf(1.5))
                    / ((v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32).powi(2)
                        + (v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34).powi(2)))
            * (self.n_phi
                * (v13 * v25 / v32 + v29 * v31 / v32 + v26 * v33 / v32)
                    .atan2(v30 * v13 / v34 + v6 * v29 / v34 + v17 * v26 / v34))
            .sin()
            * (self.n_phi * self.phi0).cos();
    }
}

/// Value of the torsional dihedral
fn phi(dihedral: &TorsionalDihedral, coordinates: &[Point]) -> f64 {
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

/// Improper dihedral defining inversion about a centre
/// ```text
///                    j
///                  /
///        i  ---  c    ω
///                 \
///                  \
///                    k
/// ```
///
/// where the angle is between the c-k axis and the plane
/// formed by atoms c,i,j. It may also be defined by the
/// c,i,j normal and the axis: ω = γ - π
pub struct InversionDihedral {
    pub(crate) c: usize,
    pub(crate) i: usize,
    pub(crate) j: usize,
    pub(crate) k: usize,

    pub(crate) c0: f64,
    pub(crate) c1: f64,
    pub(crate) c2: f64,
    pub(crate) k_cijk: f64,
}

impl InversionDihedral {
    fn e_gamma(&self, gamma: f64) -> f64 {
        self.k_cijk * (self.c0 + self.c1 * gamma.sin() + self.c2 * (2. * gamma).cos())
    }
}

impl EnergyFunction for InversionDihedral {
    fn involves_idxs(&self, idxs: Vec<usize>) -> bool {
        idxs.len() == 4
            && (ImproperDihedral {
                c: self.c,
                i: self.i,
                j: self.j,
                k: self.k,
            } == ImproperDihedral {
                c: idxs[0],
                i: idxs[1],
                j: idxs[2],
                k: idxs[3],
            })
    }

    fn force_constant(&self) -> f64 {
        self.k_cijk
    }

    fn energy(&self, coordinates: &[Point]) -> f64 {
        (self.e_gamma(gamma(self.c, self.i, self.j, self.k, coordinates))
            + self.e_gamma(gamma(self.c, self.k, self.i, self.j, coordinates))
            + self.e_gamma(gamma(self.c, self.j, self.k, self.i, coordinates)))
            / 3.
    }

    fn add_gradient(&self, coordinates: &[Point], gradient: &mut Vec<Vector3D>) {
        let x_i = coordinates[self.i].x;
        let y_i = coordinates[self.i].y;
        let z_i = coordinates[self.i].z;

        let x_j = coordinates[self.j].x;
        let y_j = coordinates[self.j].y;
        let z_j = coordinates[self.j].z;

        let x_k = coordinates[self.k].x;
        let y_k = coordinates[self.k].y;
        let z_k = coordinates[self.k].z;

        let x_c = coordinates[self.c].x;
        let y_c = coordinates[self.c].y;
        let z_c = coordinates[self.c].z;

        let v0 = (-x_c + x_i).powi(2);
        let v1 = (-y_c + y_i).powi(2);
        let v2 = -z_c + z_i;
        let v3 = v2.powi(2);
        let v4 = -x_c + x_j;
        let v5 = -x_c + x_k;
        let v6 = -y_c + y_j;
        let v7 = -z_c + z_j;
        let v8 = -z_c + z_k;
        let v9 = -y_c + y_k;
        let v10 = v0 + v1 + v3;
        let v11 = v4 * v9 - v5 * v6;
        let v12 = -x_c + x_i;
        let v13 = -y_c + y_i;
        let v14 = -v4 * v8 + v5 * v7;
        let v15 = v14.powi(2);
        let v16 = v11.powi(2);
        let v17 = v6 * v8 - v9 * v7;
        let v18 = v17.powi(2);
        let v19 = v16 + v15 + v18;
        let v20 = v12 * v17 + v13 * v14 + v2 * v11;
        let v21 = v10 * v19;
        let v22 = v20.powi(2);
        let v23 = 2. * v10 * v19;
        let v24 = v19.sqrt();
        let v25 = (-v22 / v21 + 1.).sqrt();
        gradient[self.c].x += self.k_cijk
            * (self.c1
                * (-(-2. * x_c + 2. * x_i) * v22 / (2. * v10.powi(2) * v19)
                    - (-(2. * y_j - 2. * y_k) * v11 - (-2. * z_j + 2. * z_k) * v14) * v22
                        / (2. * v10 * v19.powi(2))
                    - v20
                        * (2. * v13 * (-z_j + z_k) - 2. * v6 * v8
                            + 2. * v9 * v7
                            + 2. * (y_j - y_k) * v2)
                        / v23)
                / v25
                + 2. * self.c2
                    * (v12 * v20 / (v10.powf(1.5) * v24)
                        + (-(2. * y_j - 2. * y_k) * v11 / 2. - (-2. * z_j + 2. * z_k) * v14 / 2.)
                            * v20
                            / (v10.sqrt() * v19.powf(1.5))
                        + (v13 * (-z_j + z_k) - v6 * v8 + v9 * v7 + (y_j - y_k) * v2)
                            / (v10.sqrt() * v24))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-2. * x_c + 2. * x_j)
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2)).powi(2)
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        - (-(-2. * y_i + 2. * y_k) * (-v12 * v9 + v5 * v13)
                            - (2. * z_i - 2. * z_k) * (v12 * v8 - v5 * v2))
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                .powi(2)
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .powi(2))
                        - (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            * (2. * v13 * v8 + 2. * v6 * (z_i - z_k) - 2. * v9 * v2
                                + 2. * (-y_i + y_k) * v7)
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * (v4
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).powf(1.5)
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .sqrt())
                            + (-(-2. * y_i + 2. * y_k) * (-v12 * v9 + v5 * v13) / 2.
                                - (2. * z_i - 2. * z_k) * (v12 * v8 - v5 * v2) / 2.)
                                * (v4 * (-v13 * v8 + v9 * v2)
                                    + v6 * (v12 * v8 - v5 * v2)
                                    + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .powf(1.5))
                            + (v13 * v8 + v6 * (z_i - z_k) - v9 * v2 + (-y_i + y_k) * v7)
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-2. * x_c + 2. * x_k)
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2)).powi(2)
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        - (-(2. * y_i - 2. * y_j) * (v12 * v6 - v4 * v13)
                            - (-2. * z_i + 2. * z_j) * (-v12 * v7 + v4 * v2))
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                .powi(2)
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .powi(2))
                        - (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            * (-2. * v13 * v7
                                + 2. * v6 * v2
                                + 2. * v9 * (-z_i + z_j)
                                + 2. * (y_i - y_j) * v8)
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * (v5
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).powf(1.5)
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .sqrt())
                            + (-(2. * y_i - 2. * y_j) * (v12 * v6 - v4 * v13) / 2.
                                - (-2. * z_i + 2. * z_j) * (-v12 * v7 + v4 * v2) / 2.)
                                * (v5 * (v13 * v7 - v6 * v2)
                                    + v9 * (-v12 * v7 + v4 * v2)
                                    + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .powf(1.5))
                            + (-v13 * v7 + v6 * v2 + v9 * (-z_i + z_j) + (y_i - y_j) * v8)
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
        gradient[self.c].y += self.k_cijk
            * (self.c1
                * (-(-2. * y_c + 2. * y_i) * v22 / (2. * v10.powi(2) * v19)
                    - (-(-2. * x_j + 2. * x_k) * v11 - (2. * z_j - 2. * z_k) * v17) * v22
                        / (2. * v10 * v19.powi(2))
                    - v20
                        * (2. * v12 * (z_j - z_k) + 2. * v4 * v8 - 2. * v5 * v7
                            + 2. * (-x_j + x_k) * v2)
                        / v23)
                / v25
                + 2. * self.c2
                    * (v13 * v20 / (v10.powf(1.5) * v24)
                        + (-(-2. * x_j + 2. * x_k) * v11 / 2. - (2. * z_j - 2. * z_k) * v17 / 2.)
                            * v20
                            / (v10.sqrt() * v19.powf(1.5))
                        + (v12 * (z_j - z_k) + v4 * v8 - v5 * v7 + (-x_j + x_k) * v2)
                            / (v10.sqrt() * v24))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-2. * y_c + 2. * y_j)
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2)).powi(2)
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        - (-(2. * x_i - 2. * x_k) * (-v12 * v9 + v5 * v13)
                            - (-2. * z_i + 2. * z_k) * (-v13 * v8 + v9 * v2))
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                .powi(2)
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .powi(2))
                        - (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            * (-2. * v12 * v8
                                + 2. * v4 * (-z_i + z_k)
                                + 2. * v5 * v2
                                + 2. * (x_i - x_k) * v7)
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * (v6
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).powf(1.5)
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .sqrt())
                            + (-(2. * x_i - 2. * x_k) * (-v12 * v9 + v5 * v13) / 2.
                                - (-2. * z_i + 2. * z_k) * (-v13 * v8 + v9 * v2) / 2.)
                                * (v4 * (-v13 * v8 + v9 * v2)
                                    + v6 * (v12 * v8 - v5 * v2)
                                    + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .powf(1.5))
                            + (-v12 * v8 + v4 * (-z_i + z_k) + v5 * v2 + (x_i - x_k) * v7)
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-2. * y_c + 2. * y_k)
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2)).powi(2)
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        - (-(-2. * x_i + 2. * x_j) * (v12 * v6 - v4 * v13)
                            - (2. * z_i - 2. * z_j) * (v13 * v7 - v6 * v2))
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                .powi(2)
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .powi(2))
                        - (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            * (2. * v12 * v7 - 2. * v4 * v2
                                + 2. * v5 * (z_i - z_j)
                                + 2. * (-x_i + x_j) * v8)
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * (v9
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).powf(1.5)
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .sqrt())
                            + (-(-2. * x_i + 2. * x_j) * (v12 * v6 - v4 * v13) / 2.
                                - (2. * z_i - 2. * z_j) * (v13 * v7 - v6 * v2) / 2.)
                                * (v5 * (v13 * v7 - v6 * v2)
                                    + v9 * (-v12 * v7 + v4 * v2)
                                    + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .powf(1.5))
                            + (v12 * v7 - v4 * v2 + v5 * (z_i - z_j) + (-x_i + x_j) * v8)
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
        gradient[self.c].z += self.k_cijk
            * (self.c1
                * (-(-2. * z_c + 2. * z_i) * v22 / (2. * v10.powi(2) * v19)
                    - (-(2. * x_j - 2. * x_k) * v14 - (-2. * y_j + 2. * y_k) * v17) * v22
                        / (2. * v10 * v19.powi(2))
                    - v20
                        * (2. * v12 * (-y_j + y_k) - 2. * v4 * v9
                            + 2. * v5 * v6
                            + 2. * (x_j - x_k) * v13)
                        / v23)
                / v25
                + 2. * self.c2
                    * (v2 * v20 / (v10.powf(1.5) * v24)
                        + (-(2. * x_j - 2. * x_k) * v14 / 2. - (-2. * y_j + 2. * y_k) * v17 / 2.)
                            * v20
                            / (v10.sqrt() * v19.powf(1.5))
                        + (v12 * (-y_j + y_k) - v4 * v9 + v5 * v6 + (x_j - x_k) * v13)
                            / (v10.sqrt() * v24))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-2. * z_c + 2. * z_j)
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2)).powi(2)
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        - (-(-2. * x_i + 2. * x_k) * (v12 * v8 - v5 * v2)
                            - (2. * y_i - 2. * y_k) * (-v13 * v8 + v9 * v2))
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                .powi(2)
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .powi(2))
                        - (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            * (2. * v12 * v9 + 2. * v4 * (y_i - y_k) - 2. * v5 * v13
                                + 2. * (-x_i + x_k) * v6)
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * (v7
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).powf(1.5)
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .sqrt())
                            + (-(-2. * x_i + 2. * x_k) * (v12 * v8 - v5 * v2) / 2.
                                - (2. * y_i - 2. * y_k) * (-v13 * v8 + v9 * v2) / 2.)
                                * (v4 * (-v13 * v8 + v9 * v2)
                                    + v6 * (v12 * v8 - v5 * v2)
                                    + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .powf(1.5))
                            + (v12 * v9 + v4 * (y_i - y_k) - v5 * v13 + (-x_i + x_k) * v6)
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-2. * z_c + 2. * z_k)
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2)).powi(2)
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        - (-(2. * x_i - 2. * x_j) * (-v12 * v7 + v4 * v2)
                            - (-2. * y_i + 2. * y_j) * (v13 * v7 - v6 * v2))
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                .powi(2)
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .powi(2))
                        - (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            * (-2. * v12 * v6
                                + 2. * v4 * v13
                                + 2. * v5 * (-y_i + y_j)
                                + 2. * (x_i - x_j) * v9)
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * (v8
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).powf(1.5)
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .sqrt())
                            + (-(2. * x_i - 2. * x_j) * (-v12 * v7 + v4 * v2) / 2.
                                - (-2. * y_i + 2. * y_j) * (v13 * v7 - v6 * v2) / 2.)
                                * (v5 * (v13 * v7 - v6 * v2)
                                    + v9 * (-v12 * v7 + v4 * v2)
                                    + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .powf(1.5))
                            + (-v12 * v6 + v4 * v13 + v5 * (-y_i + y_j) + (x_i - x_j) * v9)
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
        gradient[self.i].x += self.k_cijk
            * (self.c1
                * (-(2. * x_c - 2. * x_i) * v22 / (2. * v10.powi(2) * v19)
                    - (2. * v6 * v8 - 2. * v9 * v7) * v20 / v23)
                / v25
                + 2. * self.c2
                    * ((x_c - x_i) * v20 / (v10.powf(1.5) * v24) + v17 / (v10.sqrt() * v24))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * v6 * v8 + 2. * (y_c - y_k) * v7)
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        - (-(2. * y_c - 2. * y_k) * (-v12 * v9 + v5 * v13)
                            - (-2. * z_c + 2. * z_k) * (v12 * v8 - v5 * v2))
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                .powi(2)
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .powi(2)))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((v6 * v8 + (y_c - y_k) * v7)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .sqrt())
                            + (-(2. * y_c - 2. * y_k) * (-v12 * v9 + v5 * v13) / 2.
                                - (-2. * z_c + 2. * z_k) * (v12 * v8 - v5 * v2) / 2.)
                                * (v4 * (-v13 * v8 + v9 * v2)
                                    + v6 * (v12 * v8 - v5 * v2)
                                    + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .powf(1.5)))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-(-2. * y_c + 2. * y_j) * (v12 * v6 - v4 * v13)
                        - (2. * z_c - 2. * z_j) * (-v12 * v7 + v4 * v2))
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2))
                            .powi(2))
                        - (2. * v6 * v8 + 2. * v9 * (z_c - z_j))
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((-(-2. * y_c + 2. * y_j) * (v12 * v6 - v4 * v13) / 2.
                            - (2. * z_c - 2. * z_j) * (-v12 * v7 + v4 * v2) / 2.)
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .powf(1.5))
                            + (v6 * v8 + v9 * (z_c - z_j))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
        gradient[self.i].y += self.k_cijk
            * (self.c1
                * (-(2. * y_c - 2. * y_i) * v22 / (2. * v10.powi(2) * v19)
                    - (-2. * v4 * v8 + 2. * v5 * v7) * v20 / v23)
                / v25
                + 2. * self.c2
                    * ((y_c - y_i) * v20 / (v10.powf(1.5) * v24) + v14 / (v10.sqrt() * v24))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-(-2. * x_c + 2. * x_k) * (-v12 * v9 + v5 * v13)
                        - (2. * z_c - 2. * z_k) * (-v13 * v8 + v9 * v2))
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2))
                            .powi(2))
                        - (2. * v4 * (z_c - z_k) + 2. * v5 * v7)
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((-(-2. * x_c + 2. * x_k) * (-v12 * v9 + v5 * v13) / 2.
                            - (2. * z_c - 2. * z_k) * (-v13 * v8 + v9 * v2) / 2.)
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .powf(1.5))
                            + (v4 * (z_c - z_k) + v5 * v7)
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * v5 * v7 + 2. * (x_c - x_j) * v8)
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        - (-(2. * x_c - 2. * x_j) * (v12 * v6 - v4 * v13)
                            - (-2. * z_c + 2. * z_j) * (v13 * v7 - v6 * v2))
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                .powi(2)
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .powi(2)))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((v5 * v7 + (x_c - x_j) * v8)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .sqrt())
                            + (-(2. * x_c - 2. * x_j) * (v12 * v6 - v4 * v13) / 2.
                                - (-2. * z_c + 2. * z_j) * (v13 * v7 - v6 * v2) / 2.)
                                * (v5 * (v13 * v7 - v6 * v2)
                                    + v9 * (-v12 * v7 + v4 * v2)
                                    + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .powf(1.5)))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
        gradient[self.i].z += self.k_cijk
            * (self.c1
                * (-(2. * z_c - 2. * z_i) * v22 / (2. * v10.powi(2) * v19)
                    - (2. * v4 * v9 - 2. * v5 * v6) * v20 / v23)
                / v25
                + 2. * self.c2
                    * ((z_c - z_i) * v20 / (v10.powf(1.5) * v24) + v11 / (v10.sqrt() * v24))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * v4 * v9 + 2. * (x_c - x_k) * v6)
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        - (-(2. * x_c - 2. * x_k) * (v12 * v8 - v5 * v2)
                            - (-2. * y_c + 2. * y_k) * (-v13 * v8 + v9 * v2))
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                .powi(2)
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .powi(2)))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((v4 * v9 + (x_c - x_k) * v6)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .sqrt())
                            + (-(2. * x_c - 2. * x_k) * (v12 * v8 - v5 * v2) / 2.
                                - (-2. * y_c + 2. * y_k) * (-v13 * v8 + v9 * v2) / 2.)
                                * (v4 * (-v13 * v8 + v9 * v2)
                                    + v6 * (v12 * v8 - v5 * v2)
                                    + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .powf(1.5)))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-(-2. * x_c + 2. * x_j) * (-v12 * v7 + v4 * v2)
                        - (2. * y_c - 2. * y_j) * (v13 * v7 - v6 * v2))
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2))
                            .powi(2))
                        - (2. * v4 * v9 + 2. * v5 * (y_c - y_j))
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((-(-2. * x_c + 2. * x_j) * (-v12 * v7 + v4 * v2) / 2.
                            - (2. * y_c - 2. * y_j) * (v13 * v7 - v6 * v2) / 2.)
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .powf(1.5))
                            + (v4 * v9 + v5 * (y_c - y_j))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
        gradient[self.j].x += self.k_cijk
            * (self.c1
                * (-(-(-2. * y_c + 2. * y_k) * v11 - (2. * z_c - 2. * z_k) * v14) * v22
                    / (2. * v10 * v19.powi(2))
                    - (2. * v13 * (z_c - z_k) + 2. * v9 * v2) * v20 / v23)
                / v25
                + 2. * self.c2
                    * ((-(-2. * y_c + 2. * y_k) * v11 / 2. - (2. * z_c - 2. * z_k) * v14 / 2.)
                        * v20
                        / (v10.sqrt() * v19.powf(1.5))
                        + (v13 * (z_c - z_k) + v9 * v2) / (v10.sqrt() * v24))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * x_c - 2. * x_j)
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2)).powi(2)
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        - (-2. * v13 * v8 + 2. * v9 * v2)
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((x_c - x_j)
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).powf(1.5)
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .sqrt())
                            + (-v13 * v8 + v9 * v2)
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * v9 * v2 + 2. * (y_c - y_i) * v8)
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        - (-(2. * y_c - 2. * y_i) * (v12 * v6 - v4 * v13)
                            - (-2. * z_c + 2. * z_i) * (-v12 * v7 + v4 * v2))
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                .powi(2)
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .powi(2)))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((v9 * v2 + (y_c - y_i) * v8)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .sqrt())
                            + (-(2. * y_c - 2. * y_i) * (v12 * v6 - v4 * v13) / 2.
                                - (-2. * z_c + 2. * z_i) * (-v12 * v7 + v4 * v2) / 2.)
                                * (v5 * (v13 * v7 - v6 * v2)
                                    + v9 * (-v12 * v7 + v4 * v2)
                                    + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .powf(1.5)))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
        gradient[self.j].y += self.k_cijk
            * (self.c1
                * (-(2. * v12 * v8 + 2. * (x_c - x_k) * v2) * v20 / v23
                    - (-(2. * x_c - 2. * x_k) * v11 - (-2. * z_c + 2. * z_k) * v17) * v22
                        / (2. * v10 * v19.powi(2)))
                / v25
                + 2. * self.c2
                    * ((v12 * v8 + (x_c - x_k) * v2) / (v10.sqrt() * v24)
                        + (-(2. * x_c - 2. * x_k) * v11 / 2. - (-2. * z_c + 2. * z_k) * v17 / 2.)
                            * v20
                            / (v10.sqrt() * v19.powf(1.5)))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * y_c - 2. * y_j)
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2)).powi(2)
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        - (2. * v12 * v8 - 2. * v5 * v2)
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((y_c - y_j)
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).powf(1.5)
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .sqrt())
                            + (v12 * v8 - v5 * v2)
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-(-2. * x_c + 2. * x_i) * (v12 * v6 - v4 * v13)
                        - (2. * z_c - 2. * z_i) * (v13 * v7 - v6 * v2))
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2))
                            .powi(2))
                        - (2. * v12 * v8 + 2. * v5 * (z_c - z_i))
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((-(-2. * x_c + 2. * x_i) * (v12 * v6 - v4 * v13) / 2.
                            - (2. * z_c - 2. * z_i) * (v13 * v7 - v6 * v2) / 2.)
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .powf(1.5))
                            + (v12 * v8 + v5 * (z_c - z_i))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
        gradient[self.j].z += self.k_cijk
            * (self.c1
                * (-(-(-2. * x_c + 2. * x_k) * v14 - (2. * y_c - 2. * y_k) * v17) * v22
                    / (2. * v10 * v19.powi(2))
                    - (2. * v12 * (y_c - y_k) + 2. * v5 * v13) * v20 / v23)
                / v25
                + 2. * self.c2
                    * ((-(-2. * x_c + 2. * x_k) * v14 / 2. - (2. * y_c - 2. * y_k) * v17 / 2.)
                        * v20
                        / (v10.sqrt() * v19.powf(1.5))
                        + (v12 * (y_c - y_k) + v5 * v13) / (v10.sqrt() * v24))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * z_c - 2. * z_j)
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2)).powi(2)
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        - (-2. * v12 * v9 + 2. * v5 * v13)
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((z_c - z_j)
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).powf(1.5)
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .sqrt())
                            + (-v12 * v9 + v5 * v13)
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * v5 * v13 + 2. * (x_c - x_i) * v9)
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        - (-(2. * x_c - 2. * x_i) * (-v12 * v7 + v4 * v2)
                            - (-2. * y_c + 2. * y_i) * (v13 * v7 - v6 * v2))
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                .powi(2)
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .powi(2)))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((v5 * v13 + (x_c - x_i) * v9)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .sqrt())
                            + (-(2. * x_c - 2. * x_i) * (-v12 * v7 + v4 * v2) / 2.
                                - (-2. * y_c + 2. * y_i) * (v13 * v7 - v6 * v2) / 2.)
                                * (v5 * (v13 * v7 - v6 * v2)
                                    + v9 * (-v12 * v7 + v4 * v2)
                                    + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .powf(1.5)))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
        gradient[self.k].x += self.k_cijk
            * (self.c1
                * (-(2. * v13 * v7 + 2. * (y_c - y_j) * v2) * v20 / v23
                    - (-(2. * y_c - 2. * y_j) * v11 - (-2. * z_c + 2. * z_j) * v14) * v22
                        / (2. * v10 * v19.powi(2)))
                / v25
                + 2. * self.c2
                    * ((v13 * v7 + (y_c - y_j) * v2) / (v10.sqrt() * v24)
                        + (-(2. * y_c - 2. * y_j) * v11 / 2. - (-2. * z_c + 2. * z_j) * v14 / 2.)
                            * v20
                            / (v10.sqrt() * v19.powf(1.5)))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-(-2. * y_c + 2. * y_i) * (-v12 * v9 + v5 * v13)
                        - (2. * z_c - 2. * z_i) * (v12 * v8 - v5 * v2))
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2))
                            .powi(2))
                        - (2. * v13 * v7 + 2. * v6 * (z_c - z_i))
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((-(-2. * y_c + 2. * y_i) * (-v12 * v9 + v5 * v13) / 2.
                            - (2. * z_c - 2. * z_i) * (v12 * v8 - v5 * v2) / 2.)
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .powf(1.5))
                            + (v13 * v7 + v6 * (z_c - z_i))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * x_c - 2. * x_k)
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2)).powi(2)
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        - (2. * v13 * v7 - 2. * v6 * v2)
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((x_c - x_k)
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).powf(1.5)
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .sqrt())
                            + (v13 * v7 - v6 * v2)
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
        gradient[self.k].y += self.k_cijk
            * (self.c1
                * (-(-(-2. * x_c + 2. * x_j) * v11 - (2. * z_c - 2. * z_j) * v17) * v22
                    / (2. * v10 * v19.powi(2))
                    - (2. * v12 * (z_c - z_j) + 2. * v4 * v2) * v20 / v23)
                / v25
                + 2. * self.c2
                    * ((-(-2. * x_c + 2. * x_j) * v11 / 2. - (2. * z_c - 2. * z_j) * v17 / 2.)
                        * v20
                        / (v10.sqrt() * v19.powf(1.5))
                        + (v12 * (z_c - z_j) + v4 * v2) / (v10.sqrt() * v24))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * v4 * v2 + 2. * (x_c - x_i) * v7)
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        - (-(2. * x_c - 2. * x_i) * (-v12 * v9 + v5 * v13)
                            - (-2. * z_c + 2. * z_i) * (-v13 * v8 + v9 * v2))
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                .powi(2)
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .powi(2)))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((v4 * v2 + (x_c - x_i) * v7)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .sqrt())
                            + (-(2. * x_c - 2. * x_i) * (-v12 * v9 + v5 * v13) / 2.
                                - (-2. * z_c + 2. * z_i) * (-v13 * v8 + v9 * v2) / 2.)
                                * (v4 * (-v13 * v8 + v9 * v2)
                                    + v6 * (v12 * v8 - v5 * v2)
                                    + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .powf(1.5)))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * y_c - 2. * y_k)
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2)).powi(2)
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        - (-2. * v12 * v7 + 2. * v4 * v2)
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((y_c - y_k)
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).powf(1.5)
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .sqrt())
                            + (-v12 * v7 + v4 * v2)
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
        gradient[self.k].z += self.k_cijk
            * (self.c1
                * (-(2. * v12 * v6 + 2. * (x_c - x_j) * v13) * v20 / v23
                    - (-(2. * x_c - 2. * x_j) * v14 - (-2. * y_c + 2. * y_j) * v17) * v22
                        / (2. * v10 * v19.powi(2)))
                / v25
                + 2. * self.c2
                    * ((v12 * v6 + (x_c - x_j) * v13) / (v10.sqrt() * v24)
                        + (-(2. * x_c - 2. * x_j) * v14 / 2. - (-2. * y_c + 2. * y_j) * v17 / 2.)
                            * v20
                            / (v10.sqrt() * v19.powf(1.5)))
                    * (2. * (v20 / (v10.sqrt() * v24)).acos()).sin()
                    / v25)
            / 3.
            + self.k_cijk
                * (self.c1
                    * (-(-(-2. * x_c + 2. * x_i) * (v12 * v8 - v5 * v2)
                        - (2. * y_c - 2. * y_i) * (-v13 * v8 + v9 * v2))
                        * (v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                        / (2.
                            * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2))
                            .powi(2))
                        - (2. * v12 * v6 + 2. * v4 * (y_c - y_i))
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / (2.
                                * (v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))))
                    / (-(v4 * (-v13 * v8 + v9 * v2)
                        + v6 * (v12 * v8 - v5 * v2)
                        + v7 * (-v12 * v9 + v5 * v13))
                        .powi(2)
                        / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                            * ((-v12 * v9 + v5 * v13).powi(2)
                                + (v12 * v8 - v5 * v2).powi(2)
                                + (-v13 * v8 + v9 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((-(-2. * x_c + 2. * x_i) * (v12 * v8 - v5 * v2) / 2.
                            - (2. * y_c - 2. * y_i) * (-v13 * v8 + v9 * v2) / 2.)
                            * (v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2))
                                .powf(1.5))
                            + (v12 * v6 + v4 * (y_c - y_i))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v4 * (-v13 * v8 + v9 * v2)
                                + v6 * (v12 * v8 - v5 * v2)
                                + v7 * (-v12 * v9 + v5 * v13))
                                / ((v4.powi(2) + v6.powi(2) + v7.powi(2)).sqrt()
                                    * ((-v12 * v9 + v5 * v13).powi(2)
                                        + (v12 * v8 - v5 * v2).powi(2)
                                        + (-v13 * v8 + v9 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v4 * (-v13 * v8 + v9 * v2)
                            + v6 * (v12 * v8 - v5 * v2)
                            + v7 * (-v12 * v9 + v5 * v13))
                            .powi(2)
                            / ((v4.powi(2) + v6.powi(2) + v7.powi(2))
                                * ((-v12 * v9 + v5 * v13).powi(2)
                                    + (v12 * v8 - v5 * v2).powi(2)
                                    + (-v13 * v8 + v9 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.
            + self.k_cijk
                * (self.c1
                    * (-(2. * z_c - 2. * z_k)
                        * (v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                        / (2.
                            * (v5.powi(2) + v9.powi(2) + v8.powi(2)).powi(2)
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        - (2. * v12 * v6 - 2. * v4 * v13)
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / (2.
                                * (v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))))
                    / (-(v5 * (v13 * v7 - v6 * v2)
                        + v9 * (-v12 * v7 + v4 * v2)
                        + v8 * (v12 * v6 - v4 * v13))
                        .powi(2)
                        / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                            * ((v12 * v6 - v4 * v13).powi(2)
                                + (-v12 * v7 + v4 * v2).powi(2)
                                + (v13 * v7 - v6 * v2).powi(2)))
                        + 1.)
                        .sqrt()
                    + 2. * self.c2
                        * ((z_c - z_k)
                            * (v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).powf(1.5)
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2))
                                .sqrt())
                            + (v12 * v6 - v4 * v13)
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                        * (2.
                            * ((v5 * (v13 * v7 - v6 * v2)
                                + v9 * (-v12 * v7 + v4 * v2)
                                + v8 * (v12 * v6 - v4 * v13))
                                / ((v5.powi(2) + v9.powi(2) + v8.powi(2)).sqrt()
                                    * ((v12 * v6 - v4 * v13).powi(2)
                                        + (-v12 * v7 + v4 * v2).powi(2)
                                        + (v13 * v7 - v6 * v2).powi(2))
                                    .sqrt()))
                            .acos())
                        .sin()
                        / (-(v5 * (v13 * v7 - v6 * v2)
                            + v9 * (-v12 * v7 + v4 * v2)
                            + v8 * (v12 * v6 - v4 * v13))
                            .powi(2)
                            / ((v5.powi(2) + v9.powi(2) + v8.powi(2))
                                * ((v12 * v6 - v4 * v13).powi(2)
                                    + (-v12 * v7 + v4 * v2).powi(2)
                                    + (v13 * v7 - v6 * v2).powi(2)))
                            + 1.)
                            .sqrt())
                / 3.;
    }
}

/// Inversion angle γ
#[inline(always)]
fn gamma(c: usize, i: usize, j: usize, k: usize, coordinates: &[Point]) -> f64 {
    let v0 = &coordinates[i] - &coordinates[c];
    let v1 = &coordinates[j] - &coordinates[c];

    let v2 = v0.cross(&v1);
    let v3 = &coordinates[k] - &coordinates[c];

    (v2.dot(&v3) / (v2.length() * v3.length())).acos()
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
mod tests {

    use super::*;
    use crate::utils::*;

    /// Ensure a torsion can be calculated correctly. Compared to an Avogadro reference
    #[test]
    fn test_torsional_dihedral_calculation() {
        // Coordinates for a H2O2 molecule
        let x: Vec<Point> = Vec::from([
            Point {
                x: -3.80272,
                y: 0.11331,
                z: -0.29090,
            },
            Point {
                x: -3.75673,
                y: 1.06948,
                z: -0.03303,
            },
            Point {
                x: -2.46930,
                y: 1.33955,
                z: 0.03004,
            },
            Point {
                x: -2.25240,
                y: 1.31540,
                z: 0.99712,
            },
        ]);
        let dihedral = TorsionalDihedral {
            i: 0,
            j: 1,
            k: 2,
            l: 3,
            n_phi: 0.,
            phi0: 0.,
            v_phi: 0.,
        };
        assert!(is_very_close(dihedral.force_constant(), 0.0));

        let expected_phi = -1.7593; // -100.8 degrees

        assert!(is_close(phi(&dihedral, &x), expected_phi, 1E-3));
    }
}
