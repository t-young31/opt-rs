use std::collections::HashSet;
use crate::coordinates::Point;
use crate::ff::forcefield::EnergyFunction;

/// Dihedral angle
///
///       i  ---  j
///                 \  φ
///                  \
///                    k ---- l
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

impl TorsionalDihedral {

    #[inline(always)]
    fn half_v_phi(&self) -> f64{ self.v_phi/2.0 }
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
        self.half_v_phi() * (1. - (self.n_phi * self.phi0).cos()*(self.n_phi * phi).cos())
    }

    fn add_gradient(&self,
                    coordinates:      &Vec<Point>,
                    gradient: &mut Vec<Point>) {

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
        let v6 = -v2*v3 + v0*v5;
        let v7 = v0.powi(2);
        let v8 = -y_j + y_k;
        let v9 = v8.powi(2);
        let v10 = v3.powi(2);
        let v11 = -y_k + y_l;
        let v12 = x_j - x_k;
        let v13 = v12*v11 - v4*v1;
        let v14 = y_i - y_j;
        let v15 = (v7 + v9 + v10).sqrt();
        let v16 = v6.powi(2);
        let v17 = v14*v3 - v8*v5;
        let v18 = v17.powi(2);
        let v19 = (v2*v8 - v0*v14).powi(2);
        let v20 = (v19 + v16 + v18).sqrt();
        let v21 = v13.powi(2);
        let v22 = v15*v20;
        let v23 = -z_k + z_l;
        let v24 = z_j - z_k;
        let v25 = -v0*v6/v22 + v8*v17/v22;
        let v26 = v1*v23 - v11*v24;
        let v27 = v26.powi(2);
        let v28 = (-v12*v23 + v4*v24).powi(2);
        let v29 = -v12*v23 + v4*v24;
        let v30 = v2*v8 - v0*v14;
        let v31 = v0*v30/v22 - v3*v17/v22;
        let v32 = (v21 + v28 + v27).sqrt();
        let v33 = -v8*v30/v22 + v3*v6/v22;
        let v34 = v20*v32;
        let v35 = -v13*v25/v32 - v29*v31/v32 - v26*v33/v32;
        gradient[self.i].x += 0.5*self.n_phi*self.v_phi*(v35*(v8*v13/v34 + v24*v29/v34 + v30*v13*(-(-2.*y_j + 2.*y_k)*v30/2. - (2.*z_j - 2.*z_k)*v6/2.)/((v19 + v16 + v18).powf(1.5)*v32) + v6*v29*(-(-2.*y_j + 2.*y_k)*v30/2. - (2.*z_j - 2.*z_k)*v6/2.)/((v19 + v16 + v18).powf(1.5)*v32) + v17*(-(-2.*y_j + 2.*y_k)*v30/2. - (2.*z_j - 2.*z_k)*v6/2.)*v26/((v19 + v16 + v18).powf(1.5)*v32))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v13*(-v0*v24/v22 - v0*v6*(-(-2.*y_j + 2.*y_k)*v30/2. - (2.*z_j - 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v8*v17*(-(-2.*y_j + 2.*y_k)*v30/2. - (2.*z_j - 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)))/v32 + v29*(v0*v8/v22 + v0*v30*(-(-2.*y_j + 2.*y_k)*v30/2. - (2.*z_j - 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)) - v3*v17*(-(-2.*y_j + 2.*y_k)*v30/2. - (2.*z_j - 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)))/v32 + v26*(-v9/v22 - v8*v30*(-(-2.*y_j + 2.*y_k)*v30/2. - (2.*z_j - 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v3*v24/v22 + v3*v6*(-(-2.*y_j + 2.*y_k)*v30/2. - (2.*z_j - 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)))/v32)*(v30*v13/v34 + v6*v29/v34 + v17*v26/v34)/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
        gradient[self.i].y += 0.5*self.n_phi*self.v_phi*(v35*(v12*v13/v34 + v3*v26/v34 + v30*v13*(-(2.*x_j - 2.*x_k)*v30/2. - (-2.*z_j + 2.*z_k)*v17/2.)/((v19 + v16 + v18).powf(1.5)*v32) + v6*v29*(-(2.*x_j - 2.*x_k)*v30/2. - (-2.*z_j + 2.*z_k)*v17/2.)/((v19 + v16 + v18).powf(1.5)*v32) + (-(2.*x_j - 2.*x_k)*v30/2. - (-2.*z_j + 2.*z_k)*v17/2.)*v17*v26/((v19 + v16 + v18).powf(1.5)*v32))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v13*(-v0*v6*(-(2.*x_j - 2.*x_k)*v30/2. - (-2.*z_j + 2.*z_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v8*v3/v22 + v8*(-(2.*x_j - 2.*x_k)*v30/2. - (-2.*z_j + 2.*z_k)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)))/v32 + v29*(v0*v12/v22 + v0*v30*(-(2.*x_j - 2.*x_k)*v30/2. - (-2.*z_j + 2.*z_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) - v10/v22 - v3*(-(2.*x_j - 2.*x_k)*v30/2. - (-2.*z_j + 2.*z_k)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)))/v32 + v26*(-v12*v8/v22 - v8*v30*(-(2.*x_j - 2.*x_k)*v30/2. - (-2.*z_j + 2.*z_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v3*v6*(-(2.*x_j - 2.*x_k)*v30/2. - (-2.*z_j + 2.*z_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)))/v32)*(v30*v13/v34 + v6*v29/v34 + v17*v26/v34)/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
        gradient[self.i].z += 0.5*self.n_phi*self.v_phi*(v35*(v0*v29/v34 + v1*v26/v34 + v30*(-(-2.*x_j + 2.*x_k)*v6/2. - (2.*y_j - 2.*y_k)*v17/2.)*v13/((v19 + v16 + v18).powf(1.5)*v32) + v6*(-(-2.*x_j + 2.*x_k)*v6/2. - (2.*y_j - 2.*y_k)*v17/2.)*v29/((v19 + v16 + v18).powf(1.5)*v32) + (-(-2.*x_j + 2.*x_k)*v6/2. - (2.*y_j - 2.*y_k)*v17/2.)*v17*v26/((v19 + v16 + v18).powf(1.5)*v32))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v13*(-v7/v22 - v0*v6*(-(-2.*x_j + 2.*x_k)*v6/2. - (2.*y_j - 2.*y_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v8*v1/v22 + v8*(-(-2.*x_j + 2.*x_k)*v6/2. - (2.*y_j - 2.*y_k)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)))/v32 + v29*(v0*v30*(-(-2.*x_j + 2.*x_k)*v6/2. - (2.*y_j - 2.*y_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) - v1*v3/v22 - v3*(-(-2.*x_j + 2.*x_k)*v6/2. - (2.*y_j - 2.*y_k)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)))/v32 + v26*(v0*v3/v22 - v8*v30*(-(-2.*x_j + 2.*x_k)*v6/2. - (2.*y_j - 2.*y_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v3*v6*(-(-2.*x_j + 2.*x_k)*v6/2. - (2.*y_j - 2.*y_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)))/v32)*(v30*v13/v34 + v6*v29/v34 + v17*v26/v34)/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
        gradient[self.j].x += 0.5*self.n_phi*self.v_phi*(v35*((y_i - y_k)*v13/v34 + v11*v30/v34 + (-z_i + z_k)*v29/v34 + (z_k - z_l)*v6/v34 + v30*v13*(-(2.*y_i - 2.*y_k)*v30/2. - (-2.*z_i + 2.*z_k)*v6/2.)/((v19 + v16 + v18).powf(1.5)*v32) + v30*v13*(-(-2.*y_k + 2.*y_l)*v13/2. - (2.*z_k - 2.*z_l)*v29/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + v6*v29*(-(2.*y_i - 2.*y_k)*v30/2. - (-2.*z_i + 2.*z_k)*v6/2.)/((v19 + v16 + v18).powf(1.5)*v32) + v6*v29*(-(-2.*y_k + 2.*y_l)*v13/2. - (2.*z_k - 2.*z_l)*v29/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + v17*(-(2.*y_i - 2.*y_k)*v30/2. - (-2.*z_i + 2.*z_k)*v6/2.)*v26/((v19 + v16 + v18).powf(1.5)*v32) + v17*v26*(-(-2.*y_k + 2.*y_l)*v13/2. - (2.*z_k - 2.*z_l)*v29/2.)/(v20*(v21 + v28 + v27).powf(1.5)))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34)*(v11*v25/v32 + (z_k - z_l)*v31/v32 + v13*(-(-2.*y_k + 2.*y_l)*v13/2. - (2.*z_k - 2.*z_l)*v29/2.)*v25/(v21 + v28 + v27).powf(1.5) + v13*(-v7*v6/((v7 + v9 + v10).powf(1.5)*v20) + v0*v8*v17/((v7 + v9 + v10).powf(1.5)*v20) - v0*(-z_i + z_k)/v22 - v0*v6*(-(2.*y_i - 2.*y_k)*v30/2. - (-2.*z_i + 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v8*v17*(-(2.*y_i - 2.*y_k)*v30/2. - (-2.*z_i + 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v6/v22)/v32 + v29*(-(-2.*y_k + 2.*y_l)*v13/2. - (2.*z_k - 2.*z_l)*v29/2.)*v31/(v21 + v28 + v27).powf(1.5) + v29*(v7*v30/((v7 + v9 + v10).powf(1.5)*v20) + v0*(y_i - y_k)/v22 - v0*v3*v17/((v7 + v9 + v10).powf(1.5)*v20) + v0*v30*(-(2.*y_i - 2.*y_k)*v30/2. - (-2.*z_i + 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)) - v3*v17*(-(2.*y_i - 2.*y_k)*v30/2. - (-2.*z_i + 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)) - v30/v22)/v32 + v26*(-(-2.*y_k + 2.*y_l)*v13/2. - (2.*z_k - 2.*z_l)*v29/2.)*v33/(v21 + v28 + v27).powf(1.5) + v26*(-v0*v8*v30/((v7 + v9 + v10).powf(1.5)*v20) + v0*v3*v6/((v7 + v9 + v10).powf(1.5)*v20) - (y_i - y_k)*v8/v22 - v8*v30*(-(2.*y_i - 2.*y_k)*v30/2. - (-2.*z_i + 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + (-z_i + z_k)*v3/v22 + v3*v6*(-(2.*y_i - 2.*y_k)*v30/2. - (-2.*z_i + 2.*z_k)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)))/v32)/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
        gradient[self.j].y += 0.5*self.n_phi*self.v_phi*(v35*((-x_i + x_k)*v13/v34 + (x_k - x_l)*v30/v34 + (z_i - z_k)*v26/v34 + v23*v17/v34 + (-(-2.*x_i + 2.*x_k)*v30/2. - (2.*z_i - 2.*z_k)*v17/2.)*v30*v13/((v19 + v16 + v18).powf(1.5)*v32) + (-(-2.*x_i + 2.*x_k)*v30/2. - (2.*z_i - 2.*z_k)*v17/2.)*v6*v29/((v19 + v16 + v18).powf(1.5)*v32) + (-(-2.*x_i + 2.*x_k)*v30/2. - (2.*z_i - 2.*z_k)*v17/2.)*v17*v26/((v19 + v16 + v18).powf(1.5)*v32) + v30*v13*(-(2.*x_k - 2.*x_l)*v13/2. - (-2.*z_k + 2.*z_l)*v26/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + v6*v29*(-(2.*x_k - 2.*x_l)*v13/2. - (-2.*z_k + 2.*z_l)*v26/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + (-(2.*x_k - 2.*x_l)*v13/2. - (-2.*z_k + 2.*z_l)*v26/2.)*v17*v26/(v20*(v21 + v28 + v27).powf(1.5)))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34)*((x_k - x_l)*v25/v32 + v23*v33/v32 + v13*(-(2.*x_k - 2.*x_l)*v13/2. - (-2.*z_k + 2.*z_l)*v26/2.)*v25/(v21 + v28 + v27).powf(1.5) + v13*(-v0*v8*v6/((v7 + v9 + v10).powf(1.5)*v20) - v0*(-(-2.*x_i + 2.*x_k)*v30/2. - (2.*z_i - 2.*z_k)*v17/2.)*v6/(v15*(v19 + v16 + v18).powf(1.5)) + v9*v17/((v7 + v9 + v10).powf(1.5)*v20) + v8*(z_i - z_k)/v22 + v8*(-(-2.*x_i + 2.*x_k)*v30/2. - (2.*z_i - 2.*z_k)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)) - v17/v22)/v32 + v29*(-(2.*x_k - 2.*x_l)*v13/2. - (-2.*z_k + 2.*z_l)*v26/2.)*v31/(v21 + v28 + v27).powf(1.5) + v29*((-x_i + x_k)*v0/v22 + v0*v8*v30/((v7 + v9 + v10).powf(1.5)*v20) + v0*(-(-2.*x_i + 2.*x_k)*v30/2. - (2.*z_i - 2.*z_k)*v17/2.)*v30/(v15*(v19 + v16 + v18).powf(1.5)) - v8*v3*v17/((v7 + v9 + v10).powf(1.5)*v20) - (z_i - z_k)*v3/v22 - v3*(-(-2.*x_i + 2.*x_k)*v30/2. - (2.*z_i - 2.*z_k)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)))/v32 + (-(2.*x_k - 2.*x_l)*v13/2. - (-2.*z_k + 2.*z_l)*v26/2.)*v26*v33/(v21 + v28 + v27).powf(1.5) + v26*(-(-x_i + x_k)*v8/v22 - v9*v30/((v7 + v9 + v10).powf(1.5)*v20) + v8*v3*v6/((v7 + v9 + v10).powf(1.5)*v20) - v8*(-(-2.*x_i + 2.*x_k)*v30/2. - (2.*z_i - 2.*z_k)*v17/2.)*v30/(v15*(v19 + v16 + v18).powf(1.5)) + v3*(-(-2.*x_i + 2.*x_k)*v30/2. - (2.*z_i - 2.*z_k)*v17/2.)*v6/(v15*(v19 + v16 + v18).powf(1.5)) + v30/v22)/v32)/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
        gradient[self.j].z += 0.5*self.n_phi*self.v_phi*(v35*((x_i - x_k)*v29/v34 + v4*v6/v34 + (-y_i + y_k)*v26/v34 + (y_k - y_l)*v17/v34 + v30*(-(2.*x_i - 2.*x_k)*v6/2. - (-2.*y_i + 2.*y_k)*v17/2.)*v13/((v19 + v16 + v18).powf(1.5)*v32) + v30*v13*(-(-2.*x_k + 2.*x_l)*v29/2. - (2.*y_k - 2.*y_l)*v26/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + v6*(-(2.*x_i - 2.*x_k)*v6/2. - (-2.*y_i + 2.*y_k)*v17/2.)*v29/((v19 + v16 + v18).powf(1.5)*v32) + v6*v29*(-(-2.*x_k + 2.*x_l)*v29/2. - (2.*y_k - 2.*y_l)*v26/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + (-(2.*x_i - 2.*x_k)*v6/2. - (-2.*y_i + 2.*y_k)*v17/2.)*v17*v26/((v19 + v16 + v18).powf(1.5)*v32) + (-(-2.*x_k + 2.*x_l)*v29/2. - (2.*y_k - 2.*y_l)*v26/2.)*v17*v26/(v20*(v21 + v28 + v27).powf(1.5)))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34)*(v4*v31/v32 + (y_k - y_l)*v33/v32 + v13*(-(-2.*x_k + 2.*x_l)*v29/2. - (2.*y_k - 2.*y_l)*v26/2.)*v25/(v21 + v28 + v27).powf(1.5) + v13*(-(x_i - x_k)*v0/v22 - v0*v3*v6/((v7 + v9 + v10).powf(1.5)*v20) - v0*v6*(-(2.*x_i - 2.*x_k)*v6/2. - (-2.*y_i + 2.*y_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + (-y_i + y_k)*v8/v22 + v8*v3*v17/((v7 + v9 + v10).powf(1.5)*v20) + v8*(-(2.*x_i - 2.*x_k)*v6/2. - (-2.*y_i + 2.*y_k)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)))/v32 + v29*(-(-2.*x_k + 2.*x_l)*v29/2. - (2.*y_k - 2.*y_l)*v26/2.)*v31/(v21 + v28 + v27).powf(1.5) + v29*(v0*v3*v30/((v7 + v9 + v10).powf(1.5)*v20) + v0*v30*(-(2.*x_i - 2.*x_k)*v6/2. - (-2.*y_i + 2.*y_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) - (-y_i + y_k)*v3/v22 - v10*v17/((v7 + v9 + v10).powf(1.5)*v20) - v3*(-(2.*x_i - 2.*x_k)*v6/2. - (-2.*y_i + 2.*y_k)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)) + v17/v22)/v32 + (-(-2.*x_k + 2.*x_l)*v29/2. - (2.*y_k - 2.*y_l)*v26/2.)*v26*v33/(v21 + v28 + v27).powf(1.5) + v26*((x_i - x_k)*v3/v22 - v8*v3*v30/((v7 + v9 + v10).powf(1.5)*v20) - v8*v30*(-(2.*x_i - 2.*x_k)*v6/2. - (-2.*y_i + 2.*y_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v10*v6/((v7 + v9 + v10).powf(1.5)*v20) + v3*v6*(-(2.*x_i - 2.*x_k)*v6/2. - (-2.*y_i + 2.*y_k)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) - v6/v22)/v32)/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
        gradient[self.k].x += 0.5*self.n_phi*self.v_phi*(v35*((-y_i + y_j)*v13/v34 + (y_j - y_l)*v30/v34 + v5*v29/v34 + (-z_j + z_l)*v6/v34 + v30*v13*(-(-2.*y_i + 2.*y_j)*v30/2. - (2.*z_i - 2.*z_j)*v6/2.)/((v19 + v16 + v18).powf(1.5)*v32) + v30*v13*(-(2.*y_j - 2.*y_l)*v13/2. - (-2.*z_j + 2.*z_l)*v29/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + v6*v29*(-(-2.*y_i + 2.*y_j)*v30/2. - (2.*z_i - 2.*z_j)*v6/2.)/((v19 + v16 + v18).powf(1.5)*v32) + v6*v29*(-(2.*y_j - 2.*y_l)*v13/2. - (-2.*z_j + 2.*z_l)*v29/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + (-(-2.*y_i + 2.*y_j)*v30/2. - (2.*z_i - 2.*z_j)*v6/2.)*v17*v26/((v19 + v16 + v18).powf(1.5)*v32) + v17*v26*(-(2.*y_j - 2.*y_l)*v13/2. - (-2.*z_j + 2.*z_l)*v29/2.)/(v20*(v21 + v28 + v27).powf(1.5)))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34)*((y_j - y_l)*v25/v32 + (-z_j + z_l)*v31/v32 + v13*(-(2.*y_j - 2.*y_l)*v13/2. - (-2.*z_j + 2.*z_l)*v29/2.)*v25/(v21 + v28 + v27).powf(1.5) + v13*(-v0*v12*v6/((v7 + v9 + v10).powf(1.5)*v20) - v0*v5/v22 - v0*v6*(-(-2.*y_i + 2.*y_j)*v30/2. - (2.*z_i - 2.*z_j)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v12*v8*v17/((v7 + v9 + v10).powf(1.5)*v20) + v8*(-(-2.*y_i + 2.*y_j)*v30/2. - (2.*z_i - 2.*z_j)*v6/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)) - v6/v22)/v32 + v29*(-(2.*y_j - 2.*y_l)*v13/2. - (-2.*z_j + 2.*z_l)*v29/2.)*v31/(v21 + v28 + v27).powf(1.5) + v29*(v0*v12*v30/((v7 + v9 + v10).powf(1.5)*v20) + v0*(-y_i + y_j)/v22 + v0*v30*(-(-2.*y_i + 2.*y_j)*v30/2. - (2.*z_i - 2.*z_j)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)) - v12*v3*v17/((v7 + v9 + v10).powf(1.5)*v20) - v3*(-(-2.*y_i + 2.*y_j)*v30/2. - (2.*z_i - 2.*z_j)*v6/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)) + v30/v22)/v32 + v26*(-(2.*y_j - 2.*y_l)*v13/2. - (-2.*z_j + 2.*z_l)*v29/2.)*v33/(v21 + v28 + v27).powf(1.5) + v26*(-v12*v8*v30/((v7 + v9 + v10).powf(1.5)*v20) + v12*v3*v6/((v7 + v9 + v10).powf(1.5)*v20) - (-y_i + y_j)*v8/v22 - v8*v30*(-(-2.*y_i + 2.*y_j)*v30/2. - (2.*z_i - 2.*z_j)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v5*v3/v22 + v3*v6*(-(-2.*y_i + 2.*y_j)*v30/2. - (2.*z_i - 2.*z_j)*v6/2.)/(v15*(v19 + v16 + v18).powf(1.5)))/v32)/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
        gradient[self.k].y += 0.5*self.n_phi*self.v_phi*(v35*(v2*v13/v34 + (-x_j + x_l)*v30/v34 + (-z_i + z_j)*v26/v34 + (z_j - z_l)*v17/v34 + v30*(-(2.*x_i - 2.*x_j)*v30/2. - (-2.*z_i + 2.*z_j)*v17/2.)*v13/((v19 + v16 + v18).powf(1.5)*v32) + v30*(-(-2.*x_j + 2.*x_l)*v13/2. - (2.*z_j - 2.*z_l)*v26/2.)*v13/(v20*(v21 + v28 + v27).powf(1.5)) + v6*(-(2.*x_i - 2.*x_j)*v30/2. - (-2.*z_i + 2.*z_j)*v17/2.)*v29/((v19 + v16 + v18).powf(1.5)*v32) + v6*(-(-2.*x_j + 2.*x_l)*v13/2. - (2.*z_j - 2.*z_l)*v26/2.)*v29/(v20*(v21 + v28 + v27).powf(1.5)) + (-(2.*x_i - 2.*x_j)*v30/2. - (-2.*z_i + 2.*z_j)*v17/2.)*v17*v26/((v19 + v16 + v18).powf(1.5)*v32) + (-(-2.*x_j + 2.*x_l)*v13/2. - (2.*z_j - 2.*z_l)*v26/2.)*v17*v26/(v20*(v21 + v28 + v27).powf(1.5)))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34)*((-x_j + x_l)*v25/v32 + (z_j - z_l)*v33/v32 + (-(-2.*x_j + 2.*x_l)*v13/2. - (2.*z_j - 2.*z_l)*v26/2.)*v13*v25/(v21 + v28 + v27).powf(1.5) + (-(-2.*x_j + 2.*x_l)*v13/2. - (2.*z_j - 2.*z_l)*v26/2.)*v29*v31/(v21 + v28 + v27).powf(1.5) + (-(-2.*x_j + 2.*x_l)*v13/2. - (2.*z_j - 2.*z_l)*v26/2.)*v26*v33/(v21 + v28 + v27).powf(1.5) + v13*(-v0*v1*v6/((v7 + v9 + v10).powf(1.5)*v20) - v0*v6*(-(2.*x_i - 2.*x_j)*v30/2. - (-2.*z_i + 2.*z_j)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v8*v1*v17/((v7 + v9 + v10).powf(1.5)*v20) + v8*(-z_i + z_j)/v22 + v8*(-(2.*x_i - 2.*x_j)*v30/2. - (-2.*z_i + 2.*z_j)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)) + v17/v22)/v32 + v29*(v2*v0/v22 + v0*v1*v30/((v7 + v9 + v10).powf(1.5)*v20) + v0*v30*(-(2.*x_i - 2.*x_j)*v30/2. - (-2.*z_i + 2.*z_j)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) - v1*v3*v17/((v7 + v9 + v10).powf(1.5)*v20) - (-z_i + z_j)*v3/v22 - v3*(-(2.*x_i - 2.*x_j)*v30/2. - (-2.*z_i + 2.*z_j)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)))/v32 + v26*(-v2*v8/v22 - v8*v1*v30/((v7 + v9 + v10).powf(1.5)*v20) - v8*v30*(-(2.*x_i - 2.*x_j)*v30/2. - (-2.*z_i + 2.*z_j)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) + v1*v3*v6/((v7 + v9 + v10).powf(1.5)*v20) + v3*v6*(-(2.*x_i - 2.*x_j)*v30/2. - (-2.*z_i + 2.*z_j)*v17/2.)/(v15*(v19 + v16 + v18).powf(1.5)) - v30/v22)/v32)/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
        gradient[self.k].z += 0.5*self.n_phi*self.v_phi*(v35*((-x_i + x_j)*v29/v34 + (x_j - x_l)*v6/v34 + v14*v26/v34 + (-y_j + y_l)*v17/v34 + (-(-2.*x_i + 2.*x_j)*v6/2. - (2.*y_i - 2.*y_j)*v17/2.)*v30*v13/((v19 + v16 + v18).powf(1.5)*v32) + (-(-2.*x_i + 2.*x_j)*v6/2. - (2.*y_i - 2.*y_j)*v17/2.)*v6*v29/((v19 + v16 + v18).powf(1.5)*v32) + (-(-2.*x_i + 2.*x_j)*v6/2. - (2.*y_i - 2.*y_j)*v17/2.)*v17*v26/((v19 + v16 + v18).powf(1.5)*v32) + v30*v13*(-(2.*x_j - 2.*x_l)*v29/2. - (-2.*y_j + 2.*y_l)*v26/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + v6*v29*(-(2.*x_j - 2.*x_l)*v29/2. - (-2.*y_j + 2.*y_l)*v26/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + (-(2.*x_j - 2.*x_l)*v29/2. - (-2.*y_j + 2.*y_l)*v26/2.)*v17*v26/(v20*(v21 + v28 + v27).powf(1.5)))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34)*((x_j - x_l)*v31/v32 + (-y_j + y_l)*v33/v32 + v13*(-(2.*x_j - 2.*x_l)*v29/2. - (-2.*y_j + 2.*y_l)*v26/2.)*v25/(v21 + v28 + v27).powf(1.5) + v13*(-(-x_i + x_j)*v0/v22 - v0*v24*v6/((v7 + v9 + v10).powf(1.5)*v20) - v0*(-(-2.*x_i + 2.*x_j)*v6/2. - (2.*y_i - 2.*y_j)*v17/2.)*v6/(v15*(v19 + v16 + v18).powf(1.5)) + v14*v8/v22 + v8*v24*v17/((v7 + v9 + v10).powf(1.5)*v20) + v8*(-(-2.*x_i + 2.*x_j)*v6/2. - (2.*y_i - 2.*y_j)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)))/v32 + v29*(-(2.*x_j - 2.*x_l)*v29/2. - (-2.*y_j + 2.*y_l)*v26/2.)*v31/(v21 + v28 + v27).powf(1.5) + v29*(v0*v24*v30/((v7 + v9 + v10).powf(1.5)*v20) + v0*(-(-2.*x_i + 2.*x_j)*v6/2. - (2.*y_i - 2.*y_j)*v17/2.)*v30/(v15*(v19 + v16 + v18).powf(1.5)) - v14*v3/v22 - v3*v24*v17/((v7 + v9 + v10).powf(1.5)*v20) - v3*(-(-2.*x_i + 2.*x_j)*v6/2. - (2.*y_i - 2.*y_j)*v17/2.)*v17/(v15*(v19 + v16 + v18).powf(1.5)) - v17/v22)/v32 + (-(2.*x_j - 2.*x_l)*v29/2. - (-2.*y_j + 2.*y_l)*v26/2.)*v26*v33/(v21 + v28 + v27).powf(1.5) + v26*((-x_i + x_j)*v3/v22 - v8*v24*v30/((v7 + v9 + v10).powf(1.5)*v20) - v8*(-(-2.*x_i + 2.*x_j)*v6/2. - (2.*y_i - 2.*y_j)*v17/2.)*v30/(v15*(v19 + v16 + v18).powf(1.5)) + v3*v24*v6/((v7 + v9 + v10).powf(1.5)*v20) + v3*(-(-2.*x_i + 2.*x_j)*v6/2. - (2.*y_i - 2.*y_j)*v17/2.)*v6/(v15*(v19 + v16 + v18).powf(1.5)) + v6/v22)/v32)/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
        gradient[self.l].x += 0.5*self.n_phi*self.v_phi*(v35*(v8*v30/v34 + v24*v6/v34 + v30*v13*(-(-2.*y_j + 2.*y_k)*v13/2. - (2.*z_j - 2.*z_k)*v29/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + v6*v29*(-(-2.*y_j + 2.*y_k)*v13/2. - (2.*z_j - 2.*z_k)*v29/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + v17*(-(-2.*y_j + 2.*y_k)*v13/2. - (2.*z_j - 2.*z_k)*v29/2.)*v26/(v20*(v21 + v28 + v27).powf(1.5)))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34)*(v8*v25/v32 + v24*v31/v32 + v13*(-(-2.*y_j + 2.*y_k)*v13/2. - (2.*z_j - 2.*z_k)*v29/2.)*v25/(v21 + v28 + v27).powf(1.5) + v29*(-(-2.*y_j + 2.*y_k)*v13/2. - (2.*z_j - 2.*z_k)*v29/2.)*v31/(v21 + v28 + v27).powf(1.5) + (-(-2.*y_j + 2.*y_k)*v13/2. - (2.*z_j - 2.*z_k)*v29/2.)*v26*v33/(v21 + v28 + v27).powf(1.5))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
        gradient[self.l].y += 0.5*self.n_phi*self.v_phi*(v35*(v12*v30/v34 + v3*v17/v34 + v30*v13*(-(2.*x_j - 2.*x_k)*v13/2. - (-2.*z_j + 2.*z_k)*v26/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + v6*v29*(-(2.*x_j - 2.*x_k)*v13/2. - (-2.*z_j + 2.*z_k)*v26/2.)/(v20*(v21 + v28 + v27).powf(1.5)) + (-(2.*x_j - 2.*x_k)*v13/2. - (-2.*z_j + 2.*z_k)*v26/2.)*v17*v26/(v20*(v21 + v28 + v27).powf(1.5)))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34)*(v12*v25/v32 + v3*v33/v32 + v13*(-(2.*x_j - 2.*x_k)*v13/2. - (-2.*z_j + 2.*z_k)*v26/2.)*v25/(v21 + v28 + v27).powf(1.5) + v29*(-(2.*x_j - 2.*x_k)*v13/2. - (-2.*z_j + 2.*z_k)*v26/2.)*v31/(v21 + v28 + v27).powf(1.5) + (-(2.*x_j - 2.*x_k)*v13/2. - (-2.*z_j + 2.*z_k)*v26/2.)*v26*v33/(v21 + v28 + v27).powf(1.5))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
        gradient[self.l].z += 0.5*self.n_phi*self.v_phi*(v35*(v0*v6/v34 + v1*v17/v34 + v30*(-(-2.*x_j + 2.*x_k)*v29/2. - (2.*y_j - 2.*y_k)*v26/2.)*v13/(v20*(v21 + v28 + v27).powf(1.5)) + v6*(-(-2.*x_j + 2.*x_k)*v29/2. - (2.*y_j - 2.*y_k)*v26/2.)*v29/(v20*(v21 + v28 + v27).powf(1.5)) + (-(-2.*x_j + 2.*x_k)*v29/2. - (2.*y_j - 2.*y_k)*v26/2.)*v17*v26/(v20*(v21 + v28 + v27).powf(1.5)))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34)*(v0*v31/v32 + v1*v33/v32 + (-(-2.*x_j + 2.*x_k)*v29/2. - (2.*y_j - 2.*y_k)*v26/2.)*v13*v25/(v21 + v28 + v27).powf(1.5) + (-(-2.*x_j + 2.*x_k)*v29/2. - (2.*y_j - 2.*y_k)*v26/2.)*v29*v31/(v21 + v28 + v27).powf(1.5) + (-(-2.*x_j + 2.*x_k)*v29/2. - (2.*y_j - 2.*y_k)*v26/2.)*v26*v33/(v21 + v28 + v27).powf(1.5))/((v13*v25/v32 + v29*v31/v32 + v26*v33/v32).powi(2) + (v30*v13/v34 + v6*v29/v34 + v17*v26/v34).powi(2)))*(self.n_phi*(v13*v25/v32 + v29*v31/v32 + v26*v33/v32).atan2( v30*v13/v34 + v6*v29/v34 + v17*v26/v34)).sin()*(self.n_phi*self.phi0).cos();
    }
}


/// Value of the torsional dihedral
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


/// Improper dihedral defining inversion about a centre
///
///                    j
///                  /
///        i  ---  c    ω
///                 \
///                  \
///                    k
///
/// where the angle is between the c-k axis and the plane
/// formed by atoms c,i,j. It may also be defined by the
/// c,i,j normal and the axis: ω = γ - π
pub struct InversionDihedral{

    pub(crate) c:  usize,
    pub(crate) i:  usize,
    pub(crate) j:  usize,
    pub(crate) k:  usize,

    pub(crate) gamma:  f64,
    pub(crate) k_cijk: f64,
}

impl EnergyFunction for InversionDihedral {
    fn involves_idxs(&self, idxs: Vec<usize>) -> bool {
        idxs.clone().sort() == Vec::from([self.c, self.i, self.j, self.k]).sort()
    }

    fn force_constant(&self) -> f64 { self.k_cijk }

    fn energy(&self, coordinates: &Vec<Point>) -> f64 {
        todo!()
    }

    fn add_gradient(&self, coordinates: &Vec<Point>, gradient: &mut Vec<Point>) {
        todo!()
    }
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
        assert!(is_very_close(dihedral.force_constant(), 0.0));

        let expected_phi = -1.7593;  // -100.8 degrees

        assert!(is_close(phi(&dihedral, &x), expected_phi, 1E-3));
    }

}
