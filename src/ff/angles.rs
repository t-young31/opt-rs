use crate::coordinates::{angle_value, Point, Vector3D};
use crate::ff::forcefield::EnergyFunction;
use std::collections::HashSet;

trait HarmonicAngle {
    fn i(&self) -> usize;
    fn j(&self) -> usize;
    fn k(&self) -> usize;
}

pub struct HarmonicAngleTypeA {
    pub(crate) i: usize,
    pub(crate) j: usize,
    pub(crate) k: usize,
    pub(crate) k_ijk: f64,
    pub(crate) n: f64,
}

impl EnergyFunction for HarmonicAngleTypeA {
    fn involves_idxs(&self, idxs: Vec<usize>) -> bool {
        involves_idxs(self, idxs)
    }

    fn force_constant(&self) -> f64 {
        self.k_ijk
    }

    /// Energy: k/n^2 (1-cos(nθ))
    fn energy(&self, coordinates: &[Point]) -> f64 {
        let theta = theta(self, coordinates);

        (self.k_ijk / self.n.powi(2)) * (1.0 - (self.n * theta).cos())
    }

    /// Add the gradient for this term
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

        let v0 = -x_j + x_k;
        let v1 = (z_i - z_j).powi(2);
        let v2 = (y_i - y_j).powi(2);
        let v3 = (x_i - x_j).powi(2);
        let v4 = v3 + v2 + v1;
        let v5 = -z_j + z_k;
        let v6 = v5.powi(2);
        let v7 = z_i - z_j;
        let v8 = -y_j + y_k;
        let v9 = v8.powi(2);
        let v10 = v0.powi(2);
        let v11 = y_i - y_j;
        let v12 = x_i - x_j;
        let v13 = v10 + v9 + v6;
        let v14 = v4 * v13;
        let v15 = v13.sqrt();
        let v16 = v4.sqrt();
        let v17 = v16 * v15;
        let v18 = v12 * v0 + v11 * v8 + v7 * v5;
        let v19 = v18.powi(2);
        let v20 = (-v19 / v14 + 1.).sqrt();
        let v21 = (v18 / v17).acos();
        let v22 = (self.n * v21).sin();
        gradient[self.i].x +=
            -self.k_ijk * ((-x_i + x_j) * v18 / (v4.powf(1.5) * v15) + v0 / v17) * v22
                / (self.n * v20);
        gradient[self.i].y +=
            -self.k_ijk * ((-y_i + y_j) * v18 / (v4.powf(1.5) * v15) + v8 / v17) * v22
                / (self.n * v20);
        gradient[self.i].z +=
            -self.k_ijk * ((-z_i + z_j) * v18 / (v4.powf(1.5) * v15) + v5 / v17) * v22
                / (self.n * v20);
        gradient[self.j].x += -self.k_ijk
            * (v12 * v18 / (v4.powf(1.5) * v15)
                + v0 * v18 / (v16 * v13.powf(1.5))
                + (-x_i + 2. * x_j - x_k) / v17)
            * v22
            / (self.n * v20);
        gradient[self.j].y += -self.k_ijk
            * (v11 * v18 / (v4.powf(1.5) * v15)
                + v8 * v18 / (v16 * v13.powf(1.5))
                + (-y_i + 2. * y_j - y_k) / v17)
            * v22
            / (self.n * v20);
        gradient[self.j].z += -self.k_ijk
            * (v7 * v18 / (v4.powf(1.5) * v15)
                + v5 * v18 / (v16 * v13.powf(1.5))
                + (-z_i + 2. * z_j - z_k) / v17)
            * v22
            / (self.n * v20);
        gradient[self.k].x +=
            -self.k_ijk * (v12 / v17 + (x_j - x_k) * v18 / (v16 * v13.powf(1.5))) * v22
                / (self.n * v20);
        gradient[self.k].y +=
            -self.k_ijk * (v11 / v17 + (y_j - y_k) * v18 / (v16 * v13.powf(1.5))) * v22
                / (self.n * v20);
        gradient[self.k].z +=
            -self.k_ijk * (v7 / v17 + (z_j - z_k) * v18 / (v16 * v13.powf(1.5))) * v22
                / (self.n * v20);
    }
}
impl HarmonicAngle for HarmonicAngleTypeA {
    fn i(&self) -> usize {
        self.i
    }
    fn j(&self) -> usize {
        self.j
    }
    fn k(&self) -> usize {
        self.k
    }
}

pub struct HarmonicAngleTypeB {
    pub(crate) i: usize,
    pub(crate) j: usize,
    pub(crate) k: usize,
    pub(crate) k_ijk: f64,
    pub(crate) c0: f64,
    pub(crate) c1: f64,
    pub(crate) c2: f64,
}

impl EnergyFunction for HarmonicAngleTypeB {
    fn involves_idxs(&self, idxs: Vec<usize>) -> bool {
        involves_idxs(self, idxs)
    }

    fn force_constant(&self) -> f64 {
        self.k_ijk
    }

    /// Energy: k(c0 + c1 cos(θ) + c2 cos(2θ))
    fn energy(&self, coordinates: &[Point]) -> f64 {
        let theta = theta(self, coordinates);

        self.k_ijk * (self.c0 + self.c1 * theta.cos() + self.c2 * (2. * theta).cos())
    }

    /// Add the gradient for this term
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

        let v0 = -x_j + x_k;
        let v1 = (z_i - z_j).powi(2);
        let v2 = (y_i - y_j).powi(2);
        let v3 = (x_i - x_j).powi(2);
        let v4 = v3 + v2 + v1;
        let v5 = -z_j + z_k;
        let v6 = v5.powi(2);
        let v7 = z_i - z_j;
        let v8 = -y_j + y_k;
        let v9 = v8.powi(2);
        let v10 = v0.powi(2);
        let v11 = y_i - y_j;
        let v12 = x_i - x_j;
        let v13 = v10 + v9 + v6;
        let v14 = v4 * v13;
        let v15 = v13.sqrt();
        let v16 = v4.sqrt();
        let v17 = v16 * v15;
        let v18 = v12 * v0 + v11 * v8 + v7 * v5;
        let v19 = v18.powi(2);
        let v20 = (v18 / v17).acos();
        let v21 = (-v19 / v14 + 1.).sqrt();
        let v22 = (2. * v20).sin();
        gradient[self.i].x += self.k_ijk
            * (self.c1 * (-x_i + x_j) * v18 / (v4.powf(1.5) * v15)
                + self.c1 * v0 / v17
                + 2. * self.c2 * ((-x_i + x_j) * v18 / (v4.powf(1.5) * v15) + v0 / v17) * v22
                    / v21);
        gradient[self.i].y += self.k_ijk
            * (self.c1 * (-y_i + y_j) * v18 / (v4.powf(1.5) * v15)
                + self.c1 * v8 / v17
                + 2. * self.c2 * ((-y_i + y_j) * v18 / (v4.powf(1.5) * v15) + v8 / v17) * v22
                    / v21);
        gradient[self.i].z += self.k_ijk
            * (self.c1 * (-z_i + z_j) * v18 / (v4.powf(1.5) * v15)
                + self.c1 * v5 / v17
                + 2. * self.c2 * ((-z_i + z_j) * v18 / (v4.powf(1.5) * v15) + v5 / v17) * v22
                    / v21);
        gradient[self.j].x += self.k_ijk
            * (self.c1 * v12 * v18 / (v4.powf(1.5) * v15)
                + self.c1 * v0 * v18 / (v16 * v13.powf(1.5))
                + self.c1 * (-x_i + 2. * x_j - x_k) / v17
                + 2. * self.c2
                    * (v12 * v18 / (v4.powf(1.5) * v15)
                        + v0 * v18 / (v16 * v13.powf(1.5))
                        + (-x_i + 2. * x_j - x_k) / v17)
                    * v22
                    / v21);
        gradient[self.j].y += self.k_ijk
            * (self.c1 * v11 * v18 / (v4.powf(1.5) * v15)
                + self.c1 * v8 * v18 / (v16 * v13.powf(1.5))
                + self.c1 * (-y_i + 2. * y_j - y_k) / v17
                + 2. * self.c2
                    * (v11 * v18 / (v4.powf(1.5) * v15)
                        + v8 * v18 / (v16 * v13.powf(1.5))
                        + (-y_i + 2. * y_j - y_k) / v17)
                    * v22
                    / v21);
        gradient[self.j].z += self.k_ijk
            * (self.c1 * v7 * v18 / (v4.powf(1.5) * v15)
                + self.c1 * v5 * v18 / (v16 * v13.powf(1.5))
                + self.c1 * (-z_i + 2. * z_j - z_k) / v17
                + 2. * self.c2
                    * (v7 * v18 / (v4.powf(1.5) * v15)
                        + v5 * v18 / (v16 * v13.powf(1.5))
                        + (-z_i + 2. * z_j - z_k) / v17)
                    * v22
                    / v21);
        gradient[self.k].x += self.k_ijk
            * (self.c1 * v12 / v17
                + self.c1 * (x_j - x_k) * v18 / (v16 * v13.powf(1.5))
                + 2. * self.c2 * (v12 / v17 + (x_j - x_k) * v18 / (v16 * v13.powf(1.5))) * v22
                    / v21);
        gradient[self.k].y += self.k_ijk
            * (self.c1 * v11 / v17
                + self.c1 * (y_j - y_k) * v18 / (v16 * v13.powf(1.5))
                + 2. * self.c2 * (v11 / v17 + (y_j - y_k) * v18 / (v16 * v13.powf(1.5))) * v22
                    / v21);
        gradient[self.k].z += self.k_ijk
            * (self.c1 * v7 / v17
                + self.c1 * (z_j - z_k) * v18 / (v16 * v13.powf(1.5))
                + 2. * self.c2 * (v7 / v17 + (z_j - z_k) * v18 / (v16 * v13.powf(1.5))) * v22
                    / v21);
    }
}
impl HarmonicAngle for HarmonicAngleTypeB {
    fn i(&self) -> usize {
        self.i
    }
    fn j(&self) -> usize {
        self.j
    }
    fn k(&self) -> usize {
        self.k
    }
}

fn involves_idxs(angle: &dyn HarmonicAngle, idxs: Vec<usize>) -> bool {
    idxs.len() == 3 && HashSet::from([angle.i(), angle.j(), angle.k()]) == HashSet::from_iter(idxs)
}

#[inline(always)]
fn theta(angle: &dyn HarmonicAngle, coordinates: &[Point]) -> f64 {
    angle_value(angle.i(), angle.j(), angle.k(), coordinates)
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

    struct TestHarmonicAngle {
        i: usize,
        j: usize,
        k: usize,
    }

    impl HarmonicAngle for TestHarmonicAngle {
        fn i(&self) -> usize {
            self.i
        }
        fn j(&self) -> usize {
            self.j
        }
        fn k(&self) -> usize {
            self.k
        }
    }

    /// Ensure a bend angle can be calculated correctly. Compared to an Avogadro reference
    #[test]
    fn test_angle_calculation() {
        // Coordinates for a H2O molecule
        let x: Vec<Point> = Vec::from([
            Point {
                x: -2.17305,
                y: -1.64172,
                z: 0.01027,
            },
            Point {
                x: -1.18373,
                y: -1.60261,
                z: -0.00725,
            },
            Point {
                x: -2.46242,
                y: -0.77749,
                z: -0.37700,
            },
        ]);

        let angle = theta(&TestHarmonicAngle { i: 1, j: 0, k: 2 }, &x);

        assert!(is_close(angle, 1.8239, 1E-3));
    }
}
