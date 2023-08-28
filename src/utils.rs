use crate::ff::forcefield::Forcefield;
use crate::molecule::Molecule;

pub fn is_very_close(x: f64, y: f64) -> bool {
    is_close(x, y, 1E-8)
}

/// Are two numbers close to within an absolute tolerance?
#[inline(always)]
pub fn is_close(x: f64, y: f64, atol: f64) -> bool {
    if (x - y).abs() <= atol {
        return true;
    }

    // println!("\nleft = {}\nright = {}", x, y);
    false
}

pub trait IsVeryClose {
    fn is_very_close(&self, other: &f64) -> bool;
}

impl IsVeryClose for f64 {
    fn is_very_close(&self, other: &f64) -> bool {
        (*self - *other).abs() < 1E-8
    }
}

/// Is a analytical gradient close to a numerical one?
#[cfg(test)]
pub(crate) fn num_and_anal_gradient_are_close(mol: &mut Molecule, ff: &mut dyn Forcefield) -> bool {
    let grad = mol.gradient(ff);
    let num_grad = mol.numerical_gradient(ff);

    for i in 0..grad.len() {
        if !grad[i].is_close_to(&num_grad[i], 1E-6) {
            println!("{:?} not close to {:?}", grad[i], num_grad[i]);
            return false;
        }
    }
    true
}

#[cfg(test)]
pub trait IsPerfectSquare {
    fn is_perfect_square(&self) -> bool;
}

#[cfg(test)]
impl IsPerfectSquare for usize {
    fn is_perfect_square(&self) -> bool {
        *self == ((*self as f64).sqrt() as usize).pow(2)
    }
}

#[cfg(test)]
pub(crate) fn remove_file_or_panic(filename: &str) {
    std::fs::remove_file(filename).expect("Failed to remove file")
}
#[cfg(test)]
pub(crate) fn print_methane_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "5\n\n\
                        C     0.00000   0.00000   0.00000\n\
                        H    -0.65860  -0.85220  -0.30120\n\
                        H    -0.45940   0.97110  -0.28590\n\
                        H     0.08440  -0.02940   1.10060\n\
                        H     1.02910  -0.10990  -0.41250\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_dihydrogen_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "2\n\n\
                        H     0.00000   0.00000   0.00000\n\
                        H     0.77000   0.00000   0.00000\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_water_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "3\n\n\
                        O    -0.00110   0.36310  -0.00000\n\
                        H    -0.82500  -0.18190  -0.00000\n\
                        H     0.82610  -0.18120   0.00000\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_ethene_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "6\n\n\
C         -4.22149        2.30283        0.00000\n\
C         -2.90167        2.47646        0.00000\n\
H         -4.69362        1.69525        0.76527\n\
H         -4.83470        2.76763       -0.76527\n\
H         -2.28846        2.01166        0.76527\n\
H         -2.42954        3.08405       -0.76527\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_h2coh_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "5\n\n\
C         -4.15592        2.32790       -0.02306\n\
O         -2.94565        2.49183       -0.02968\n\
H         -4.61267        1.71790        0.74844\n\
H         -4.78306        2.78690       -0.78277\n\
H         -2.52046        3.03273       -0.71007\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_pdcl2nh3h_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "11\n\n
Pd        -1.47702       -0.03376        0.00040\n\
Cl        -1.30219       -2.34599       -0.00127\n\
Cl        -1.65184        2.27864       -0.00187\n\
N         -3.48953       -0.16197        0.00064\n\
N          0.53554        0.09400        0.00017\n\
H         -3.86878        0.39258        0.80030\n\
H         -3.79567       -1.15526        0.09806\n\
H         -3.85655        0.22609       -0.89684\n\
H          0.84241        1.09170        0.00023\n\
H          0.90826       -0.38020       -0.85254\n\
H          0.90845       -0.38034        0.85272\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_aume2_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "9\n\n
Au        -0.63866       -1.10981        0.00000\n\
C         -2.64047       -0.85699       -0.00000\n\
H         -3.11649       -1.67438        0.57968\n\
H         -3.01553       -0.87670       -1.04387\n\
H         -2.88947        0.11933        0.46418\n\
C          1.36324       -1.36193        0.00000\n\
H          1.59955       -2.44588        0.00001\n\
H          1.79519       -0.88910        0.90584\n\
H          1.79519       -0.88911       -0.90583\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_h2o2_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "4\n\n\
H         -3.80272        0.11331       -0.29090\n\
O         -3.75673        1.06948       -0.03303\n\
O         -2.46930        1.33955        0.03004\n\
H         -2.46930        1.31540        0.99712\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_ph3_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "4\n\n\
  P    0.0056  -0.0028   0.5805\n\
  H    1.0584  -0.5395  -0.2058\n\
  H   -0.9999  -0.6465  -0.1872\n\
  H   -0.0642   1.1889  -0.1874\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_benzene_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "12\n\n
C         -0.50959        1.40925        0.00000\n\
C         -0.50172       -0.01013        0.00000\n\
C          0.71485        2.13236        0.00000\n\
C          1.91595        1.39203       -0.00000\n\
C          1.91230       -0.00722       -0.00000\n\
C          0.71648       -0.73162        0.00000\n\
H          0.75399       -1.83203        0.00000\n\
H         -1.42864       -0.58129        0.00000\n\
H         -1.47287        1.94140        0.00000\n\
H          0.78724        3.23589        0.00000\n\
H          2.86816        1.90649       -0.00000\n\
H          2.85791       -0.54513       -0.00000\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_c2h2_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "4\n\n
C         -3.99338        0.82887        0.00000\n\
C         -2.79056        0.90688        0.00000\n\
H         -5.04954        0.76038       -0.00000\n\
H         -1.73439        0.97537        0.00000\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_ph_ene_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "14\n\n
C         -2.45353       -0.18354        0.00485\n\
C         -2.40963       -1.58150       -0.02148\n\
C         -1.26583        0.55451        0.01877\n\
C         -0.03251       -0.10452        0.00639\n\
C          0.01489       -1.50241       -0.01994\n\
C         -1.17447       -2.24441       -0.03394\n\
C         -1.12784       -3.68294       -0.06108\n\
C         -1.08938       -4.88772       -0.08386\n\
H         -1.05597       -5.94527       -0.10388\n\
H          0.97559       -2.00320       -0.02935\n\
H         -1.30118        1.63608        0.03914\n\
H          0.88620        0.46764        0.01719\n\
H         -3.40763        0.32742        0.01445\n\
H         -3.33559       -2.14389       -0.03210\n",
    )
    .expect(filename)
}
#[cfg(test)]
pub(crate) fn print_distorted_water_xyz_file(filename: &str) {
    std::fs::write(
        filename,
        "3\n\n
O     2.208695   1.357633   2.375361\n\
H     2.434634   2.077241   2.871151\n\
H     1.697517   0.615572   2.323087\n",
    )
    .expect(filename)
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

    #[test]
    fn test_is_perfect_square() {
        assert!((9 as usize).is_perfect_square());
        assert!(!(10 as usize).is_perfect_square());
        assert!((1 as usize).is_perfect_square());
        assert!((144 as usize).is_perfect_square());
    }
}
