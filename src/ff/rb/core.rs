use crate::coordinates::{Point, Vector3D};
use crate::ff::bonds::HarmonicBond;
use crate::ff::forcefield::EnergyFunction;
use crate::ff::forcefield::Forcefield;
use crate::ff::nonbonded::RepulsiveInverseDistance;
use crate::molecule::Molecule;

#[derive(Default)]
pub struct RB {
    energy: f64,
    gradient: Vec<Vector3D>,
    energy_functions: Vec<Box<dyn EnergyFunction>>,
    exponent: RepulsiveExponent,
}

/// Exponent in the repulsive component of the potential
#[derive(Clone)]
pub struct RepulsiveExponent {
    pub value: i32,
}

impl Default for RepulsiveExponent {
    fn default() -> Self {
        RepulsiveExponent { value: 2 }
    }
}

impl RB {
    /// Add harmonic bonds between atoms
    fn add_bond_stretches(&mut self, molecule: &Molecule) {
        let atoms = molecule.atoms();

        for bond in molecule.bonds() {
            let i = bond.pair.i;
            let j = bond.pair.j;
            let r0 = atoms[i].covalent_radius() + atoms[j].covalent_radius();

            // yes, k is hard coded
            self.energy_functions.push(Box::new(HarmonicBond {
                i,
                j,
                r0,
                k_ij: 1000.,
            }));
        }
    }

    /// Add repulsion between atoms
    fn add_repulsive_terms(&mut self, molecule: &Molecule) {
        for nb_pair in molecule.non_bonded_pairs.iter() {
            self.energy_functions
                .push(Box::new(RepulsiveInverseDistance {
                    i: nb_pair.pair.i,
                    j: nb_pair.pair.j,
                    c: 10., // hard coded :-(
                    exponent: self.exponent.clone(),
                }));
        }
    }
}

impl Forcefield for RB {
    /// Create a new, bespoke, RB forcefield for a molecule
    fn new(molecule: &Molecule) -> Self {
        let mut ff = RB::default();

        ff.gradient = molecule.blank_gradient();
        ff.add_bond_stretches(molecule);
        ff.add_repulsive_terms(molecule);

        ff
    }

    fn energy(&mut self, coordinates: &[Point]) -> f64 {
        self.energy = self
            .energy_functions
            .iter()
            .map(|f| f.energy(coordinates))
            .sum();
        self.energy
    }

    fn gradient(&mut self, coordinates: &[Point]) -> &Vec<Vector3D> {
        self.gradient.iter_mut().for_each(|v| v.zero());
        self.energy_functions
            .iter()
            .for_each(|f| f.add_gradient(coordinates, &mut self.gradient));
        &self.gradient
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
mod tests {

    use super::*;
    use crate::utils::*;

    #[test]
    fn test_numerical_vs_analytic() {
        let filename = "ph_ene_tnva.xyz";
        print_ph_ene_xyz_file(filename);
        let mut mol = Molecule::from_xyz_file(filename);
        let mut ff = RB::new(&mol);
        remove_file_or_panic(filename);

        assert!(num_and_anal_gradient_are_close(&mut mol, &mut ff));
    }
}
