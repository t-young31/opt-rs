use std::cmp::Ordering;
use crate::{Forcefield, Molecule};
use crate::atoms::CartesianCoordinate;
use crate::ff::bonds::HarmonicBond;
use crate::ff::forcefield::EnergyFunction;
use crate::ff::uff::atom_typing::UFFAtomType;
use crate::ff::uff::atom_types::ATOM_TYPES;


#[derive(Default)]
pub(crate) struct UFF{

    atom_types:       Vec<UFFAtomType>,
    energy_functions: Vec<Box<dyn EnergyFunction>>,

    energy:           f64,
    gradient:         Vec<CartesianCoordinate>,
}

impl UFF {

    /// Set a zero gradient for all atoms that this force field applies to
    fn set_zero_gradient(&mut self, molecule: &Molecule){

        for _ in 0..molecule.num_atoms(){
            self.gradient.push(CartesianCoordinate::default());
        }
    }

    /// Zero an existing gradient vector
    fn zero_gradient(&mut self){

        for vec in self.gradient.iter_mut(){
            vec.x = 0.0;
            vec.y = 0.0;
            vec.z = 0.0;
        }
    }

    /// Add a bond stretching term to the FF for all bonds present in a molecule
    fn add_bond_stretch(&mut self, molecule: &Molecule){

        for bond in molecule.bonds(){

            let i = bond.pair.i;
            let j = bond.pair.j;
            let r0 = self.r0(i, j, bond.order.value());
            let k = self.k_bond(i, j, r0);

            self.energy_functions.push(Box::new(HarmonicBond{i, j, r0, k}));
        }
    }

    /// Equilibrium bond distance (Å) between two atoms indexed with i and j  [eqn. 2]
    fn r0(&self, i: usize, j: usize, bo: f64) -> f64{
        self.atom_types[i].r + self.atom_types[j].r + self.r_bo(i, j, bo) + self.r_en(i, j)
    }

    /// Bond order correction to the bond distance between two atoms  [eqn. 3]
    fn r_bo(&self, i: usize, j: usize, bo: f64) -> f64{
        -0.1332 * (self.atom_types[i].r + self.atom_types[j].r) * bo.ln()
    }

    /// Electronegativity correction to the bond distance between two atoms  [eqn. 4]
    fn r_en(&self, i: usize, j: usize) -> f64{
        let r_i = self.atom_types[i].r;
        let r_j = self.atom_types[j].r;

        let chi_i = self.atom_types[i].gmp_electronegativity();
        let chi_j = self.atom_types[j].gmp_electronegativity();

        r_i * r_j * ((chi_i.sqrt() - chi_j.sqrt()).powi(2)
                     / (chi_i*r_i + chi_j*r_j))
    }

    /// Force constant for a harmonic bond [eqn. 6]
    fn k_bond(&self, i: usize, j: usize, r0: f64) -> f64{
        664.12 * (self.atom_types[i].z_eff * self.atom_types[j].z_eff) / r0.powi(3)
    }

}


impl Forcefield for UFF {

    /// Create a new, bespoke, forcefield for a molecule by setting atom types which
    fn new(molecule: &Molecule) -> Self{

        let mut ff = UFF::default();
        ff.set_atom_types(molecule);
        ff.set_zero_gradient(molecule);
        ff.add_bond_stretch(molecule);     // E_R

        ff
    }

    /// Generate atom types for each atom in a molecule, defined by the connectivity and charge on
    /// each atom. As each atom may have multiple matches this function should set the best match
    /// to each atom
    fn set_atom_types(&mut self, molecule: &Molecule) {

        let mut match_qualities: [i64; ATOM_TYPES.len()] = [0; ATOM_TYPES.len()];

        for atom in molecule.atoms().iter(){

            for (i, atom_type) in ATOM_TYPES.iter().enumerate(){
                match_qualities[i] = atom_type.match_quality(atom, molecule);
            }

            let best_match = match_qualities
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                .map(|(index, _)| index)
                .expect("Failed to find a best match");

            self.atom_types.push(ATOM_TYPES[best_match].clone());
        }
    }

    /// Evaluate the total energy of the sustem
    fn energy(&mut self, coordinates: &Vec<CartesianCoordinate>) -> f64{

        self.energy = 0.0;

        for function in self.energy_functions.iter(){
            self.energy += function.energy(coordinates);
        }

        self.energy
    }

    /// Evaluate the gradient {dE/dX_ik, ...} for atom i and Cartesian component k
    fn gradient(&mut self, coordinates: &Vec<CartesianCoordinate>) -> &Vec<CartesianCoordinate>{

        self.zero_gradient();

        for function in self.energy_functions.iter(){
            function.add_gradient(coordinates, &mut self.gradient);
        }

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
mod tests{
    use super::*;
    use crate::utils::*;

    /// Generate a valid H2 molecule
    fn h2() -> Molecule{

        let mut mol = Molecule::from_atomic_symbols(&["H", "H"]);
        mol.coordinates[0].x = 0.77;
        mol.add_bonds();

        mol
    }

    /// Given a H atom then the optimal atom type should be selected
    #[test]
    fn test_h_atom_ff(){

        let mol = Molecule::from_atomic_symbols(&["H"]);

        let uff = UFF::new(&mol);
        assert_eq!(uff.atom_types.len(), 1);

        let atom_type = uff.atom_types[0].clone();
        assert_eq!(&atom_type, ATOM_TYPES.get(0).expect("Failed to get H atom type"));
    }

    /// Given a forcefield initialised for a new molecule
    #[test]
    fn test_initial_zero_gradient(){

        let uff = UFF::new(&h2());

        for i in 0..2{
            assert!(is_very_close(uff.gradient[i].x, 0.0));
            assert!(is_very_close(uff.gradient[i].y, 0.0));
            assert!(is_very_close(uff.gradient[i].z, 0.0));
        }
    }

    /// Given a H2 molecule then the atom types should be H_ and an energy be non-zero
    #[test]
    fn test_h2_ff(){

        let h2 = h2();

        let mut uff = UFF::new(&h2);
        assert_eq!(uff.atom_types.len(), 2);

        for atom_type in uff.atom_types.iter(){
            assert_eq!(atom_type.name, "H_");
        }

        // Should have a single bond
        println!("{:?}", h2.bonds());
        assert_eq!(h2.bonds().iter().count(), 1);

        // Energy is non-zero
        assert!(!is_close(uff.energy(&h2.coordinates), 0.0, 1E-5))
    }

    /// Ensure the electronegativity correction on a bond length is close to that in the paper
    #[test]
    fn test_r_en_correction(){

        let mut uff = UFF::new(&Molecule::blank());

        for atom_type_str in ["Si3", "O_3_z"]{

            uff.atom_types.push(ATOM_TYPES.iter()
                .find(|a| a.name == atom_type_str)
                .unwrap()
                .clone())
        }

        assert!(is_close(uff.r_en(0, 1), 0.0533, 1E-2));
    }

    /// Test numerical vs analytical gradient evaluation for H2
    #[test]
    fn test_num_vs_anal_grad_h2(){

        let mut h2 = h2();
        let mut uff = UFF::new(&h2);

        let grad = h2.gradient(&mut uff);
        let num_grad = h2.numerical_gradient(&mut uff);


        for i in 0..h2.num_atoms(){
            for k in 0..3{
                assert!(is_close(grad[i][k], num_grad[i][k], 1E-6));
            }
        }
    }

    /// Test that numerical gradients doesn't shift atoms
    #[test]
    fn test_numerical_gradient_leaves_coordinates_unchanged(){

        let mut h2 = h2();
        let mut uff = UFF::new(&h2);
        let init_energy = h2.energy(&mut uff);

        let _ = h2.numerical_gradient(&mut uff);
        assert!(is_very_close(h2.energy(&mut uff), init_energy));
    }

}