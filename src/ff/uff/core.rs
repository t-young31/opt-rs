use crate::coordinates::{Point, Vector3D};
use crate::ff::angles::{HarmonicAngleTypeA, HarmonicAngleTypeB};
use crate::ff::bonds::HarmonicBond;
use crate::ff::dihedrals::{InversionDihedral, TorsionalDihedral};
use crate::ff::forcefield::EnergyFunction;
use crate::ff::forcefield::Forcefield;
use crate::ff::nonbonded::LennardJones12x6;
use crate::ff::uff::atom_types::ATOM_TYPES;
use crate::ff::uff::atom_typing::UFFAtomType;
use crate::ff::uff::dihedral_bond::DihedralBond;
use crate::ff::uff::inversion_centers::INVERSION_CENTERS;
use crate::molecule::Molecule;
use crate::pairs::AtomPair;
use std::cmp::Ordering;
use std::collections::HashMap;

#[derive(Default)]
pub struct UFF {
    atom_types: Vec<UFFAtomType>,
    energy_functions: Vec<Box<dyn EnergyFunction>>,

    r0_cache: HashMap<AtomPair, f64>, // Cache of equilibrium bond distances (r0)
    energy: f64,
    gradient: Vec<Vector3D>,
}

impl UFF {
    /// Generate atom types for each atom in a molecule, defined by the connectivity and charge on
    /// each atom. As each atom may have multiple matches this function should set the best match
    /// to each atom
    fn set_atom_types(&mut self, molecule: &Molecule) {
        let mut match_qualities: [f64; ATOM_TYPES.len()] = [0.; ATOM_TYPES.len()];

        for atom in molecule.atoms().iter() {
            for (i, atom_type) in ATOM_TYPES.iter().enumerate() {
                match_qualities[i] = atom_type.match_quality(atom, molecule);
            }

            let best_match = match_qualities
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                .map(|(index, _)| index)
                .expect("Failed to find a best match");

            let mut atom_type = ATOM_TYPES[best_match].clone();
            atom_type.set_coordination_environment(atom);
            // println!("{:?}", atom_type);

            self.atom_types.push(atom_type);
        }
    }

    /// Set a zero gradient for all atoms that this force field applies to
    fn set_zero_gradient(&mut self, molecule: &Molecule) {
        self.gradient = molecule.blank_gradient();
    }

    /// Add a bond stretching term to the FF for all bonds present in a molecule
    fn add_bond_stretches(&mut self, molecule: &Molecule) {
        for bond in molecule.bonds() {
            let i = bond.pair.i;
            let j = bond.pair.j;
            let r0 = self.r0(i, j, bond.order.value());
            self.r0_cache.insert(bond.pair.clone(), r0);

            let k_ij = self.k_ij(i, j, r0);

            self.energy_functions
                .push(Box::new(HarmonicBond { i, j, r0, k_ij }));
        }
    }

    /// Equilibrium bond distance (Å) between two atoms indexed with i and j  [eqn. 2]
    fn r0(&self, i: usize, j: usize, bo: f64) -> f64 {
        self.atom_types[i].r + self.atom_types[j].r + self.r_bo(i, j, bo) + self.r_en(i, j)
    }

    /// Bond order correction to the bond distance between two atoms  [eqn. 3]
    fn r_bo(&self, i: usize, j: usize, bo: f64) -> f64 {
        -0.1332 * (self.atom_types[i].r + self.atom_types[j].r) * bo.ln()
    }

    /// Electronegativity correction to the bond distance between two atoms  [eqn. 4]
    fn r_en(&self, i: usize, j: usize) -> f64 {
        let r_i = self.atom_types[i].r;
        let r_j = self.atom_types[j].r;

        let chi_i = self.atom_types[i].gmp_electronegativity();
        let chi_j = self.atom_types[j].gmp_electronegativity();

        r_i * r_j * ((chi_i.sqrt() - chi_j.sqrt()).powi(2) / (chi_i * r_i + chi_j * r_j))
    }

    /// Force constant for a harmonic bond [eqn. 6]
    fn k_ij(&self, i: usize, j: usize, r0: f64) -> f64 {
        664.12 * (self.atom_types[i].z_eff * self.atom_types[j].z_eff) / r0.powi(3)
    }

    /// Add a term for an angle bend
    fn add_angle_bends(&mut self, molecule: &Molecule) {
        for angle in molecule.angles() {
            let i = angle.i;
            let j = angle.j;
            let k = angle.k;

            let atom_type = &self.atom_types[j];
            let k_ijk = self.k_ijk(i, j, k);

            match atom_type.bend_type() {
                'A' => {
                    self.energy_functions.push(Box::new(HarmonicAngleTypeA {
                        i,
                        j,
                        k,
                        k_ijk,
                        n: atom_type.bend_n(),
                    }));
                }
                'B' => {
                    let c2 = 1. / (4. * atom_type.theta.sin().powi(2));
                    let c1 = -4. * c2 * atom_type.theta.cos();
                    let c0 = c2 * (2. * atom_type.theta.cos().powi(2) + 1.);

                    self.energy_functions.push(Box::new(HarmonicAngleTypeB {
                        i,
                        j,
                        k,
                        k_ijk,
                        c0,
                        c1,
                        c2,
                    }));
                }
                _ => panic!("Cannot match on unsupported atom type"),
            }
        }
    }

    /// Bend force constant
    fn k_ijk(&self, i: usize, j: usize, k: usize) -> f64 {
        let theta0 = self.atom_types[j].theta;

        let r0_ij = self
            .r0_cache
            .get(&AtomPair { i, j })
            .expect("Failed to get r_ij");
        let r0_jk = self
            .r0_cache
            .get(&AtomPair { i: j, j: k })
            .expect("Failed to get r_jk");

        // Cosine rule to calculate the equilibrium distance between the (probably) non-bonded atoms
        let r0_ik = (r0_ij.powi(2) + r0_jk.powi(2) - 2. * r0_ij * r0_jk * theta0.cos()).sqrt();

        let z_i = self.atom_types[i].z_eff;
        let z_k = self.atom_types[k].z_eff;

        let beta = 664.12 / (r0_ij * r0_jk);

        beta * ((z_i * z_k) / r0_ik.powi(5))
            * (r0_ij * r0_jk * (1. - theta0.cos().powi(2)) - r0_ik.powi(2) * theta0.cos())
    }

    /// Add dihedral terms between quadruples bonded atoms in a sequence
    fn add_dihedral_torsions(&mut self, molecule: &Molecule) {
        for dihedral in molecule.proper_dihedrals().iter() {
            let i = dihedral.i;
            let j = dihedral.j;
            let k = dihedral.k;
            let l = dihedral.l;

            let central_bond =
                DihedralBond::from_atom_types(&self.atom_types[j], &self.atom_types[k]);

            if !central_bond.contains_only_main_group_elements() {
                // Only main group atoms have non-zero torsional potentials
                continue;
            }

            if molecule.angle_is_close_to_linear(i, j, k)
                || molecule.angle_is_close_to_linear(j, k, l)
            {
                // Bad things happen when dihedrals with close to linear angles are included
                continue;
            }

            self.energy_functions.push(Box::new(TorsionalDihedral {
                i,
                j,
                k,
                l,
                phi0: central_bond.phi0,
                n_phi: central_bond.n,
                v_phi: central_bond.v,
            }));
        }
    }

    /// Add terms for distortions from (trigonal) pyramidal geometries
    fn add_dihedral_inversions(&mut self, molecule: &Molecule) {
        for improper in molecule.improper_dihedrals().iter() {
            let mut dihedral = InversionDihedral {
                c: improper.c,
                i: improper.i,
                j: improper.j,
                k: improper.k,
                c0: 0.,
                c1: 0.,
                c2: 0.,
                k_cijk: 0.,
            };

            match self.atom_types[improper.c].name {
                "C_2" | "C_R" => {
                    dihedral.c0 = 1.;
                    dihedral.c1 = -1.;
                    dihedral.c2 = 0.;

                    match self.is_bonded_to_o_2(improper.c, molecule) {
                        true => {
                            dihedral.k_cijk = 50.;
                        }
                        false => {
                            dihedral.k_cijk = 6.;
                        }
                    }
                }

                "O_2" | "S_2" => {
                    todo!();
                }

                x => {
                    let centre = INVERSION_CENTERS.iter().find(|y| x.contains(y.name));
                    match centre {
                        Some(y) => {
                            dihedral.c0 = y.c0;
                            dihedral.c1 = y.c1;
                            dihedral.c2 = y.c2;
                            dihedral.k_cijk = y.k;
                        }
                        None => {
                            continue;
                        }
                    }
                }
            }

            self.energy_functions.push(Box::new(dihedral));
        }
    }

    /// Is an atom bonded to an O_2 atom type?
    fn is_bonded_to_o_2(&self, atom_idx: usize, molecule: &Molecule) -> bool {
        let neighbours = &molecule.atoms()[atom_idx].bonded_neighbours;
        neighbours
            .iter()
            .any(|i| self.atom_types[*i].name == "O_2" && i != &atom_idx)
    }

    /// Add van der Waals terms for repulsion and dispersive attraction
    fn add_vdw(&mut self, molecule: &Molecule) {
        // TODO: Should 1-3 exclusions be included? Not included in the original paper
        // but are included in: https://doi.org/10.1016/j.fluid.2006.07.014

        for nb_pair in molecule.non_bonded_pairs.iter() {
            let type_i = &self.atom_types[nb_pair.pair.i];
            let type_j = &self.atom_types[nb_pair.pair.j];

            self.energy_functions.push(Box::new(LennardJones12x6 {
                i: nb_pair.pair.i,
                j: nb_pair.pair.j,
                d: (type_i.d * type_j.d).sqrt(),
                sigma: (type_i.r * type_j.r).sqrt(),
            }));
        }
    }
}

impl Forcefield for UFF {
    /// Create a new, bespoke, UFF forcefield for a molecule
    fn new(molecule: &Molecule) -> Self {
        let mut ff = UFF::default();

        ff.set_zero_gradient(molecule);
        ff.set_atom_types(molecule);
        ff.add_bond_stretches(molecule);
        ff.add_angle_bends(molecule);
        ff.add_dihedral_torsions(molecule);
        ff.add_dihedral_inversions(molecule);
        ff.add_vdw(molecule);

        ff
    }

    /// Evaluate the total energy of the sustem
    fn energy(&mut self, coordinates: &[Point]) -> f64 {
        self.energy = 0.0;

        for function in self.energy_functions.iter() {
            self.energy += function.energy(coordinates);
        }

        self.energy
    }

    /// Evaluate the gradient {dE/dX_ik, ...} for atom i and Cartesian component k
    fn gradient(&mut self, coordinates: &[Point]) -> &Vec<Vector3D> {
        self.gradient.iter_mut().for_each(|v| v.zero());

        for function in self.energy_functions.iter() {
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
mod tests {
    use super::*;
    use crate::connectivity::bonds::{Bond, BondOrder};
    use crate::ff::uff::atom_typing::Hybridisation;
    use crate::pairs::distance;
    use crate::utils::*;
    use std::f64::consts::PI;

    /// Generate a valid H2 molecule
    fn h2() -> Molecule {
        let mut mol = Molecule::from_atomic_symbols(&["H", "H"]);
        mol.coordinates[0].x = 0.77;
        mol.add_bonds();

        mol
    }

    /// Find a UFFAtomType from a string label of it
    fn atom_type(name: &str) -> UFFAtomType {
        ATOM_TYPES
            .iter()
            .filter(|t| t.name == name)
            .next()
            .unwrap()
            .clone()
    }

    /// Given a H atom then the optimal atom type should be selected
    #[test]
    fn test_h_atom_ff() {
        let mol = Molecule::from_atomic_symbols(&["H"]);

        let uff = UFF::new(&mol);
        assert_eq!(uff.atom_types.len(), 1);

        let atom_type = uff.atom_types[0].clone();
        assert_eq!(
            &atom_type,
            ATOM_TYPES.get(0).expect("Failed to get H atom type")
        );
    }

    /// Given a forcefield initialised for a new molecule
    #[test]
    fn test_initial_zero_gradient() {
        let uff = UFF::new(&h2());

        for i in 0..2 {
            assert!(is_very_close(uff.gradient[i].x, 0.0));
            assert!(is_very_close(uff.gradient[i].y, 0.0));
            assert!(is_very_close(uff.gradient[i].z, 0.0));
        }
    }

    /// Given a H2 molecule then the atom types should be H_ and an energy be non-zero
    #[test]
    fn test_h2_ff() {
        let h2 = h2();

        let mut uff = UFF::new(&h2);
        assert_eq!(uff.atom_types.len(), 2);

        for atom_type in uff.atom_types.iter() {
            assert_eq!(atom_type.name, "H_");
        }

        // Should have a single bond
        assert_eq!(h2.bonds().iter().count(), 1);

        // Energy is non-zero
        assert!(!is_close(uff.energy(&h2.coordinates), 0.0, 1E-5))
    }

    /// Ensure the electronegativity correction on a bond length is close to that in the paper
    #[test]
    fn test_r_en_correction() {
        let mut uff = UFF::new(&Molecule::blank());

        for atom_type_str in ["Si3", "O_3_z"] {
            uff.atom_types.push(
                ATOM_TYPES
                    .iter()
                    .find(|a| a.name == atom_type_str)
                    .unwrap()
                    .clone(),
            )
        }

        assert!(is_close(uff.r_en(0, 1), 0.0533, 1E-2));
    }

    /// Test numerical vs analytical gradient evaluation for H2
    #[test]
    fn test_num_vs_anal_grad_h2() {
        let mut h2 = h2();
        let mut uff = UFF::new(&h2);

        assert!(num_and_anal_gradient_are_close(&mut h2, &mut uff));
    }

    /// Test that numerical gradients doesn't shift atoms
    #[test]
    fn test_numerical_gradient_leaves_coordinates_unchanged() {
        let mut h2 = h2();
        let mut uff = UFF::new(&h2);
        let init_energy = h2.energy(&mut uff);

        let _ = h2.numerical_gradient(&mut uff);
        assert!(is_very_close(h2.energy(&mut uff), init_energy));
    }

    /// The UFF force field should have a minimum in the energy of H2 around 0.741 Å
    #[test]
    fn test_minimum_in_h2_is_near_true_bond_length() {
        let mut min_energy = f64::MAX;
        let mut min_r = f64::MAX;

        let mut h2 = h2();
        let mut uff = UFF::new(&h2);

        assert_eq!(
            h2.bonds()
                .get(&Bond::from_atom_indices(0, 1))
                .unwrap()
                .order,
            BondOrder::Single
        );

        // Ensure the minimum energy bond length is the same as the one obtained from avogadro
        assert!(is_close(
            uff.r0(0, 1, BondOrder::Single.value()),
            0.708,
            1E-3
        ));

        for i in 0..100 {
            let r = 0.5 + (i as f64 * 0.01); // 0.5 + x Å

            h2.coordinates[0].x = r;

            assert!(is_very_close(distance(0, 1, &h2.coordinates), r));

            let energy = uff.energy(&h2.coordinates);

            if energy < min_energy {
                min_energy = energy;
                min_r = r;
            }
        }

        assert!(is_close(min_r, 0.741, 1E-1))
    }

    /// Ensure that the equilibrium angle π/n is undefined for a general bent atom
    #[test]
    fn test_bend_n_for_bent_angle_is_undefined() {
        assert!(is_very_close(atom_type("O_2").bend_n(), 0.));
    }

    #[test]
    fn test_angle_bend_amide() {
        let mut ff = UFF::default();

        ff.atom_types.push(atom_type("C_R"));
        ff.atom_types.push(atom_type("N_R"));
        ff.atom_types.push(atom_type("C_3"));

        let mut mol = Molecule::from_atomic_symbols(&["C", "N", "C"]);
        mol.coordinates[0].x = -1.0;
        mol.coordinates[2].x = 1.0;
        mol.add_bonds();
        mol.add_angles();

        assert_eq!(mol.bonds().len(), 2); // Two C-N bonds
        assert_eq!(mol.angles().len(), 1); // Single amide angle

        ff.add_bond_stretches(&mol);
        ff.add_angle_bends(&mol);

        let bend = ff
            .energy_functions
            .iter()
            .find(|f| f.as_ref().involves_idxs(Vec::from([0, 1, 2])))
            .unwrap()
            .as_ref();

        // Ensure the force constant is close that in the UFF paper (kcal mol-1 rad-2)
        assert!(is_close(bend.force_constant(), 105.5, 15.));
    }

    /// Test numerical vs analytical gradient evaluation
    #[test]
    fn test_num_vs_anal_grad_h2o() {
        let filename = "water_tnvagh.xyz";
        print_water_xyz_file(filename);
        let mut h2o = Molecule::from_xyz_file(filename);
        let mut uff = UFF::new(&h2o);

        assert!(num_and_anal_gradient_are_close(&mut h2o, &mut uff));

        remove_file_or_panic(filename);
    }
    #[test]
    fn test_num_vs_anal_grad_linear_metal_complex() {
        let filename = "water_tnvaglmc.xyz";
        print_aume2_xyz_file(filename);
        let mut au_me2 = Molecule::from_xyz_file(filename);
        let mut uff = UFF::new(&au_me2);

        assert!(num_and_anal_gradient_are_close(&mut au_me2, &mut uff));

        remove_file_or_panic(filename);
    }
    #[test]
    fn test_num_vs_anal_grad_h2o2() {
        let filename = "h2o2_tnvagh.xyz";
        print_h2o2_xyz_file(filename);
        let mut h2o2 = Molecule::from_xyz_file(filename);

        let o_o_bond = Bond {
            pair: AtomPair { i: 1, j: 2 },
            order: BondOrder::Single,
        };
        assert!(h2o2.bonds().contains(&o_o_bond));
        assert_eq!(h2o2.bonds().len(), 3);

        assert_eq!(h2o2.proper_dihedrals().len(), 1); // Should have a single dihedral

        let mut uff = UFF::new(&h2o2);

        let n_dihedral_terms = uff
            .energy_functions
            .iter()
            .filter(|f| f.involves_idxs(Vec::from([0, 1, 2, 3])))
            .count();

        assert_eq!(n_dihedral_terms, 1);
        assert!(num_and_anal_gradient_are_close(&mut h2o2, &mut uff));

        remove_file_or_panic(filename);
    }
    #[test]
    fn test_num_vs_anal_grad_ethene() {
        let filename = "ethene_tnvage.xyz";
        print_ethene_xyz_file(filename);
        let mut ethene = Molecule::from_xyz_file(filename);
        let mut uff = UFF::new(&ethene);

        assert!(num_and_anal_gradient_are_close(&mut ethene, &mut uff));

        remove_file_or_panic(filename);
    }
    #[test]
    fn test_num_vs_anal_grad_ph3() {
        let filename = "ph3_tnvagp.xyz";
        print_ph3_xyz_file(filename);
        let mut mol = Molecule::from_xyz_file(filename);
        assert_eq!(mol.improper_dihedrals().len(), 1);

        let mut uff = UFF::new(&mol);
        assert_eq!(uff.atom_types[0].name, "P_3+3");

        assert_eq!(
            uff.energy_functions
                .iter()
                .filter(|f| f.involves_idxs(Vec::from([0, 1, 2, 3])))
                .count(),
            1
        );

        assert!(num_and_anal_gradient_are_close(&mut mol, &mut uff));

        remove_file_or_panic(filename);
    }

    /// Ensure that two inversion centers with non-zero force constants have been added
    /// to hinder pyramidalisation of the carbon atoms
    #[test]
    fn test_two_inversion_centres_added_for_ethene() {
        let filename = "ethene_tticafe.xyz";
        print_ethene_xyz_file(filename);
        let uff = UFF::new(&Molecule::from_xyz_file(filename));

        let functions: Vec<&Box<dyn EnergyFunction>> = uff
            .energy_functions
            .iter()
            .filter(|f| f.involves_idxs(Vec::from([0, 1, 2, 3])))
            .collect();

        assert_eq!(functions.len(), 1);
        assert!(is_very_close(functions[0].force_constant(), 6.));

        remove_file_or_panic(filename);
    }

    #[test]
    fn test_valency_assigned_for_carbon_in_methane() {
        let filename = "methane_tvafcim.xyz";
        print_methane_xyz_file(filename);

        let mol = Molecule::from_xyz_file(filename);
        let ff = UFF::new(&mol);

        assert_eq!(ff.atom_types[0].atomic_symbol, "C");
        assert_eq!(ff.atom_types[0].hybridisation(), Hybridisation::SP3);

        // Hydrogen don't have any hybridisation
        for i in 1..mol.num_atoms() {
            assert_eq!(ff.atom_types[i].hybridisation(), Hybridisation::None);
        }

        remove_file_or_panic(filename);
    }

    #[test]
    fn test_benzene_has_aromatic_carbons() {
        let filename = "benzene_tbhac.xyz";
        print_benzene_xyz_file(filename);

        let ff = UFF::new(&Molecule::from_xyz_file(filename));

        for atom_type in ff.atom_types.iter() {
            if atom_type.atomic_symbol == "C" {
                assert_eq!(atom_type.name, "C_R");
            }
        }

        remove_file_or_panic(filename);
    }

    #[test]
    fn test_phenyl_acetylene_has_correct_topology() {
        let filename = "phenylacetylene_tpahct.xyz";
        print_ph_ene_xyz_file(filename);
        let mol = Molecule::from_xyz_file(filename);
        remove_file_or_panic(filename);

        let num_triple_bonds = mol
            .bonds()
            .iter()
            .filter(|b| b.order == BondOrder::Triple)
            .count();
        assert_eq!(num_triple_bonds, 1);

        let ff = UFF::new(&mol);

        let num_linear_carbons = ff.atom_types.iter().filter(|x| x.name == "C_1").count();
        assert_eq!(num_linear_carbons, 2);

        let num_linear_carbon_bends = ff
            .atom_types
            .iter()
            .filter(|x| x.atomic_symbol == "C" && is_close(x.theta, PI, 1E-3))
            .count();
        assert_eq!(num_linear_carbon_bends, 2);
    }

    #[test]
    fn test_atom_typing_water() {
        let filename = "water_tatw.xyz";
        print_distorted_water_xyz_file(filename);
        let mol = Molecule::from_xyz_file(filename);
        remove_file_or_panic(filename);

        let ff = UFF::new(&mol);
        let o_type = ff
            .atom_types
            .iter()
            .filter(|a| a.atomic_symbol == "O")
            .next()
            .unwrap();

        assert_eq!(o_type.name, "O_3");
    }
}
