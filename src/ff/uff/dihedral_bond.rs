use std::collections::HashSet;
use std::f64::consts::PI;
use crate::atoms::Atom;
use crate::ff::uff::atom_typing::{UFFAtomType, Hybridisation};

/// The central bond of a dihedral e.g.
///                  l
///                 /
///        j ---- k
///      /
///    i
///
/// then this is the j-k bond
pub struct DihedralBond{
    type_j: UFFAtomType,
    type_k: UFFAtomType,

    n:    f64,    // n_ϕ
    phi0: f64,    // ϕ_0
    v:    f64     // V_ϕ
}


impl DihedralBond{

    /// Construct a central bond within a dihedral from two atoms types
    pub(crate) fn from_atom_types(type_j: &UFFAtomType, type_k: &UFFAtomType) -> Self{
        let mut bond = DihedralBond{type_j: type_j.clone(),
                       type_k: type_k.clone(),
                       n: 0.,
                       phi0: 0.,
                       v: 0.};

        bond.set_n_phi0_and_v();
        bond
    }

    /// Does this bond only contain main group elements?
    pub fn contains_only_main_group_elements(&self) -> bool{
        return self.type_j.is_main_group() && self.type_k.is_main_group()
    }

    /// Does this bond only contain main group elements?
    pub fn contains_only_group_16_elements(&self) -> bool{

        fn is_group_16(atom_type: &UFFAtomType) -> bool{
            Atom::from_atomic_symbol(atom_type.atomic_symbol).group() == 16
        }

        return is_group_16(&self.type_j) && is_group_16(&self.type_k)
    }

    /// Joint hybridisation of both atoms within this bond
    pub fn hybridisation(&self) -> JointHybridisation{
        JointHybridisation::from(&self.type_j.hybridisation(), &self.type_k.hybridisation())
    }

    /// Set the multiplicity (n), equilibrium torsional angle (ϕ) and force constant (V)
    fn set_n_phi0_and_v(&mut self){

        match self.hybridisation() {
            JointHybridisation::SP3SP3 => {
                match self.contains_only_group_16_elements() {
                    true => {self.n = 2.; self.phi0 = PI/3.},
                    false => {self.n =3.; self.phi0 = PI}
                }

                self.v = (self.type_j.v_phi * self.type_k.v_phi).sqrt();
            },
            JointHybridisation::SP2SP3 => {
                todo!()
            },
            JointHybridisation::SP2SP2 => {
                todo!()
            },
            JointHybridisation::None => {}
        }
    }
}


pub enum JointHybridisation{
    SP3SP3,
    SP2SP3,
    SP2SP2,
    None
}

impl JointHybridisation {

    /// Create a joint hybridisation from a pair of hybridisations for two atoms
    pub fn from(x: &Hybridisation, y: &Hybridisation) -> JointHybridisation{

        match (x, y) {
            (Hybridisation::SP3, Hybridisation::SP3) => JointHybridisation::SP3SP3,
            (Hybridisation::SP3, Hybridisation::SP2) => JointHybridisation::SP2SP3,
            (Hybridisation::SP2, Hybridisation::SP3) => JointHybridisation::SP2SP3,
            (Hybridisation::SP2, Hybridisation::SP2) => JointHybridisation::SP2SP2,
            _ => JointHybridisation::None
        }
    }

}

