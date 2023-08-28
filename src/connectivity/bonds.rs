use crate::atoms::Atom;
use crate::connectivity::traits::OrderedAtomIndexes;
use crate::pairs::AtomPair;
use crate::utils::IsVeryClose;
use std::hash::{Hash, Hasher};
use std::ops::Index;
use std::slice::Iter;

#[derive(Default, Debug, Clone)]
pub struct Bond {
    pub(crate) pair: AtomPair,
    pub(crate) order: BondOrder,
}

impl Bond {
    /// Construct a bond between two atoms
    pub fn from_atoms(atom_i: &Atom, atom_j: &Atom) -> Self {
        Bond::from_atom_indices(atom_i.idx, atom_j.idx)
    }

    /// Construct a bond between two atom indices. These must be distinct
    pub fn from_atom_indices(i: usize, j: usize) -> Self {
        if i == j {
            panic!("Cannot create a bonded pair between identical atoms")
        }

        Bond {
            pair: AtomPair { i, j },
            order: BondOrder::Single,
        }
    }

    /// Does this bond contain a particular atom?
    pub fn contains(&self, atom: &Atom) -> bool {
        atom.idx == self.pair.i || atom.idx == self.pair.j
    }

    /// Does this bond contain a particular atomic index?
    pub fn contains_index(&self, idx: usize) -> bool {
        idx == self.pair.i || idx == self.pair.j
    }

    /// Given the index of an atom that may, or may not be present in this bond, return the other
    /// atom index in the bond, if present.
    pub fn other(&self, idx: usize) -> Option<usize> {
        if idx == self.pair.i {
            Some(self.pair.j)
        } else if idx == self.pair.j {
            Some(self.pair.i)
        } else {
            None
        }
    }

    /// Set a best guess of a bond order given the connectivity
    pub fn set_possible_bond_order(&mut self, atoms: &[Atom]) {
        self.order = BondOrder::Single; // Default to a single bond

        let atom_i = &atoms[self.pair.i];
        let atom_j = &atoms[self.pair.j];

        if !(atom_i.can_form_multiple_bonds() && atom_j.can_form_multiple_bonds()) {
            return;
        }

        let n = atom_i.num_possible_unpaired_electrons();
        let m = atom_j.num_possible_unpaired_electrons();

        match n.max(m) {
            1 => {
                self.order = BondOrder::Double;
            }
            2 => {
                self.order = BondOrder::Triple;
            }
            _ => {
                self.order = BondOrder::Single;
            }
        }
    }
}

impl OrderedAtomIndexes for Bond {
    fn ordered(&self) -> Vec<usize> {
        self.pair.ordered()
    }
}

impl Hash for Bond {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.pair.hash(state);
        self.ordered().hash(state);
    }
}

impl PartialEq for Bond {
    fn eq(&self, other: &Self) -> bool {
        self.ordered() == other.ordered()
    }
}

impl Eq for Bond {}

impl Index<usize> for Bond {
    type Output = usize;

    /// Index of the atoms that comprise this bond
    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.pair.i,
            1 => &self.pair.j,
            n => panic!("Invalid index for bond: {}. Must be [0, 1]", n),
        }
    }
}

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub enum BondOrder {
    Single,
    Aromatic,
    Double,
    Triple,
    Quadruple,
}

impl Default for BondOrder {
    fn default() -> Self {
        BondOrder::Single
    }
}

impl BondOrder {
    fn iterator() -> Iter<'static, BondOrder> {
        [
            BondOrder::Single,
            BondOrder::Aromatic,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Quadruple,
        ]
        .iter()
    }

    /// Determine the fractional bond order as a floa
    pub(crate) fn value(&self) -> f64 {
        match self {
            BondOrder::Single => 1.0,
            BondOrder::Aromatic => 1.5,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Quadruple => 4.0,
        }
    }

    /// Create a bond order given the value e.g. 1.0 -> BondOrder::Single
    pub(crate) fn from_value(value: &f64) -> Self {
        for bond_order in BondOrder::iterator() {
            if bond_order.value().is_very_close(value) {
                return bond_order.clone();
            }
        }

        panic!("Failed to create a bond order. Unsupported value");
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

    /// Check that bonds are index-able
    #[test]
    fn test_valid_bond_indexing() {
        let bond = Bond::from_atom_indices(3, 5);
        assert_eq!(bond[0], 3);
        assert_eq!(bond[1], 5);
    }

    #[test]
    #[should_panic]
    fn test_invalid_bond_indexing() {
        let bond = Bond::from_atom_indices(3, 5);
        let _ = bond[2]; // No 3rd index in a bond
    }

    /// Check that bond orders default to single bonds
    #[test]
    fn test_default_bond_order() {
        assert_eq!(BondOrder::default(), BondOrder::Single);
    }

    /// Given bond with indices in either order then they should be the same
    #[test]
    fn test_bond_equality() {
        assert_eq!(Bond::from_atom_indices(0, 1), Bond::from_atom_indices(1, 0));
    }

    /// Given bond with indices then it contains one of the atoms
    #[test]
    fn test_bond_contains() {
        let atom_a = Atom::default();
        let mut atom_b = Atom::default();
        atom_b.idx = 1;

        assert!(Bond::from_atoms(&atom_a, &atom_b).contains(&atom_a));
        assert!(Bond::from_atoms(&atom_a, &atom_b).contains(&atom_b));

        let mut atom_c = atom_b.clone();
        atom_c.idx = 3;
        assert!(!Bond::from_atoms(&atom_a, &atom_b).contains(&atom_c));
    }

    /// Test that the other atom index within a bond can be easily accessed
    #[test]
    fn test_bond_other_idx() {
        assert_eq!(Bond::from_atom_indices(0, 1).other(0).unwrap(), 1);
        assert_eq!(Bond::from_atom_indices(0, 1).other(1).unwrap(), 0);
        assert!(Bond::from_atom_indices(0, 1).other(2).is_none());
    }

    /// Test bond order values are ordered correctly
    #[test]
    fn test_bond_order_value() {
        assert!(BondOrder::Single.value() < BondOrder::Aromatic.value());
        assert!(BondOrder::Aromatic.value() < BondOrder::Double.value());
        assert!(BondOrder::Double.value() < BondOrder::Triple.value());
        assert!(BondOrder::Triple.value() < BondOrder::Quadruple.value());
    }

    /// Test that an ordered bond returns atom indices low->high
    #[test]
    fn test_bond_ordering() {
        assert_eq!(Bond::from_atom_indices(0, 1).ordered(), vec![0, 1]);
    }
}
