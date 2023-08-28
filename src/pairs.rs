use crate::atoms::Atom;
use crate::connectivity::traits::OrderedAtomIndexes;
use crate::coordinates::Point;
use std::hash::{Hash, Hasher};

#[inline(always)]
pub fn distance(i: usize, j: usize, x: &[Point]) -> f64 {
    ((x[i].x - x[j].x).powi(2) + (x[i].y - x[j].y).powi(2) + (x[i].z - x[j].z).powi(2)).sqrt()
}

#[derive(Default, Debug, Clone)]
pub(crate) struct AtomPair {
    pub i: usize,
    pub j: usize,
}

impl OrderedAtomIndexes for AtomPair {
    fn ordered(&self) -> Vec<usize> {
        if self.i < self.j {
            vec![self.i, self.j]
        } else {
            vec![self.j, self.i]
        }
    }
}

impl PartialEq for AtomPair {
    fn eq(&self, other: &Self) -> bool {
        self.ordered() == other.ordered()
    }
}

impl Eq for AtomPair {}

impl Hash for AtomPair {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ordered().hash(state);
    }
}

#[derive(Default, Debug, Clone)]
pub(crate) struct NBPair {
    pub(crate) pair: AtomPair,
}

impl NBPair {
    /// Construct a non-bonded pair from two atoms
    pub fn from_atoms(atom_i: &Atom, atom_j: &Atom) -> Self {
        if atom_i == atom_j {
            panic!("Cannot create a non-bonded pair between identical atoms")
        }

        NBPair {
            pair: AtomPair {
                i: atom_i.idx,
                j: atom_j.idx,
            },
        }
    }
}

impl PartialEq for NBPair {
    fn eq(&self, other: &Self) -> bool {
        self.pair == other.pair
    }
}

impl Eq for NBPair {}
