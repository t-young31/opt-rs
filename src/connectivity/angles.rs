use crate::connectivity::traits::OrderedAtomIndexes;
use std::hash::{Hash, Hasher};

#[derive(Default, Debug)]
pub struct Angle {
    pub i: usize,
    pub j: usize,
    pub k: usize,
}

impl OrderedAtomIndexes for Angle {
    fn ordered(&self) -> Vec<usize> {
        if self.i < self.k {
            vec![self.i, self.j, self.k]
        } else {
            vec![self.k, self.j, self.i]
        }
    }
}

impl PartialEq for Angle {
    fn eq(&self, other: &Self) -> bool {
        self.ordered() == other.ordered()
    }
}

impl Hash for Angle {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ordered().hash(state);
    }
}

impl Eq for Angle {}
