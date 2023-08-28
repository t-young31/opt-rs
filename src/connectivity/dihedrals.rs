use crate::connectivity::traits::OrderedAtomIndexes;
use std::hash::{Hash, Hasher};

#[derive(Debug)]
pub struct ProperDihedral {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
}

impl OrderedAtomIndexes for ProperDihedral {
    fn ordered(&self) -> Vec<usize> {
        if self.i < self.l {
            vec![self.i, self.j, self.k, self.l]
        } else {
            vec![self.l, self.k, self.j, self.i]
        }
    }
}

impl PartialEq for ProperDihedral {
    fn eq(&self, other: &Self) -> bool {
        self.ordered() == other.ordered()
    }
}

impl Hash for ProperDihedral {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ordered().hash(state);
    }
}

impl Eq for ProperDihedral {}

#[derive(Debug)]
pub struct ImproperDihedral {
    pub c: usize, // Central atom index
    pub i: usize,
    pub j: usize,
    pub k: usize,
}

impl OrderedAtomIndexes for ImproperDihedral {
    fn ordered(&self) -> Vec<usize> {
        let mut vec = vec![self.i, self.j, self.k];
        vec.sort();
        vec.insert(0, self.c);
        vec
    }
}

impl PartialEq for ImproperDihedral {
    fn eq(&self, other: &Self) -> bool {
        self.ordered() == other.ordered()
    }
}

impl Hash for ImproperDihedral {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ordered().hash(state);
    }
}

impl Eq for ImproperDihedral {}
