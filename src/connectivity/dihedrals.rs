use std::hash::{Hash, Hasher};


#[derive(Debug)]
pub struct ProperDihedral{

    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize
}

impl PartialEq for ProperDihedral{
    fn eq(&self, other: &Self) -> bool {
        (self.i == other.i && self.j == other.j && self.k == other.k && self.l == other.l)
            || (self.i == other.l && self.j == other.k && self.k == other.j && self.l == other.i)
    }
}

impl Eq for ProperDihedral {}

impl Hash for ProperDihedral {
    fn hash<H: Hasher>(&self, state: &mut H) {
        Vec::from([self.i, self.j, self.k, self.l]).sort().hash(state);
    }
}

#[derive(Debug)]
pub struct ImproperDihedral {

    pub c: usize,  // Central atom index
    pub i: usize,
    pub j: usize,
    pub k: usize
}

impl PartialEq for ImproperDihedral {
    fn eq(&self, other: &Self) -> bool {
        self.c == other.c
        && [self.i, self.j, self.k].to_vec().sort() == [other.i, other.j, other.k].to_vec().sort()
    }
}

impl Eq for ImproperDihedral {}

impl Hash for ImproperDihedral {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.c.hash(state);
        Vec::from([self.i, self.j, self.k]).sort().hash(state);
    }
}

