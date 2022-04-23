use std::hash::{Hash, Hasher};

#[derive(Default, Debug)]
pub struct Dihedral{
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize
}

impl PartialEq for Dihedral{
    fn eq(&self, other: &Self) -> bool {
        (self.i == other.i && self.j == other.j && self.k == other.k && self.l == other.l)
            || (self.i == other.l && self.j == other.k && self.k == other.j && self.l == other.i)
    }
}

impl Eq for Dihedral {}

impl Hash for Dihedral {
    fn hash<H: Hasher>(&self, state: &mut H) {
        Vec::from([self.i, self.j, self.k, self.l]).sort().hash(state);
    }
}
