use std::hash::{Hash, Hasher};

#[derive(Default, Debug)]
pub struct Angle {
    pub i: usize,
    pub j: usize,
    pub k: usize,
}

impl PartialEq for Angle {
    fn eq(&self, other: &Self) -> bool {
        self.j == other.j
            && (self.i == other.i && self.k == other.k || self.i == other.k && self.k == other.i)
    }
}

impl Hash for Angle {
    fn hash<H: Hasher>(&self, state: &mut H) {
        Vec::from([self.i, self.j, self.k]).sort().hash(state);
    }
}

impl Eq for Angle {}
