#[derive(Default, Debug, Hash)]
pub struct Angle{
    pub i: usize,
    pub j: usize,
    pub k: usize
}

impl PartialEq for Angle {
    fn eq(&self, other: &Self) -> bool {
        self.j == other.j
            && (self.i == other.i && self.k == other.k || self.i == other.k && self.k == other.i)
    }
}

impl Eq for Angle {}