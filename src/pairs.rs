use crate::atoms::Atom;
use crate::coordinates::CartesianCoordinate;


#[inline(always)]
pub fn distance(i: usize, j: usize, x: &Vec<CartesianCoordinate>) -> f64{
    ((x[i].x - x[j].x).powi(2) + (x[i].y - x[j].y).powi(2) + (x[i].z - x[j].z).powi(2)).sqrt()
}


#[derive(Default, Hash, Debug, Clone)]
pub(crate) struct AtomPair{
    pub i: usize,
    pub j: usize
}

impl PartialEq for AtomPair {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i && self.j == other.j || self.i == other.j && self.j == other.i
    }
}

impl Eq for AtomPair {}


#[derive(Default, Hash, Debug, Clone)]
pub(crate) struct NBPair{
    pub(crate) pair: AtomPair
}

impl NBPair {

    /// Construct a non-bonded pair from teo atoms
    pub fn from_atoms(atom_i: &Atom, atom_j: &Atom) -> Self{

        if atom_i == atom_j{
            panic!("Cannot create a non-bonded pair between identical atoms")
        }

        NBPair{pair: AtomPair{i: atom_i.idx, j: atom_j.idx}}
    }
}

impl PartialEq for NBPair {
    fn eq(&self, other: &Self) -> bool { self.pair == other.pair }
}

impl Eq for NBPair {}