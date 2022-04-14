use crate::atoms::Atom;
use crate::pairs::AtomPair;


#[derive(Default, Debug, Hash)]
pub struct Bond{
    pub(crate) pair: AtomPair
}

impl Bond {

    /// Construct a bond between two atoms
    pub fn from_atoms(atom_i: &Atom, atom_j: &Atom) -> Self{
        Bond::from_atom_indices(atom_i.idx, atom_j.idx)
    }

    /// Construct a bond between two atom indices. These must be distinct
    pub fn from_atom_indices(i: usize, j: usize) -> Self{

        if i == j{
            panic!("Cannot create a bonded pair between identical atoms")
        }

        Bond{pair: AtomPair{i, j}}
    }

    /// Does this bond contain a particular atom?
    pub fn contains(&self, atom: &Atom) -> bool{
        atom.idx == self.pair.i || atom.idx == self.pair.j
    }

    /// Given the index of an atom that may, or may not be present in this bond, return the other
    /// atom index in the bond, if present.
    pub fn other(&self, idx: usize) -> Option<usize>{

        if      idx == self.pair.i { Some(self.pair.j) }
        else if idx == self.pair.j { Some(self.pair.i) }
        else{                       None     }
    }
}

impl PartialEq for Bond {
    fn eq(&self, other: &Self) -> bool { self.pair == other.pair }
}

impl Eq for Bond {}
