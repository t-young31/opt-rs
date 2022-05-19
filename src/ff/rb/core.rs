use crate::ff::forcefield::Forcefield;
use crate::coordinates::{Point, Vector3D};
use crate::ff::forcefield::EnergyFunction;
use crate::Molecule;


#[derive(Default)]
pub struct RB{

    energy_functions: Vec<Box<dyn EnergyFunction>>,

    energy:           f64,
    gradient:         Vec<Vector3D>,
}


impl RB {
    // TODO
}

impl Forcefield for RB {

    fn new(molecule: &Molecule) -> Self where Self: Sized {
        todo!()
    }

    fn set_atom_types(&mut self, molecule: &Molecule) {
        todo!()
    }

    fn energy(&mut self, coordinates: &Vec<Point>) -> f64 {
        self.energy = self.energy_functions.iter().map(|f| f.energy(coordinates)).sum();
        self.energy
    }

    fn gradient(&mut self, coordinates: &Vec<Point>) -> &Vec<Vector3D> {
        self.gradient.iter_mut().for_each(|v| v.zero());
        self.energy_functions.iter().for_each(|f| f.add_gradient(coordinates, &mut self.gradient));
        &self.gradient
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
mod tests{

    use super::*;
    use crate::utils::*;

}