use std::str::FromStr;
use std::ops::{Add, AddAssign, Index, IndexMut, Sub, SubAssign};


#[derive(Default, Clone, Debug)]
pub struct CartesianCoordinate{
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl CartesianCoordinate {

    /// Create a cartesian coordinate from a set of optional strings
    pub fn from_option_strings(x: Option<&str>,
                               y: Option<&str>,
                               z: Option<&str>
    ) -> Result<Self, &'static str>{

        for k in [x, y, z]{
            if k.is_none() {
                return Err("An optional was None")
            }
            if f64::from_str(k.unwrap()).is_err(){
                return Err("Failed to parse float")
            };
        }

        let coord = CartesianCoordinate{
            x: x.unwrap().parse::<f64>().unwrap(),
            y: y.unwrap().parse::<f64>().unwrap(),
            z: z.unwrap().parse::<f64>().unwrap()
        };

        Ok(coord)
    }
}

impl Add for CartesianCoordinate {        // +  operator
type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z}
    }
}

impl Sub for CartesianCoordinate {        // -  operator
type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z}
    }
}

impl AddAssign for CartesianCoordinate {   // +=  operator
fn add_assign(&mut self, rhs: Self) {
    self.x += rhs.x;
    self.y += rhs.y;
    self.z += rhs.z;
}
}

impl SubAssign for CartesianCoordinate {
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl Index<usize> for CartesianCoordinate {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            n => panic!("Invalid Vector3d index: {}", n)
        }
    }
}

impl IndexMut<usize> for CartesianCoordinate {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            n => panic!("Invalid Vector3d index: {}", n)
        }
    }
}

