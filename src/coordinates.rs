use std::str::FromStr;
use std::ops::{Add, Sub, Index, IndexMut};


/// Point in 3D space
#[derive(Default, Clone, Debug)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl Point {

    /// Create a cartesian coordinate (point in 3D space) from a set of optional strings
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

        let coord = Point {
            x: x.unwrap().parse::<f64>().unwrap(),
            y: y.unwrap().parse::<f64>().unwrap(),
            z: z.unwrap().parse::<f64>().unwrap()
        };

        Ok(coord)
    }
}

pub struct Vector3D {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl Vector3D {

    /// Length of this vector
    #[inline(always)]
    pub fn length(&self) -> f64{
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    #[inline(always)]
    pub fn dot(&self, other: &Vector3D) -> f64{
        self.x * other.x + self.y * other.y + self.z * other.z
    }

}

impl<'a, 'b> Add<&'b Point> for &'a Point {        // +  operator
type Output = Vector3D;

    fn add(self, other: &'b Point) -> Vector3D {
        Vector3D{x: self.x + other.x,
                 y: self.y + other.y,
                 z: self.z + other.z}
    }
}
impl<'a, 'b> Sub<&'b Point> for &'a Point {        // -  operator
    type Output = Vector3D;

    fn sub(self, other: &'b Point) -> Vector3D {
        Vector3D{x: self.x - other.x,
                 y: self.y - other.y,
                 z: self.z - other.z}
    }
}

impl Index<usize> for Point {
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

impl IndexMut<usize> for Point {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            n => panic!("Invalid Vector3d index: {}", n)
        }
    }
}
