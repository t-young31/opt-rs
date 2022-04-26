use std::str::FromStr;
use std::ops::{Add, Sub, Index, IndexMut, Neg, DivAssign};
use std::process::Output;


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


#[derive(Debug, Clone)]
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

    /// Dot (scalar) product
    #[inline(always)]
    pub fn dot(&self, other: &Vector3D) -> f64{
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product
    #[inline(always)]
    pub fn cross(&self, other: &Vector3D) -> Self{
        Vector3D{x: self.y*other.z - self.z*other.y,
                 y: self.z*other.x - self.x*other.z,
                 z: self.x*other.y - self.y*other.x}
    }

    /// Divide this vector by a scalar
    pub fn divide_by(&mut self, value: f64){
        self.x /= value;
        self.y /= value;
        self.z /= value;
    }
}

impl Neg for Vector3D {
    type Output = Vector3D;

    fn neg(self) -> Self::Output {
        Vector3D{x: -self.x, y: -self.y, z: -self.z}
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



#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::*;

    #[test]
    fn test_point_addition() {

        let vec = &Point{x: 1., y: -2., z: -3.} + &Point{x: -1., y: 2., z: 3.};

        assert!(is_very_close(vec.x, 0.));
        assert!(is_very_close(vec.y, 0.));
        assert!(is_very_close(vec.z, 0.));
    }

    #[test]
    fn test_point_subtraction() {

        let vec = &Point{x: 1., y: -2., z: -3.} - &Point{x: -1., y: 2., z: 3.};

        assert!(is_very_close(vec.x, 2.));
        assert!(is_very_close(vec.y, -4.));
        assert!(is_very_close(vec.z, -6.));
    }

    #[test]
    fn test_vector_dot(){

        let v0 = Vector3D{x: 1., y: 1., z: 1.};
        let v1 = Vector3D{x: -1., y: 2., z: -1.};

        assert!(is_very_close(v0.dot(&v1), 0.));
    }

    #[test]
    fn test_cross_product(){

        let v0 = Vector3D{x: 1., y: 1., z: 1.};
        let v1 = Vector3D{x: -1., y: 2., z: -1.};
        let v2 = v0.cross(&v1);

        assert!(is_very_close(v0.dot(&v2), 0.));
        assert!(is_very_close(v1.cross(&v1).length(), 0.));

    }
}