
pub fn is_very_close(x: f64, y: f64) -> bool{
    is_close(x, y, 1E-8)
}


pub fn is_close(x: f64, y: f64, atol: f64) -> bool{
    // Are two numbers close to within an absolute tolerance?

    if (x - y).abs() <= atol{
        return true
    }

    println!("\nleft = {}\nright = {}", x, y);
    return false
}


pub(crate) fn remove_file_or_panic(filename: &str){
    std::fs::remove_file(filename).expect("Failed to remove file")
}

pub(crate) fn print_methane_xyz_file(filename: &str){

    std::fs::write(filename,
                   "5\n\n\
                        C     0.00000   0.00000   0.00000\n\
                        H    -0.65860  -0.85220  -0.30120\n\
                        H    -0.45940   0.97110  -0.28590\n\
                        H     0.08440  -0.02940   1.10060\n\
                        H     1.02910  -0.10990  -0.41250\n")
        .expect(filename)
}

pub(crate) fn print_dihydrogen_xyz_file(filename: &str){

    std::fs::write(filename,
                   "2\n\n\
                        H     0.00000   0.00000   0.00000\n\
                        H     0.77000   0.00000   0.00000\n")
        .expect(filename)
}

pub(crate) fn print_water_xyz_file(filename: &str){

    std::fs::write(filename,
                   "3\n\n\
                        O    -0.00110   0.36310  -0.00000\n\
                        H    -0.82500  -0.18190  -0.00000\n\
                        H     0.82610  -0.18120   0.00000\n")
        .expect(filename)
}


pub(crate) fn print_ethene_xyz_file(filename: &str){

    std::fs::write(filename,
                   "6\n\n\
C         -4.22149        2.30283        0.00000\n\
C         -2.90167        2.47646        0.00000\n\
H         -4.69362        1.69525        0.76527\n\
H         -4.83470        2.76763       -0.76527\n\
H         -2.28846        2.01166        0.76527\n\
H         -2.42954        3.08405       -0.76527\n")
        .expect(filename)
}

pub(crate) fn print_h2coh_xyz_file(filename: &str){

    std::fs::write(filename,
                   "5\n\n\
C         -4.15592        2.32790       -0.02306\n\
O         -2.94565        2.49183       -0.02968\n\
H         -4.61267        1.71790        0.74844\n\
H         -4.78306        2.78690       -0.78277\n\
H         -2.52046        3.03273       -0.71007\n")
        .expect(filename)
}
