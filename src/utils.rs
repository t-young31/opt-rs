
pub fn is_very_close(x: f64, y: f64) -> bool{
    is_close(x, y, 1E-8)
}


pub fn is_close(x: f64, y: f64, atol: f64) -> bool{
    // Are two numbers close to within an absolute tolerance?
    println!("\nleft = {}\nright = {}", x, y);
    (x - y).abs() <= atol
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
