
fn is_very_close(x: f64, y: f64) -> bool{
    is_close(x, y, 1E-8)
}


fn is_close(x: f64, y: f64, atol: f64) -> bool{
    // Are two numbers close to within an absolute tolerance?
    println!("\nleft = {}\nright = {}", x, y);
    (x - y).abs() <= atol
}
