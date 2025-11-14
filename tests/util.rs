/// checks almost equality
///
/// used in testing other functions
pub fn almost_equal(x: f32, y: f32, tol: f32) -> bool {
    eprintln!("almost_equal? x:{} y:{}", x, y);
    return f32::abs(x - y) < tol;
}



