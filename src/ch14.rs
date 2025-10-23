// statistical description of data

/// variance
///
/// xm:mean
pub fn var(x:&Vec<f32>, mean:f32) -> f32 {
    let n = x.len();
    let mut s:f32=0.0;
    for j in 0..n {
	s+= x[j] * x[j] 
    }
    return 1.0/(n-1) as f32 * ( s - n as f32 * mean * mean); 
}


#[test]
fn var_test() {
    let x = vec![1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0];
    let xm = x.iter().sum::<f32>()/x.len() as f32;
    let v = var(&x,xm);
    assert!(f32::abs(v-9.166667) < 1e-6);
}
