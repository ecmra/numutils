mod util;

use numutils::ch5;

fn t0_plus_t1(x:f32) -> f32 {
    return 1.0 + x;
}

#[test]
fn chebft_test() {
    let (a,b) = (-1.0f32, 1.0f32);
    let n = 2u8;
    let func = t0_plus_t1;
    let mut c:[f32;2]=[0.0,0.0];
    ch5::chebft(a,b,n,func,&mut c);
    println!("c={:?}", c);
    assert!(util::almost_equal(c[0], 2.0, 1e-7) && util::almost_equal(c[1], 1.0, 1e-7) );
}
