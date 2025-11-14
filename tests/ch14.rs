mod util;

use numutils::ch14;

#[test]
fn var_test() {
    let x = vec![1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0];
    let xm = x.iter().sum::<f32>()/x.len() as f32;
    let v = ch14::var(&x,xm);
    assert!(util::almost_equal(v, 9.166667 , 1e-6));
}

#[test]
fn avevar_test() {
    let x = vec![1.0, 2.0, 3.0, 4.0];
    let avevar = ch14::avevar(&x);
    assert!(util::almost_equal(avevar.0, 2.5, 1e-6) && 
	    util::almost_equal(avevar.1,1.666666667 , 1e-6) );
}


#[test]
fn tutest_test() {
    let x1 = vec![1.0, 2.0, 3.0];
    let x2 = vec![4.0, 5.0, 6.0];
    let res = ch14::tutest(&x1, &x2);
    assert!(util::almost_equal(res.0, -3.67423, 1e-5) &&
	    util::almost_equal(res.1, 0.0213116, 1e-6));
}

    
