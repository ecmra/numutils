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

    

#[test]
fn padjust_test() {
    // Reference result: p.adjust(c(0.01, 0.04, 0.03, 0.002), method = "BH")
    let pvals = vec![0.01, 0.04, 0.03, 0.002];
    let adjusted = ch14::padjust(&pvals);
    let expected = [0.02, 0.04, 0.04, 0.008];

    assert_eq!(adjusted.len(), expected.len());
    for (actual, expected) in adjusted.iter().zip(expected) {
        assert!(util::almost_equal(*actual, expected, 1e-6));
    }
}

#[test]
fn padjust_handles_ties_and_caps_values_at_one() {
    let pvals = vec![0.01, 0.01, 0.05, 0.7, 0.8, 0.9];
    let adjusted = ch14::padjust(&pvals);
    let expected = [0.03, 0.03, 0.1, 0.9, 0.9, 0.9];

    assert_eq!(adjusted.len(), expected.len());
    for (actual, expected) in adjusted.iter().zip(expected) {
        assert!(util::almost_equal(*actual, expected, 1e-6));
        assert!(*actual <= 1.0);
    }
}
