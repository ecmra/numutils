mod util;
use numutils::ch6;

#[test]
fn gammaln_test() {
    assert!(f32::abs(ch6::gammln(0.3) - 1.095797995) < 1e-6);
    assert!(f32::abs(ch6::gammln(3.0) - 0.6931471806) < 1e-6);
    assert!(f32::abs(ch6::gammln(30f32) - 71.25703897) < 1e-6);
}
#[test]
fn lnfact_test() {
    assert!(util::almost_equal(ch6::lnfact(12), 19.98721450, 1e-6));
}
#[test]
/// Mathematica has a different definition of incomplete gamma.
/// Hence I use ig[a_, z_] := 1 - Gamma[a, z]/Gamma[a]
/// N[ig[1, 3], 10] = 0.9502129316
/// N[ig[1, 3/10], 10] = 0.2591817793
fn gammq_test() {
    let res1 = ch6::gammq(1f32, 3f32);
    println!("gammq1 {}", res1);
    let res2 = ch6::gammq(1f32, 0.3);
    println!("gammq2 {}", res2);
    assert!(util::almost_equal(res1, 1.0 - 0.9502129316, 1e-6));
    assert!(util::almost_equal(res2, 1.0 - 0.2591817793, 1e-6));
}
#[test]
fn gammp_test() {
    let res1 = ch6::gammp(1f32, 3f32);
    println!("gammp1 {}", res1);
    let res2 = ch6::gammp(1f32, 0.3);
    println!("gammp2 {}", res2);
    assert!(util::almost_equal(res1, 0.9502129316, 1e-6));
    assert!(util::almost_equal(res2, 0.2591817793, 1e-6));
}

#[test]
fn chi2_test() {
    // In R this is pchisq(2.0,2)
    assert!(util::almost_equal(ch6::chi2(2.0, 2), 0.6321205588, 1e-6));
    assert!(util::almost_equal(ch6::chi2(2.0, 3), 0.4275932955, 1e-6));
}

#[test]
// BetaRegularized[0.5, 1, 3] 
fn betai_test(){
    let (a,b,x)= (1.0, 3.0, 0.5);
    let res = ch6::betai(a,b,x);
    assert!(util::almost_equal(res, 0.875, 1e-6));
}
    
#[test]
// Log[Beta[2.0, 3.0]] = -2.48491
fn betaln_test(){
    let (a,b)= (2.0, 3.0);
    let res = ch6::betaln(a,b);
    assert!(util::almost_equal(res, -2.48491, 1e-5));
}

