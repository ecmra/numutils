/// return the value of ln(gamma(xx)) xx>0
///
/// from numerical recipes in C, second edition
/// http://s3.amazonaws.com/nrbook.com/book_C210.html
/// or ch6 p214
pub fn gammln(xx: f32) -> f32 {
    let (x, mut y, mut tmp, mut ser): (f64, f64, f64, f64);
    let cof: [f64; 6] = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5,
    ];
    x = xx as f64;
    y = x;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * f64::ln(tmp);
    ser = 1.000000000190015;
    for j in 0..6 {
        y += 1.0;
        ser += cof[j] / y;
    }
    return (-tmp + f64::ln(2.5066282746310005 * ser / x)) as f32;
}


pub fn lnfact(n:u32) -> f32{
    return gammln((n+1) as f32);
}


/// return log beta(z,w)
pub fn betaln(z: f32, w: f32) -> f32 {
    return gammln(z) + gammln(w) - gammln(z + w);
}


/// return incomplete gamma function Q evaluated
/// by continued fraction representation
/// NR 6.2 p219
pub fn gammq(a: f32, x: f32) -> f32 {
    const EPS: f32 = 1e-7;
    const FPMIN: f32 = 1.0e-30;
    const MAXIT: i32 = 100;
    let mut b: f32 = x + 1.0 - a;
    let mut c: f32 = 1.0 / FPMIN;
    let mut d: f32 = 1.0 / b;
    let mut an: f32;
    let mut delta: f32;
    let mut h = d;
    let mut i = 0;
    loop {
        an = -i as f32 * (i as f32 - a);
        b += 2.0;
        d = b + an * d;
        if f32::abs(d) < FPMIN {
            d = FPMIN
        }
        c = b + an / c;
        if f32::abs(c) < FPMIN {
            c = FPMIN
        }
        d = 1.0 / d;
        delta = c * d;
        h *= delta;
        if f32::abs(delta - 1.0) < EPS {
            break;
        }
        i += 1;
        if i > MAXIT {
            break;
        }
    }
    return f32::exp(-x + a * f32::ln(x) - gammln(a)) * h;
}

/// incomplete gamma function
///
/// NR in C, 2nd edition, p217ff
/// here I'm fixing the number of terms in the series,
/// should check convergence instead.
pub fn gammp(a: f32, x: f32) -> f32 {
    let (mut num, mut den, mut a1) = (1.0, a, a);
    let maxit = 100;
    let mut s = num / den;
    let mut nit = 0;
    if x < (a + 1.0) {
        loop {
            num = num * x;
            a1 += 1.0;
            den = den * a1;
            s += num / den;
            if nit > maxit {
                break;
            }
            nit += 1;
        }
        return f32::exp(-x + a * f32::ln(x) - gammln(a)) * s;
    } else {
        return 1.0 - gammq(a, x);
    }
}

/// chi squared (cumulative distribution function)
///
/// NR in C, 2nd edition, p221ff
pub fn chi2(obs: f32, df: u32) -> f32 {
    return gammp(df as f32 / 2.0, obs / 2.0);
}

/// log binomial coefficient
///
/// NR in C, 2nd edition, p215ff
pub fn logbico(n: u32, k: u32) -> f32 {
    return lnfact(n) - lnfact(k) - lnfact(n - k);
}




#[cfg(test)]
mod tests {
    use super::*;
    /// checks almost equality
    ///
    /// used in testing other functions
    fn almost_equal(x: f32, y: f32, tol: f32) -> bool {
        eprintln!("almost_equal? x:{} y:{}", x, y);
        return f32::abs(x - y) < tol;
    }
    #[test]
    fn gammaln_test() {
        assert!(f32::abs(gammln(0.3) - 1.095797995) < 1e-6);
        assert!(f32::abs(gammln(3.0) - 0.6931471806) < 1e-6);
        assert!(f32::abs(gammln(30f32) - 71.25703897) < 1e-6);
    }
    #[test]
    fn lnfact_test() {
        assert!(almost_equal(lnfact(12), 19.98721450, 1e-6));
    }
    #[test]
    /// Mathematica has a different definition of incomplete gamma.
    /// Hence I use ig[a_, z_] := 1 - Gamma[a, z]/Gamma[a]
    /// N[ig[1, 3], 10] = 0.9502129316
    /// N[ig[1, 3/10], 10] = 0.2591817793
    fn gammq_test() {
        let res1 = gammq(1f32, 3f32);
        println!("gammq1 {}", res1);
        let res2 = gammq(1f32, 0.3);
        println!("gammq2 {}", res2);
        assert!(almost_equal(res1, 1.0 - 0.9502129316, 1e-6));
        assert!(almost_equal(res2, 1.0 - 0.2591817793, 1e-6));
    }
    #[test]
    fn gammp_test() {
        let res1 = gammp(1f32, 3f32);
        println!("gammp1 {}", res1);
        let res2 = gammp(1f32, 0.3);
        println!("gammp2 {}", res2);
        assert!(almost_equal(res1, 0.9502129316, 1e-6));
        assert!(almost_equal(res2, 0.2591817793, 1e-6));
    }

    #[test]
    fn chi2_test() {
        // In R this is pchisq(2.0,2)
        assert!(almost_equal(chi2(2.0, 2), 0.6321205588, 1e-6));
        assert!(almost_equal(chi2(2.0, 3), 0.4275932955, 1e-6));
    }
}
