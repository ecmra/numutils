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


/// incomplete beta
///
/// p227
pub fn betai(a:f32, b:f32, x:f32) -> f32 {
    let num0ln = betaln(a + 1.0 ,1.0) + f32::ln(x);
    let den0ln = betaln(a+b, 1.0);
    let (mut numln,  mut denln) = (num0ln, den0ln);
    let maxit:u32 = 100;
    let mut s:f32 = f32::exp(num0ln - den0ln);
    for n in 1..=maxit {
	numln = numln + f32::ln(n as f32) - f32::ln(a + n as f32 + 1.0) + f32::ln(x) ;
	denln = denln + f32::ln(n as f32) - f32::ln( a + b + n as f32);
	s += f32::exp(numln - denln);
    }
    let f1 = f32::powf(x, a) * f32::powf(1.0 - x, b) /(a * f32::exp(betaln(a,b))); 
    return f1 * (1.0 + s);
}




