// statistical description of data

use crate::ch6;

/// variance
pub fn var(x:&Vec<f32>, mean:f32) -> f32 {
    let n = x.len();
    let mut s:f32=0.0;
    for j in 0..n {
	s+= x[j] * x[j] 
    }
    return 1.0/(n-1) as f32 * ( s - n as f32 * mean * mean); 
}

/// average, variance
pub fn avevar(x:&Vec<f32>) -> (f32, f32) {
    let n = x.len();
    let ave = x.iter().sum::<f32>()/(n as f32);
    let (mut var, mut ep) = (0f32,0f32); 
    for j in 0..n {
	let s  = x[j] - ave;
	var += s*s;
	ep += s;
    }
    (ave, (var - ep*ep/(n as f32))/( n as f32 - 1.0) )
}

macro_rules! sqr {
    ( $( $x:expr )?) => {
	$(
	    $x * $x
	)?
    };
}


/// t-test, unequal variances
///
/// return t statistics, 1-pval
pub fn tutest(x1:&Vec<f32>, x2:&Vec<f32>) -> (f32, f32) {
    let (n1, n2)= (x1.len() as f32 , x2.len() as f32);
    let (ave1, var1) = avevar(x1);
    let(ave2, var2) = avevar(x2);
    let t = (ave1 - ave2)/f32::sqrt(var1/n1 + var2/n2 );
    let df = sqr!(var1/n1+var2/n2)/(sqr!(var1/n1)/(n1-1f32)+sqr!(var2/n2)/(n2-1f32));
    let pval = ch6::betai(0.5*df,0.5,df/(df+sqr!(t)));
    (t, pval)
}
 
