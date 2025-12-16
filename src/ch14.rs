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
 
/// return index which sort a vector
///
/// eg order(c(8,1,3)) [1] 1 2 0
fn order(v: &Vec<f32>, dec: bool) -> Vec<usize> {
    let mut idx: Vec<usize> = (1..v.len() + 1).into_iter().collect();
    //println!("{:?}", idx);
    if dec {
        idx.sort_by(|i, j| v[j - 1].total_cmp(&v[i - 1]));
    } else {
        idx.sort_by(|i, j| v[i - 1].total_cmp(&v[j - 1]));
    }
    idx
}

/// cumulative minimum
///
/// cummin(v)[j] = min(v[0..j])
fn cummin(v: &Vec<f32>) -> Vec<f32> {
    let mut res: Vec<f32> = vec![0.0; v.len()];
    res[0] = v[0];
    for i in 1..v.len() {
        if v[i] < res[i - 1] {
            res[i] = v[i]
        } else {
            res[i] = res[i - 1]
        }
    }
    res
}

/// parallel min of v1 and v2
///
/// create a vector which stores at each position
/// the minimum of v1 and v2 at the same index
fn pmin(v1: &Vec<f32>, v2: &Vec<f32>) -> Vec<f32> {
    let mut res: Vec<f32> = vec![0.0; v1.len()];
    for i in 0..v1.len() {
        if v1[i] < v2[i] {
            res[i] = v1[i];
        } else {
            res[i] = v2[i];
        }
    }
    res
}

/// adjust pvalues according to BH method
///
/// translation of the R function p.adjust
/// R help for the p.ajust function says:
/// n: number of comparisons, must be at least ‘length(p)’; only set
///          this (to non-default) when you know what you are doing!
///R code for p.adjust (BH method)
/// BH = {        i <- lp:1L
/// o <- order(p, decreasing = TRUE)
///ro <- order(o)
///     pmin(1, cummin(n/i * p[o]))[ro]
///}
pub fn padjust(pvals: &Vec<f32>) -> Vec<f32> {
    let lp = pvals.len();
    let o = order(pvals, true);
    //println!("order of pvals:{:?}", o);
    let of32: Vec<f32> = o.iter().map(|x| *x as f32).collect();
    let ro = order(&of32, false);
    //println!("order of order of pvals:{:?}", ro);
    let ones = vec![1.0; lp];
    let tmpv: Vec<f32> = (1..=lp)
        .rev()
        .map(|i| (lp as f32) / (i as f32) * pvals[o[lp - i] - 1])
        .collect();
    //println!("tmpv:{:?}", tmpv);
    //println!("cummin(tmpv):{:?}", &cummin(&tmpv));
    let tmpv2 = pmin(&ones, &cummin(&tmpv));
    (0..ro.len()).map(|i| tmpv2[ro[i] - 1]).collect()
}
