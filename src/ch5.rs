/// compute coefficients of Chepyshev polynomials (ch5.8 p190ff)
/// 
/// func is defined on [a,b].
/// since Chebyshev polynomials are defined on [-1,1]
/// here we approximate g, a function defined on [-1,1]
/// so that func[x] = g[x1]  where x \in [a,b]
/// and x1 = (x-0.5(b+a))/(0.5*(b-a)) \in [-1,1].
/// When computing an approximation to func[x] at point x0
/// we use the polynomial for g[x1(x0)]
pub fn chebft(a:f32, b:f32, n:u8, func:fn(f32) -> f32, c:&mut [f32]){
    let pi = std::f32::consts::PI;
    // compute cj
    for j in 0..n as usize {
	c[j]=0.0;
	for k in 1..(n+1) as usize {
	    let x1 = f32::cos(pi*(k as f32 - 0.5)/n as f32);
	    let x = x1 * (b-a) / 2.0 + 0.5*(b+a);
	    c[j] += func(x)*f32::cos(pi * j as f32 * (k as f32 - 0.5)/ n as f32)
	}
	c[j] = 2.0*c[j]/ n as f32;
    }
}

