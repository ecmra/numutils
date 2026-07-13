use numutils::ch7;

fn main (){
    let mut seed = 3u32;
    for _i in 0..1_000_000{
	let r = seed as f32/u64::pow(2,32) as f32;
	println!("{}", r);
	//println!("{}", seed as f32/u32::pow(2,32) as f32);
	seed = ch7::ran_u32(seed);
    }
}
