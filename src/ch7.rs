
/**
- keep calling with the result of the
  previous call
- "an even quicker generator", NR 7.1
*/
pub fn ran_u32(seed:u32) -> u32 {
    return 1664525u32 * seed + 1013904223u32
}
