//! Ice Ih census through the safe seams crate.
use seams::System;
use std::env;

fn main() {
    let path = env::args().nth(1).expect("fixture path");
    let mut system = System::new().expect("system");
    system.read(&path, 1, 1).expect("read");
    let n = system.count();
    let (_q4, q6) = system.q4q6(3.5).expect("q6");
    let mean = q6.iter().sum::<f64>() / q6.len() as f64;
    let (_labels, counts) = system
        .chill_plus(
            3.5,
            seams_sys::SEAMS_GRAPH_CUTOFF,
            seams_sys::SEAMS_CHILL_PUBLISHED,
        )
        .expect("chill+");
    let (six, _fraction) = system.six_ring_fraction(3.5).expect("rings");
    println!(
        "rust ice nop {} q6 {:.6} hexagonal {} rings {}",
        n, mean, counts[1], six
    );
}
