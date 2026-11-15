// Copyright (c) d-SEAMS developers
//! Ice Ih census through the safe seams crate.
use seams::{ChillOptions, Label, System};
use std::env;

fn main() {
    let path = env::args().nth(1).expect("fixture path");
    let mut system = System::new().expect("system");
    system.read(&path, 1, 1).expect("read");
    let n = system.count();
    let order = system.q4q6(3.5).expect("q6");
    let mean = order.q6.iter().sum::<f64>() / order.q6.len() as f64;
    let ice = system.chill_plus(&ChillOptions::default()).expect("chill+");
    let (six, _fraction) = system.six_ring_fraction(3.5).expect("rings");
    println!(
        "rust ice nop {} q6 {:.6} hexagonal {} rings {}",
        n,
        mean,
        ice.counts[Label::Hexagonal as usize],
        six
    );
}
