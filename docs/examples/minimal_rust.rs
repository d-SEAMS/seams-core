//! Read the GenIce sI cell through the safe seams crate.
use seams::{System, num_threads, set_num_threads, version};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn read_dump(path: &str) -> (Vec<f64>, [f64; 3]) {
    let file = File::open(path).expect("open fixture");
    let mut xyz = Vec::new();
    let mut edges = [0.0; 3];
    let mut tilts = [0.0; 3];
    let mut lows = [0.0; 3];
    let mut axis = 0;
    let mut reading = false;
    let mut bounds = false;
    for line in BufReader::new(file).lines() {
        let line = line.expect("line");
        if line.starts_with("ITEM: BOX BOUNDS") {
            bounds = true;
            axis = 0;
            reading = false;
        } else if bounds {
            let mut parts = line.split_whitespace();
            let lo: f64 = parts.next().unwrap().parse().unwrap();
            let hi: f64 = parts.next().unwrap().parse().unwrap();
            let tilt: f64 = parts.next().and_then(|item| item.parse().ok()).unwrap_or(0.0);
            // Recovered after all three rows. Store the bound span and
            // the tilt, then convert once the third row arrives.
            edges[axis] = hi - lo;
            tilts[axis] = tilt;
            lows[axis] = lo;
            axis += 1;
            if axis == 3 {
                let xy = tilts[0];
                let xz = tilts[1];
                let yz = tilts[2];
                let xmin = xy.min(xz).min(xy + xz).min(0.0);
                let xmax = xy.max(xz).max(xy + xz).max(0.0);
                let ymin = yz.min(0.0);
                let ymax = yz.max(0.0);
                edges[0] = edges[0] - xmax + xmin;
                edges[1] = edges[1] - ymax + ymin;
                bounds = false;
                let _ = lows;
            }
        } else if line.starts_with("ITEM: ATOMS") {
            reading = true;
        } else if reading {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 5 {
                break;
            }
            xyz.extend(parts[2..5].iter().map(|v| v.parse::<f64>().unwrap()));
        }
    }
    (xyz, edges)
}

fn main() {
    let path = env::args().nth(1).expect("fixture path");
    let (xyz, edges) = read_dump(&path);
    assert_eq!(xyz.len() / 3, 46);
    set_num_threads(1).expect("threads");
    assert_eq!(num_threads(), 1);
    let mut system = System::new().expect("system");
    system
        .set_ortho(&xyz, edges[0], edges[1], edges[2])
        .expect("frame");
    assert_eq!(system.count(), 46);
    let degree = system.neighbour_degree(3.5).expect("neighbours");
    assert_eq!(degree.len(), 46);
    let (_q4, q6) = system.q4q6(3.5).expect("q6");
    let mean = q6.iter().sum::<f64>() / q6.len() as f64;
    let (_labels, counts) = system
        .chill_plus(3.5, seams_sys::SEAMS_GRAPH_CUTOFF, seams_sys::SEAMS_CHILL_PUBLISHED)
        .expect("chill+");
    let (six, _fraction) = system.six_ring_fraction(3.5).expect("rings");
    assert_eq!(counts[4], 46);
    assert_eq!(six, 54);
    assert!(mean > 0.5);
    println!("rust example ok {} q6 {:.3}", version(), mean);
}
