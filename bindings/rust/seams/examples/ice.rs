/* Copyright (c) d-SEAMS developers */
fn diamond(nrep: usize) -> (Vec<f64>, f64) {
    let a = 2.75 * 4.0 / 3.0_f64.sqrt();
    let fcc = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let basis = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]];
    let mut xyz = Vec::new();
    for i in 0..nrep {
        for j in 0..nrep {
            for k in 0..nrep {
                for shift in fcc {
                    for site in basis {
                        xyz.push((i as f64 + shift[0] + site[0]) * a);
                        xyz.push((j as f64 + shift[1] + site[1]) * a);
                        xyz.push((k as f64 + shift[2] + site[2]) * a);
                    }
                }
            }
        }
    }
    (xyz, nrep as f64 * a)
}

fn main() {
    let (xyz, length) = diamond(2);
    let mut sys = seams::System::new().expect("system");
    sys.set_ortho(&xyz, length, length, length).expect("frame");
    let ice = sys
        .chill_plus(&seams::ChillOptions::default())
        .expect("chill+");
    println!(
        "dseams {} cubic {} of {}",
        seams::version(),
        ice.counts[0],
        ice.labels.len()
    );
}

#[test]
fn example_runs() {
    main();
}
