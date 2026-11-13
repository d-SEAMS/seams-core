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

#[test]
fn cubic_ice_through_the_c_api() {
    seams::set_num_threads(1).unwrap();
    assert_eq!(seams::num_threads(), 1);
    assert!(!seams::version().is_empty());
    let (xyz, length) = diamond(2);
    let mut sys = seams::System::new().unwrap();
    sys.set_ortho(&xyz, length, length, length).unwrap();
    let ice = sys.chill_plus(&seams::ChillOptions::default()).unwrap();
    assert!(ice.labels.iter().all(|label| *label == 0));
    assert_eq!(ice.counts[0] as usize, ice.labels.len());
    assert_eq!(sys.count() as usize, ice.labels.len());
    let mutual = sys
        .chill_plus(&seams::ChillOptions {
            graph: seams::Graph::Mutual,
            rule: seams::ChillRule::Ovito,
            ..seams::ChillOptions::default()
        })
        .unwrap();
    assert_eq!(mutual.counts[0] as usize, ice.labels.len());
    let orders = sys.q4q6(3.5).unwrap();
    let mean = orders.q6.iter().sum::<f64>() / orders.q6.len() as f64;
    assert!(mean > 0.3);
    assert!(orders.q4.iter().all(|v| v.is_finite()));
    let degree = sys.neighbour_degree(3.5).unwrap();
    assert!(degree.iter().all(|d| *d == 4));
    let (six, fraction) = sys.six_ring_fraction(3.5).unwrap();
    assert!(six > 0);
    assert!(fraction > 0.5);
    let mut tilted = [0.0; 9];
    tilted[3] = length;
    tilted[4] = length;
    tilted[5] = length;
    tilted[6] = 0.3;
    sys.set_frame(&xyz, &tilted).unwrap();
    sys.chill_plus(&seams::ChillOptions::default()).unwrap();
}
