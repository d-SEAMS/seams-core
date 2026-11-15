/* Copyright (c) d-SEAMS developers */
//! Compare the safe wrappers with the C++ engine on ice Ih, sI, XYZ, CON,
//! and the triclinic sH dump.
use std::collections::HashMap;
use std::env;
use std::process::Command;

fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().sum::<f64>() / values.len() as f64
}

fn load() -> HashMap<String, Vec<f64>> {
    let oracle = env::var("SEAMS_CPP_ORACLE").expect("SEAMS_CPP_ORACLE");
    let root = env::var("SEAMS_SOURCE_ROOT").expect("SEAMS_SOURCE_ROOT");
    let mut command = Command::new(oracle);
    command.arg(&root).env("OMP_NUM_THREADS", "1");
    let output = command.output().expect("run the C++ oracle");
    assert!(
        output.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut rows = HashMap::new();
    for line in String::from_utf8(output.stdout).unwrap().lines() {
        let mut parts = line.split_whitespace();
        let key = format!("{} {}", parts.next().unwrap(), parts.next().unwrap());
        let values = parts.map(|token| token.parse::<f64>().unwrap()).collect();
        rows.insert(key, values);
    }
    rows
}

fn close(got: f64, want: f64) {
    let scale = got.abs().max(want.abs()).max(1.0);
    assert!((got - want).abs() <= 1e-8 * scale, "{got} != {want}");
}

fn read(root: &str, path: &str, filter: i32) -> seams::System {
    let mut sys = seams::System::new().unwrap();
    sys.read(&format!("{root}/{path}"), 1, filter).unwrap();
    sys
}

#[test]
fn analyses_match_the_cpp_path() {
    seams::set_num_threads(1).unwrap();
    let root = env::var("SEAMS_SOURCE_ROOT").unwrap();
    let want = load();

    let mut ih = read(&root, "input/traj/ice1h.lammpstrj", -1);
    assert_eq!(ih.count() as f64, want["ih n"][0]);
    let ice = ih
        .chill(&seams::ChillOptions {
            cutoff: 3.5,
            species: 0,
            graph: seams::Graph::Cutoff,
            rule: seams::ChillRule::Published,
        })
        .unwrap();
    for (got, want_i) in ice.counts.iter().zip(want["ih chill0"].iter()) {
        assert_eq!(*got as f64, *want_i);
    }
    println!("rustcall hexagonal {}", ice.counts[1]);
    let q3 = ih.steinhardt(&order(3, false)).unwrap();
    let q4 = ih.steinhardt(&order(4, false)).unwrap();
    let q6 = ih.steinhardt(&order(6, true)).unwrap();
    let q8 = ih.steinhardt(&order(8, false)).unwrap();
    let q12 = ih.steinhardt(&order(12, false)).unwrap();
    close(mean(&q3.local), want["ih q3"][0]);
    close(mean(&q4.local), want["ih q4"][0]);
    close(mean(&q6.local), want["ih q6"][0]);
    close(mean(q6.average.as_ref().unwrap()), want["ih q6bar"][0]);
    println!("rustcall q6bar {:.12}", mean(q6.average.as_ref().unwrap()));
    close(mean(&q8.local), want["ih q8"][0]);
    close(mean(&q12.local), want["ih q12"][0]);
    let rings = ih
        .rings(&seams::RingOptions {
            cutoff: 3.5,
            max_size: 6,
            ..seams::RingOptions::default()
        })
        .unwrap();
    let first = rings.vertices(0).unwrap();
    assert_eq!(first.len() as f64, want["ih vert"][0]);
    assert_eq!(
        first.first().copied().unwrap_or(-1) as f64,
        want["ih vert"][1]
    );
    println!(
        "rustcall vert {} {}",
        first.len(),
        first.first().copied().unwrap_or(-1)
    );
    let cages = ih
        .cages(&seams::CageOptions {
            cutoff: 3.5,
            graph: seams::Graph::Seeded,
            k: 4,
            complete: false,
            ..seams::CageOptions::default()
        })
        .unwrap();
    let counts = cages.counts().unwrap();
    println!(
        "rustcall cages {} {} {} {}",
        counts.hexagonal, counts.cubic, counts.water, counts.both
    );
    let cage_got = [counts.hexagonal, counts.cubic, counts.water, counts.both];
    for (got, want_i) in cage_got.iter().zip(want["ih cages"].iter()) {
        assert_eq!(*got as f64, *want_i);
    }
    let knn = ih
        .knn(&seams::KnnOptions {
            k: 4,
            cutoff: 3.5,
            mutual: true,
            ..seams::KnnOptions::default()
        })
        .unwrap();
    assert_eq!(knn.nnz().unwrap() as f64, want["ih knn"][0]);
    println!("rustcall knn {}", knn.nnz().unwrap());
    let rdf = ih
        .rdf(&seams::DistributionOptions {
            rmax: 6.0,
            bins: 8,
            ..seams::DistributionOptions::default()
        })
        .unwrap();
    assert_eq!(rdf.radius.len() as f64, want["ih rdf"][0]);
    close(rdf.volume, want["ih rdf"][1]);
    close(*rdf.g.first().unwrap_or(&0.0), want["ih rdf"][2]);
    let running = ih
        .running_cn(&seams::DistributionOptions {
            rmax: 6.0,
            bins: 8,
            ..seams::DistributionOptions::default()
        })
        .unwrap();
    close(
        *running.coordination.last().unwrap_or(&0.0),
        want["ih rcn"][0],
    );
    close(
        ih.coordination(&seams::DistributionOptions {
            rmax: 3.5,
            bins: 8,
            ..seams::DistributionOptions::default()
        })
        .unwrap(),
        want["ih cn"][0],
    );
    assert_eq!(
        ih.hbonds(&seams::BondOptions {
            cutoff: 3.5,
            oxygen: 1,
            hydrogen: 2,
            distance: 2.42,
            angle: 30.0,
        })
        .unwrap() as f64,
        want["ih hbond"][0]
    );
    let density = ih
        .density(&seams::DensityOptions {
            species: 0,
            bins: 4,
            axis: seams::Axis::Z,
        })
        .unwrap();
    assert_eq!(density.density.len() as f64, want["ih rho"][0]);
    close(*density.density.first().unwrap_or(&0.0), want["ih rho"][1]);
    let fp = ih
        .fingerprint(&seams::FingerprintOptions {
            cutoff: 3.5,
            hops: 2,
            ..seams::FingerprintOptions::default()
        })
        .unwrap();
    let ring_sum: i32 = fp.rings().unwrap().iter().sum();
    assert_eq!(fp.key().unwrap().len() as f64, want["ih fp"][0]);
    assert_eq!(ring_sum as f64, want["ih fp"][1]);
    assert_eq!(
        ih.site_table("1=waterO").unwrap()[0] as f64,
        want["ih site"][0]
    );
    let contacts = ih.contact_pairs("1=cationHead").unwrap();
    assert_eq!(contacts.pairs as f64, want["ih pairs"][0]);
    assert_eq!(contacts.cations as f64, want["ih pairs"][1]);
    assert_eq!(contacts.anions as f64, want["ih pairs"][2]);
    let domain = ih
        .domains(&seams::DomainOptions {
            cutoff: 3.5,
            spec: "1=waterO".into(),
            subset: seams::DomainSubset::Apolar,
        })
        .unwrap();
    assert_eq!(domain.subset as f64, want["ih domain"][0]);
    assert_eq!(domain.largest as f64, want["ih domain"][1]);
    close(domain.percolation, want["ih domain"][2]);
    let prisms = ih.prisms(&seams::PrismOptions::default()).unwrap();
    assert_eq!(prisms.perfect[6] as f64, want["ih prism"][0]);
    let sheet = ih
        .monolayer(&seams::SheetOptions {
            cutoff: 3.5,
            max_size: 6,
            sheet_area: 100.0,
            ..seams::SheetOptions::default()
        })
        .unwrap();
    assert_eq!(sheet[6] as f64, want["ih sheet"][0]);
    println!("rustcall sheet {}", sheet[6]);

    let mut guest = read(&root, "input/traj/genice_sI.lammpstrj", 1);
    assert_eq!(guest.count() as f64, want["si n"][0]);
    let found = guest
        .signature_cages(&seams::SignatureOptions {
            cutoff: 3.5,
            graph: seams::Graph::Cutoff,
            spec: "512".into(),
            ..seams::SignatureOptions::default()
        })
        .unwrap();
    assert_eq!(found.cages as f64, want["si sig"][0]);
    assert_eq!(found.atoms as f64, want["si sig"][1]);
    println!("rustcall sig {} {}", found.cages, found.atoms);

    let mut xyz = read(&root, "input/xyz/configuration.xyz", -1);
    assert_eq!(xyz.count() as f64, want["xyz n"][0]);
    close(
        xyz.coordination(&seams::DistributionOptions {
            rmax: 6.0,
            bins: 8,
            ..seams::DistributionOptions::default()
        })
        .unwrap(),
        want["xyz cn"][0],
    );

    let mut con = read(&root, "input/con/tiny_multi_cuh2.con", -1);
    assert_eq!(con.count() as f64, want["con n"][0]);
    assert_eq!(
        con.hbonds(&seams::BondOptions {
            cutoff: 3.5,
            oxygen: 0,
            hydrogen: 2,
            distance: 2.42,
            angle: 30.0,
        })
        .unwrap() as f64,
        want["con hbond"][0]
    );

    let mut tri = read(&root, "input/traj/genice_sH.lammpstrj", 1);
    assert_eq!(tri.count() as f64, want["tri n"][0]);
    let q6 = tri.steinhardt(&order(6, false)).unwrap();
    close(mean(&q6.local), want["tri q6"][0]);
}

fn order(degree: i32, average: bool) -> seams::SteinhardtOptions {
    seams::SteinhardtOptions {
        cutoff: 3.5,
        order: degree,
        average,
        ..seams::SteinhardtOptions::default()
    }
}
