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
        let key = format!(
            "{} {}",
            parts.next().unwrap(),
            parts.next().unwrap()
        );
        let values = parts.map(|token| token.parse::<f64>().unwrap()).collect();
        rows.insert(key, values);
    }
    rows
}

fn close(got: f64, want: f64) {
    let scale = got.abs().max(want.abs()).max(1.0);
    assert!(
        (got - want).abs() <= 1e-8 * scale,
        "{got} != {want}"
    );
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
    let (_labels, counts) = ih
        .chill(3.5, 0, seams_sys::SEAMS_GRAPH_CUTOFF, seams_sys::SEAMS_CHILL_PUBLISHED)
        .unwrap();
    for (got, want_i) in counts.iter().zip(want["ih chill0"].iter()) {
        assert_eq!(*got as f64, *want_i);
    }
    let (q3, _) = ih.steinhardt(3.5, 3, false).unwrap();
    let (q4, _) = ih.steinhardt(3.5, 4, false).unwrap();
    let (q6, bar) = ih.steinhardt(3.5, 6, true).unwrap();
    let (q8, _) = ih.steinhardt(3.5, 8, false).unwrap();
    let (q12, _) = ih.steinhardt(3.5, 12, false).unwrap();
    close(mean(&q3), want["ih q3"][0]);
    close(mean(&q4), want["ih q4"][0]);
    close(mean(&q6), want["ih q6"][0]);
    close(mean(bar.as_ref().unwrap()), want["ih q6bar"][0]);
    close(mean(&q8), want["ih q8"][0]);
    close(mean(&q12), want["ih q12"][0]);
    let rings = ih.rings(3.5, 6).unwrap();
    let first = rings.vertices(0).unwrap();
    assert_eq!(first.len() as f64, want["ih vert"][0]);
    assert_eq!(first.first().copied().unwrap_or(-1) as f64, want["ih vert"][1]);
    let cages = ih
        .cages(3.5, seams_sys::SEAMS_GRAPH_SEEDED, 4, false)
        .unwrap();
    let (n_hex, n_cubic, n_water, n_both) = cages.counts().unwrap();
    let cage_got = [n_hex, n_cubic, n_water, n_both];
    for (got, want_i) in cage_got.iter().zip(want["ih cages"].iter()) {
        assert_eq!(*got as f64, *want_i);
    }
    let knn = ih.knn(4, 3.5, true).unwrap();
    assert_eq!(knn.nnz().unwrap() as f64, want["ih knn"][0]);
    let (r, g, volume) = ih.rdf(6.0, 8).unwrap();
    assert_eq!(r.len() as f64, want["ih rdf"][0]);
    close(volume, want["ih rdf"][1]);
    close(*g.first().unwrap_or(&0.0), want["ih rdf"][2]);
    let (_rr, cn) = ih.running_cn(6.0, 8).unwrap();
    close(*cn.last().unwrap_or(&0.0), want["ih rcn"][0]);
    close(ih.coordination(3.5, 8).unwrap(), want["ih cn"][0]);
    assert_eq!(
        ih.hbonds(3.5, 1, 2, 2.42, 30.0).unwrap() as f64,
        want["ih hbond"][0]
    );
    let (_z, rho) = ih.density_axis(0, 4, 2).unwrap();
    assert_eq!(rho.len() as f64, want["ih rho"][0]);
    close(*rho.first().unwrap_or(&0.0), want["ih rho"][1]);
    let fp = ih.fingerprint(3.5, 2).unwrap();
    let ring_sum: i32 = fp.rings().unwrap().iter().sum();
    assert_eq!(fp.key().unwrap().len() as f64, want["ih fp"][0]);
    assert_eq!(ring_sum as f64, want["ih fp"][1]);
    assert_eq!(ih.site_table("1=waterO").unwrap()[0] as f64, want["ih site"][0]);
    let (n_pairs, n_cation, n_anion) = ih.contact_pairs("1=cationHead").unwrap();
    assert_eq!(n_pairs as f64, want["ih pairs"][0]);
    assert_eq!(n_cation as f64, want["ih pairs"][1]);
    assert_eq!(n_anion as f64, want["ih pairs"][2]);
    let (n_subset, largest, percolation) = ih.domains(3.5, "1=waterO", 0).unwrap();
    assert_eq!(n_subset as f64, want["ih domain"][0]);
    assert_eq!(largest as f64, want["ih domain"][1]);
    close(percolation, want["ih domain"][2]);
    let (perfect, _imperfect) = ih.prisms(3.5).unwrap();
    assert_eq!(perfect[6] as f64, want["ih prism"][0]);
    let sheet = ih.monolayer(3.5, 6, 100.0).unwrap();
    assert_eq!(sheet[6] as f64, want["ih sheet"][0]);

    let mut guest = read(&root, "input/traj/genice_sI.lammpstrj", 1);
    assert_eq!(guest.count() as f64, want["si n"][0]);
    let (n_cages, n_atoms) = guest
        .signature_cages(3.5, seams_sys::SEAMS_GRAPH_CUTOFF, "512")
        .unwrap();
    assert_eq!(n_cages as f64, want["si sig"][0]);
    assert_eq!(n_atoms as f64, want["si sig"][1]);

    let mut xyz = read(&root, "input/xyz/configuration.xyz", -1);
    assert_eq!(xyz.count() as f64, want["xyz n"][0]);
    close(xyz.coordination(6.0, 8).unwrap(), want["xyz cn"][0]);

    let mut con = read(&root, "input/con/tiny_multi_cuh2.con", -1);
    assert_eq!(con.count() as f64, want["con n"][0]);
    assert_eq!(
        con.hbonds(3.5, 0, 2, 2.42, 30.0).unwrap() as f64,
        want["con hbond"][0]
    );

    let mut tri = read(&root, "input/traj/genice_sH.lammpstrj", 1);
    assert_eq!(tri.count() as f64, want["tri n"][0]);
    let (q6, _) = tri.steinhardt(3.5, 6, false).unwrap();
    close(mean(&q6), want["tri q6"][0]);
}
