/* Copyright (c) d-SEAMS developers
 * Times the safe Rust wrapper. One process, commands on stdin.
 */
use seams::System;
use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::time::Instant;

fn retune(threads: i32) -> Result<(), Box<dyn std::error::Error>> {
    let text = threads.to_string();
    std::env::set_var("OMP_NUM_THREADS", &text);
    std::env::set_var("PLUMED_NUM_THREADS", &text);
    seams::set_num_threads(threads)?;
    Ok(())
}

fn main() {
    let mut frames: HashMap<String, System> = HashMap::new();
    let mut loaded: HashMap<String, bool> = HashMap::new();
    let mut sink = 0.0_f64;
    let stdin = std::io::stdin();
    let mut out = std::io::stdout();
    for line in stdin.lock().lines() {
        let line = line.unwrap_or_default();
        if line.is_empty() {
            continue;
        }
        if line.starts_with("quit") {
            break;
        }
        let mut parts = line.split_whitespace();
        let threads: i32 = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let op = parts.next().unwrap_or("").to_string();
        let path = parts.next().unwrap_or("").to_string();
        if threads < 1 || path.is_empty() {
            eprintln!("error command");
            std::process::exit(1);
        }
        if let Err(err) = retune(threads) {
            eprintln!("error {err}");
            std::process::exit(1);
        }
        let sys = frames.entry(path.clone()).or_insert_with(|| {
            System::new().unwrap_or_else(|err| {
                eprintln!("error {err}");
                std::process::exit(1);
            })
        });
        let ready = loaded.get(&path).copied().unwrap_or(false);
        if op != "read" && !ready {
            if let Err(err) = sys.read(&path, 1, -1) {
                eprintln!("error {err}");
                std::process::exit(1);
            }
            loaded.insert(path.clone(), true);
        }
        let t0 = Instant::now();
        let result = match op.as_str() {
            "read" => sys.read(&path, 1, -1).map(|_| {
                loaded.insert(path.clone(), true);
                sink += f64::from(sys.count());
            }),
            "knn" => sys.knn(4, 3.5, true).and_then(|nb| {
                sink += nb.nnz()? as f64;
                Ok(())
            }),
            "steinhardt" => sys.steinhardt(3.5, 6, true).map(|(ql, bar)| {
                let extra = bar.as_ref().map(|v| v[0]).unwrap_or(0.0);
                sink += ql[0] + extra;
            }),
            "chill" => sys.chill(3.5, 0, 0, 0).map(|(_labels, counts)| {
                sink += f64::from(counts[1]);
            }),
            "rings" => sys.rings(3.5, 6).map(|rings| {
                sink += f64::from(rings.count());
            }),
            _ => {
                eprintln!("error {op}");
                std::process::exit(1);
            }
        };
        let ms = t0.elapsed().as_secs_f64() * 1e3;
        if let Err(err) = result {
            eprintln!("error {err}");
            std::process::exit(1);
        }
        writeln!(out, "{ms:.9}").ok();
        out.flush().ok();
    }
    eprintln!("sink {sink:.6}");
}
