use std::env;
use std::fs;
use std::path::{Path, PathBuf};

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=sources.txt");
    println!("cargo:rerun-if-env-changed=SEAMS_LIB_DIR");
    println!("cargo:rerun-if-env-changed=EIGEN3_INCLUDE_DIR");

    // A prebuilt libseams from a meson build tree. The in-tree tests use this.
    if let Ok(dir) = env::var("SEAMS_LIB_DIR") {
        println!("cargo:rustc-link-search=native={dir}");
        println!("cargo:rustc-link-lib=dylib=seams");
        println!("cargo:rustc-link-arg=-Wl,-rpath,{dir}");
        return;
    }

    // rustdoc does not link, and docs.rs has no business compiling Eigen.
    if env::var_os("DOCS_RS").is_some() {
        return;
    }

    if env::var_os("CARGO_FEATURE_SYSTEM").is_some() {
        let lib = pkg_config::Config::new()
            .atleast_version(&env::var("CARGO_PKG_VERSION").unwrap())
            .probe("seams")
            .unwrap_or_else(|e| panic!("feature `system` needs libseams through pkg-config: {e}"));
        for dir in &lib.include_paths {
            println!("cargo:include={}", dir.display());
        }
        return;
    }

    build_vendored();
}

fn source_root(manifest: &Path) -> PathBuf {
    let vendored = manifest.join("vendor");
    if vendored.join("src/seams_api.cpp").is_file() {
        return vendored;
    }
    let checkout = manifest.join("../../..");
    if checkout.join("src/seams_api.cpp").is_file() {
        return checkout;
    }
    panic!(
        "no C++ sources: run `python3 bindings/rust/seams-sys/vendor.py`, \
         set SEAMS_LIB_DIR to a directory holding libseams, or enable the `system` feature"
    );
}

fn eigen_include(root: &Path) -> Vec<PathBuf> {
    let vendored = root.join("eigen");
    if vendored.join("Eigen/Core").is_file() {
        return vec![vendored];
    }
    if let Ok(dir) = env::var("EIGEN3_INCLUDE_DIR") {
        return vec![PathBuf::from(dir)];
    }
    match pkg_config::Config::new()
        .atleast_version("3.4.0")
        .cargo_metadata(false)
        .probe("eigen3")
    {
        Ok(lib) => lib.include_paths,
        Err(e) => panic!("Eigen >= 3.4.0 not found (vendor/eigen, EIGEN3_INCLUDE_DIR, pkg-config): {e}"),
    }
}

fn build_vendored() {
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    if target_arch != "x86_64" {
        panic!(
            "the bundled libseams needs x86_64 (shell_classify.cpp uses AVX2 intrinsics); \
             on {target_arch} install libseams and enable the `system` feature"
        );
    }
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let root = source_root(&manifest);
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());

    let template = fs::read_to_string(root.join("include/seams_version.h.in")).unwrap();
    let header = template
        .replace("@MAJOR@", &env::var("CARGO_PKG_VERSION_MAJOR").unwrap())
        .replace("@MINOR@", &env::var("CARGO_PKG_VERSION_MINOR").unwrap())
        .replace("@PATCH@", &env::var("CARGO_PKG_VERSION_PATCH").unwrap())
        .replace(
            "@STRING@",
            &format!(
                "{}.{}.{}",
                env::var("CARGO_PKG_VERSION_MAJOR").unwrap(),
                env::var("CARGO_PKG_VERSION_MINOR").unwrap(),
                env::var("CARGO_PKG_VERSION_PATCH").unwrap()
            ),
        );
    let generated = out.join("include");
    fs::create_dir_all(&generated).unwrap();
    fs::write(generated.join("seams_version.h"), header).unwrap();

    let mut core = Vec::new();
    let mut avx2 = Vec::new();
    let mut avx512 = Vec::new();
    for line in fs::read_to_string(manifest.join("sources.txt")).unwrap().lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (bucket, rel) = match line.split_once(':') {
            Some(("avx2", rel)) => (&mut avx2, rel.trim()),
            Some(("avx512", rel)) => (&mut avx512, rel.trim()),
            Some((tag, _)) => panic!("sources.txt: unknown tag {tag}"),
            None => (&mut core, line),
        };
        let path = root.join(rel);
        println!("cargo:rerun-if-changed={}", path.display());
        bucket.push(path);
    }

    let eigen = eigen_include(&root);
    let base = || {
        let mut build = cc::Build::new();
        build
            .cpp(true)
            .std("c++20")
            .opt_level(2)
            .warnings(false)
            .flag("-w")
            .flag("-fopenmp")
            .define("SEAMS_HAS_OPENMP", "1")
            .define("EIGEN_MPL2_ONLY", None)
            .include(root.join("src/include/internal"))
            .include(root.join("include"))
            .include(&generated)
            .includes(&eigen);
        if env::var_os("CARGO_CFG_DEBUG_ASSERTIONS").is_none() {
            build.define("NDEBUG", None);
        }
        build
    };

    let mut objects = Vec::new();
    let mut simd = base();
    simd.flag("-mavx2").flag("-mfma").files(&avx2);
    objects.extend(simd.compile_intermediates());
    let mut wide = base();
    wide.flag("-mavx512f").flag("-mavx512dq").flag("-mfma").files(&avx512);
    objects.extend(wide.compile_intermediates());

    let mut lib = base();
    lib.define("SEAMS_YLM_AVX512", "1").files(&core).objects(&objects);
    let compiler = lib.get_compiler();
    lib.compile("seams_bundled");

    if compiler.is_like_clang() {
        println!("cargo:rustc-link-lib=omp");
    } else {
        println!("cargo:rustc-link-lib=gomp");
    }
    println!("cargo:include={}", root.join("include").display());
    println!("cargo:include_generated={}", generated.display());
}
