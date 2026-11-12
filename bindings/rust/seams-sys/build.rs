fn main() {
    let dir = std::env::var("SEAMS_LIB_DIR").expect("SEAMS_LIB_DIR points at libseams");
    println!("cargo:rustc-link-search=native={dir}");
    println!("cargo:rustc-link-lib=dylib=seams");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{dir}");
    println!("cargo:rerun-if-env-changed=SEAMS_LIB_DIR");
}
