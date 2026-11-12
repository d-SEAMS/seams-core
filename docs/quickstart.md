# Quickstart

Build the C library and the API test:

```bash
meson setup build --wrap-mode=nofallback -Dwith_tests=false
meson compile -C build test_seams_api
./build/test_seams_api
```

`libseams.so` exports `seams.h`. The language fronts and the engine hooks link
that library and do not include the C++ headers.

```bash
export SEAMS_LIB_DIR=$PWD/build/src
export LD_LIBRARY_PATH=$SEAMS_LIB_DIR:$LD_LIBRARY_PATH
export SEAMS_LIB=$SEAMS_LIB_DIR/libseams.so
```

Then run one front:

```bash
python -m dseams
cargo run --manifest-path bindings/rust/Cargo.toml --example ice
gfortran -Ibindings/fortran bindings/fortran/test_seams.f90 -L$SEAMS_LIB_DIR -lseams
julia --project=bindings/julia/Seams bindings/julia/Seams/examples/ice.jl
```

A frame is a coordinate array and a cell. Orthorhombic boxes pass three
lengths. Triclinic boxes pass the nine LAMMPS numbers. CHILL+ accepts the
published rule and the OVITO rule. q4 and q6 are the local Steinhardt
parameters. Rings are primitive rings with a symmetry label on each ring.
