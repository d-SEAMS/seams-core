# WASM classify

The C ABI is `seams_chill_plus` in `src/include/internal/seams_c_api.h`.
It is the same CHILL+ path as `seams chill-plus` (mutual 4-NN).

Emscripten, when the host has `emcc`:

```
meson setup bwasm --cross-file wasm/emscripten.ini \
  -Dwith_tests=false -Dwith_cli=false -Dwith_python=false \
  -Dwith_openmp_offload=disabled -Dwith_mpi=disabled
meson compile -C bwasm
```

The paper viewer (`dseams2_paper/viewer/index.html`) is the chemiscope
widget Jayant can open in a browser. It loads precomputed frames today.
The WASM module plugs the same `labels[]` integers into that JSON when
the emscripten build exists.
