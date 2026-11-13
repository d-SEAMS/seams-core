Clang builds record the directory of clang's `libomp.so` as a runpath when it lies outside `/usr`, so the libraries, tests and libFuzzer harnesses start without `LD_LIBRARY_PATH`.
