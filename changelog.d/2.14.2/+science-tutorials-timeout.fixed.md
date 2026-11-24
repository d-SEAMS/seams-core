`docs_science_tutorials` keeps a 180 s timeout. An AddressSanitizer build runs the test serially instead of allowing 900 s: under contention the tutorial took 212 s, and alone it took 22.6 s.
