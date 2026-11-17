The vendored fmt 6 in chemfiles 0.10.4 counts format strings without `std::char_traits`, so the macOS build compiles under libc++ 19.
