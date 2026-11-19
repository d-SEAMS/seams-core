<!-- https://docs.dseams.info/lang/tasks.md -->

# The same three tasks {#the-same-three-tasks}

This repository carries two of the language fronts: C
(`include/seams.h`) and C++ (`bindings/cxx/seams_lang.h`). Each page
below reads the shipped Ice Ih frame `input/traj/ice1h.lammpstrj`
(every atom, cutoff 3.5 Å, the cutoff graph) and prints one line: the
directed neighbour pairs, the primitive ring histogram for sizes 0
through 6, and the nine CHILL+ label counts.

The Python, R, Julia, Rust, Fortran and Lua bindings run the same three
tasks in their own repositories, against the seams-core tag with the
same version. [Language APIs](index.md) lists them.

[C](task-c.md)[C++](task-cpp.md)
