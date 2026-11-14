The same three tasks
====================

The tree has eight language fronts: C (``include/seams.h``), C++
(``bindings/cxx/seams_lang.h``), Python (``bindings/python``), Rust
(``bindings/rust``), Julia (``bindings/julia/Dseams``), Lua
(``bindings/lua``), Fortran (``bindings/fortran``), and R
(``bindings/r``). Each page below reads the shipped Ice Ih frame
``input/traj/ice1h.lammpstrj`` (every atom, cutoff 3.5 Å, the cutoff
graph) and prints one line: the directed neighbour pairs, the primitive
ring histogram for sizes 0 through 6, and the nine CHILL+ label counts.

A missing compiler or an unbuilt binding exits 77, which Meson records
as a skip. The census of the languages that did run is generated from
those runs.

.. include:: parity.rst

.. toctree::
   :maxdepth: 1

   task-c
   task-cpp
   task-python
   task-rust
   task-julia
   task-lua
   task-fortran
   task-r
