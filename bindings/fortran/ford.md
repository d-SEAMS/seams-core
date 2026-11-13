---
project: d-SEAMS
summary: Fortran bindings for the d-SEAMS C API
project_github: https://github.com/d-SEAMS/seams-core
project_website: https://docs.dseams.info/
author: Rohit Goswami
author_description: d-SEAMS developers
license: MIT
version: 2.14.0
src_dir: ./src
exclude_dir: ./example
             ./test
exclude: test_seams.f90
         seams_support.f90
output_dir: ../../docs/ford
preprocess: false
display: public
         protected
source: true
warn: true
graph: false
search: true
coloured_edges: true
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
---

`seams` is the module a program uses. `seams_c` is the `iso_c_binding`
image of `include/seams.h`. `seams_kinds` holds the `iso_fortran_env`
kinds. `seams_support` is the reference-counted handle box and is not
part of the public API.

Coordinates are `real(seams_real)` with shape `(3, n)`. A missing
`stat` and a missing `errmsg` make a failed call stop with the library
message. The last finalizer of a handle destroys the C object.
Assignment shares that object.

Build with fpm from this directory, or pass the sources to gfortran in
the order `seams_kinds.f90`, `seams_c.f90`, `seams_support.f90`,
`seams.f90`, `seams_methods.f90`. Link `-lseams`.
