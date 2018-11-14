# structureFactor
This C++ program reads in lammps trajectory files
or XYZ files to calculate RDF, in-plane RDF and the structure factor.

\brief Post-processing code for RDF, in-plane RDF, structure factor

# Compilation

## Dependency Management
Unless you are sure that all the libraries necessary are available on your
system, it is expedient to ensure `conan` is available for managing
dependencies. A `Pipfile` is also provided for ease of usage. This is preferred
since it will ensure a reproducable environment which does not interfere with
the system `python`.

``` bash
# This will get you conan in a local environment
pipenv shell
```

## Build Things
To generate the executable *yodaStruct* in `bin/` use:

```bash
# Always prefer an out of tree build
mkdir build
cd build
cmake ..
make
```

# Running
To run the sample inputs, simply move the binary to the project root, or to a
directory where `input/` is a child directory.

```bash
# Assuming you are in the project root
cp bin/yodaStruct .
./yodaStruct
``` 

# Details

* For example, if you want to find the in-plane RDF in the
XY plane, you should use the Rdf2D class. The equation used for 2D-RDF for the \f$n^{th}\f$ layer is:
\f[
 g^n(r) = \frac{1}{(\rho^n)^2 A \delta z} \Sigma_{i \neq j} \delta(r - r_{ij}) \left[ \Theta\left( \frac{\delta z}{2}-|z_i-z^2| \right) \times \Theta\left( \frac{\delta z}{2}-|z_j-z^n| \right) \right] 
\f]
  For detailed instructions, see the Rdf2D class documentation
 
# Acknowledgements
The following libraries and tools are used in this project:
- [CMake](https://cmake.org/) for compilation ([cmake-init](https://github.com/cginternals/cmake-init) was used as a reference)
- [Doxygen](https://www.doxygen.org) for the developer API
