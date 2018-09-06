# structureFactor
This C++ program reads in lammps trajectory files
or XYZ files to calculate RDF, in-plane RDF and the structure factor.

\brief Post-processing code for RDF, in-plane RDF, structure factor

# Compilation

To generate the executable *runme* use:

```bash
make 
```

An executable *runme* will be formed in the base directory. 

# Running

```bash
./runme
``` 

# Details

* For example, if you want to find the in-plane RDF in the
XY plane, you should use the Rdf2D class. The equation used for 2D-RDF for the \f$n^{th}\f$ layer is:
\f[
 g^n(r) = \frac{1}{(\rho^n)^2 A \delta z} \Sigma_{i \neq j} \delta(r - r_{ij}) \left[ \Theta\left( \frac{\delta z}{2}-|z_i-z^2| \right) \times \Theta\left( \frac{\delta z}{2}-|z_j-z^n| \right) \right] 
\f]
  For detailed instructions, see the Rdf2D class documentation
 
