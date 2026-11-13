# shape #

## About 

This code was originally developed by Prof. David Avnir and Prof. Mark Pinsky at the Hebrew University of Jerusalem. 
The code calculates the degree of a general reference shape in a studied shape, and particularly the degree of a given ideal polyhedral shape, in a distorted polyhedron. The continuous shape measure coincides with the continuous symmetry measure for the Platonic polyhedra which are composed of the minimal number of vertices which defines it. 

* An online calculator is available at: https://csm.ouproj.org.il. 
* A related code developed by Prof. Santiago Alvarez and coworkers is available <a href="http://www.ee.ub.edu/index.php?option=com_jdownloads&Itemid=529&view=viewcategory&catid=4">here</a>. 


## Citations
Pinsky M., Avnir D., <a href="https://pubs.acs.org/doi/10.1021/ic9804925"> "Continuous symmetry measures. 5. The classical polyhedral"</a>, Inorganic Chemistry 37: 5575-5582 (1998).

## Usage

shape input_file case_flag [ref_shape if case_flag is 0] permute_flag [perms_file if permute_flag is 2]  [allow_improper]

### Input Parameters:
* **input_file** - The molecular structure  
* **case_flag** - The shape to compare to (0 for a user defined reference shape)  
* **ref_shape** (optional) - The reference shape - only if case_flag = 0  
* **permute_flag** -   
  * 0 for enumeration on all permutations  
  * 1 for identity perm only  
  * 2 for enumeration on permutations in a given permutation file  
* **perms_file** (optional) - A permutations file, only valid if permute_flag = 2.  
The file is of the format: num_permutations + a list of permutations  
The permutations are performed on the reference shape (ref_shape)  
* **allow_improper** (optional) -   
  * 0 - doesn't allow improper rotations (default)   
  * 1 - allow improper rotations  

### Reference shapes:
1. Tetrahedron 
2. Tetrahedron with Center 
3. Bipyramide 
4. Bipyramide with Center
5. Octahedron
6. Octahedron with Center
7. Cube
8. Cube with Center
9. Icosahedron
10. Icosahedron with Center
11. Dodecahedron
12. Dodecahedron with Center
13. Buckyball (Fullerene, C60)
14. Buckyball with Center (Fullerene, C60)
15. Square with Center
16. Bypiramode (Equidistance)
17. Bypiramode with Center (Equidistance)
18. Trigonal Prism (D3H symmetry) 
19. Trigonal Prism with Center (D3H symmetry)
20. Trigonal Equlateral Prism
21. Trigonal Equlateral Prism with Center

## Contact ##

For questions about the code feel free to use the CoSyM website users group at: https://groups.google.com/g/csm-openu. 

## License ##
This project is provided under the GNU-GPL v2 license. Look at `LICENSE.MD` for more information.

