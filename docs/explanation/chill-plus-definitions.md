# Two CHILL+ definitions

`seams_chill_plus` takes a definition. The default is the paper. The
other is the OVITO Chill+ modifier. The numbers -0.8, -0.35, and 0.25
stay. The paper closes those edges and averages \(q_{lm}\) over the four
closest neighbours. The modifier opens the edges and averages \(q_{lm}\)
over every neighbour inside 3.5 A.

## Paper

Nguyen and Molinero, J. Phys. Chem. B 2015, 119, 9369, equation (1),
average the spherical harmonic over the four closest neighbours:

    q_lm(i) = (1/4) sum_{j=1}^{4} Y_lm(r_ij)

The text calls those "the four closest neighbours". Table 2 then closes
both ends of the windows: an eclipsed bond is
\(0.25 \ge c(i,j) \ge -0.35\), and a staggered bond is \(c(i,j) \le -0.8\).
Footnote a of Table 1 is the neighbour test for interfacial ice. Pass
`SEAMS_CHILL_PUBLISHED`, or `rule="paper"` in Python. That is the default.

## OVITO

The OVITO User Manual 3.17.2, Chill+ modifier, defines the neighbours
by distance rather than by a count of four. It says the algorithm
considers "two neighboring water molecules (O-O distance less than
3.5 A)", and that "the structure is determined by whether the
correlation is between 0.25 and -0.35 (eclipsed) or less than -0.8
(staggered)". Interfacial ice is "4 bonds of which 2 or 3 staggered",
with hexagonal preferred when the bonds are 3 staggered and 1 eclipsed.
There is no Table 1 neighbour test.

Pass `SEAMS_CHILL_OVITO`, or `rule="ovito"`. \(q_{lm}\) is then the sum
over every neighbour inside the cutoff. On a crystal where that shell
already holds four neighbours, the two definitions see the same bonds.
They split when a neighbour of a four-coordinated oxygen has a fifth
oxygen inside 3.5 A: the paper keeps the four closest, and the modifier
folds the extra oxygen into \(q_{lm}\).

```c
seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts);
seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, labels, counts);
```

```python
labels, counts = dseams.chill_plus(positions, box, rule="paper")
labels, counts = dseams.chill_plus(positions, box, rule="ovito")
```
