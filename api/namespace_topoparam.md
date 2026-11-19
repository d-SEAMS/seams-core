<!-- https://docs.dseams.info/api/namespace_topoparam.md -->

# namespace topoparam {#namespace-topoparam}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Functions {#global-functions}

Calculates the height%, an average measure of filled volume. The average height of a prism can be taken to be 2.75-2.85 Angstrom. (Koga et. al., 2001)

Calculates the coverage area%, an area-based measure of relative proportion of monolayer ices.

Calculates the projected area on the XY, YZ and XZ planes.

Per-atom Rodger F4 = mean of cos(3 phi) over H-O…O-H dihedrals of this oxygen with its oxygen neighbours. Hydrogens share molID with their oxygen. An oxygen with no hydrogens (mW) is quiet_NaN. nList is by atom ID with the leading self entry, as neighListO.

Mean of the finite per-atom F4 values; quiet_NaN when none are finite.

Jump-rotor tau_90 from two frames. Each water is the H-H vector of the two hydrogens that share the oxygen molID. Returns dt when any molecule rotates by at least 90 degrees, quiet_NaN otherwise.

Literature I_sd reference: CHILL+ molecules binned by coordinate. TUM/rings stacking: only HC-basal and DDC-equatorial six-rings vote. Prismatic/peripheral rings and non-plane cages (clathrate 5^12) do not write a letter. Ions stay off the ring graph the caller passes. phiC is N_equatorial / (N_basal + N_equatorial) over plane rings.
