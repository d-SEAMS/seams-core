#!/bin/sh
# Copyright (c) d-SEAMS developers
# Run a few frozen mdrun steps and classify the trajectory.
set -e
root=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
tool=$1
libdir=$2
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
awk 'NR==1{print; next} NR==2{n=$1; print; left=n; next} left>0{print; left--; next} {print; exit}' \
  "$root/ice.gro" > "$work/ice.gro"
cp "$root/topol.top" "$root/md.mdp" "$work/"
(
  cd "$work"
  env -u LD_LIBRARY_PATH -u LIBRARY_PATH \
    gmx grompp -f md.mdp -c ice.gro -p topol.top -o ice.tpr -maxwarn 2
  env -u LD_LIBRARY_PATH -u LIBRARY_PATH \
    gmx mdrun -s ice.tpr -deffnm ice -nsteps 4
)
env -u LIBRARY_PATH LD_LIBRARY_PATH="$libdir" \
  "$tool" -f "$work/ice.xtc" -cutoff 3.5
