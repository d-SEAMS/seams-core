#!/bin/sh
# Copyright (c) d-SEAMS developers
# Run a few frozen mdrun steps and classify the trajectory.
set -e
case $0 in
  */*) root=${0%/*} ;;
  *) root=. ;;
esac
root=$(CDPATH= cd -- "$root" && pwd)

tool=$1
libdir=$2
# mktemp and gmx abort when libasan is preloaded. The classifier below
# is not instrumented; libseams brings the runtime in when it is mapped.
unset LD_PRELOAD
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
awk 'NR==1{print; next} NR==2{n=$1; print; left=n; next} left>0{print; left--; next} {print; exit}' \
  "$root/ice.gro" > "$work/ice.gro"
cp "$root/topol.top" "$root/md.mdp" "$work/"
(
  cd "$work"
  env -u LD_LIBRARY_PATH -u LIBRARY_PATH \
    OMP_NUM_THREADS=1 GMX_MAXBACKUP=-1 \
    gmx grompp -f md.mdp -c ice.gro -p topol.top -o ice.tpr -maxwarn 2
  env -u LD_LIBRARY_PATH -u LIBRARY_PATH \
    OMP_NUM_THREADS=1 GMX_MAXBACKUP=-1 \
    gmx mdrun -s ice.tpr -deffnm ice -nsteps 4 -ntmpi 1
)
# libseams is instrumented and libgromacs is not, so libasan is not the
# first library this process maps. The tool itself is built without the
# sanitizer: an instrumented main makes gmx_fio_open jump to address 0.
env -u LIBRARY_PATH \
  ASAN_OPTIONS="verify_asan_link_order=0${ASAN_OPTIONS:+:$ASAN_OPTIONS}" \
  LD_LIBRARY_PATH="$libdir${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  "$tool" -f "$work/ice.xtc" -cutoff 3.5
