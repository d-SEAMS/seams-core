#!/bin/sh
# Release comparison on one core. Twelve calls on the 110592-molecule frame
# and six on the 995328-molecule frame. Writes the JSON record and the TeX
# the paper inputs.
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OVITO_THREAD_COUNT=1
export SEAMS_BASE="${SEAMS_BASE:-seams-base/bbdir-release}"
export SEAMS_TIP="${SEAMS_TIP:-bbdir-release}"
export SEAMS_1E6="${SEAMS_1E6:-v37-frames/ice-1e6.lammpstrj}"
PY="${PYTHON:-python3}"
"$PY" paper/tex/repro/shipped_campaign.py
"$PY" paper/tex/repro/check_shipped_numbers.py
