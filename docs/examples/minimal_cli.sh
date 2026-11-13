#!/bin/sh
# Read the GenIce sI cell with the seams CLI.
# Usage: minimal_cli.sh SEAMS FIXTURE
set -eu
SEAMS=$1
FIXTURE=$2
NO_COLOR=1 "$SEAMS" read "$FIXTURE" | grep -q "nop 46"
NO_COLOR=1 "$SEAMS" cn "$FIXTURE" --cutoff 3.5 | grep -q "nI 46"
NO_COLOR=1 "$SEAMS" steinhardt "$FIXTURE" --cutoff 3.5 --l 6 | grep -q "q6 0.595"
NO_COLOR=1 "$SEAMS" chill-plus "$FIXTURE" --cutoff 3.5 | grep -q "clathrate 46"
NO_COLOR=1 "$SEAMS" rings "$FIXTURE" --cutoff 3.5 | grep -q "^6 54"
echo "cli example ok"
