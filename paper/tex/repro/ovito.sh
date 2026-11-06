H2H_CORE=$CORE taskset -c $CORE python scripts/h2h_ovito_rings.py \
  build/tests/bench_paper liquid.lammpstrj dump-240-square.lammpstrj \
  $D/h2h/ovito-rings.tsv
python paper/tex/figures/chill_neighbours.py build/tests/bench_paper \
  liquid.lammpstrj > $D/h2h/chill-neighbours.txt
