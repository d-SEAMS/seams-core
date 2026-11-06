python3 paper/tex/figures/worked_example.py \
  build/bindings/libseams_lang.so worked
cp worked/worked_example.tsv $D/worked_example.tsv
python3 paper/tex/figures/worked_update.py $D/worked_update.tsv
for f in ice-Ih ice-Ic; do
  echo "# seams cages --signature hc+ddc $f.lammpstrj"
  build/seams cages --signature hc+ddc worked/$f.lammpstrj
done > $D/ideal-cages.txt
for f in input/traj/genice_sI.lammpstrj input/traj/genice_sII.lammpstrj; do
  echo "# seams cages --signature 512+51262+51264 $f"
  build/seams cages --signature 512+51262+51264 $f
done >> $D/ideal-cages.txt
