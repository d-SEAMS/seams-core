# Timed calls of a case. Both scales use at least six calls, with no best of them.
REPEATS=6
REPEATS_LARGEST=6
# Seconds before a case of a classifier or of the CHILL+ library, a ring case
# and one run of RINGS stop.
CASE_TIMEOUT_S=900
RING_CASE_TIMEOUT_S=180
RINGS_RUN_TIMEOUT_S=180
python3 -m venv h2h-venv
. h2h-venv/bin/activate
python -m pip install numpy==2.5.4 scipy==1.18.1 ovito==3.17.2 \
  freud-analysis==3.6.1 pyscal3==4.1.0 matscipy==1.3.1 \
  GenIce2==2.2.13.3 cycless==0.7
git clone https://github.com/palmergroupUH/Order_Parameters_Library ext/palmer
git -C ext/palmer checkout ee8a3cb05e3dcbc975cca1ad8820570871973d12
git -C ext/palmer apply "$PWD/$D/h2h/palmer-chill-gnu-fortran.patch"
make -C ext/palmer
PALMER=$PWD/ext/palmer/order_parameter/lib/lib_CHILL.so
tar -xzf rings-code-v1.3.5.tar.gz -C ext
(cd ext/rings-code-v1.3.5 && ./configure FC=gfortran && make -j1)
RINGS=$PWD/ext/rings-code-v1.3.5/src/rings
python scripts/h2h_campaign.py prepare --data h2h-data \
  --liquid liquid.lammpstrj --confined dump-240-square.lammpstrj
for system in ic ih liquid confined; do
  for scale in 1e3 1e4 1e5 1e6; do
    repeats=$REPEATS; test $scale = 1e6 && repeats=$REPEATS_LARGEST
    H2H_CORE=$CORE H2H_REPEATS=$repeats python scripts/h2h_campaign.py \
      case $system $scale --data h2h-data --output h2h-out \
      --bench build/tests/bench_paper --palmer $PALMER \
      --timeout $CASE_TIMEOUT_S --ring-timeout $RING_CASE_TIMEOUT_S \
      --reference-timeout $CASE_TIMEOUT_S
  done
done
for system in ic ih liquid confined; do
  for scale in 1e3 1e4 1e5 1e6; do
    H2H_CORE=$CORE H2H_REPEATS=$REPEATS python scripts/h2h_campaign.py \
      rings $system $scale --data h2h-data --output h2h-out \
      --bench build/tests/bench_paper --rings $RINGS \
      --timeout $RINGS_RUN_TIMEOUT_S
  done
done
python scripts/h2h_campaign.py summarize --output h2h-out
python scripts/h2h_campaign.py environment --output h2h-out \
  --bench build/tests/bench_paper --library build/src/libyodaLib.so \
  --source . --palmer-source ext/palmer --palmer-library $PALMER \
  --rings-archive rings-code-v1.3.5.tar.gz --rings $RINGS
