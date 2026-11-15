# The brine front PDF is not in git. paper/tex/repro/brine_front.sh
# runs rule brine_front_pdf.
meson setup build -Dbuildtype=release -Dwith_tests=true -Dwith_cli=true \
  -Dwith_linkcell=enabled -Dwith_minimage=enabled -Dwith_gpulite=disabled \
  -Dwith_ira=disabled -Dwith_sphericart=disabled
meson compile -C build
meson test -C build
# Every timed process runs pinned to this core, with one thread.
CORE=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export RAYON_NUM_THREADS=1 OVITO_THREAD_COUNT=1
D=paper/tex/figures/data
