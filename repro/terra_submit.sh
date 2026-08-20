#!/usr/bin/env bash
# Not the paper campaign. Paper timings are Elja
# (repro/elja_submit.sh, scripts/elja_paper_benches.sh, repro/elja_gpu_job.sh).
# Terra Slurm cluster helper only (partition cpu).
#
#   repro/terra_submit.sh prep
#   repro/terra_submit.sh submit
#   repro/terra_submit.sh run      # sbatch body
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(dirname "$SCRIPT_DIR")
BASE_SHA=${BASE_SHA:-$(grep '^base_sha' "$SCRIPT_DIR/config.yaml" | sed 's/.*"\(.*\)".*/\1/')}
BASE_TREE=$ROOT/../seams-base-repro
SOURCE_ROOT=${DSEAMS_SOURCE_ROOT:-$ROOT/../dseams-repro-sources}
HQ_VERSION=${HQ_VERSION:-v0.19.0}
export PATH=$HOME/.pixi/bin:$ROOT/repro/bin:$PATH
export SLURM_CONF=${SLURM_CONF:-/etc/slurm-llnl/slurm.conf}
export DSEAMS_SOURCE_ROOT=$SOURCE_ROOT

case "${1:-}" in
prep)
  cd "$ROOT"
  pixi install -e repro
  pixi run -e repro -- python repro/scripts/ecosystem_sources.py fetch \
    --root "$SOURCE_ROOT"
  pixi run -e repro -- python repro/scripts/ecosystem_sources.py wire \
    --root "$SOURCE_ROOT" --core "$ROOT"
  pixi run -e repro -- meson subprojects download
  git rev-parse HEAD > .tip_sha
  if [ ! -d "$BASE_TREE" ]; then
    git worktree add --detach "$BASE_TREE" "$BASE_SHA"
  fi
  (cd "$BASE_TREE" &&
   pixi run -e repro --manifest-path "$ROOT/pixi.toml" -- \
     meson subprojects download)
  pixi run -e repro -- python repro/scripts/figshare_demos.py fetch repro/figshare
  if ! command -v hq > /dev/null; then
    mkdir -p repro/bin
    curl -sL "https://github.com/It4innovations/hyperqueue/releases/download/${HQ_VERSION}/hq-${HQ_VERSION}-linux-x64.tar.gz" |
      tar -xz -C repro/bin
  fi
  echo "prep done: tip $(cat .tip_sha), sources $SOURCE_ROOT, base tree $BASE_TREE, hq $(hq --version)"
  ;;
submit)
  cd "$ROOT"
  sbatch --partition="${TERRA_PARTITION:-cpu}" --exclusive \
    --ntasks=1 --cpus-per-task=8 --hint=nomultithread --time=04:00:00 \
    --mem=32G --job-name=seams-repro \
    --output=repro/results/terra-%j.out --wrap "repro/terra_submit.sh run"
  ;;
run)
  cd "$ROOT"
  mkdir -p repro/results
  export HQ_SERVER_DIR=$ROOT/repro/results/hq-server
  mkdir -p "$HQ_SERVER_DIR"
  hq server start > repro/results/hq-server.log 2>&1 &
  HQ_SERVER_PID=$!
  sleep 3
  pixi run -e repro -- hq worker start --cpus "${SLURM_CPUS_PER_TASK:-8}" \
    > repro/results/hq-worker.log 2>&1 &
  HQ_WORKER_PID=$!
  sleep 2
  trap 'hq server stop > /dev/null 2>&1 || true; kill $HQ_WORKER_PID $HQ_SERVER_PID 2>/dev/null || true' EXIT
  export SEAMS_BUILD_ROOT=/tmp/seams-repro-${SLURM_JOB_ID:-manual}
  mkdir -p "$SEAMS_BUILD_ROOT"
  echo "loadavg_before_workflow: $(cut -d' ' -f1 /proc/loadavg)"
  pixi run -e repro -- snakemake -s repro/Snakefile --cores all --printshellcmds
  echo "manifest: repro/results/paper_manifest.json"
  ;;
*)
  echo "usage: $0 {prep|submit|run}" >&2
  exit 2
  ;;
esac
