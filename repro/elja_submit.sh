#!/usr/bin/env bash
# Reproducibility campaign driver for an IRHPC-style Slurm cluster.
#
#   repro/elja_submit.sh prep     # login node: env, wraps, baseline, hq
#   repro/elja_submit.sh submit   # sbatch the exclusive campaign
#   repro/elja_submit.sh run      # the sbatch body (hq server + snakemake)
#
# Environment:
#   ELJA_ACCOUNT    Slurm account for submit (required for submit)
#   ELJA_PARTITION  partition (default 64cpu_256mem)
#   BASE_SHA        baseline commit (default: repro/config.yaml value)
#   HQ_VERSION      HyperQueue release to fetch (default v0.19.0)
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(dirname "$SCRIPT_DIR")
BASE_SHA=${BASE_SHA:-$(grep '^base_sha' "$SCRIPT_DIR/config.yaml" | sed 's/.*"\(.*\)".*/\1/')}
BASE_TREE=$ROOT/../seams-base-repro
HQ_VERSION=${HQ_VERSION:-v0.19.0}
export PATH=$HOME/.pixi/bin:$ROOT/repro/bin:$PATH

case "${1:-}" in
prep)
  cd "$ROOT"
  pixi install -e repro
  pixi run -e repro -- python -m pip install --no-deps pydseamslib==2.6.0
  pixi run -e repro -- meson subprojects download || true
  git rev-parse HEAD > .tip_sha
  if [ ! -d "$BASE_TREE" ]; then
    git worktree add --detach "$BASE_TREE" "$BASE_SHA"
  fi
  (cd "$BASE_TREE" &&
   pixi run -e repro --manifest-path "$ROOT/pixi.toml" -- \
     meson subprojects download || true)
  # Compute nodes are offline; the figshare deposits download here
  pixi run -e repro -- python repro/scripts/figshare_demos.py fetch repro/figshare
  # Lua module lives in yodaStruct 2.x. A sibling named yodaStruct may
  # be a symlink to the v1 baseline tree (example_lua only).
  YODA=${YODASTRUCT_ROOT:-$ROOT/../yodaStruct}
  if [ ! -d "$YODA/lua" ]; then
    YODA=$ROOT/../yodaStruct-2.6
    if [ ! -d "$YODA/lua" ]; then
      git clone --depth 1 --branch v2.6.0 https://github.com/d-SEAMS/yodaStruct.git "$YODA"
    fi
  fi
  mkdir -p "$YODA/subprojects"
  if [ -e "$YODA/subprojects/seams-core" ] && [ ! -L "$YODA/subprojects/seams-core" ]; then
    rm -rf "$YODA/subprojects/seams-core"
  fi
  ln -sfn "$ROOT" "$YODA/subprojects/seams-core"
  export YODASTRUCT_ROOT=$YODA
  # HyperQueue is a single static binary
  if ! command -v hq > /dev/null; then
    mkdir -p repro/bin
    curl -sL "https://github.com/It4innovations/hyperqueue/releases/download/${HQ_VERSION}/hq-${HQ_VERSION}-linux-x64.tar.gz" |
      tar -xz -C repro/bin
  fi
  echo "prep done: tip $(cat .tip_sha), base tree $BASE_TREE, hq $(hq --version)"
  ;;
submit)
  : "${ELJA_ACCOUNT:=chem-ui}"
  cd "$ROOT"
  sbatch --partition="${ELJA_PARTITION:-64cpu_256mem}" --exclusive \
    --ntasks=1 --cpus-per-task=32 --hint=nomultithread --time=04:00:00 \
    --mem=32G --account="$ELJA_ACCOUNT" --job-name=seams-repro \
    --output=repro-%j.out --wrap "repro/elja_submit.sh run"
  ;;
run)
  cd "$ROOT"
  mkdir -p repro/results
  if [ -d "$ROOT/../yodaStruct/lua" ]; then
    export YODASTRUCT_ROOT=$ROOT/../yodaStruct
  elif [ -d "$ROOT/../yodaStruct-2.6/lua" ]; then
    export YODASTRUCT_ROOT=$ROOT/../yodaStruct-2.6
  fi
  # One HyperQueue server and one worker own the allocation; Snakemake
  # submits every heavy rule through it
  export HQ_SERVER_DIR=$ROOT/repro/results/hq-server
  mkdir -p "$HQ_SERVER_DIR"
  hq server start > repro/results/hq-server.log 2>&1 &
  HQ_SERVER_PID=$!
  sleep 3
  # Tasks inherit the worker environment, so the worker starts inside the
  # repro pixi environment; meson, ninja and python resolve there
  pixi run -e repro -- hq worker start --cpus "${SLURM_CPUS_PER_TASK:-32}" \
    > repro/results/hq-worker.log 2>&1 &
  HQ_WORKER_PID=$!
  sleep 2
  trap 'hq server stop > /dev/null 2>&1 || true; kill $HQ_WORKER_PID $HQ_SERVER_PID 2> /dev/null || true' EXIT

  # Node-local build root: the cluster NFS clock skews against the nodes,
  # which meson refuses at configure time
  export SEAMS_BUILD_ROOT=/tmp/seams-repro-${SLURM_JOB_ID:-manual}
  mkdir -p "$SEAMS_BUILD_ROOT"
  LOAD=$(cut -d' ' -f1 /proc/loadavg)
  echo "loadavg_before_workflow: $LOAD"
  pixi run -e repro -- snakemake -s repro/Snakefile --cores all --printshellcmds
  echo "manifest: repro/results/paper_manifest.json"
  ;;
*)
  echo "usage: $0 {prep|submit|run}" >&2
  exit 2
  ;;
esac
