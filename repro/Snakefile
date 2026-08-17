# Reproducibility DAG for the d-SEAMS 2.0 paper numbers.
#
# Run from the repository root inside the repro pixi environment, on one
# exclusively allocated node (repro/elja_submit.sh does both):
#
#   pixi run -e repro -- snakemake -s repro/Snakefile --cores all
#
# Heavy rules execute through HyperQueue when an `hq` server is reachable
# (the submit script starts one per allocation); set SEAMS_NO_HQ=1 or run
# without hq on PATH to execute the same commands directly. The baseline
# tree and subproject wraps are prepared on the login node by
# `repro/elja_submit.sh prep`, since compute nodes carry neither git nor
# network access.

import os
import shutil

configfile: "repro/config.yaml"

R = "repro/results"
BASE_TREE = config["base_tree"]
# Build directories default to the working tree but move to node-local
# storage (SEAMS_BUILD_ROOT) on clusters whose NFS clocks skew against
# the nodes, which meson refuses
BUILD = os.path.abspath(os.environ.get("SEAMS_BUILD_ROOT", "."))
TIP_BUILD = os.path.join(BUILD, "build-repro")
BASE_BUILD = os.path.join(BUILD, "build-base")
def _yoda_root():
    env = os.environ.get("YODASTRUCT_ROOT")
    if env:
        return os.path.abspath(env)
    for cand in (os.path.join("..", "yodaStruct"), os.path.join("..", "yodaStruct-2.6")):
        if os.path.isdir(os.path.join(cand, "lua")):
            return os.path.abspath(cand)
    return os.path.abspath(os.path.join("..", "yodaStruct"))


YODA_ROOT = _yoda_root()
YODA_BUILD = os.path.join(BUILD, "build-yoda")


def hq(cpus):
    if os.environ.get("SEAMS_NO_HQ") or shutil.which("hq") is None:
        return ""
    return ("hq submit --wait --cpus {} "
            "--stdout repro/results/hq-%{{JOB_ID}}.stdout "
            "--stderr repro/results/hq-%{{JOB_ID}}.stderr --").format(cpus)


FIGSHARE_DIR = config.get("figshare_dir", "repro/figshare")
FIGSHARE_FILES = [
    "nucleation.lammpstrj",
    "mW_cubic.lammpstrj",
    "dump-240-square.lammpstrj",
    "dump-6-320-310.lammpstrj",
    "dump-320.lammpstrj",
]


rule all:
    input:
        R + "/paper_manifest.json",


rule conditions:
    output:
        R + "/conditions.txt",
    params:
        base=config["base_sha"],
    shell:
        r"""
        {{
          echo "hostname: $(hostname)"
          echo "jobid: ${{SLURM_JOB_ID:-none}}"
          echo "partition: ${{SLURM_JOB_PARTITION:-none}}"
          echo "tip_sha: $(cat .tip_sha 2>/dev/null || echo unknown)"
          echo "base_sha: {params.base}"
          echo "hq: $(hq --version 2>/dev/null || echo absent)"
          lscpu | grep -E 'Model name|Socket|Core|Thread' || true
          echo "loadavg_at_start: $(cut -d' ' -f1 /proc/loadavg)"
        }} > {output}
        """


rule setup_tip:
    output:
        touch(R + "/tip-setup.done"),
    params:
        log=R + "/tip-setup.log",
        bdir=TIP_BUILD,
    shell:
        "meson setup {params.bdir} --prefix=$CONDA_PREFIX --buildtype=release "
        "-Dwith_python=false -Dwith_tests=true > {params.log} 2>&1"


rule build_tip:
    input:
        R + "/tip-setup.done",
    output:
        touch(R + "/tip-build.done"),
    params:
        hq=lambda wc: hq(8),
        bdir=TIP_BUILD,
    shell:
        "{params.hq} meson compile -C {params.bdir}"


rule identity_gate:
    # A red suite aborts the DAG; no table row is emitted from a red tree
    input:
        R + "/tip-build.done",
    output:
        R + "/tip-test.log",
    params:
        hq=lambda wc: hq(8),
        bdir=TIP_BUILD,
    shell:
        "{params.hq} bash -c 'meson test -C {params.bdir} > {output} 2>&1'"


rule install_python:
    # Bindings are pydseamslib on PyPI, not this tree.
    input:
        R + "/tip-test.log",
    output:
        touch(R + "/py-install.done"),
    shell:
        "python -c 'import pydseams as ds; assert ds.__version__ == \"2.6.0\", ds.__version__'"


rule build_base:
    output:
        touch(R + "/base-build.done"),
    params:
        tree=BASE_TREE,
        bbuild=BASE_BUILD,
        hq=lambda wc: hq(8),
    shell:
        r"""
        test -d {params.tree} || {{ echo "baseline tree missing; run repro/elja_submit.sh prep" >&2; exit 1; }}
        meson setup {params.bbuild} {params.tree} --buildtype=release \
          -Dwith_python=false -Dwith_tests=false > repro/results/base-setup.log 2>&1
        {params.hq} meson compile -C {params.bbuild}
        """


rule base_drivers:
    # The baseline predates the bench targets; its drivers compile directly
    # against the baseline library
    input:
        R + "/base-build.done",
    output:
        scaling=R + "/base_bench_scaling",
        cages=R + "/base_bench_cages",
    params:
        tree=BASE_TREE,
        bbuild=BASE_BUILD,
    shell:
        r"""
        EIGEN=$CONDA_PREFIX/include/eigen3
        INC="-I{params.tree}/src/include/internal -I{params.tree}/src/include/external -I$EIGEN"
        LINK="-L{params.bbuild}/src -lyodaLib -Wl,-rpath,{params.bbuild}/src -Wl,--enable-new-dtags -Wl,--allow-shlib-undefined"
        g++ -std=c++20 -O2 -o {output.scaling} tests/bench_scaling.cpp $INC $LINK
        g++ -std=c++20 -O2 -o {output.cages} tests/bench_cages_base.cpp $INC $LINK
        """


rule tip_scaling:
    input:
        R + "/tip-test.log",
    output:
        R + "/tip-scaling.txt",
    params:
        hq=lambda wc: hq(1),
        tbuild=TIP_BUILD,
        sizes=" ".join(str(s) for s in config["scaling_sizes"]),
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "{params.tbuild}/tests/bench_scaling {params.sizes} > ../{output}'"


rule base_scaling:
    input:
        R + "/base_bench_scaling",
    output:
        R + "/base-scaling.txt",
    params:
        hq=lambda wc: hq(1),
        sizes=" ".join(str(s) for s in config["scaling_sizes"]),
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "../{input} {params.sizes} > ../{output}'"


rule tip_cages:
    input:
        R + "/tip-test.log",
    output:
        R + "/tip-cages.txt",
    params:
        hq=lambda wc: hq(1),
        tbuild=TIP_BUILD,
        traj=config["trajectory"],
        reps=config["cage_reps"],
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "{params.tbuild}/tests/bench_cages {params.traj} 1 {params.reps} "
        "> ../{output}'"


rule base_cages:
    input:
        R + "/base_bench_cages",
    output:
        R + "/base-cages.txt",
    params:
        hq=lambda wc: hq(1),
        traj=config["trajectory"],
        reps=config["cage_reps"],
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "../{input} {params.traj} 1 {params.reps} > ../{output}'"


rule tip_overhead:
    input:
        R + "/tip-test.log",
    output:
        R + "/tip-overhead.txt",
    params:
        hq=lambda wc: hq(1),
        tbuild=TIP_BUILD,
        sizes=" ".join(str(s) for s in config["overhead_sizes"]),
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "{params.tbuild}/tests/bench_overhead {params.sizes} > ../{output}'"


rule tip_strong:
    input:
        R + "/tip-test.log",
    output:
        R + "/tip-strong-t{threads}.txt",
    params:
        hq=lambda wc: hq(int(wc.threads)),
        tbuild=TIP_BUILD,
        atoms=config["strong_atoms"],
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS={wildcards.threads} "
        "{params.tbuild}/tests/bench_strong {params.atoms} > ../{output}'"


rule tip_pipeline:
    input:
        R + "/tip-test.log",
    output:
        R + "/tip-pipeline.txt",
    params:
        hq=lambda wc: hq(1),
        tbuild=TIP_BUILD,
        sizes=" ".join(str(s) for s in config["scaling_sizes"]),
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "{params.tbuild}/tests/bench_pipeline {params.sizes} > ../{output}'"


rule tip_incremental:
    input:
        R + "/tip-test.log",
    output:
        R + "/tip-incremental-{n}.txt",
    params:
        hq=lambda wc: hq(1),
        tbuild=TIP_BUILD,
        frames=config["incremental_frames"],
        sites=config["incremental_sites"],
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "{params.tbuild}/tests/bench_rings_incremental {wildcards.n} "
        "{params.frames} {params.sites} > ../{output}'"


rule tip_stages_cubic:
    input:
        R + "/tip-test.log",
    output:
        R + "/tip-stages-cubic.txt",
    params:
        hq=lambda wc: hq(1),
        tbuild=TIP_BUILD,
        traj=config["trajectory"],
        reps=config["cage_reps"],
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "{params.tbuild}/tests/bench_stages {params.traj} 1 1 {params.reps} "
        "> ../{output}'"


rule tip_stages_nucleation:
    input:
        gate=R + "/tip-test.log",
        traj=FIGSHARE_DIR + "/nucleation.lammpstrj",
    output:
        R + "/tip-stages-nucleation.txt",
    params:
        hq=lambda wc: hq(1),
        tbuild=TIP_BUILD,
        reps=config["cage_reps"],
    shell:
        "{params.hq} bash -c 'OMP_NUM_THREADS=1 "
        "{params.tbuild}/tests/bench_stages {input.traj} 1 1 {params.reps} "
        "> {output}'"


rule trajectory_incremental:
    # Exactness referee and steady-state timings for the incremental rings
    # and the incremental cage affiliation; nonzero exit on any inequality
    input:
        R + "/tip-test.log",
    output:
        R + "/trajectory-incremental.txt",
    params:
        hq=lambda wc: hq(1),
        tbuild=TIP_BUILD,
        traj=config["trajectory"],
        frames=config["trajectory_frames"],
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "{params.tbuild}/tests/bench_trajectory_incremental {params.traj} "
        "{params.frames} > ../{output}'"


rule ql_compare_dseams:
    input:
        R + "/tip-test.log",
    output:
        R + "/ql-dseams.txt",
    params:
        hq=lambda wc: hq(1),
        tbuild=TIP_BUILD,
    shell:
        "{params.hq} bash -c 'OMP_NUM_THREADS=1 "
        "{params.tbuild}/tests/compare_structure_desc > {output}'"


rule ql_compare_python:
    input:
        R + "/py-install.done",
    output:
        R + "/ql-python.txt",
    params:
        hq=lambda wc: hq(1),
    shell:
        "{params.hq} bash -c "
        "'LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH "
        "python scripts/compare_ql_literature.py > {output}'"


rule figshare_fetch:
    # Network access; on clusters whose compute nodes are offline this
    # runs during `repro/elja_submit.sh prep` and verifies here as a no-op
    output:
        expand(FIGSHARE_DIR + "/{f}", f=FIGSHARE_FILES),
    shell:
        "python repro/scripts/figshare_demos.py fetch " + FIGSHARE_DIR


rule build_yoda:
    # yodaStruct is require("dseams"), not an in-tree seams-core binary.
    # Point its seams-core wrap at this tip so the Lua module links the
    # same engine the benches just gated.
    input:
        R + "/tip-test.log",
    output:
        touch(R + "/yoda-build.done"),
    params:
        hq=lambda wc: hq(8),
        yoda=YODA_ROOT,
        bdir=YODA_BUILD,
        tip=os.path.abspath("."),
    shell:
        r"""
        test -d {params.yoda}/lua || {{ echo "YODASTRUCT_ROOT missing: {params.yoda}" >&2; exit 1; }}
        mkdir -p {params.yoda}/subprojects
        if [ -e {params.yoda}/subprojects/seams-core ] && [ ! -L {params.yoda}/subprojects/seams-core ]; then
          rm -rf {params.yoda}/subprojects/seams-core
        fi
        ln -sfn {params.tip} {params.yoda}/subprojects/seams-core
        meson setup {params.bdir} {params.yoda} --buildtype=release --wrap-mode=nodownload \
          > repro/results/yoda-setup.log 2>&1
        {params.hq} meson compile -C {params.bdir}
        """


rule figshare_demos:
    # The five 1.0 figshare deposits through require("dseams"); nonzero
    # exit on any failed demo
    input:
        yoda=R + "/yoda-build.done",
        traj=expand(FIGSHARE_DIR + "/{f}", f=FIGSHARE_FILES),
    output:
        R + "/figshare-demos/figshare-demos.json",
    params:
        hq=lambda wc: hq(4),
        yoda=YODA_ROOT,
        bdir=YODA_BUILD,
    shell:
        "{params.hq} bash -c 'OMP_NUM_THREADS=4 "
        "python repro/scripts/figshare_demos.py demos "
        "{params.yoda} {params.bdir} " + FIGSHARE_DIR + " "
        + R + "/figshare-demos'"


NOTEBOOKS = [
    "01_chillplus_ic",
    "02_nucleation_cages",
    "03_nanotube_prisms",
    "04_monolayer_rings",
    "05_rdf2d",
]


rule figshare_notebook:
    # Percent-format sources under repro/notebooks/; jupytext converts
    # and papermill executes. Execution is the test (each notebook
    # asserts its own headline numbers) and the executed .ipynb is the
    # artifact. The nucleation notebook also writes
    # figshare-incremental.json.
    input:
        py=R + "/py-install.done",
        nb="repro/notebooks/{nb}.py",
        runner="repro/scripts/execute_notebook.sh",
        traj=expand(FIGSHARE_DIR + "/{f}", f=FIGSHARE_FILES),
    output:
        R + "/notebooks/{nb}.ipynb",
    params:
        hq=lambda wc: hq(2),
    shell:
        "{params.hq} bash -c "
        "'LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH "
        "OMP_NUM_THREADS=2 repro/scripts/execute_notebook.sh {input.nb} {output}'"


rule figshare_incremental:
    # Produced by the nucleation notebook: per-frame incremental rings
    # (refereed against batch every frame), incremental affiliation, and
    # seeded classification across the whole deposit
    input:
        R + "/notebooks/02_nucleation_cages.ipynb",
    output:
        R + "/figshare-incremental.json",
    shell:
        "test -s {output}"


rule aggregate:
    input:
        conditions=R + "/conditions.txt",
        gate=R + "/tip-test.log",
        tip_scaling=R + "/tip-scaling.txt",
        base_scaling=R + "/base-scaling.txt",
        tip_cages=R + "/tip-cages.txt",
        base_cages=R + "/base-cages.txt",
        overhead=R + "/tip-overhead.txt",
        strong=expand(
            R + "/tip-strong-t{t}.txt", t=config["strong_threads"]
        ),
        trajectory=R + "/trajectory-incremental.txt",
        ql_dseams=R + "/ql-dseams.txt",
        ql_python=R + "/ql-python.txt",
        figshare_demos=R + "/figshare-demos/figshare-demos.json",
        figshare_notebooks=expand(
            R + "/notebooks/{nb}.ipynb", nb=NOTEBOOKS
        ),
        figshare_incremental=R + "/figshare-incremental.json",
        pipeline=R + "/tip-pipeline.txt",
        incremental=expand(
            R + "/tip-incremental-{n}.txt", n=config["scaling_sizes"]
        ),
        stages_cubic=R + "/tip-stages-cubic.txt",
        stages_nucleation=R + "/tip-stages-nucleation.txt",
    output:
        R + "/paper_manifest.json",
    params:
        results=R,
    shell:
        "python repro/scripts/aggregate.py {params.results} > {output}"
