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

import json
import os
import shutil

configfile: "repro/config.yaml"

R = "repro/results"
LOCK = "repro/ecosystem-lock.json"
BASE_TREE = config["base_tree"]
# Build directories default to the working tree but move to node-local
# storage (SEAMS_BUILD_ROOT) on clusters whose NFS clocks skew against
# the nodes, which meson refuses
BUILD = os.path.abspath(os.environ.get("SEAMS_BUILD_ROOT", "."))
TIP_BUILD = os.path.join(BUILD, "build-repro")
BASE_BUILD = os.path.join(BUILD, "build-base")
SOURCE_ROOT = os.path.abspath(
    os.environ.get("DSEAMS_SOURCE_ROOT", os.path.join("..", "dseams-repro-sources"))
)
with open(LOCK, encoding="utf-8") as lock_file:
    SOURCE_LOCK = json.load(lock_file)
PYTHON_ROOT = os.path.join(
    SOURCE_ROOT, SOURCE_LOCK["components"]["pydseamslib"]["directory"]
)
YODA_ROOT = os.path.join(
    SOURCE_ROOT, SOURCE_LOCK["components"]["yodastruct"]["directory"]
)
PYTHON_SITE = os.path.join(BUILD, "python-site")
YODA_BUILD = os.path.join(BUILD, "build-yoda")
LUA_EXE = shutil.which("lua") or "lua"


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


rule source_manifest:
    input:
        lock=LOCK,
    output:
        R + "/source-manifest.json",
    params:
        root=SOURCE_ROOT,
        core=os.path.abspath("."),
    shell:
        "python repro/scripts/ecosystem_sources.py --lock {input.lock} manifest "
        "--root {params.root} --core {params.core} --output {output}"


rule setup_tip:
    input:
        R + "/source-manifest.json",
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
    input:
        gate=R + "/tip-test.log",
        source=R + "/source-manifest.json",
    output:
        touch(R + "/py-install.done"),
    params:
        source=PYTHON_ROOT,
        site=PYTHON_SITE,
    shell:
        "python -m pip install --no-build-isolation --no-deps --upgrade "
        "--target {params.site} {params.source} && "
        "PYTHONPATH={params.site}:${{PYTHONPATH:-}} "
        "python -c 'import pydseams; print(pydseams.__file__)'"


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
        site=PYTHON_SITE,
    shell:
        "{params.hq} bash -c "
        "'LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH "
        "PYTHONPATH={params.site}:${{PYTHONPATH:-}} "
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
    # The source manifest verifies that its wrap points at this engine.
    input:
        gate=R + "/tip-test.log",
        source=R + "/source-manifest.json",
    output:
        touch(R + "/yoda-build.done"),
    params:
        hq=lambda wc: hq(8),
        yoda=YODA_ROOT,
        bdir=YODA_BUILD,
    shell:
        r"""
        test -d {params.yoda}/lua || {{ echo "locked yodaStruct source missing: {params.yoda}" >&2; exit 1; }}
        meson setup {params.bdir} {params.yoda} --buildtype=release --wrap-mode=nodownload \
          > repro/results/yoda-setup.log 2>&1
        {params.hq} meson compile -C {params.bdir}
        """


rule yoda_identity_gate:
    input:
        R + "/yoda-build.done",
    output:
        R + "/yoda-test.log",
    params:
        hq=lambda wc: hq(8),
        bdir=YODA_BUILD,
    shell:
        "{params.hq} bash -c 'meson test -C {params.bdir} "
        "--print-errorlogs > {output} 2>&1'"


rule workflow_parity:
    input:
        gate=R + "/tip-test.log",
        python=R + "/py-install.done",
        yoda=R + "/yoda-test.log",
        source=R + "/source-manifest.json",
        water="input/traj/exampleTraj.lammpstrj",
        ice="input/traj/mW_cubic.lammpstrj",
        ions="repro/fixtures/tiny-ions.lammpstrj",
    output:
        R + "/workflow-parity.json",
    params:
        cli=os.path.join(TIP_BUILD, "seams"),
        lua=LUA_EXE,
        yoda=YODA_ROOT,
        yoda_build=YODA_BUILD,
        site=PYTHON_SITE,
    shell:
        "PYTHONPATH={params.site}:${{PYTHONPATH:-}} "
        "python repro/scripts/workflow_parity.py "
        "--seams {params.cli} --lua {params.lua} "
        "--lua-source {params.yoda} --lua-build {params.yoda_build} "
        "--water {input.water} --ice {input.ice} --ions {input.ions} "
        "--output {output}"


rule figshare_demos:
    # The five 1.0 figshare deposits through require("dseams"); nonzero
    # exit on any failed demo
    input:
        yoda=R + "/yoda-test.log",
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

NOTEBOOK_INPUTS = dict(
    py=R + "/py-install.done",
    runner="repro/scripts/execute_notebook.sh",
    traj=expand(FIGSHARE_DIR + "/{f}", f=FIGSHARE_FILES),
)


rule figshare_notebook:
    # Percent-format sources under repro/notebooks/; jupytext converts
    # and papermill executes. Execution is the test (each notebook
    # asserts its own headline numbers) and the executed .ipynb is the
    # artifact. 02_nucleation_cages is a separate rule because it also
    # writes figshare-incremental.json (Snakemake deletes declared
    # outputs before the job, so a later test -s stamp cannot see a
    # side-effect file).
    input:
        nb="repro/notebooks/{nb}.py",
        **NOTEBOOK_INPUTS,
    output:
        R + "/notebooks/{nb}.ipynb",
    wildcard_constraints:
        nb="(?!02_nucleation_cages).+",
    params:
        hq=lambda wc: hq(2),
        site=PYTHON_SITE,
    shell:
        "{params.hq} bash -c "
        "'LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH "
        "PYTHONPATH={params.site}:${{PYTHONPATH:-}} OMP_NUM_THREADS=2 "
        "repro/scripts/execute_notebook.sh {input.nb} {output}'"


rule figshare_nucleation_notebook:
    input:
        nb="repro/notebooks/02_nucleation_cages.py",
        **NOTEBOOK_INPUTS,
    output:
        ipynb=R + "/notebooks/02_nucleation_cages.ipynb",
        incremental=R + "/figshare-incremental.json",
    params:
        hq=lambda wc: hq(2),
        site=PYTHON_SITE,
    shell:
        "{params.hq} bash -c "
        "'LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH "
        "PYTHONPATH={params.site}:${{PYTHONPATH:-}} OMP_NUM_THREADS=2 "
        "repro/scripts/execute_notebook.sh {input.nb} {output.ipynb}'"


rule aggregate:
    input:
        conditions=R + "/conditions.txt",
        source=R + "/source-manifest.json",
        parity=R + "/workflow-parity.json",
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
