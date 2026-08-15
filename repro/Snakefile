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


def hq(cpus):
    if os.environ.get("SEAMS_NO_HQ") or shutil.which("hq") is None:
        return ""
    return "hq submit --wait --cpus {} --".format(cpus)


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
    shell:
        "meson setup build-repro --buildtype=release "
        "-Dwith_python=true -Dwith_tests=true > {params.log} 2>&1"


rule build_tip:
    input:
        R + "/tip-setup.done",
    output:
        touch(R + "/tip-build.done"),
    params:
        hq=hq(8),
    shell:
        "{params.hq} meson compile -C build-repro"


rule identity_gate:
    # A red suite aborts the DAG; no table row is emitted from a red tree
    input:
        R + "/tip-build.done",
    output:
        R + "/tip-test.log",
    params:
        hq=hq(8),
    shell:
        "{params.hq} bash -c 'meson test -C build-repro > {output} 2>&1'"


rule install_python:
    # The descriptor comparison imports the built bindings from the prefix
    input:
        R + "/tip-test.log",
    output:
        touch(R + "/py-install.done"),
    shell:
        "meson install -C build-repro > /dev/null 2>&1"


rule build_base:
    output:
        touch(R + "/base-build.done"),
    params:
        tree=BASE_TREE,
        hq=hq(8),
    shell:
        r"""
        test -d {params.tree} || {{ echo "baseline tree missing; run repro/elja_submit.sh prep" >&2; exit 1; }}
        (cd {params.tree} && meson setup build-base --buildtype=release \
          -Dwith_python=false -Dwith_tests=false > /dev/null 2>&1 || true)
        {params.hq} meson compile -C {params.tree}/build-base
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
    shell:
        r"""
        EIGEN=$CONDA_PREFIX/include/eigen3
        INC="-I{params.tree}/src/include/internal -I{params.tree}/src/include/external -I$EIGEN"
        LINK="-L{params.tree}/build-base/src -lyodaLib -Wl,-rpath,$(readlink -f {params.tree})/build-base/src -Wl,--enable-new-dtags -Wl,--allow-shlib-undefined"
        g++ -std=c++20 -O2 -o {output.scaling} tests/bench_scaling.cpp $INC $LINK
        g++ -std=c++20 -O2 -o {output.cages} tests/bench_cages_base.cpp $INC $LINK
        """


rule tip_scaling:
    input:
        R + "/tip-test.log",
    output:
        R + "/tip-scaling.txt",
    params:
        hq=hq(1),
        sizes=" ".join(str(s) for s in config["scaling_sizes"]),
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "../build-repro/tests/bench_scaling {params.sizes} > ../{output}'"


rule base_scaling:
    input:
        R + "/base_bench_scaling",
    output:
        R + "/base-scaling.txt",
    params:
        hq=hq(1),
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
        hq=hq(1),
        traj=config["trajectory"],
        reps=config["cage_reps"],
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "../build-repro/tests/bench_cages {params.traj} 1 {params.reps} "
        "> ../{output}'"


rule base_cages:
    input:
        R + "/base_bench_cages",
    output:
        R + "/base-cages.txt",
    params:
        hq=hq(1),
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
        hq=hq(1),
        sizes=" ".join(str(s) for s in config["overhead_sizes"]),
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "../build-repro/tests/bench_overhead {params.sizes} > ../{output}'"


rule tip_strong:
    input:
        R + "/tip-test.log",
    output:
        R + "/tip-strong-t{threads}.txt",
    params:
        hq=hq(4),
        atoms=config["strong_atoms"],
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS={wildcards.threads} "
        "../build-repro/tests/bench_strong {params.atoms} > ../{output}'"


rule trajectory_incremental:
    # Exactness referee and steady-state timings for the incremental rings
    # and the incremental cage affiliation; nonzero exit on any inequality
    input:
        R + "/tip-test.log",
    output:
        R + "/trajectory-incremental.txt",
    params:
        hq=hq(1),
        traj=config["trajectory"],
        frames=config["trajectory_frames"],
    shell:
        "{params.hq} bash -c 'cd input && OMP_NUM_THREADS=1 "
        "../build-repro/tests/bench_trajectory_incremental {params.traj} "
        "{params.frames} > ../{output}'"


rule ql_compare_dseams:
    input:
        R + "/tip-test.log",
    output:
        R + "/ql-dseams.txt",
    params:
        hq=hq(1),
    shell:
        "{params.hq} bash -c 'OMP_NUM_THREADS=1 "
        "build-repro/tests/compare_structure_desc > {output}'"


rule ql_compare_python:
    input:
        R + "/py-install.done",
    output:
        R + "/ql-python.txt",
    params:
        hq=hq(1),
    shell:
        "{params.hq} bash -c "
        "'LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH "
        "python scripts/compare_ql_literature.py > {output}'"


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
    output:
        R + "/paper_manifest.json",
    params:
        results=R,
    shell:
        "python repro/scripts/aggregate.py {params.results} > {output}"
