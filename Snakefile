# Engine records this repository regenerates. The d-SEAMS 2.0 paper
# campaign, its figure and table rules, lives in the paper repository.

# Engine end-to-end records. Each engine may need its own toolchain, so
# E2E_BUILD_<ENGINE> names the d-SEAMS build directory (default build) and
# E2E_SHELL_<ENGINE> a command prefix that sets up the engine, for example
# one that sources a module environment. <ENGINE> is upper case with - as _.
# An engine that is not found is a skip row in the table.
import os

ENGINES_E2E = ["lammps", "lammps-np2", "plumed", "gromacs", "gromacs-plumed", "openmm"]


def _e2e_env(engine, kind, default):
    return os.environ.get(f"E2E_{kind}_{engine.upper().replace('-', '_')}", default)


rule engine_e2e:
    input:
        "scripts/engines_e2e.py",
        "engines/e2e_common.py",
    output:
        "data/engines-e2e/{engine}.json",
    wildcard_constraints:
        engine="|".join(ENGINES_E2E),
    params:
        build=lambda w: _e2e_env(w.engine, "BUILD", "build"),
        shell=lambda w: _e2e_env(w.engine, "SHELL", ""),
    shell:
        "{params.shell} python3 scripts/engines_e2e.py run {wildcards.engine} "
        "{params.build} {output}"


rule engines_e2e:
    input:
        expand("data/engines-e2e/{engine}.json", engine=ENGINES_E2E),
    output:
        "data/engines-e2e.csv",
    shell:
        "python3 scripts/engines_e2e.py table {input} --out {output}"


# Ideal ice Ic, cutoff 3.5 A, King's shortest-path rings and q_6.
# One core, twenty rounds, no pre-warm. FAS_SEAMS_PYTHON, FAS220_PYTHON
# and FAS230_PYTHON select the interpreters.
rule fas230_samples:
    input:
        "scripts/h2h_fas230.py",
        "scripts/h2h_fas230_summary.py",
        "scripts/h2h_fastatomstruct.py",
        "scripts/h2h_watch.py",
        "scripts/h2h_protocol.py",
        "scripts/release_build.py",
    output:
        "benchmarks/fas230-samples.csv",
        "benchmarks/fas230-lscpu.txt",
    shell:
        "\"${{FAS_SEAMS_PYTHON:-python3}}\" scripts/h2h_fas230.py "
        "--samples {output[0]} --lscpu {output[1]}"


rule fas230_summary:
    input:
        samples="benchmarks/fas230-samples.csv",
        script="scripts/h2h_fas230_summary.py",
    output:
        "benchmarks/fas230-summary.csv",
    shell:
        "python3 {input.script} --samples {input.samples} --summary {output}"
