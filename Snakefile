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
