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


# Head-to-head campaign: d-SEAMS against freud, OVITO and the other tools
# scripts/h2h_campaign.py covers, for neighbour lists and ring statistics on
# water and ice from 136 to 110592 molecules, and the timing lanes of the
# asv suite in benchmarks/. Every rule runs one subcommand of that script and
# writes a done marker under bbdir-h2h/<stamp>/. The samples go to
# data/h2h-campaign/<stamp>.jsonl. That record is append-only and is not a
# Snakemake output, so Snakemake never deletes it. A rule that is interrupted
# and run again resumes its run and skips the samples already recorded.
#
# The rules form one chain so that nothing else shares the machine, and each
# claims every core Snakemake has. Every single-core rule pins to the same
# h2h_core (default 0) with no warm-up call. The arms and the cells
# interleave round by round and the odd rounds run in reverse (n>=12).
# d-SEAMS is measured twice in every cell, as its own A/A band. asv runs one
# invocation per sample.
#
#   snakemake --cores 1 h2h_campaign --config h2h_stamp=2026-11-18-r2
#
# Config: h2h_stamp names the run and the record (default today),
# h2h_core the pinned CPU id for campaign and single-core lane runs (default 0),
# h2h_rounds the rounds of the head-to-heads and lanes (default 12, at least 12),
# h2h_asv_rounds the rounds of the asv cross-check of the head-to-heads
# (default 12). h2h_attempt is a suffix such as -2 for the run ids, the raw logs
# and the done markers of h2h_lanes, h2h_asv and h2h_inventory.
# H2H_OVITO_LD_LIBRARY_PATH, when set, is added to the library path of the
# OVITO arms only.
import datetime

H2H_STAMP = str(config.get("h2h_stamp", datetime.date.today().isoformat()))
H2H_CORE = str(config.get("h2h_core", "0"))
H2H_ROUNDS = int(config.get("h2h_rounds", 12))
H2H_ASV_ROUNDS = int(config.get("h2h_asv_rounds", 12))
H2H_ATTEMPT = str(config.get("h2h_attempt", ""))
H2H_DIR = f"bbdir-h2h/{H2H_STAMP}"
H2H_RECORD = f"data/h2h-campaign/{H2H_STAMP}.jsonl"
H2H_COMMON = (
    f"--stamp {H2H_STAMP} --core {H2H_CORE} --record {H2H_RECORD} --resume --verbose"
)


def _h2h_args(name, attempt=""):
    tag = f"{name}{attempt}"
    return (
        f"{H2H_COMMON} --run-id {H2H_STAMP}-{tag} "
        f"--raw {H2H_DIR}/raw/{tag} --work {H2H_DIR}/work/{tag}"
    )


rule h2h_release_build:
    output:
        f"{H2H_DIR}/release-build.done",
    threads: 1000
    shell:
        "python3 scripts/release_build.py && touch {output}"


rule h2h_neighbours:
    input:
        rules.h2h_release_build.output,
    output:
        f"{H2H_DIR}/neighbours.done",
    threads: 1000
    params:
        args=_h2h_args("neighbours"),
    shell:
        "python3 scripts/h2h_campaign.py run --head-to-head neighbours "
        f"--rounds {H2H_ROUNDS} {{params.args}} --marker {{output}}"


rule h2h_rings:
    input:
        rules.h2h_neighbours.output,
    output:
        f"{H2H_DIR}/rings.done",
    threads: 1000
    params:
        args=_h2h_args("rings"),
    shell:
        "python3 scripts/h2h_campaign.py run --head-to-head rings "
        f"--rounds {H2H_ROUNDS} {{params.args}} --marker {{output}}"


rule h2h_lanes:
    input:
        rules.h2h_rings.output,
    output:
        f"{H2H_DIR}/lanes{H2H_ATTEMPT}.done",
    threads: 1000
    params:
        args=_h2h_args("lanes", H2H_ATTEMPT),
    shell:
        "python3 scripts/h2h_campaign.py lanes --lanes kernel,rings,scaling "
        f"--rounds {H2H_ROUNDS} {{params.args}} --marker {{output}}"


rule h2h_asv:
    input:
        rules.h2h_lanes.output,
    output:
        f"{H2H_DIR}/asv{H2H_ATTEMPT}.done",
    threads: 1000
    params:
        args=_h2h_args("asv", H2H_ATTEMPT),
    shell:
        "python3 scripts/h2h_campaign.py lanes --lanes h2h "
        f"--rounds {H2H_ASV_ROUNDS} {{params.args}} --marker {{output}}"


rule h2h_inventory:
    input:
        rules.h2h_asv.output,
    output:
        f"{H2H_DIR}/inventory{H2H_ATTEMPT}.done",
    threads: 1000
    params:
        args=_h2h_args("inventory", H2H_ATTEMPT),
    shell:
        "python3 scripts/h2h_campaign.py lanes --lanes inventory --rounds 1 "
        "--min-rounds 1 --timeout 300 {params.args} --marker {output}"


rule h2h_summary:
    input:
        rules.h2h_inventory.output,
    output:
        f"data/h2h-campaign/{H2H_STAMP}.md",
    params:
        record=H2H_RECORD,
    shell:
        "python3 scripts/h2h_campaign.py summary --record {params.record} > {output}"


rule h2h_campaign:
    input:
        rules.h2h_summary.output,
