# Figure products the paper reads and git does not store.
# python3 paper/tex/figures/stats.py writes the plot tables and, when
# ring-head2head.json is present, the ring table and its macros.

PLOTS = [
    "plot-chill-dseams-chill.dat",
    "plot-chill-ovito-chill.dat",
    "plot-chill-palmer-chill.dat",
    "plot-churn.dat",
    "plot-nucleation.dat",
    "plot-ovrings-confined-dseams-rings.dat",
    "plot-ovrings-confined-ovito-rings.dat",
    "plot-ovrings-ic-dseams-rings.dat",
    "plot-ovrings-ic-ovito-rings.dat",
    "plot-ovrings-ih-dseams-rings.dat",
    "plot-ovrings-ih-ovito-rings.dat",
    "plot-ovrings-liquid-dseams-rings.dat",
    "plot-ovrings-liquid-ovito-rings.dat",
    "plot-q6-dseams-q6.dat",
    "plot-q6-freud-q6.dat",
    "plot-q6-pyscal-q6.dat",
    "plot-rings-dseams-rings.dat",
    "plot-rings-matscipy-rings-kernel.dat",
    "plot-seeded-dseams-seeded.dat",
    "plot-seeded-ovito-cna.dat",
    "plot-seeded-ovito-diamond.dat",
    "plot-seeded-ovito-ptm.dat",
    "plot-seeded-pyscal-diamond.dat",
    "plot-versions-one-ddc.dat",
    "plot-versions-one-hc.dat",
    "plot-versions-one-path.dat",
    "plot-versions-one-ring.dat",
    "plot-versions-two-ddc.dat",
    "plot-versions-two-hc.dat",
    "plot-versions-two-path.dat",
    "plot-versions-two-ring.dat",
]


rule all:
    input:
        expand("paper/tex/figures/data/{name}", name=PLOTS),
        "paper/tex/figures/data/plot-serial-rings.dat",
        "paper/tex/figures/data/ringperf_macros.tex",
        "paper/tex/figures/data/tab-ring-perf.tex",
        "paper/tex/figures/data/tab-neigh-perf.tex",
        "paper/tex/figures/data/h2h/input-manifest.json",
        "paper/tex/figures/data/nucleation-host.txt",
        "paper/tex/figures/data/.records-ready",
        "paper/tex/repro/inputs.sha256",
        "paper/tex/figures/xyzrender_xyz/si_hydrate.xyz",
        "paper/tex/figures/xyzrender_xyz/sii_hydrate.xyz",
        "paper/tex/figures/data/neigh-head2head.json",
        "bbdir/neigh-h2h.log",
        "paper/tex/figures/fig-cages.png",
        "paper/tex/figures/systems.png",
        "paper/tex/figures/interface.png",
        "paper/tex/figures/si_hydrate_xyzrender.png",
        "paper/tex/figures/sii_hydrate_xyzrender.png",


rule machine_records:
    output:
        "paper/tex/figures/data/host.txt",
        "paper/tex/figures/data/churn-host.txt",
        "paper/tex/figures/data/h2h/environment.json",
    shell:
        "python3 scripts/write_machine_records.py"


rule nucleation_host:
    output:
        "paper/tex/figures/data/nucleation-host.txt",
    shell:
        "python3 paper/tex/figures/generate_figure_data.py host"


rule input_manifest:
    output:
        "paper/tex/figures/data/h2h/input-manifest.json",
    shell:
        "python3 paper/tex/figures/generate_figure_data.py manifest"


rule hydrate_png:
    input:
        "paper/tex/figures/xyzrender_xyz/si_hydrate.xyz",
        "paper/tex/figures/xyzrender_xyz/sii_hydrate.xyz",
        "paper/tex/figures/render_hydrate.py",
    output:
        "paper/tex/figures/si_hydrate_xyzrender.png",
        "paper/tex/figures/sii_hydrate_xyzrender.png",
    shell:
        ".pixi/envs/default/bin/python paper/tex/figures/render_hydrate.py"


rule hydrate_xyz:
    input:
        "input/traj/genice_sI.lammpstrj",
        "input/traj/genice_sII.lammpstrj",
    output:
        "paper/tex/figures/xyzrender_xyz/si_hydrate.xyz",
        "paper/tex/figures/xyzrender_xyz/sii_hydrate.xyz",
    shell:
        "python3 paper/tex/figures/generate_figure_data.py hydrate"


rule plots:
    input:
        "paper/tex/figures/stats.py",
        "paper/tex/figures/data/.records-ready",
        "paper/tex/figures/data/nucleation-host.txt",
        "paper/tex/figures/data/host.txt",
        "paper/tex/figures/data/churn-host.txt",
        "paper/tex/figures/data/h2h/environment.json",
    output:
        expand("paper/tex/figures/data/{name}", name=PLOTS),
    shell:
        "python3 paper/tex/figures/stats.py"


rule ring_panel:
    input:
        "paper/tex/figures/data/ring-head2head.json",
        "paper/tex/figures/stats.py",
        "paper/tex/figures/data/.records-ready",
    output:
        "paper/tex/figures/data/plot-serial-rings.dat",
        "paper/tex/figures/data/ringperf_macros.tex",
        "paper/tex/figures/data/tab-ring-perf.tex",
        "paper/tex/figures/data/macros.tex",
    shell:
        "python3 paper/tex/figures/stats.py"


rule recorded_inputs:
    output:
        "paper/tex/figures/data/.records-ready",
        "paper/tex/repro/inputs.sha256",
    shell:
        "python3 paper/tex/figures/generate_figure_data.py records"


rule neigh_table:
    input:
        "paper/tex/figures/data/neigh-head2head.json",
        "paper/tex/figures/data/.records-ready",
        "paper/tex/figures/stats.py",
        "paper/tex/figures/data/ring-head2head.json",
    output:
        "paper/tex/figures/data/tab-neigh-perf.tex",
    shell:
        "python3 paper/tex/figures/stats.py"


rule neigh_record:
    input:
        "tests/bench_h2h_neighbours.py",
        "bbdir-eessi/tests/bench_neigh_only",
        "scripts/neigh_h2h.sh",
    output:
        "paper/tex/figures/data/neigh-head2head.json",
        "bbdir/neigh-h2h.log",
    shell:
        "bash scripts/neigh_h2h.sh"


rule scale_logs:
    input:
        "bbdir-eessi/tests/bench_scale",
        "scripts/scale_campaign.py",
    output:
        "bbdir/scale/campaign.done",
    shell:
        "python3 scripts/scale_campaign.py all && touch bbdir/scale/campaign.done"


rule scaling:
    input:
        "bbdir/scale/campaign.done",
        "scripts/scale_reduce.py",
    output:
        "paper/tex/figures/data/scale_macros.tex",
        "paper/tex/figures/data/scale-strong.dat",
        "paper/tex/figures/data/scale-weak.dat",
        "paper/tex/figures/data/tab-scale.tex",
        "paper/tex/figures/data/tab-scale-full.tex",
        "paper/tex/figures/data/tab-phase.tex",
        "paper/tex/figures/data/scale-number-log.tsv",
        "paper/tex/figures/data/scale-commit.txt",
        "paper/tex/figures/scale-strong.png",
        "paper/tex/figures/scale-weak.png",
    shell:
        "python3 scripts/scale_reduce.py"


rule figure_pngs:
    input:
        "paper/tex/figures/data/macros.tex",
        "paper/tex/figures/data/ring-head2head.json",
        "paper/tex/figures/rasterize_tikz.py",
        "paper/tex/figures/render_finish.py",
    output:
        "paper/tex/figures/fig-cages.png",
        "paper/tex/figures/fig-update.png",
        "paper/tex/figures/fig-versions.png",
        "paper/tex/figures/fig-codes.png",
        "paper/tex/figures/fig-churn.png",
        "paper/tex/figures/fig-nucleation.png",
        "paper/tex/figures/systems.png",
        "paper/tex/figures/classify.png",
        "paper/tex/figures/ringperf.png",
        "paper/tex/figures/interface.png",
    shell:
        "python3 paper/tex/figures/render_finish.py && python3 paper/tex/figures/rasterize_tikz.py"


rule related_watch_v5:
    input:
        "scripts/h2h_protocol.py",
        "scripts/h2h_v5.py",
        "scripts/bench_symmetry_ice.cpp",
    output:
        "benchmarks/h2h-watch-v5.csv",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_v5.py "
        "--python \"$PY\" "
        "--out benchmarks/h2h-watch-v5.csv"


rule related_watch_v4:
    input:
        "scripts/h2h_protocol.py",
        "scripts/h2h_v4.py",
        "scripts/release_build.py",
    output:
        "benchmarks/h2h-watch-v4.csv",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_v4.py "
        "--python \"$PY\" "
        "--out benchmarks/h2h-watch-v4.csv"


rule related_watch:
    input:
        "scripts/h2h_protocol.py",
        "scripts/h2h_watch.py",
        "scripts/h2h_mdapy.py",
        "scripts/h2h_ovito.py",
        "scripts/h2h_fastatomstruct.py",
        "scripts/release_build.py",
        "benchmarks/h2h_watch.py",
    output:
        "benchmarks/h2h-watch.csv",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_watch.py "
        "--python \"$PY\" "
        "--out benchmarks/h2h-watch.csv"


rule shipped_campaign:
    input:
        "paper/tex/repro/shipped_campaign.py",
        "paper/tex/repro/shipped.sh",
        "benchmarks/shipped.py",
        "tests/bench_h2h_rings.py",
    output:
        "paper/tex/figures/data/shipped-campaign.json",
        "paper/tex/figures/data/shipped_macros.tex",
        "paper/tex/figures/data/tab-shipped.tex",
    shell:
        "bash paper/tex/repro/shipped.sh"
