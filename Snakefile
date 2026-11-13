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


rule related_watch_v11:
    input:
        "scripts/h2h_protocol.py",
        "scripts/h2h_v11_clavis.py",
        "scripts/release_build.py",
    output:
        "benchmarks/h2h-watch-v11.csv",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_v11_clavis.py "
        "--python \"$PY\" "
        "--out benchmarks/h2h-watch-v11.csv"


rule related_watch_v12:
    input:
        "scripts/h2h_protocol.py",
        "scripts/h2h_v12_rings.py",
        "scripts/release_build.py",
    output:
        "benchmarks/h2h-watch-v12.csv",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_v12_rings.py "
        "--python \"$PY\" "
        "--out benchmarks/h2h-watch-v12.csv"


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


rule ql_eight:
    input:
        "scripts/bench_ql_threads.py",
        "tests/bench_ql.cpp",
    output:
        "discovery/ql-eight.csv",
    shell:
        "PY=\"${{QL_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/bench_ql_threads.py --out discovery/ql-eight.csv"


rule related_watch:
    input:
        "scripts/h2h_protocol.py",
        "scripts/h2h_watch.py",
        "scripts/h2h_mdapy.py",
        "scripts/h2h_mdapy_q6.py",
        "scripts/h2h_mdapy_rings.py",
        "scripts/h2h_ovito.py",
        "scripts/h2h_fastatomstruct.py",
        "scripts/h2h_vesin.py",
        "scripts/release_build.py",
        "benchmarks/h2h_watch.py",
    output:
        "benchmarks/h2h-watch.csv",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_watch.py "
        "--python \"$PY\" "
        "--out benchmarks/h2h-watch.csv"


# One sitting. Ideal ice Ic at 110592 and 884736 molecules, one thread
# and eight, six timed calls after a warm-up. GenIce cages and the King
# ladder are in the same run. benchmarks/h2h_watch.py is the one-call
# Airspeed Velocity arm.
rule rc3_sitting:
    input:
        "scripts/h2h_protocol.py",
        "scripts/h2h_rc3_sitting.py",
        "scripts/h2h_literature.py",
        "scripts/release_build.py",
        "bbdir-release/bindings/libseams_lang.so",
    output:
        "benchmarks/rc3-sitting.csv",
        "benchmarks/rc3-sitting.txt",
        "benchmarks/rc3-vesin-gap.txt",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_rc3_sitting.py "
        "--out {output[0]} --log {output[1]}"


rule linkcell_scale:
    input:
        "scripts/rank_scale.py",
    output:
        "benchmarks/rank-scale.csv",
    shell:
        "python3 scripts/rank_scale.py"


rule literature_h2h:
    input:
        "scripts/h2h_protocol.py",
        "scripts/h2h_literature.py",
        "scripts/gen_ideal_ice.py",
        "benchmarks/h2h_literature.py",
    output:
        "benchmarks/h2h-literature.csv",
        "benchmarks/h2h-literature.txt",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_literature.py "
        "--out benchmarks/h2h-literature.csv "
        "--log benchmarks/h2h-literature.txt"


rule hop_stride:
    input:
        expand("benchmarks/hop-stride-{stride}.txt", stride=[10, 100, 1000]),
        "benchmarks/hop-stride-crossover.txt",
        expand("benchmarks/hop-stride-tip4p-{stride}.txt", stride=[5, 50, 500]),


rule hop_stride_generate:
    input:
        "input/hop/regenerate.py",
        "input/hop/in.equilibrate",
        "input/hop/in.stride",
        "input/mw215/mW.sw",
    output:
        "data/hop/stride-10.lammpstrj",
        "data/hop/stride-100.lammpstrj",
        "data/hop/stride-1000.lammpstrj",
        "data/hop/nucleation.lammpstrj",
        "data/hop/manifest.txt",
    shell:
        "python3 input/hop/regenerate.py --frame 104 "
        "--strides 10,100,1000 --transitions 200 "
        "--nequil 2000 --seed 15780000 --threads 4"


rule hop_stride_time:
    input:
        "build/tests/bench_hop_stride",
        "data/hop/manifest.txt",
        "data/hop/stride-{stride}.lammpstrj",
    output:
        "benchmarks/hop-stride-{stride}.txt",
    shell:
        "cat data/hop/manifest.txt > {output}; "
        "./build/tests/bench_hop_stride time "
        "data/hop/stride-{stride}.lammpstrj "
        "--graph both --reps 8 --warmup 1 --stride {stride} >> {output}"


rule hop_tip4p_generate:
    input:
        "input/hop/regenerate_tip4p.py",
        "input/hop/in.tip4p.equilibrate",
        "input/hop/in.tip4p.stride",
        "input/brine/regenerate.py",
        "input/brine/grow.py",
    output:
        "data/hop/tip4p-stride-5.lammpstrj",
        "data/hop/tip4p-stride-50.lammpstrj",
        "data/hop/tip4p-stride-500.lammpstrj",
        "data/hop/tip4p-manifest.txt",
    shell:
        "python3 input/hop/regenerate_tip4p.py "
        "--strides 5,50,500 --transitions 200 "
        "--nmelt 10000 --ncool 5000 --nequil 10000 "
        "--seed 2522005 --temp 252.00 --threads 8"


rule hop_tip4p_time:
    input:
        "build/tests/bench_hop_stride",
        "data/hop/tip4p-manifest.txt",
        "data/hop/tip4p-stride-{stride}.lammpstrj",
    output:
        "benchmarks/hop-stride-tip4p-{stride}.txt",
    shell:
        "cat data/hop/tip4p-manifest.txt > {output}; "
        "./build/tests/bench_hop_stride time "
        "data/hop/tip4p-stride-{stride}.lammpstrj "
        "--graph both --reps 8 --warmup 1 --stride {stride} "
        "--type 1 >> {output}"


rule hop_stride_crossover:
    input:
        "build/tests/bench_hop_stride",
        "data/hop/stride-1000.lammpstrj",
    output:
        "benchmarks/hop-stride-crossover.txt",
    shell:
        "./build/tests/bench_hop_stride crossover "
        "data/hop/stride-1000.lammpstrj "
        "--graph both --frame-a 1 --frame-b 2 "
        "--reps 8 --warmup 1 > {output}"


rule rings_one_thread:
    input:
        "bbdir-release/tests/bench_rings_h2h",
        "scripts/rings_one_thread.py",
        "input/traj/genice_sH.lammpstrj",
        "input/traj/mW_cubic.lammpstrj",
    output:
        "markdown/rings-one-thread.json",
    shell:
        "PY=\"${{RINGS_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/rings_one_thread.py "
        "--bench {input[0]} --python \"$PY\" --json {output}"


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


rule ddc_regression:
    input:
        "bbdir/tests/bench_ddc_ladder",
        "scripts/ddc_ladder.sh",
        "paper/tex/figures/ladder_frames.py",
        "paper/tex/figures/worked_example.py",
    output:
        "bbdir/ddc-regression/ic-216.log",
        "bbdir/ddc-regression/ic-512.log",
        "bbdir/ddc-regression/ic-1000.log",
        "bbdir/ddc-regression/ic-1728.log",
        "bbdir/ddc-regression/ic-2744.log",
        "bbdir/ddc-regression/ic-4096.log",
    shell:
        "bash scripts/ddc_ladder.sh bbdir/tests/bench_ddc_ladder bbdir/ddc-regression"

# Basal, primary prism and secondary prism at 220 to 265 K.
# Frames stay in $SI4_DATA (default si4-data/) and are not inputs here.
rule classify_carry:
    input:
        "scripts/bench_classify_carry.py",
        "scripts/bench_classify_carry_worker.py",
        "bbdir-release/tests/bench_ql",
        "bbdir-release/tests/bench_trajectory_incremental",
        "bbdir-release/src/libyodaLib.so",
        "input/traj/genice_sI.lammpstrj",
        "input/traj/genice_sII.lammpstrj",
    output:
        "benchmarks/classify-carry-sitting.txt",
    shell:
        "python3 scripts/bench_classify_carry.py "
        "--lib bbdir-release/src/libyodaLib.so "
        "--ql bbdir-release/tests/bench_ql "
        "--traj \"$CLASSIFY_TRAJ\" "
        "--out {output}"


# Paper cells that were unmeasured. One rule per producer. The stride
# dumps come from hop_stride_generate and stay out of git. Run
# `snakemake -j1 paper_cells`. The hop rule times carried cages against
# a full rebuild at 1 and 8 threads in one sitting and writes
# cage_speedup on each row of tab-hop-stride.csv.
rule paper_hop_stride:
    input:
        "scripts/paper_hop_stride.py",
        "scripts/paper_stats.py",
        "bbdir-release/tests/bench_trajectory_incremental",
        "data/hop/stride-10.lammpstrj",
        "data/hop/stride-100.lammpstrj",
        "data/hop/stride-1000.lammpstrj",
        "data/hop/nucleation.lammpstrj",
    output:
        "benchmarks/paper-cells/tab-hop-stride.csv",
        "benchmarks/paper-cells/raw-hop-stride.txt",
    shell:
        "python3 scripts/paper_hop_stride.py --out benchmarks/paper-cells"


rule paper_ptm:
    input:
        "scripts/paper_ptm.py",
        "scripts/paper_stats.py",
        "scripts/bench_neighbours_v2.py",
        "bbdir-release/tests/bench_v2",
    output:
        "benchmarks/paper-cells/tab-sym-eq.csv",
        "benchmarks/paper-cells/frag-cli-ptm.csv",
        "benchmarks/paper-cells/raw-ptm.txt",
    shell:
        "python3 scripts/paper_ptm.py --out benchmarks/paper-cells"


rule paper_ptm_templates:
    input:
        "scripts/paper_ptm_templates.py",
        "scripts/paper_stats.py",
        "bbdir-release/src/libyodaLib.so",
    output:
        "benchmarks/paper-cells/ptm-default.csv",
        "benchmarks/paper-cells/ptm-ice.csv",
        "benchmarks/paper-cells/ptm-agreement.csv",
        "benchmarks/paper-cells/raw-ptm-templates.txt",
    shell:
        "python3 scripts/paper_ptm_templates.py --out benchmarks/paper-cells"


rule paper_symmetry_vs_ptm:
    input:
        "scripts/paper_symmetry_sitting.py",
        "scripts/paper_ptm_templates.py",
        "scripts/paper_stats.py",
        "bbdir-release/tests/bench_v2",
        "bbdir-release/src/libyodaLib.so",
    output:
        "benchmarks/paper-cells/symmetry-vs-ptm.csv",
        "benchmarks/paper-cells/raw-symmetry-vs-ptm.txt",
    shell:
        "python3 scripts/paper_symmetry_sitting.py --timing --out benchmarks/paper-cells"


rule paper_q6:
    input:
        "scripts/thread_scale_time.py",
        "scripts/paper_stats.py",
        "bbdir-release/bindings/libseams_lang.so",
        "bbdir-release/src/libyodaLib.so",
    output:
        "benchmarks/paper-cells/tab-ql-eq.csv",
        "benchmarks/paper-cells/frag-cli-q.csv",
        "benchmarks/paper-cells/frag-lit-q6.csv",
        "benchmarks/paper-cells/raw-q6.txt",
    shell:
        "python3 scripts/thread_scale_time.py --paper-cells benchmarks/paper-cells "
        "--repeats 8 --threads 8"


rule paper_rings:
    input:
        "scripts/paper_rings.py",
        "scripts/paper_stats.py",
        "bbdir-release/bindings/libseams_lang.so",
    output:
        "benchmarks/paper-cells/tab-rings-214.csv",
        "benchmarks/paper-cells/frag-cli-rsm.csv",
        "benchmarks/paper-cells/raw-rings.txt",
    shell:
        "python3 scripts/paper_rings.py --out benchmarks/paper-cells"


rule paper_cages:
    input:
        "scripts/paper_cages.py",
        "scripts/paper_stats.py",
        "scripts/rerun_main_text.py",
        "bbdir-release/tests/bench_cage_h2h",
        "bbdir-release/seams",
        "input/traj/genice_sI.lammpstrj",
        "input/traj/genice_sII.lammpstrj",
    output:
        "benchmarks/paper-cells/frag-lit-cages.csv",
        "benchmarks/paper-cells/raw-cages.txt",
        "benchmarks/paper-cells/rerun-main/report.txt",
    shell:
        "python3 scripts/paper_cages.py --out benchmarks/paper-cells"


rule paper_cells:
    input:
        "scripts/paper_cells_join.py",
        # tab-hop-stride.csv carries cage_speedup from the 1-thread and
        # 8-thread sitting (paper_hop_stride).
        "benchmarks/paper-cells/tab-hop-stride.csv",
        "benchmarks/paper-cells/tab-sym-eq.csv",
        "benchmarks/paper-cells/tab-ql-eq.csv",
        "benchmarks/paper-cells/tab-rings-214.csv",
        "benchmarks/paper-cells/frag-cli-ptm.csv",
        "benchmarks/paper-cells/frag-cli-q.csv",
        "benchmarks/paper-cells/frag-cli-rsm.csv",
        "benchmarks/paper-cells/frag-lit-q6.csv",
        "benchmarks/paper-cells/frag-lit-cages.csv",
    output:
        "benchmarks/paper-cells/tab-cli.csv",
        "benchmarks/paper-cells/tab-lit.csv",
    shell:
        "python3 scripts/paper_cells_join.py benchmarks/paper-cells"


rule growth_curve:
    input:
        "input/interface/growth_curve.py",
        "input/interface/espinosa_faces.py",
        "input/interface/espinosa_protocol.py",
        "input/traj/ice1h.lammpstrj",
        "input/mw215/mW.sw",
        "scripts/growth_curve_rates.py",
        "scripts/espinosa_candidates.py",
        "scripts/q6bar_baseline.py",
        "scripts/q6bar_timing.py",
        "bbdir/tests/interface_order",
    output:
        "docs/orgmode/explanation/growth-curve.tsv",
        "docs/orgmode/explanation/growth-curve-summary.csv",
        "docs/orgmode/explanation/growth-curve-fit.csv",
        "docs/orgmode/explanation/growth-curve-compare.csv",
        "docs/orgmode/explanation/growth-curve-timing.csv",
    shell:
        "python3 input/interface/growth_curve.py && "
        "python3 scripts/growth_curve_rates.py --bin bbdir/tests/interface_order && "
        "python3 scripts/q6bar_timing.py --bin bbdir/tests/interface_order"
