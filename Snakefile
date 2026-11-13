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


rule related_watch_v13:
    input:
        "scripts/h2h_protocol.py",
        "scripts/h2h_v13_rings.py",
        "scripts/release_build.py",
    output:
        "benchmarks/h2h-watch-v13.csv",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_v13_rings.py "
        "--python \"$PY\" "
        "--out benchmarks/h2h-watch-v13.csv"
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


rule rings_ih:
    output:
        "benchmarks/rings-ih.txt",
    shell:
        "OMP_DYNAMIC=FALSE "
        "bbdir-release/tests/bench_rings_ih --session --reps 8 --warmup 1 "
        "> {output}"


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


# Ice Ih basal face against liquid TIP4P/Ice, and the mW nucleation census.
# The threaded kernels are science-front (CHILL+, q6, shell label).
# Dumps larger than 1 MB stay out of git. dump_ps in the meta file is the
# spacing passed to the interface rule. SCIENCE_FRONT overrides the binary.

rule science_tip4p_basal:
    output:
        "data/science-frames/tip4p-basal.lammpstrj",
        "paper/tex/figures/data/science/interface-meta.csv",
        "paper/tex/figures/data/science/large-frames.csv",
    shell:
        "${{SCIENCE_PYTHON:-python3}} scripts/make_tip4p_basal.py"


rule science_nucleation_dump:
    output:
        "data/science-frames/nucleation.lammpstrj",
    shell:
        "python3 scripts/fetch_mw_nucleation.py"


rule science_interface:
    input:
        "data/science-frames/tip4p-basal.lammpstrj",
        "paper/tex/figures/data/science/interface-meta.csv",
    output:
        "paper/tex/figures/data/science/interface-layers.csv",
        "paper/tex/figures/data/science/interface-fronts.csv",
        "paper/tex/figures/data/science/interface-velocity.csv",
        "paper/tex/figures/data/science/interface-atoms.csv",
        "paper/tex/figures/data/science/interface-causes.csv",
    shell:
        "ps=$(awk '$1==\"dump_ps\"{{print $2}}' paper/tex/figures/data/science/interface-meta.csv); "
        "OMP_NUM_THREADS=${{SCIENCE_THREADS:-4}} "
        "${{SCIENCE_FRONT:-bbdir/science-front}} interface {input[0]} "
        "--ps \"$ps\" --bin 1 --out-dir paper/tex/figures/data/science "
        "--threads ${{SCIENCE_THREADS:-4}}"


rule science_nucleation:
    input:
        "data/science-frames/nucleation.lammpstrj",
    output:
        "paper/tex/figures/data/science/nucleation-clusters.csv",
        "paper/tex/figures/data/science/nucleation-summary.csv",
        "paper/tex/figures/data/science/nucleus-atoms.csv",
        "paper/tex/figures/data/science/nucleation-causes.csv",
    shell:
        "OMP_NUM_THREADS=${{SCIENCE_THREADS:-4}} "
        "${{SCIENCE_FRONT:-bbdir/science-front}} nucleation {input} "
        "--out-dir paper/tex/figures/data/science "
        "--threads ${{SCIENCE_THREADS:-4}}"


rule science_validate:
    input:
        "paper/tex/figures/data/science/interface-atoms.csv",
        "paper/tex/figures/data/science/nucleus-atoms.csv",
        "paper/tex/figures/data/science/interface-causes.csv",
        "paper/tex/figures/data/science/nucleation-causes.csv",
    output:
        "paper/tex/figures/data/science/validation.csv",
    shell:
        "OVITO_THREAD_COUNT=${{SCIENCE_THREADS:-4}} "
        "${{SCIENCE_PYTHON:-python3}} scripts/science_validate.py"


rule science_figures:
    input:
        "paper/tex/figures/data/science/interface-layers.csv",
        "paper/tex/figures/data/science/interface-fronts.csv",
        "paper/tex/figures/data/science/interface-velocity.csv",
        "paper/tex/figures/data/science/interface-atoms.csv",
        "paper/tex/figures/data/science/nucleation-clusters.csv",
        "paper/tex/figures/data/science/nucleation-summary.csv",
        "paper/tex/figures/data/science/nucleus-atoms.csv",
        "paper/tex/figures/data/science/validation.csv",
    output:
        "paper/tex/figures/science-pdf/interface-snapshot.pdf",
        "paper/tex/figures/science-pdf/interface-profiles.pdf",
        "paper/tex/figures/science-pdf/cluster-size.pdf",
        "paper/tex/figures/science-pdf/nucleus.pdf",
    shell:
        "${{SCIENCE_PYTHON:-python3}} scripts/science_figures.py "
        "--out paper/tex/figures/science-pdf"


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
# Hydrate signature search at 1, 2, 4 and 8 threads. Replicates larger
# than the GenIce unit cell are built in the benchmark and are not
# committed. One sitting, serial.
# One sitting of the neighbour-fed kernels at 1, 2, 4 and 8 threads.
# Ice Ic at 4096, 32768, 110592 and 884736 molecules is built in the
# benchmark and is not written out.
rule pool_audit:
    input:
        "scripts/pool_audit_sitting.py",
        "scripts/paper_stats.py",
        "bbdir-release/tests/bench_pool_audit",
    output:
        "benchmarks/pool-audit.csv",
        "benchmarks/pool-audit.txt",
    shell:
        "python3 scripts/pool_audit_sitting.py --out benchmarks"


rule trace_threads:
    input:
        "scripts/trace_threads_sitting.py",
        "scripts/paper_stats.py",
        "bbdir-release/tests/bench_trace_threads",
        "input/traj/genice_sI.lammpstrj",
        "input/traj/genice_sII.lammpstrj",
    output:
        "benchmarks/trace-threads.csv",
        "benchmarks/trace-threads.txt",
    shell:
        "python3 scripts/trace_threads_sitting.py --out benchmarks"


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
        "benchmarks/paper-cells/tab-sym-eq.csv",


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
        "benchmarks/paper-cells/symmetry-vs-ptm.csv",


# One sitting. Every quoted timing is a row of sitting.csv.
rule paper_sitting:
    input:
        "scripts/paper_sitting.py",
        "scripts/paper_stats.py",
        "scripts/h2h_protocol.py",
        "scripts/h2h_literature.py",
        "bbdir-release/bindings/libseams_lang.so",
        "bbdir-release/src/libyodaLib.so",
        "bbdir-release/seams",
        "bbdir-release/tests/bench_v2",
        "input/traj/genice_sI.lammpstrj",
        "input/traj/genice_sII.lammpstrj",
    output:
        "benchmarks/paper-cells/sitting.csv",
        "benchmarks/paper-cells/raw-sitting.txt",
    shell:
        "python3 scripts/paper_sitting.py --out benchmarks/paper-cells"


# Same ice operations as sitting.csv, on an 8x ice Ic cell, at 1, 2, 4 and 8 threads.
rule paper_sitting_large:
    input:
        "scripts/paper_sitting_large.py",
        "scripts/paper_sitting.py",
        "scripts/paper_stats.py",
        "scripts/h2h_protocol.py",
        "benchmarks/paper-cells/tip4pice-liquid-4096.txt",
        "bbdir-release/bindings/libseams_lang.so",
        "bbdir-release/src/libyodaLib.so",
        "bbdir-release/seams",
    output:
        "benchmarks/paper-cells/sitting-large.csv",
        "benchmarks/paper-cells/raw-sitting-large.txt",
    shell:
        "python3 scripts/paper_sitting_large.py --out benchmarks/paper-cells"


# One command per table. Each file is a selection of sitting.csv.
rule paper_tables:
    input:
        "scripts/paper_tables.py",
        "benchmarks/paper-cells/sitting.csv",
    output:
        "benchmarks/paper-cells/tab-rings-214.csv",
        "benchmarks/paper-cells/tab-ql-eq.csv",
        "benchmarks/paper-cells/tab-neighbours.csv",
        "benchmarks/paper-cells/tab-chill.csv",
        "benchmarks/paper-cells/tab-cages.csv",
        "benchmarks/paper-cells/frag-cli-q.csv",
        "benchmarks/paper-cells/frag-cli-rsm.csv",
        "benchmarks/paper-cells/frag-lit-q6.csv",
        "benchmarks/paper-cells/frag-lit-cages.csv",
        "benchmarks/paper-cells/symmetry-vs-ptm.csv",
        "benchmarks/paper-cells/tab-sym-eq.csv",
        "benchmarks/paper-cells/frag-cli-ptm.csv",
    shell:
        "python3 scripts/paper_tables.py benchmarks/paper-cells"


rule paper_rings:
    input:
        "benchmarks/paper-cells/tab-rings-214.csv",


rule paper_q6:
    input:
        "benchmarks/paper-cells/tab-ql-eq.csv",


rule paper_neighbours:
    input:
        "benchmarks/paper-cells/tab-neighbours.csv",


rule paper_chill:
    input:
        "benchmarks/paper-cells/tab-chill.csv",


rule cage_h2h_v8:
    input:
        "scripts/cage_h2h_v8.py",
        "scripts/paper_stats.py",
        "bbdir-release/seams",
        "bbdir-release/bindings/libseams_lang.so",
    output:
        "benchmarks/cage-h2h-v8.txt",
    shell:
        "python3 scripts/cage_h2h_v8.py --out benchmarks/cage-h2h-v8.txt"


rule paper_cages:
    input:
        "benchmarks/paper-cells/tab-cages.csv",
        "benchmarks/paper-cells/rerun-main/report.txt",


rule paper_rerun:
    input:
        "scripts/rerun_main_text.py",
        "bbdir-release/seams",
        "input/traj/mW_cubic.lammpstrj",
        "input/traj/genice_sI.lammpstrj",
        "input/traj/genice_sII.lammpstrj",
        "input/traj/genice_sH.lammpstrj",
        "input/traj/nucleation_f1.lammpstrj",
    output:
        "benchmarks/paper-cells/rerun-main/report.txt",
    shell:
        "python3 scripts/rerun_main_text.py "
        "--seams bbdir-release/seams --input input/traj "
        "--archive input/traj/genice_sI.lammpstrj --source . "
        "--out benchmarks/paper-cells/rerun-main"


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


# One process, eight samples, d-SEAMS, mdapy, pyscal3, vesin, and OVITO.
# CHILL+ is timed as paper (four closest) and ovito (every neighbour inside 3.5 A).
# Core count, no host name.
rule h2h_mdapy_pyscal3_vesin:
    output:
        "benchmarks/h2h-mdapy-pyscal3-vesin.txt",
    shell:
        "LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu "
        "python3 scripts/h2h_mdapy_pyscal3_vesin.py > {output}"


# One process, eight samples, d-SEAMS, freud 3.6, and OVITO 3.17.
# Ice cells and the mixed frames. Core count, no host name.
rule h2h_freud_ovito:
    output:
        "benchmarks/h2h-freud-ovito.txt",
    shell:
        "LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu "
        "python3 scripts/h2h_freud_ovito.py > {output}"


# One process, eight samples, d-SEAMS, freud, and OVITO. Core count, no host name.
rule h2h_session:
    output:
        "benchmarks/h2h-session.txt",
    shell:
        "OVITO_THREAD_COUNT=8 python3 scripts/h2h_session.py > {output}"


# One process, eight samples. The record names the core count.
rule omp_scale_time:
    output:
        "benchmarks/omp-scale.txt",
    shell:
        "python3 scripts/omp_scale_time.py > {output}"


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


rule binding_overhead:
    input:
        "scripts/bench_binding_overhead.py",
        "benchmarks/binding_overhead_c.c",
        "benchmarks/binding_overhead_cpp.cpp",
        "benchmarks/binding_overhead.f90",
        "benchmarks/binding_overhead.jl",
        "benchmarks/binding_overhead.R",
        "benchmarks/binding_overhead.lua",
        "benchmarks/binding_overhead.py",
        "bindings/rust/seams/examples/binding_overhead.rs",
        "include/seams.hpp",
        "build/src/libseams.so",
    output:
        "benchmarks/binding-overhead.csv",
    shell:
        "python3 scripts/bench_binding_overhead.py --out benchmarks/binding-overhead.csv --reps 8"
