<!-- https://docs.dseams.info/user/howto/index.md -->

# Tasks {#tasks}

A how-to solves a goal you already have. It is not a learning path.

| guide | goal |
| --- | --- |
| Upgrade from 2.13 | API, engines, threads, and the remaining --types alias |
| Migrate from 1.x | Lua and yodaStruct -c mapped onto seams and the C API |
| Install the engine | Get a seams binary |
| Language bindings | Python, Lua, R, Julia, Rust, and Fortran, each in its own repository |
| Benchmarks | Which driver to time, and what to record |
| Literature head-to-head | mdapy, pyscal3, vesin, IRA and SQQ on ideal ice |
| Reproduce the head-to-head | Thread control for freud, OVITO, mdapy, vesin, PLUMED, PTM, and rsm |
| Cutoff neighbour lists | Triclinic cutoff list, thread count, and the comparison with vesin, freud, matscipy, mdapy, and OVITO |
| OpenMP target offload | Build the Steinhardt kernel and the TUM ice score for an NVIDIA GPU |
| Classify a dump | Run read / chill-plus / cages on a file you have |
| Primitive rings | Ring histogram of one frame |
| Cages and TRACE | HC and DDC counts, a ring-size signature, and the TRACE comparison |
| Local q4 and q6 | seams steinhardt and seams_q4q6 |
| CHILL+ | Cutoff list, mutual list, OVITO rule |
| Choosing a cutoff and padding | Open-face padding against the bulk-liquid check and the interface layers |
| Symmetry classifier | seams symmetry and seams_frame_symmetry |
| Linkcell and carried rings | The cutoff CSR and the hop-ball update |
| Sheared cells | Triclinic dump: bound spans are not edge lengths; recover H |
| Ionic-liquid dumps | Keep one species with --type; unlike-pair seams rdf --pair-types I,J; confined \rho(z) |
| ASE and solvis | to_ase and to_solvis. The ice Ih calls are on Python Frame |
| Confined ice | Nanotube and monolayer classification |
| Trajectories, dipoles, and time-correlation functions | Multi-frame TrajectorySequence, triclinic unwrapping, water dipoles, axial/lateral MSD, C_{\text{HB}}/S_{\text{HB}}, and C_1/C_2 |
| Basal profiles and nucleation | Layer profile on one frame, and the critical-size census |
| From a trajectory to the ice front | CHILL+, rings, and the front on a dump in input/traj/ |
| GROMACS | Build gmx-seams and read the CHILL+ census |
| LAMMPS | Build the fix hook and read the CHILL+ census |
| OpenMM | Build the reporter and read the CHILL+ census |
| PLUMED plugin v0.4.1 | PLUMED_NUM_THREADS, then q6, ice, and rings |
| In-situ CHILL+ | Labels and the ring histogram from libseams_lang |
| FAQ | Short answers |
| Troubleshooting | Symptom, cause, fix |
