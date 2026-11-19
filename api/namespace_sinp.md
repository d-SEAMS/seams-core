<!-- https://docs.dseams.info/api/namespace_sinp.md -->

# namespace sinp {#namespace-sinp}

## Overview {#overview}

Functions for the d-SEAMS readers. [More…](#details-namespacesinp)

## Detailed Documentation {#detailed-documentation}

Functions for the d-SEAMS readers.

This namespace contains functions that are used for reading in the formats supported by d-SEAMS. LAMMPS trajectory files and XYZ files are currently supported, though it is recommended to use LAMMPS trajectory files, since the simulation box size etc. are not normally present in XYZ files, and many analyses may fail without the correct box dimensions.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Dec 26, 2019
- Rohit Goswami \[[rog32@hi.is](mailto:rog32@hi.is)\]; date modified: Mar 20, 2021

#### Global Functions {#global-functions}

Get file list inside the input folder.

Number of ITEM: TIMESTEP frames in a LAMMPS dump, or 0 if the file cannot be read. The first full count per path walks the file once and later reads seek, matching chemfiles Trajectory::nsteps and readcon’s frame-offset table. Sequential load_frame walks reuse a live read position the way LAMMPS ReaderNative does on rerun.

Drop a cached dump session (tests that rewrite a path in place).

Text of the exception currently being handled. Call from a catch clause. A std::exception yields what(); anything else is unnamed.

Call fn(frame, cloud) for each frame from max(first, 1) to last (1-based, inclusive) that falls to part of parts. last <= 0 means every ITEM: TIMESTEP. Frames are dealt in rounds of parts, in serpentine order (part p takes slot p of even rounds and slot parts - 1 - p of odd ones), so a cost that drifts or alternates along the run reaches every part. typeFilter <= 0 keeps every atom. nThreads <= 0 uses the OpenMP default; 1 is serial. Each worker opens its own handle and seeks the shared offset table. Incremental RingUpdater / AffiliationUpdater state cannot be shared across workers: with nThreads 1 the frames arrive in order on one thread and may carry it, otherwise use the batch classifiers.

Function for reading in a specified frame (frame number and not timestep value) One id/type/x/y/z frame of the default oxygen selection: type 2 when any atom has that type, otherwise type 1. Returns false for any other column layout, leaving the cloud untouched by a successful read.

One LAMMPS dump frame, one atom type. The type argument is any LAMMPS type (the O in the name is historical). If isSlice is true, each kept point gets inSlice set; atoms outside the box stay in the cloud. nop is the type-filtered count.

One LAMMPS dump frame, one atom type, dropping atoms outside the slice when isSlice is true. nop is the kept count. An axis with lo == hi is unconstrained.

Function for reading in atom coordinates from an XYZ file.

Reads bonds into a vector of vectors from a file with a specific format.

Number of steps chemfiles can read, or -1 when chemfiles is off or the file is not a trajectory it knows. seams_lang calls this with or without chemfiles.

Postprocessing entry. CON stays on readcon. Every other file goes through chemfiles when that library is linked. typeI follows the CLI: negative keeps every atom, positive keeps one species, zero guesses oxygen.

True when each component lies in \[lo, hi\], or that axis has lo == hi (unconstrained). (\[0,0,0\], \[50,0,0\]) is x in \[0, 50\], y and z open.
