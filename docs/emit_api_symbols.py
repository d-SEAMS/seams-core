#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Write the per-symbol blocks of the Sphinx API reference."""

import re
from pathlib import Path

import check_c_symbols

ROOT = Path(__file__).resolve().parents[1]
API = ROOT / "docs" / "api"
BEGIN = "<!-- API_SYMBOLS_BEGIN -->"
END = "<!-- API_SYMBOLS_END -->"
PAGES = {
    "c": "c.md",
    "cxx": "cpp.md",
    "python": "python.md",
    "rust": "rust.md",
    "fortran": "fortran.md",
    "julia": "julia.md",
    "r": "r.md",
}
EXAMPLES = {
    "c": "docs/examples/minimal_c.c",
    "cxx": "docs/examples/minimal_cpp.cpp",
    "python": "docs/examples/minimal_py.py",
    "rust": "docs/examples/minimal_rust.rs",
    "fortran": "docs/examples/minimal_fortran.f90",
    "julia": "docs/examples/minimal_julia.jl",
    "r": "docs/examples/minimal_r.R",
}

# name, purpose, parameters, units
RAW = {
    "c": [
        ("SEAMS_ABI", "ABI number compiled into the header.", "none", "dimensionless integer 1"),
        ("SEAMS_CHILL_OVITO", "CHILL+ rule used by OVITO.", "none", "dimensionless code 1"),
        ("SEAMS_CHILL_PUBLISHED", "CHILL+ rule from Nguyen and Molinero, 2015.", "none", "dimensionless code 0"),
        ("SEAMS_ERR_ARG", "A caller argument is the wrong shape.", "none", "dimensionless status code 4"),
        ("SEAMS_ERR_BOX", "The cell is missing or has no volume.", "none", "dimensionless status code 3"),
        ("SEAMS_ERR_INTERNAL", "The library hit an internal failure.", "none", "dimensionless status code 7"),
        ("SEAMS_ERR_NOMEM", "An allocation failed.", "none", "dimensionless status code 6"),
        ("SEAMS_ERR_NULL", "A required pointer is null.", "none", "dimensionless status code 1"),
        ("SEAMS_ERR_RANGE", "A number is outside the accepted range.", "none", "dimensionless status code 2"),
        ("SEAMS_ERR_STATE", "The handle has no frame yet.", "none", "dimensionless status code 5"),
        ("SEAMS_GRAPH_CUTOFF", "Neighbour graph of every pair inside the cutoff.", "none", "dimensionless code 0"),
        ("SEAMS_GRAPH_MUTUAL", "Mutual four-nearest neighbour graph.", "none", "dimensionless code 1"),
        ("SEAMS_LABEL_CLATHRATE", "CHILL+ clathrate label.", "none", "dimensionless code 4"),
        ("SEAMS_LABEL_CUBIC", "CHILL+ cubic-ice label.", "none", "dimensionless code 0"),
        ("SEAMS_LABEL_HEXAGONAL", "CHILL+ hexagonal-ice label.", "none", "dimensionless code 1"),
        ("SEAMS_LABEL_INTERCLATHRATE", "CHILL+ interfacial-clathrate label.", "none", "dimensionless code 5"),
        ("SEAMS_LABEL_INTERFACIAL", "CHILL+ interfacial-ice label.", "none", "dimensionless code 3"),
        ("SEAMS_LABEL_RECUBIC", "CHILL+ reclassified cubic label.", "none", "dimensionless code 7"),
        ("SEAMS_LABEL_REHEX", "CHILL+ reclassified hexagonal label.", "none", "dimensionless code 8"),
        ("SEAMS_LABEL_UNCLASSIFIED", "CHILL+ label for an atom that matched no class.", "none", "dimensionless code 6"),
        ("SEAMS_LABEL_WATER", "CHILL+ liquid-water label.", "none", "dimensionless code 2"),
        ("SEAMS_OK", "The call succeeded.", "none", "dimensionless status code 0"),
        ("SEAMS_RING_BOTH_BASAL", "Six-ring that is a basal face of both cage types.", "none", "dimensionless code 4"),
        ("SEAMS_RING_BOTH_PRISMATIC", "Six-ring that is a prismatic face of both cage types.", "none", "dimensionless code 5"),
        ("SEAMS_RING_DDC", "Six-ring of a double-diamond cage.", "none", "dimensionless code 1"),
        ("SEAMS_RING_DEFORMED_PRISM", "Deformed prismatic six-ring.", "none", "dimensionless code 7"),
        ("SEAMS_RING_HC_BASAL", "Basal six-ring of a hexagonal cage.", "none", "dimensionless code 2"),
        ("SEAMS_RING_HC_PRISMATIC", "Prismatic six-ring of a hexagonal cage.", "none", "dimensionless code 3"),
        ("SEAMS_RING_HEXAGON", "Six-ring that is not a cage face.", "none", "dimensionless code 10"),
        ("SEAMS_RING_MIXED_PRISM", "Prismatic six-ring mixing the two cage types.", "none", "dimensionless code 8"),
        ("SEAMS_RING_OTHER", "Ring that is neither a pentagon nor a hexagon cage face.", "none", "dimensionless code 11"),
        ("SEAMS_RING_PENTAGON", "Five-membered ring.", "none", "dimensionless code 9"),
        ("SEAMS_RING_PRISM", "Prismatic six-ring.", "none", "dimensionless code 6"),
        ("SEAMS_RING_UNCLASSIFIED", "Ring with no symmetry class yet.", "none", "dimensionless code 0"),
        ("seams_abi", "Return the ABI number of this libseams.", "none", "dimensionless integer"),
        ("seams_chill_plus", "Write one CHILL+ label per atom.", "sys; cutoff; type; graph; rule; labels; counts9", "cutoff and coordinates in angstrom; labels and counts are dimensionless"),
        ("seams_chill_rule", "Enumerator of the published and OVITO CHILL+ rules.", "none", "dimensionless"),
        ("seams_get_num_threads", "Report the OpenMP width the next call will use.", "none", "dimensionless thread count"),
        ("seams_graph", "Enumerator of the cutoff graph and the mutual graph.", "none", "dimensionless"),
        ("seams_label", "Enumerator of the nine CHILL+ labels.", "none", "dimensionless"),
        ("seams_label_name", "Map a CHILL+ code to its name.", "label", "dimensionless code in, text out"),
        ("seams_last_error", "Return the thread-local error detail.", "none", "text"),
        ("seams_neighbours", "Opaque handle of one neighbour list.", "none", "handle"),
        ("seams_neighbours_build", "Build a cutoff or mutual neighbour list.", "sys; cutoff; type; graph; out", "cutoff in angstrom; type 0 keeps every atom"),
        ("seams_neighbours_count", "Return the number of neighbour rows.", "nb", "dimensionless count"),
        ("seams_neighbours_csr", "Return the neighbour list as CSR offsets and columns.", "nb; offsets; n_offsets; cols; nnz", "dimensionless indices; the self entry is omitted"),
        ("seams_neighbours_destroy", "Free a neighbour-list handle.", "nb", "handle"),
        ("seams_q4q6", "Write local q4 and q6 on the cutoff shell.", "sys; cutoff; type; q4; q6", "cutoff in angstrom; q4 and q6 are dimensionless"),
        ("seams_ring_symmetry", "Enumerator of primitive-ring symmetry classes.", "none", "dimensionless"),
        ("seams_ring_symmetry_name", "Map a ring-symmetry code to its name.", "symmetry", "dimensionless code in, text out"),
        ("seams_rings", "Opaque handle of one primitive-ring set.", "none", "handle"),
        ("seams_rings_count", "Return the number of primitive rings.", "rings", "dimensionless count"),
        ("seams_rings_destroy", "Free a ring-set handle.", "rings", "handle"),
        ("seams_rings_find", "Find primitive rings up to a maximum size.", "sys; cutoff; type; graph; max_depth; out", "cutoff in angstrom; max_depth is a dimensionless ring size from 3"),
        ("seams_rings_fraction", "Fraction of atoms that sit on a ring of one size.", "rings; size; fraction", "size 0 counts every ring; fraction is dimensionless"),
        ("seams_rings_histogram", "Count primitive rings by size.", "rings; counts; counts_len", "counts[k] is the dimensionless number of rings of size k"),
        ("seams_rings_labels", "Write the size and symmetry of each ring.", "rings; sizes; symmetry", "dimensionless size and symmetry code"),
        ("seams_set_num_threads", "Set the OpenMP width for the next call.", "n", "dimensionless thread count; 0 restores the library default"),
        ("seams_status", "Enumerator of the status codes.", "none", "dimensionless"),
        ("seams_status_string", "Map a status code to its name.", "status", "dimensionless code in, text out"),
        ("seams_system", "Opaque handle of one frame.", "none", "handle"),
        ("seams_system_count", "Return the number of points in the frame.", "sys", "dimensionless count"),
        ("seams_system_create", "Allocate an empty frame handle.", "none", "handle"),
        ("seams_system_destroy", "Free a frame handle.", "sys", "handle"),
        ("seams_system_set_frame", "Store points and a general LAMMPS cell.", "sys; xyz; n; box9; types; ids", "coordinates and box lengths in angstrom; n is a dimensionless count"),
        ("seams_system_set_ortho", "Store points and an orthorhombic cell.", "sys; xyz; n; ox; oy; oz; lx; ly; lz; types; ids", "origin and edge lengths in angstrom"),
        ("SEAMS_CAGE_BOTH", "Atom that belongs to both cage types.", "none", "dimensionless code 3"),
        ("SEAMS_CAGE_CUBIC", "Atom that belongs only to a double-diamond cage.", "none", "dimensionless code 2"),
        ("SEAMS_CAGE_HEX", "Atom that belongs only to a hexagonal cage.", "none", "dimensionless code 1"),
        ("SEAMS_CAGE_WATER", "Atom that belongs to neither cage.", "none", "dimensionless code 0"),
        ("SEAMS_DEPRECATED", "Mark a shipped name and name its replacement.", "message", "attribute text"),
        ("SEAMS_GRAPH_SEEDED", "Cage graph: mutual seeds and union completions.", "none", "dimensionless code 2"),
        ("SEAMS_VERSION_MAJOR", "Major component of this release.", "none", "dimensionless integer"),
        ("SEAMS_VERSION_MINOR", "Minor component of this release.", "none", "dimensionless integer"),
        ("SEAMS_VERSION_NUMBER", "Packed version, (major << 16) | (minor << 8) | patch.", "none", "dimensionless integer"),
        ("SEAMS_VERSION_PATCH", "Patch component of this release.", "none", "dimensionless integer"),
        ("SEAMS_VERSION_STRING", "Release text written by the build.", "none", "text"),
        ("seams_brine", "Opaque handle of one ice-front measurement.", "none", "handle"),
        ("seams_brine_analyse", "Measure a planar ice front from an ice assignment.", "sys; ice; n_ice; water_type; ion_type; n_bins; axis; out", "lengths in angstrom; ice flags are dimensionless"),
        ("seams_brine_destroy", "Free a brine-front handle.", "front", "handle"),
        ("seams_brine_fronts", "Report the two planar crossings of the ice slab.", "front; z_enter; z_exit; thickness; found", "lengths in angstrom"),
        ("seams_brine_pockets", "Count brine pockets and the liquid reservoir.", "front; n_pockets; reservoir_waters", "dimensionless counts"),
        ("seams_brine_profile", "Copy the concentration profile along the front axis.", "front; z_mid; molality; ice_fraction; cap; n_out", "z in angstrom; molality in mol/kg; ice fraction is dimensionless"),
        ("seams_brine_rejection", "Report the molal rejection coefficient.", "front; c_ice; c_liquid; k; defined", "molality in mol/kg; k is dimensionless"),
        ("seams_cage_label", "Enumerator of per-atom cage membership.", "none", "dimensionless"),
        ("SEAMS_SITE_ACCEPTOR", "Site kind of an acceptor atom.", "none", "dimensionless code 5"),
        ("SEAMS_SITE_ANION", "Site kind of an anion.", "none", "dimensionless code 2"),
        ("SEAMS_SITE_APOLAR", "Site kind of an apolar atom.", "none", "dimensionless code 7"),
        ("SEAMS_SITE_CATION", "Site kind of a cation head.", "none", "dimensionless code 1"),
        ("SEAMS_SITE_DONOR_H", "Site kind of a donor hydrogen.", "none", "dimensionless code 4"),
        ("SEAMS_SITE_POLAR", "Site kind of a polar atom.", "none", "dimensionless code 6"),
        ("SEAMS_SITE_SOLVENT", "Site kind of a solvent atom.", "none", "dimensionless code 10"),
        ("SEAMS_SITE_TAIL", "Site kind of a tail atom.", "none", "dimensionless code 3"),
        ("SEAMS_SITE_UNSPECIFIED", "Site kind that the table did not name.", "none", "dimensionless code 0"),
        ("SEAMS_SITE_WATER_H", "Site kind of a water hydrogen.", "none", "dimensionless code 9"),
        ("SEAMS_SITE_WATER_O", "Site kind of a water oxygen.", "none", "dimensionless code 8"),
        ("seams_cage_set", "Opaque handle of hexagonal and double-diamond membership.", "none", "handle"),
        ("seams_cages", "Build hexagonal and double-diamond cage membership.", "sys; cutoff; type; graph; k; complete; out", "cutoff in angstrom"),
        ("seams_contact_pairs", "Count mutual nearest unlike contacts of a site table.", "sys; site_spec; n_pairs; n_cation; n_anion", "dimensionless counts"),
        ("seams_coordination", "Coordination number integrated to rmax.", "sys; rmax; bins; type_i; type_j; cn", "rmax in angstrom; cn is dimensionless"),
        ("seams_density_axis", "Number density along one cell axis.", "sys; type; bins; axis; z; rho; cap; n_out", "z in angstrom; rho is number density"),
        ("seams_domains", "Largest polar or apolar domain.", "sys; cutoff; site_spec; subset; n_subset; largest; percolation", "cutoff in angstrom; percolation is dimensionless"),
        ("seams_fingerprint", "Topological fingerprint of the bonded graph.", "sys; cutoff; type; k; graph; hops; colour_types; out", "cutoff in angstrom; hops is a dimensionless depth"),
        ("seams_fingerprint_atom", "Local key of one atom.", "fp; index; buf; cap; n_out", "text"),
        ("seams_fingerprint_destroy", "Free a fingerprint handle.", "fp", "handle"),
        ("seams_fingerprint_key", "Frame key of a fingerprint.", "fp; buf; cap; n_out", "text"),
        ("seams_fingerprint_rings", "Primitive-ring census stored on a fingerprint.", "fp; counts; counts_len", "dimensionless ring counts"),
        ("seams_hbonds", "Count hydrogen bonds on one frame.", "sys; cutoff; oxygen_type; hydrogen_type; dist_cutoff; angle_cutoff; n_bonds", "distances in angstrom; angle in degrees; n_bonds is dimensionless"),
        ("seams_monolayer", "Primitive-ring histogram of a monolayer.", "sys; cutoff; type; max_depth; sheet_area; counts; counts_len", "cutoff in angstrom; sheet_area in square angstrom"),
        ("seams_neighbours_knn", "Build a k-nearest neighbour list.", "sys; k; cutoff; type; mutual; out", "cutoff in angstrom; k is a dimensionless neighbour count"),
        ("seams_prisms", "Perfect and imperfect prism counts by ring size.", "sys; cutoff; type; perfect; imperfect; len", "cutoff in angstrom; counts are dimensionless"),
        ("seams_rdf", "Partial radial distribution.", "sys; rmax; bins; type_i; type_j; r; g; count; cap; n_out; volume", "r and rmax in angstrom; g is dimensionless; volume in cubic angstrom"),
        ("seams_running_cn", "Running coordination number along the partial RDF.", "sys; rmax; bins; type_i; type_j; r; cn; cap; n_out", "r in angstrom; cn is dimensionless"),
        ("seams_signature_cages", "Count closed cages of one ring census.", "sys; cutoff; type; graph; spec; n_cages; n_atoms", "cutoff in angstrom; counts are dimensionless"),
        ("seams_site_kind", "Enumerator of site-table kinds.", "none", "dimensionless"),
        ("seams_site_table", "Per-atom site kind from a site-table spec.", "sys; spec; kinds; len", "dimensionless kind codes"),
        ("seams_steinhardt", "Local Steinhardt q_l on the cutoff shell, with the average on request.", "sys; cutoff; type; order_l; ql; ql_len; qlbar; qlbar_len", "cutoff in angstrom; q_l is dimensionless"),
        ("seams_topo", "Opaque handle of one topological fingerprint.", "none", "handle"),
        ("seams_cages_counts", "Disjoint cage population counts.", "cages; n_hex; n_cubic; n_water; n_both", "dimensionless counts"),
        ("seams_cages_destroy", "Free a cage-membership handle.", "cages", "handle"),
        ("seams_cages_find", "Score hexagonal and double-diamond cages.", "sys; cutoff; type; graph; k; complete; out", "cutoff in angstrom; k is a dimensionless neighbour count"),
        ("seams_cages_labels", "Write one cage code per atom.", "cages; per_atom; len", "dimensionless codes"),
        ("seams_cages_signature", "Count closed cages of one ring census.", "sys; cutoff; type; graph; spec; n_cages; n_atoms", "cutoff in angstrom; counts are dimensionless"),
        ("seams_chill", "Write one original CHILL label per atom.", "sys; cutoff; type; graph; rule; labels; labels_len; counts9", "cutoff in angstrom; labels are dimensionless"),
        ("seams_handle_error", "Return the detail string stored on a handle.", "handle", "text"),
        ("seams_neighbours_copy", "Copy the neighbour list into caller buffers.", "nb; offsets; offsets_cap; n_offsets; cols; cols_cap; nnz", "dimensionless indices; a null buffer queries the length"),
        ("seams_path_nframes", "Count frames in a trajectory path.", "path; nframes", "dimensionless count"),
        ("seams_q6bar", "Write Lechner-Dellago q6-bar and the ice count above a cut.", "sys; cutoff; type; ice_cut; q6bar; q6bar_len; mean; median; n_ice", "cutoff in angstrom; q6-bar is dimensionless"),
        ("seams_ql", "Write local q_l and, on request, its neighbour average.", "sys; cutoff; type; order_l; graph; ql; ql_len; qlbar; qlbar_len", "cutoff in angstrom; q_l is dimensionless"),
        ("seams_rings_vertices", "Copy the atom indices of one ring.", "rings; index; atoms; cap; n_out", "dimensionless indices"),
        ("seams_system_read", "Read one frame from a trajectory path.", "sys; path; frame; type_filter", "frame index is 1-based; coordinates in angstrom"),
        ("seams_system_set_atoms", "Store points and a cell, with size_t lengths.", "sys; xyz; n; box9; box_len; types; ids", "coordinates and box lengths in angstrom"),
        ("seams_system_symmetry", "Write one ice-shell label per atom.", "sys; shell_rmsd; labels; labels_len", "shell_rmsd is a dimensionless residual; labels are dimensionless codes"),
        ("seams_version", "Return the packed release number.", "none", "dimensionless integer, (major << 16) | (minor << 8) | patch"),
        ("seams_version_major", "Return the major version.", "none", "dimensionless integer"),
        ("seams_version_minor", "Return the minor version.", "none", "dimensionless integer"),
        ("seams_version_patch", "Return the patch version.", "none", "dimensionless integer"),
        ("seams_version_string", "Return the version as text.", "none", "text"),
    ],
}


def cxx_rows():
    cut = "cutoff in angstrom; coordinates in angstrom"
    rows = []

    def add(name, purpose, parameters, units):
        rows.append((name, purpose, parameters, units))

    add("seams_frame", "Opaque frame used by the language bindings.", "none", "handle")
    add("seams_frame_box", "Copy the nine LAMMPS cell numbers.", "frame; box9", "lengths in angstrom")
    add("seams_frame_cages", "Score hexagonal and double-diamond cage membership.", "frame; cutoff; type; k; graph; complete; per_atom; n_hex; n_cubic; n_water; n_both", cut + "; k is a dimensionless neighbour count")
    add("seams_frame_chill", "Write one original CHILL label per atom.", "frame; cutoff; type; graph; labels; counts9", cut)
    add("seams_frame_chill_plus", "Write one CHILL+ label per atom.", "frame; cutoff; type; graph; labels; counts9", cut)
    add("seams_frame_chill_plus_rule", "Write CHILL+ labels for a named rule.", "frame; cutoff; type; graph; rule; labels; counts9", cut + "; rule is text")
    add("seams_frame_chill_rings", "Run CHILL+ and primitive rings on one neighbour list.", "frame; cutoff; type; graph; max_ring; skin; labels; counts9; ring_counts; ring_counts_len; times", cut + "; skin in angstrom; phase times in seconds")
    add("seams_frame_free", "Free a language-binding frame.", "frame", "handle")
    add("seams_frame_from_arrays", "Build a frame from coordinates and a cell.", "xyz; n; box9; types; ids", "coordinates and box lengths in angstrom")
    add("seams_frame_from_species", "Build a frame and keep molecule ids.", "xyz; n; box9; types; ids; mols", "coordinates and box lengths in angstrom")
    add("seams_frame_ids", "Pointer to the atom ids stored on the frame.", "frame", "dimensionless ids")
    add("seams_frame_index", "Return the frame index that was read.", "frame", "dimensionless index, 1-based in the file")
    add("seams_frame_monolayer", "Count rings in sheets of a given area.", "frame; cutoff; type; max_depth; sheet_area; output_dir; counts; counts_len", cut + "; sheet_area in square angstrom")
    add("seams_frame_neighbour_pairs", "Return the number of directed neighbour pairs.", "frame", "dimensionless count")
    add("seams_frame_neighbours", "Build the neighbour list and discard the rows.", "frame; cutoff; type; graph", cut)
    add("seams_frame_nop", "Return the number of atoms in the frame.", "frame", "dimensionless count")
    add("seams_frame_positions", "Pointer to the coordinates stored on the frame.", "frame", "angstrom")
    add("seams_frame_prisms", "Count perfect and imperfect prism stacks.", "frame; cutoff; type; max_depth; output_dir; rows; n_rows", cut)
    add("seams_frame_q4_q6", "Write local and averaged q4 and q6.", "frame; cutoff; type; q4; q6; q4bar; q6bar", cut + "; order parameters are dimensionless")
    add("seams_frame_q4_q6_chill", "Write q4, q6, and CHILL+ from one cutoff list.", "frame; cutoff; type; q4; q6; q4bar; q6bar; labels; counts9", cut + "; order parameters are dimensionless")
    add("seams_frame_rdf", "Radial distribution of one type pair.", "frame; rmax; nbins; type_i; type_j; r; g; count; volume", "rmax and r in angstrom; g is dimensionless; volume in cubic angstrom")
    add("seams_frame_rdf2d", "Two-dimensional radial distribution in a monolayer.", "frame; cutoff; binwidth; output_dir; r; g", "cutoff, binwidth, and r in angstrom; g is dimensionless")
    add("seams_frame_rdf2d_nbin", "Number of bins in the monolayer RDF.", "cutoff; binwidth", "cutoff and binwidth in angstrom; result is a dimensionless count")
    add("seams_frame_rdf_nbins", "Number of bins in the bulk RDF.", "rmax; nbins", "rmax in angstrom; nbins is a dimensionless count")
    add("seams_frame_read", "Read one frame from a trajectory path.", "path; frame; type_filter", "frame index is dimensionless; coordinates in the file are angstrom")
    add("seams_frame_rings", "Histogram of primitive rings.", "frame; cutoff; type; graph; max_depth; counts; counts_len", cut + "; counts[k] is a dimensionless ring count")
    add("seams_frame_select", "Keep selected oxygens for later classifiers.", "frame; type; mols; nmol", "type and molecule ids are dimensionless")
    add("seams_frame_set_positions", "Replace coordinates on an array-built frame.", "frame; xyz; n; box9", "coordinates and box lengths in angstrom")
    add("seams_frame_signature", "Count closed cages whose ring census matches a spec.", "frame; cutoff; type; graph; spec; n_cages; n_atoms", cut + "; n_cages and n_atoms are dimensionless")
    add("seams_frame_signature_edge", "Same cage census with an explicit face-edge window.", "frame; cutoff; type; graph; spec; edge; n_cages; n_atoms", cut + "; edge lengths in angstrom")
    add("seams_frame_steinhardt", "Local Steinhardt q_l on the cutoff shell.", "frame; cutoff; type; order_l; ql; qlbar", cut + "; ql is dimensionless")
    add("seams_frame_steinhardt_graph", "Local Steinhardt q_l on a named graph.", "frame; cutoff; type; order_l; graph; ql; qlbar", cut + "; ql is dimensionless")
    add("seams_frame_symmetry", "Per-atom shell label for Ih, Ic, or stacking disorder.", "frame; shell_rmsd; labels", "shell_rmsd is a dimensionless residual; labels are dimensionless codes")
    add("seams_frame_types", "Pointer to the atom types stored on the frame.", "frame", "dimensionless types")
    add("seams_frame_wl", "Normalised Steinhardt w_l for degree 4 or 6.", "frame; cutoff; type; order_l; wl", cut + "; wl is dimensionless")
    add("seams_frames_jump_tau90", "Time for an H-H vector to rotate 90 degrees.", "first; second; dt; oxygen_type; hydrogen_type; tau", "dt and tau in the caller's time unit; types are dimensionless")
    add("seams_lang_version", "Version string of the language-binding library.", "none", "text")
    add("seams_last_error", "Return the thread-local error detail.", "none", "text")
    add("seams_path_nframes", "Count frames in a trajectory path.", "path", "dimensionless count")
    add("seams_phase_times", "Seconds spent in neighbours, CHILL+, and rings.", "neighbours; chill; rings; reused", "times in seconds; reused is 0 or 1")
    add("seams_prism_row", "One row of a prism stack: size, perfect, imperfect.", "size; perfect; imperfect", "dimensionless counts")
    return rows


def other_rows():
    ang = "coordinates and cutoff in angstrom"
    return {
        "python": [
            ("CUTOFF", "Graph id of the cutoff neighbour list.", "none", "dimensionless"),
            ("Frame", "One stored cell and the classifiers that read it.", "none", "object"),
            ("MUTUAL", "Graph id of the mutual neighbour list.", "none", "dimensionless"),
            ("OVITO", "CHILL+ rule id matching OVITO.", "none", "dimensionless"),
            ("PUBLISHED", "CHILL+ rule id of the 2015 paper.", "none", "dimensionless"),
            ("cages", "Score hexagonal and double-diamond cage membership.", "positions; box; or a Frame", ang + "; counts are dimensionless"),
            ("chill_plus", "Label each atom with CHILL+.", "positions; box; rule; graph; cutoff", ang),
            ("count", "Number of points stored on a Frame.", "frame", "dimensionless count"),
            ("from_ase", "Turn an ASE atoms object into positions and a LAMMPS cell.", "atoms", "positions and box lengths in angstrom"),
            ("from_mdanalysis", "Turn an MDAnalysis universe into positions and a LAMMPS cell.", "universe; select", "positions in angstrom; angles inside the universe are degrees"),
            ("neighbours", "Cutoff or mutual neighbour list as CSR.", "positions; box; cutoff; graph", ang),
            ("num_threads", "OpenMP width the next call will use.", "none", "dimensionless thread count"),
            ("q4q6", "Local q4 and q6 on the cutoff shell.", "positions; box; cutoff", ang + "; q4 and q6 are dimensionless"),
            ("rings", "Primitive-ring sizes, symmetry, histogram, and fraction.", "positions; box; cutoff; max_depth", ang + "; max_depth is a dimensionless ring size"),
            ("set_frame", "Store one cell on a Frame.", "frame; positions; box", ang),
            ("set_num_threads", "Set the OpenMP width for the next call.", "n", "dimensionless thread count; 0 restores the library default"),
            ("signature_cages", "Count closed cages matching a ring-census spec.", "positions; box; spec", ang + "; spec is text such as 512"),
            ("symmetry", "Shell label for Ih, Ic, stacking disorder, or other.", "positions; box", ang + "; labels are dimensionless codes"),
            ("version", "Library version text.", "none", "text"),
        ],
        "rust": [
            ("Axis", "Cell axis.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("BondOptions", "Knobs for a hydrogen-bond count.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Brine", "Owned planar ice-front measurement.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("BrineFronts", "The two planar crossings of an ice slab.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("BrineOptions", "Knobs for a planar ice front.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("BrineProfile", "Concentration profile along the brine axis.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("CageCounts", "Disjoint cage populations.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("CageLabel", "Per-atom cage membership.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("CageOptions", "Knobs for hexagonal and double-diamond cages.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Cages", "Owned hexagonal and double-diamond membership.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("ChillOptions", "Knobs for CHILL and CHILL+.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("ChillRule", "CHILL+ threshold table.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Classification", "CHILL or CHILL+ labels and the nine population counts.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Contacts", "Mutual nearest unlike contacts of a site table.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("DensityOptions", "Knobs for the axial number density.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("DensityProfile", "Number density along one cell axis.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Distribution", "Partial radial distribution.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("DistributionOptions", "Knobs for a partial RDF, its running coordination, or the integral.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("DomainOptions", "Knobs for a domain search.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("DomainReport", "Largest polar or apolar domain.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("DomainSubset", "Domain subset.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Error", "A failed call.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Fingerprint", "Owned topological fingerprint.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("FingerprintOptions", "Knobs for a topological fingerprint.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Graph", "Neighbour graph.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("KnnOptions", "Knobs for a k-nearest neighbour list.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Label", "CHILL and CHILL+ atom label.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("NeighbourOptions", "Knobs for a neighbour list.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("NeighbourPair", "One directed neighbour pair.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("NeighbourPairs", "Iterator of directed neighbour pairs.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Neighbours", "Owned neighbour list.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("OrderOptions", "Knobs for q_l on a named graph.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Pockets", "Brine pockets trapped in the ice.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("PrismCensus", "Perfect and imperfect prism counts by ring size.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("PrismOptions", "Knobs for prism stacks.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Q4Q6", "Local q4 and q6 on the cutoff shell.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Q6Bar", "Lechner-Dellago q6-bar and the ice count above a cut.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Q6BarOptions", "Knobs for the q6-bar ice summary.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Rejection", "Molal rejection coefficient.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("RingIter", "Iterator of primitive rings.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("RingOptions", "Knobs for primitive rings.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("RingSymmetry", "Primitive-ring symmetry class.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Rings", "Owned primitive-ring set.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("RunningCoordination", "Running coordination number along a partial RDF.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("SheetOptions", "Knobs for a monolayer ring histogram.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("SignatureMatch", "Closed cages whose ring census matches a signature.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("SignatureOptions", "Knobs for a signature cage census.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("SiteKind", "Site-table kind.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Status", "Status codes.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("Steinhardt", "Local Steinhardt parameter, and the average when it was requested.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("SteinhardtOptions", "Knobs for a local Steinhardt parameter on the cutoff shell.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("System", "One frame: coordinates, a cell, and the analyses that read them.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("abi", "ABI number compiled into this library.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("atom", "Local key of one atom.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("brine", "Planar ice front from an ice assignment.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("cages", "Hexagonal and double-diamond cage membership.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("chill", "Original CHILL labels and nine population counts.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("chill_plus", "CHILL+ labels and nine population counts.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("code", "Status integer.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("contact_pairs", "Mutual nearest unlike contacts.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("coordination", "Coordination number integrated to rmax.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("count", "Number of points in the frame.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("counts", "Disjoint population counts.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("degrees", "Degree of each atom.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("density", "Number density along one cell axis.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("diamond8", "Eight-atom cubic diamond used by the crate examples.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("domains", "Largest polar or apolar domain.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("fingerprint", "Topological fingerprint of the bonded graph.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("fraction", "Fraction of atoms that sit on a ring of one size.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("fronts", "The two planar crossings.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("hbonds", "Hydrogen-bond count.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("histogram", "Histogram of ring size.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("iter", "Walk the rings.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("key", "Frame fingerprint text.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("knn", "k-nearest neighbour list.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("label_name", "Short name of a CHILL+ label.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("labels", "One cage code per atom.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("message", "Detail text.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("monolayer", "Primitive-ring histogram of a monolayer.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("neighbour_degree", "Degree of each atom on the cutoff graph.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("neighbours", "Cutoff or mutual neighbour list.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("new", "Allocate an empty frame.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("nframes", "Count frames in a trajectory path.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("nnz", "Number of neighbour pairs, omitting self.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("num_threads", "OpenMP width the next call will use.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("pairs", "Directed pairs.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("pockets", "Pockets trapped in the ice, and the liquid reservoir.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("prisms", "Perfect and imperfect prism counts by ring size.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("profile", "Concentration profile along the front axis.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("q4q6", "Local q4 and q6 on the cutoff shell.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("q6bar", "q6-bar, its mean and median, and the ice count above the cut.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("ql", "Local q_l on a named graph.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("rdf", "Partial radial distribution and the cell volume.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("read", "Read one frame from a trajectory path.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("rejection", "Molal rejection coefficient.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("ring_symmetry_name", "Short name of a ring-symmetry code.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("rings", "Primitive-ring census stored on the fingerprint.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("running_cn", "Running coordination number along the partial RDF.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("set_atoms", "Store points, a cell, and optional species and ids.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("set_frame", "Store points and a general nine-number cell.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("set_num_threads", "Set the OpenMP width for the next call.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("set_ortho", "Store points in an orthorhombic cell whose origin is zero.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("signature_cages", "Closed cages of one ring census.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("site_table", "Per-atom site kind.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("six_ring_fraction", "Count of six-rings and the fraction of atoms that sit on one.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("size", "Number of primitive rings.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("status", "Status when the failure came from the library and the code is known.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("steinhardt", "Local q_l on the cutoff shell.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("symmetry", "Per-atom ice-shell label.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("version", "Library version text, for example 2.14.0.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
            ("vertices", "Atom indices of one ring.", "named in the rustdoc on this item", "named in the rustdoc on this item; lengths are angstroms"),
        ],
        "fortran": [
            ("SEAMS_CHILL_OVITO", "CHILL+ rule used by OVITO.", "none", "dimensionless code"),
            ("SEAMS_CHILL_PUBLISHED", "CHILL+ rule from the 2015 paper.", "none", "dimensionless code"),
            ("SEAMS_GRAPH_CUTOFF", "Neighbour graph of pairs inside the cutoff.", "none", "dimensionless code"),
            ("SEAMS_GRAPH_MUTUAL", "Mutual four-nearest neighbour graph.", "none", "dimensionless code"),
            ("SEAMS_LABEL_CUBIC", "CHILL+ cubic-ice label. Fortran arrays of labels are 1-based, so clathrate is counts(5).", "none", "dimensionless code 0"),
            ("SEAMS_OK", "The call succeeded.", "none", "dimensionless status code"),
            ("seams_chill_plus", "Write one CHILL+ label per atom.", "sys; cutoff; type; graph; rule; labels; counts", "cutoff in angstrom; counts are dimensionless"),
            ("seams_get_num_threads", "Report the OpenMP width the next call will use.", "none", "dimensionless thread count"),
            ("seams_last_error", "Return the thread-local error detail.", "none", "text"),
            ("seams_neighbours_build", "Build a cutoff or mutual neighbour list.", "sys; cutoff; type; graph; out", "cutoff in angstrom"),
            ("seams_neighbours_csr", "Return the neighbour list as CSR.", "nb; offsets; n_offsets; cols; nnz", "dimensionless indices"),
            ("seams_neighbours_destroy", "Free a neighbour-list handle.", "nb", "handle"),
            ("seams_q4q6", "Write local q4 and q6.", "sys; cutoff; type; q4; q6", "cutoff in angstrom; q4 and q6 are dimensionless"),
            ("seams_rings_destroy", "Free a ring-set handle.", "rings", "handle"),
            ("seams_rings_find", "Find primitive rings up to a maximum size.", "sys; cutoff; type; graph; max_depth; out", "cutoff in angstrom; max_depth is a dimensionless ring size"),
            ("seams_rings_fraction", "Fraction of atoms on a ring of one size.", "rings; size; fraction", "dimensionless fraction"),
            ("seams_rings_histogram", "Count primitive rings by size.", "rings; counts; counts_len", "dimensionless counts"),
            ("seams_set_num_threads", "Set the OpenMP width for the next call.", "n", "dimensionless thread count"),
            ("seams_system_create", "Allocate an empty frame handle.", "none", "handle"),
            ("seams_system_destroy", "Free a frame handle.", "sys", "handle"),
            ("seams_system_set_frame", "Store points and a nine-number cell.", "sys; xyz; n; box; types; ids", "coordinates and box lengths in angstrom"),
            ("seams_system_set_ortho", "Store points and an orthorhombic cell.", "sys; xyz; n; ox; oy; oz; lx; ly; lz; types; ids", "origin and edge lengths in angstrom"),
            ("seams_version_string", "Return the version as text.", "none", "text"),
        ],
        "julia": [
            ("Frame", "One frame read from a dump or built in memory.", "none", "object"),
            ("cages", "Hexagonal and double-diamond cage counts.", "frame; cutoff", "cutoff in angstrom; counts are dimensionless"),
            ("chill_plus", "CHILL+ label of each atom.", "frame; cutoff", "cutoff in angstrom; labels are dimensionless codes"),
            ("jump_tau90", "Time for an H-H vector to rotate 90 degrees between two frames.", "first; second; dt", "dt and the result in the caller's time unit"),
            ("nframes", "Number of frames in a trajectory path.", "path", "dimensionless count"),
            ("nop", "Number of atoms in a frame.", "frame", "dimensionless count"),
            ("primitive_rings", "Histogram of primitive rings.", "frame; cutoff; max_ring", "cutoff in angstrom; bin k is a dimensionless ring count"),
            ("q4_q6", "Local q4 and q6.", "frame; cutoff", "cutoff in angstrom; q4 and q6 are dimensionless"),
            ("read_frame", "Read one frame from a dump.", "path; frame; type", "frame index is dimensionless; coordinates are angstrom"),
            ("signature_cages", "Closed cages matching a ring-census spec.", "frame; cutoff; spec", "cutoff in angstrom; spec is text"),
            ("steinhardt", "Mean local Steinhardt q_l.", "frame; cutoff; order_l", "cutoff in angstrom; ql is dimensionless"),
            ("symmetry", "Shell label for Ih, Ic, or stacking disorder.", "frame", "labels are dimensionless codes"),
        ],
        "r": [
            ("cages", "Cage counts for a dump path.", "path; frame; cutoff", "cutoff in angstrom; frame index is dimensionless"),
            ("cages.dseams_frame", "Cage counts for a frame object.", "frame; cutoff", "cutoff in angstrom; counts are dimensionless"),
            ("chill_plus", "CHILL+ labels for a dump path.", "path; frame; cutoff", "cutoff in angstrom"),
            ("chill_plus.dseams_frame", "CHILL+ labels for a frame object.", "frame; cutoff", "cutoff in angstrom; labels are dimensionless codes"),
            ("jump_tau90", "Rotation time of H-H vectors between two frames.", "first; second; dt", "dt and the result in the caller's time unit"),
            ("nframes", "Number of frames in a trajectory path.", "path", "dimensionless count"),
            ("positions", "Coordinates of a dump frame.", "path; frame", "angstrom"),
            ("positions.dseams_frame", "Coordinates stored on a frame object.", "frame", "angstrom"),
            ("primitive_rings", "Ring histogram for a dump path.", "path; frame; cutoff", "cutoff in angstrom"),
            ("primitive_rings.dseams_frame", "Ring histogram for a frame object.", "frame; cutoff", "cutoff in angstrom; counts are dimensionless"),
            ("print.dseams_frame", "Print the atom count of a frame.", "frame", "dimensionless count"),
            ("q4_q6", "Local q4 and q6 for a dump path.", "path; frame; cutoff", "cutoff in angstrom; q4 and q6 are dimensionless"),
            ("q4_q6.dseams_frame", "Local q4 and q6 for a frame object.", "frame; cutoff", "cutoff in angstrom; q4 and q6 are dimensionless"),
            ("read_frame", "Read one frame from a dump.", "path; frame; type", "frame index is dimensionless; coordinates are angstrom"),
            ("signature_cages", "Closed cages for a dump path.", "path; frame; spec", "spec is text; cutoff used by the call is in angstrom"),
            ("signature_cages.dseams_frame", "Closed cages for a frame object.", "frame; spec", "spec is text"),
            ("steinhardt", "Mean local q_l for a dump path.", "path; frame; order_l", "order parameter is dimensionless"),
            ("steinhardt.dseams_frame", "Mean local q_l for a frame object.", "frame; order_l", "order parameter is dimensionless"),
            ("symmetry", "Shell labels for a dump path.", "path; frame", "dimensionless label codes"),
            ("symmetry.dseams_frame", "Shell labels for a frame object.", "frame", "dimensionless label codes"),
        ],
    }


def groups():
    found = {"c": list(RAW["c"]), "cxx": cxx_rows()}
    found.update(other_rows())
    return found


def render(rows, example):
    lines = ["", "## Symbol reference", ""]
    for name, purpose, parameters, units in rows:
        lines.append("### {}".format(name))
        lines.append("")
        lines.append("Purpose: {}".format(purpose))
        lines.append("")
        lines.append("Parameters: {}".format(parameters))
        lines.append("")
        lines.append("Units: {}".format(units))
        lines.append("")
        lines.append("Example: `{}`".format(example))
        lines.append("")
    return "\n".join(lines)


def splice(path, block):
    text = path.read_text()
    wrapped = BEGIN + "\n" + block + END + "\n"
    if BEGIN in text and END in text:
        text = re.sub(
            re.escape(BEGIN) + r".*?" + re.escape(END) + r"\n?",
            wrapped,
            text,
            count=1,
            flags=re.S,
        )
    else:
        if not text.endswith("\n"):
            text += "\n"
        text += "\n" + wrapped
    path.write_text(text)


def main():
    import sys

    wanted = {
        "c": check_c_symbols.c_header_symbols(),
        "cxx": check_c_symbols.cxx_symbols(),
        "python": check_c_symbols.python_symbols(),
        "rust": check_c_symbols.rust_symbols(),
        "fortran": check_c_symbols.fortran_symbols(),
        "julia": check_c_symbols.julia_symbols(),
        "r": check_c_symbols.r_symbols(),
    }
    have = groups()
    errors = []
    for lang, names in wanted.items():
        got = [row[0] for row in have[lang]]
        if got != names and sorted(got) != list(names):
            missing = sorted(set(names) - set(got))
            extra = sorted(set(got) - set(names))
            if missing or extra:
                errors.append("{} missing {} extra {}".format(lang, missing, extra))
        if len(got) != len(set(got)):
            errors.append("{} has a duplicate name".format(lang))
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    order = {lang: {name: i for i, name in enumerate(names)} for lang, names in wanted.items()}
    for lang, rows in have.items():
        rows.sort(key=lambda row: order[lang][row[0]])
        splice(API / PAGES[lang], render(rows, EXAMPLES[lang]))
    print("wrote symbol reference for {} names".format(sum(len(v) for v in wanted.values())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
