/* Copyright (c) d-SEAMS developers */
/* Generated from include/seams.h by gen_sys.py. Do not edit. */
#![allow(non_camel_case_types, non_snake_case)]

use std::os::raw::{c_char, c_int, c_void};

pub const SEAMS_OK: c_int = 0;
pub const SEAMS_ERR_NULL: c_int = 1;
pub const SEAMS_ERR_RANGE: c_int = 2;
pub const SEAMS_ERR_BOX: c_int = 3;
pub const SEAMS_ERR_ARG: c_int = 4;
pub const SEAMS_ERR_STATE: c_int = 5;
pub const SEAMS_ERR_NOMEM: c_int = 6;
pub const SEAMS_ERR_INTERNAL: c_int = 7;
pub const SEAMS_LABEL_CUBIC: c_int = 0;
pub const SEAMS_LABEL_HEXAGONAL: c_int = 1;
pub const SEAMS_LABEL_WATER: c_int = 2;
pub const SEAMS_LABEL_INTERFACIAL: c_int = 3;
pub const SEAMS_LABEL_CLATHRATE: c_int = 4;
pub const SEAMS_LABEL_INTERCLATHRATE: c_int = 5;
pub const SEAMS_LABEL_UNCLASSIFIED: c_int = 6;
pub const SEAMS_LABEL_RECUBIC: c_int = 7;
pub const SEAMS_LABEL_REHEX: c_int = 8;
pub const SEAMS_RING_UNCLASSIFIED: c_int = 0;
pub const SEAMS_RING_DDC: c_int = 1;
pub const SEAMS_RING_HC_BASAL: c_int = 2;
pub const SEAMS_RING_HC_PRISMATIC: c_int = 3;
pub const SEAMS_RING_BOTH_BASAL: c_int = 4;
pub const SEAMS_RING_BOTH_PRISMATIC: c_int = 5;
pub const SEAMS_RING_PRISM: c_int = 6;
pub const SEAMS_RING_DEFORMED_PRISM: c_int = 7;
pub const SEAMS_RING_MIXED_PRISM: c_int = 8;
pub const SEAMS_RING_PENTAGON: c_int = 9;
pub const SEAMS_RING_HEXAGON: c_int = 10;
pub const SEAMS_RING_OTHER: c_int = 11;
pub const SEAMS_CHILL_PUBLISHED: c_int = 0;
pub const SEAMS_CHILL_OVITO: c_int = 1;
pub const SEAMS_GRAPH_CUTOFF: c_int = 0;
pub const SEAMS_GRAPH_MUTUAL: c_int = 1;
pub const SEAMS_GRAPH_SEEDED: c_int = 2;
pub const SEAMS_CAGE_WATER: c_int = 0;
pub const SEAMS_CAGE_HEX: c_int = 1;
pub const SEAMS_CAGE_CUBIC: c_int = 2;
pub const SEAMS_CAGE_BOTH: c_int = 3;
pub const SEAMS_SITE_UNSPECIFIED: c_int = 0;
pub const SEAMS_SITE_CATION: c_int = 1;
pub const SEAMS_SITE_ANION: c_int = 2;
pub const SEAMS_SITE_TAIL: c_int = 3;
pub const SEAMS_SITE_DONOR_H: c_int = 4;
pub const SEAMS_SITE_ACCEPTOR: c_int = 5;
pub const SEAMS_SITE_POLAR: c_int = 6;
pub const SEAMS_SITE_APOLAR: c_int = 7;
pub const SEAMS_SITE_WATER_O: c_int = 8;
pub const SEAMS_SITE_WATER_H: c_int = 9;
pub const SEAMS_SITE_SOLVENT: c_int = 10;

#[repr(C)]
pub struct seams_system {
    _private: [u8; 0],
}
#[repr(C)]
pub struct seams_neighbours {
    _private: [u8; 0],
}
#[repr(C)]
pub struct seams_rings {
    _private: [u8; 0],
}
#[repr(C)]
pub struct seams_cage_set {
    _private: [u8; 0],
}
#[repr(C)]
pub struct seams_brine {
    _private: [u8; 0],
}
#[repr(C)]
pub struct seams_topo {
    _private: [u8; 0],
}

extern "C" {
    pub fn seams_abi() -> c_int;
    pub fn seams_version() -> c_int;
    pub fn seams_version_major() -> c_int;
    pub fn seams_version_minor() -> c_int;
    pub fn seams_version_patch() -> c_int;
    pub fn seams_version_string() -> *const c_char;
    pub fn seams_status_string(
        status: c_int,
    ) -> *const c_char;
    pub fn seams_last_error() -> *const c_char;
    pub fn seams_handle_error(
        handle: *const c_void,
    ) -> *const c_char;
    pub fn seams_label_name(
        label: c_int,
    ) -> *const c_char;
    pub fn seams_ring_symmetry_name(
        symmetry: c_int,
    ) -> *const c_char;
    pub fn seams_set_num_threads(
        n: c_int,
    ) -> c_int;
    pub fn seams_get_num_threads() -> c_int;
    pub fn seams_system_create() -> *mut seams_system;
    pub fn seams_system_destroy(
        sys: *mut seams_system,
    ) -> ();
    pub fn seams_system_set_atoms(
        sys: *mut seams_system,
        xyz: *const f64,
        n: usize,
        box9: *const f64,
        box_len: usize,
        types: *const c_int,
        ids: *const c_int,
    ) -> c_int;
    pub fn seams_system_set_frame(
        sys: *mut seams_system,
        xyz: *const f64,
        n: c_int,
        box9: *const f64,
        types: *const c_int,
        ids: *const c_int,
    ) -> c_int;
    pub fn seams_system_set_ortho(
        sys: *mut seams_system,
        xyz: *const f64,
        n: c_int,
        ox: f64,
        oy: f64,
        oz: f64,
        lx: f64,
        ly: f64,
        lz: f64,
        types: *const c_int,
        ids: *const c_int,
    ) -> c_int;
    pub fn seams_system_read(
        sys: *mut seams_system,
        path: *const c_char,
        frame: c_int,
        type_filter: c_int,
    ) -> c_int;
    pub fn seams_path_nframes(
        path: *const c_char,
        nframes: *mut usize,
    ) -> c_int;
    pub fn seams_system_count(
        sys: *const seams_system,
    ) -> c_int;
    pub fn seams_neighbours_build(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        graph: c_int,
        out: *mut *mut seams_neighbours,
    ) -> c_int;
    pub fn seams_neighbours_destroy(
        nb: *mut seams_neighbours,
    ) -> ();
    pub fn seams_neighbours_count(
        nb: *const seams_neighbours,
    ) -> c_int;
    pub fn seams_neighbours_csr(
        nb: *const seams_neighbours,
        offsets: *mut *const c_int,
        n_offsets: *mut c_int,
        cols: *mut *const c_int,
        nnz: *mut c_int,
    ) -> c_int;
    pub fn seams_neighbours_copy(
        nb: *const seams_neighbours,
        offsets: *mut c_int,
        offsets_cap: usize,
        n_offsets: *mut usize,
        cols: *mut c_int,
        cols_cap: usize,
        nnz: *mut usize,
    ) -> c_int;
    pub fn seams_chill(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        graph: c_int,
        rule: c_int,
        labels: *mut c_int,
        labels_len: usize,
        counts9: *mut c_int,
    ) -> c_int;
    pub fn seams_chill_plus(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        graph: c_int,
        rule: c_int,
        labels: *mut c_int,
        counts9: *mut c_int,
    ) -> c_int;
    pub fn seams_ql(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        order_l: c_int,
        graph: c_int,
        ql: *mut f64,
        ql_len: usize,
        qlbar: *mut f64,
        qlbar_len: usize,
    ) -> c_int;
    pub fn seams_q4q6(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        q4: *mut f64,
        q6: *mut f64,
    ) -> c_int;
    pub fn seams_q6bar(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        ice_cut: f64,
        q6bar: *mut f64,
        q6bar_len: usize,
        mean: *mut f64,
        median: *mut f64,
        n_ice: *mut c_int,
    ) -> c_int;
    pub fn seams_rings_find(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        graph: c_int,
        max_depth: c_int,
        out: *mut *mut seams_rings,
    ) -> c_int;
    pub fn seams_rings_destroy(
        rings: *mut seams_rings,
    ) -> ();
    pub fn seams_rings_count(
        rings: *const seams_rings,
    ) -> c_int;
    pub fn seams_rings_histogram(
        rings: *const seams_rings,
        counts: *mut c_int,
        counts_len: c_int,
    ) -> c_int;
    pub fn seams_rings_labels(
        rings: *const seams_rings,
        sizes: *mut c_int,
        symmetry: *mut c_int,
    ) -> c_int;
    pub fn seams_rings_vertices(
        rings: *const seams_rings,
        index: c_int,
        atoms: *mut c_int,
        cap: usize,
        n_out: *mut usize,
    ) -> c_int;
    pub fn seams_rings_fraction(
        rings: *const seams_rings,
        size: c_int,
        fraction: *mut f64,
    ) -> c_int;
    pub fn seams_cages_find(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        graph: c_int,
        k: c_int,
        complete: c_int,
        out: *mut *mut seams_cage_set,
    ) -> c_int;
    pub fn seams_cages_destroy(
        cages: *mut seams_cage_set,
    ) -> ();
    pub fn seams_cages_counts(
        cages: *const seams_cage_set,
        n_hex: *mut c_int,
        n_cubic: *mut c_int,
        n_water: *mut c_int,
        n_both: *mut c_int,
    ) -> c_int;
    pub fn seams_cages_labels(
        cages: *const seams_cage_set,
        per_atom: *mut c_int,
        len: usize,
    ) -> c_int;
    pub fn seams_cages_signature(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        graph: c_int,
        spec: *const c_char,
        n_cages: *mut usize,
        n_atoms: *mut usize,
    ) -> c_int;
    pub fn seams_system_symmetry(
        sys: *mut seams_system,
        shell_rmsd: f64,
        labels: *mut c_int,
        labels_len: usize,
    ) -> c_int;
    pub fn seams_brine_analyse(
        sys: *mut seams_system,
        ice: *const c_int,
        n_ice: usize,
        water_type: c_int,
        ion_type: c_int,
        n_bins: c_int,
        axis: c_int,
        out: *mut *mut seams_brine,
    ) -> c_int;
    pub fn seams_brine_destroy(
        front: *mut seams_brine,
    ) -> ();
    pub fn seams_brine_fronts(
        front: *const seams_brine,
        z_enter: *mut f64,
        z_exit: *mut f64,
        thickness: *mut f64,
        found: *mut c_int,
    ) -> c_int;
    pub fn seams_brine_rejection(
        front: *const seams_brine,
        c_ice: *mut f64,
        c_liquid: *mut f64,
        k: *mut f64,
        defined: *mut c_int,
    ) -> c_int;
    pub fn seams_brine_profile(
        front: *const seams_brine,
        z_mid: *mut f64,
        molality: *mut f64,
        ice_fraction: *mut f64,
        cap: usize,
        n_out: *mut usize,
    ) -> c_int;
    pub fn seams_brine_pockets(
        front: *const seams_brine,
        n_pockets: *mut c_int,
        reservoir_waters: *mut c_int,
    ) -> c_int;
    pub fn seams_cages(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        graph: c_int,
        k: c_int,
        complete: c_int,
        out: *mut *mut seams_cage_set,
    ) -> c_int;
    pub fn seams_signature_cages(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        graph: c_int,
        spec: *const c_char,
        n_cages: *mut usize,
        n_atoms: *mut usize,
    ) -> c_int;
    pub fn seams_steinhardt(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        order_l: c_int,
        ql: *mut f64,
        ql_len: usize,
        qlbar: *mut f64,
        qlbar_len: usize,
    ) -> c_int;
    pub fn seams_neighbours_knn(
        sys: *mut seams_system,
        k: c_int,
        cutoff: f64,
        r#type: c_int,
        mutual: c_int,
        out: *mut *mut seams_neighbours,
    ) -> c_int;
    pub fn seams_rdf(
        sys: *mut seams_system,
        rmax: f64,
        bins: c_int,
        type_i: c_int,
        type_j: c_int,
        r: *mut f64,
        g: *mut f64,
        count: *mut c_int,
        cap: usize,
        n_out: *mut usize,
        volume: *mut f64,
    ) -> c_int;
    pub fn seams_running_cn(
        sys: *mut seams_system,
        rmax: f64,
        bins: c_int,
        type_i: c_int,
        type_j: c_int,
        r: *mut f64,
        cn: *mut f64,
        cap: usize,
        n_out: *mut usize,
    ) -> c_int;
    pub fn seams_coordination(
        sys: *mut seams_system,
        rmax: f64,
        bins: c_int,
        type_i: c_int,
        type_j: c_int,
        cn: *mut f64,
    ) -> c_int;
    pub fn seams_hbonds(
        sys: *mut seams_system,
        cutoff: f64,
        oxygen_type: c_int,
        hydrogen_type: c_int,
        dist_cutoff: f64,
        angle_cutoff: f64,
        n_bonds: *mut c_int,
    ) -> c_int;
    pub fn seams_density_axis(
        sys: *mut seams_system,
        r#type: c_int,
        bins: c_int,
        axis: c_int,
        z: *mut f64,
        rho: *mut f64,
        cap: usize,
        n_out: *mut usize,
    ) -> c_int;
    pub fn seams_fingerprint(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        k: c_int,
        graph: c_int,
        hops: c_int,
        colour_types: c_int,
        out: *mut *mut seams_topo,
    ) -> c_int;
    pub fn seams_fingerprint_destroy(
        fp: *mut seams_topo,
    ) -> ();
    pub fn seams_fingerprint_key(
        fp: *const seams_topo,
        buf: *mut c_char,
        cap: usize,
        n_out: *mut usize,
    ) -> c_int;
    pub fn seams_fingerprint_atom(
        fp: *const seams_topo,
        index: c_int,
        buf: *mut c_char,
        cap: usize,
        n_out: *mut usize,
    ) -> c_int;
    pub fn seams_fingerprint_rings(
        fp: *const seams_topo,
        counts: *mut c_int,
        counts_len: usize,
    ) -> c_int;
    pub fn seams_contact_pairs(
        sys: *mut seams_system,
        site_spec: *const c_char,
        n_pairs: *mut usize,
        n_cation: *mut c_int,
        n_anion: *mut c_int,
    ) -> c_int;
    pub fn seams_domains(
        sys: *mut seams_system,
        cutoff: f64,
        site_spec: *const c_char,
        subset: c_int,
        n_subset: *mut c_int,
        largest: *mut c_int,
        percolation: *mut f64,
    ) -> c_int;
    pub fn seams_prisms(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        perfect: *mut c_int,
        imperfect: *mut c_int,
        len: usize,
    ) -> c_int;
    pub fn seams_monolayer(
        sys: *mut seams_system,
        cutoff: f64,
        r#type: c_int,
        max_depth: c_int,
        sheet_area: f64,
        counts: *mut c_int,
        counts_len: usize,
    ) -> c_int;
    pub fn seams_site_table(
        sys: *mut seams_system,
        spec: *const c_char,
        kinds: *mut c_int,
        len: usize,
    ) -> c_int;
}
