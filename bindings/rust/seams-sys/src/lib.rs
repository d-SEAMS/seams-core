/* Copyright (c) d-SEAMS developers */
#![allow(non_camel_case_types)]

use std::os::raw::{c_char, c_int};

pub const SEAMS_OK: c_int = 0;
pub const SEAMS_ERR_NULL: c_int = 1;
pub const SEAMS_ERR_RANGE: c_int = 2;
pub const SEAMS_ERR_BOX: c_int = 3;
pub const SEAMS_ERR_ARG: c_int = 4;
pub const SEAMS_ERR_STATE: c_int = 5;
pub const SEAMS_ERR_NOMEM: c_int = 6;
pub const SEAMS_ERR_INTERNAL: c_int = 7;

pub const SEAMS_CHILL_PUBLISHED: c_int = 0;
pub const SEAMS_CHILL_OVITO: c_int = 1;
pub const SEAMS_GRAPH_CUTOFF: c_int = 0;
pub const SEAMS_GRAPH_MUTUAL: c_int = 1;
pub const SEAMS_LABEL_CUBIC: c_int = 0;

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

extern "C" {
    pub fn seams_abi() -> c_int;
    pub fn seams_version_string() -> *const c_char;
    pub fn seams_status_string(status: c_int) -> *const c_char;
    pub fn seams_last_error() -> *const c_char;
    pub fn seams_set_num_threads(n: c_int) -> c_int;
    pub fn seams_get_num_threads() -> c_int;
    pub fn seams_system_create() -> *mut seams_system;
    pub fn seams_system_destroy(sys: *mut seams_system);
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
    pub fn seams_system_count(sys: *const seams_system) -> c_int;
    pub fn seams_neighbours_build(
        sys: *mut seams_system,
        cutoff: f64,
        type_filter: c_int,
        graph: c_int,
        out: *mut *mut seams_neighbours,
    ) -> c_int;
    pub fn seams_neighbours_destroy(nb: *mut seams_neighbours);
    pub fn seams_neighbours_csr(
        nb: *const seams_neighbours,
        offsets: *mut *const c_int,
        n_offsets: *mut c_int,
        cols: *mut *const c_int,
        nnz: *mut c_int,
    ) -> c_int;
    pub fn seams_chill_plus(
        sys: *mut seams_system,
        cutoff: f64,
        type_filter: c_int,
        graph: c_int,
        rule: c_int,
        labels: *mut c_int,
        counts9: *mut c_int,
    ) -> c_int;
    pub fn seams_q4q6(
        sys: *mut seams_system,
        cutoff: f64,
        type_filter: c_int,
        q4: *mut f64,
        q6: *mut f64,
    ) -> c_int;
    pub fn seams_rings_find(
        sys: *mut seams_system,
        cutoff: f64,
        type_filter: c_int,
        graph: c_int,
        max_depth: c_int,
        out: *mut *mut seams_rings,
    ) -> c_int;
    pub fn seams_rings_destroy(rings: *mut seams_rings);
    pub fn seams_rings_count(rings: *const seams_rings) -> c_int;
    pub fn seams_rings_histogram(rings: *const seams_rings, counts: *mut c_int, counts_len: c_int) -> c_int;
    pub fn seams_rings_labels(rings: *const seams_rings, sizes: *mut c_int, symmetry: *mut c_int) -> c_int;
    pub fn seams_rings_fraction(rings: *const seams_rings, size: c_int, fraction: *mut f64) -> c_int;
}
