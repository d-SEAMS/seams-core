/* Copyright (c) d-SEAMS developers */
//! Safe calls into the d-SEAMS C API.

use seams_sys as ffi;
use std::ffi::{CStr, CString};
use std::os::raw::{c_int, c_void};

#[derive(Debug)]
pub struct Error {
    pub code: c_int,
    pub message: String,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ({})", self.message, self.code)
    }
}

impl std::error::Error for Error {}

fn message_of(handle: *const c_void) -> String {
    unsafe {
        let ptr = ffi::seams_handle_error(handle);
        if ptr.is_null() {
            return String::new();
        }
        CStr::from_ptr(ptr).to_string_lossy().into_owned()
    }
}

fn check(code: c_int, handle: *const c_void) -> Result<(), Error> {
    if code == ffi::SEAMS_OK {
        return Ok(());
    }
    Err(Error {
        code,
        message: message_of(handle),
    })
}

fn arg_error(message: &str) -> Error {
    Error {
        code: ffi::SEAMS_ERR_ARG,
        message: message.to_string(),
    }
}

fn c_text(text: &str) -> Result<CString, Error> {
    CString::new(text).map_err(|_| arg_error("text contains a nul"))
}

pub fn version() -> String {
    unsafe { CStr::from_ptr(ffi::seams_version_string()) }
        .to_string_lossy()
        .into_owned()
}

pub fn set_num_threads(n: i32) -> Result<(), Error> {
    check(unsafe { ffi::seams_set_num_threads(n) }, std::ptr::null())
}

pub fn num_threads() -> i32 {
    unsafe { ffi::seams_get_num_threads() }
}

pub struct System {
    ptr: *mut ffi::seams_system,
}

impl Drop for System {
    fn drop(&mut self) {
        unsafe { ffi::seams_system_destroy(self.ptr) }
    }
}

pub struct Neighbours {
    ptr: *mut ffi::seams_neighbours,
}

impl Drop for Neighbours {
    fn drop(&mut self) {
        unsafe { ffi::seams_neighbours_destroy(self.ptr) }
    }
}

impl Neighbours {
    pub fn nnz(&self) -> Result<usize, Error> {
        let handle = self.ptr as *const c_void;
        let mut n_offsets = 0usize;
        let mut nnz = 0usize;
        check(
            unsafe {
                ffi::seams_neighbours_copy(
                    self.ptr,
                    std::ptr::null_mut(),
                    0,
                    &mut n_offsets,
                    std::ptr::null_mut(),
                    0,
                    &mut nnz,
                )
            },
            handle,
        )?;
        Ok(nnz)
    }
}

pub struct Rings {
    ptr: *mut ffi::seams_rings,
}

impl Drop for Rings {
    fn drop(&mut self) {
        unsafe { ffi::seams_rings_destroy(self.ptr) }
    }
}

impl Rings {
    pub fn count(&self) -> i32 {
        unsafe { ffi::seams_rings_count(self.ptr) }
    }

    pub fn vertices(&self, index: i32) -> Result<Vec<i32>, Error> {
        let handle = self.ptr as *const c_void;
        let mut n_out = 0usize;
        check(
            unsafe {
                ffi::seams_rings_vertices(self.ptr, index, std::ptr::null_mut(), 0, &mut n_out)
            },
            handle,
        )?;
        let mut atoms = vec![0; n_out];
        check(
            unsafe {
                ffi::seams_rings_vertices(
                    self.ptr,
                    index,
                    atoms.as_mut_ptr(),
                    n_out,
                    &mut n_out,
                )
            },
            handle,
        )?;
        atoms.truncate(n_out);
        Ok(atoms)
    }
}

pub struct Cages {
    ptr: *mut ffi::seams_cage_set,
}

impl Drop for Cages {
    fn drop(&mut self) {
        unsafe { ffi::seams_cages_destroy(self.ptr) }
    }
}

impl Cages {
    pub fn counts(&self) -> Result<(i32, i32, i32, i32), Error> {
        let handle = self.ptr as *const c_void;
        let mut n_hex = 0;
        let mut n_cubic = 0;
        let mut n_water = 0;
        let mut n_both = 0;
        check(
            unsafe {
                ffi::seams_cages_counts(
                    self.ptr,
                    &mut n_hex,
                    &mut n_cubic,
                    &mut n_water,
                    &mut n_both,
                )
            },
            handle,
        )?;
        Ok((n_hex, n_cubic, n_water, n_both))
    }
}

pub struct Fingerprint {
    ptr: *mut ffi::seams_topo,
}

impl Drop for Fingerprint {
    fn drop(&mut self) {
        unsafe { ffi::seams_fingerprint_destroy(self.ptr) }
    }
}

impl Fingerprint {
    pub fn key(&self) -> Result<String, Error> {
        let handle = self.ptr as *const c_void;
        let mut n_out = 0usize;
        check(
            unsafe { ffi::seams_fingerprint_key(self.ptr, std::ptr::null_mut(), 0, &mut n_out) },
            handle,
        )?;
        let mut buf = vec![0u8; n_out + 1];
        check(
            unsafe {
                ffi::seams_fingerprint_key(
                    self.ptr,
                    buf.as_mut_ptr() as *mut _,
                    buf.len(),
                    &mut n_out,
                )
            },
            handle,
        )?;
        buf.truncate(n_out);
        Ok(String::from_utf8_lossy(&buf).into_owned())
    }

    pub fn rings(&self) -> Result<Vec<i32>, Error> {
        let handle = self.ptr as *const c_void;
        let mut counts = vec![0; 13];
        check(
            unsafe { ffi::seams_fingerprint_rings(self.ptr, counts.as_mut_ptr(), counts.len()) },
            handle,
        )?;
        Ok(counts)
    }
}

impl System {
    pub fn new() -> Result<Self, Error> {
        let ptr = unsafe { ffi::seams_system_create() };
        if ptr.is_null() {
            return check(ffi::SEAMS_ERR_NOMEM, std::ptr::null()).map(|_| unreachable!());
        }
        Ok(Self { ptr })
    }

    fn handle(&self) -> *const c_void {
        self.ptr as *const c_void
    }

    pub fn set_ortho(&mut self, xyz: &[f64], lx: f64, ly: f64, lz: f64) -> Result<(), Error> {
        if xyz.len() % 3 != 0 {
            return Err(arg_error("coordinates are not a multiple of 3"));
        }
        let n = (xyz.len() / 3) as c_int;
        check(
            unsafe {
                ffi::seams_system_set_ortho(
                    self.ptr,
                    xyz.as_ptr(),
                    n,
                    0.0,
                    0.0,
                    0.0,
                    lx,
                    ly,
                    lz,
                    std::ptr::null(),
                    std::ptr::null(),
                )
            },
            self.handle(),
        )
    }

    pub fn set_frame(&mut self, xyz: &[f64], box9: &[f64; 9]) -> Result<(), Error> {
        if xyz.len() % 3 != 0 {
            return Err(arg_error("coordinates are not a multiple of 3"));
        }
        let n = (xyz.len() / 3) as c_int;
        check(
            unsafe {
                ffi::seams_system_set_frame(
                    self.ptr,
                    xyz.as_ptr(),
                    n,
                    box9.as_ptr(),
                    std::ptr::null(),
                    std::ptr::null(),
                )
            },
            self.handle(),
        )
    }

    pub fn read(&mut self, path: &str, frame: i32, type_filter: i32) -> Result<(), Error> {
        let text = c_text(path)?;
        check(
            unsafe { ffi::seams_system_read(self.ptr, text.as_ptr(), frame, type_filter) },
            self.handle(),
        )
    }

    pub fn count(&self) -> i32 {
        unsafe { ffi::seams_system_count(self.ptr) }
    }

    pub fn chill(
        &mut self,
        cutoff: f64,
        type_filter: i32,
        graph: i32,
        rule: i32,
    ) -> Result<(Vec<i32>, [i32; 9]), Error> {
        let n = self.count() as usize;
        let mut labels = vec![0; n];
        let mut counts = [0; 9];
        check(
            unsafe {
                ffi::seams_chill(
                    self.ptr,
                    cutoff,
                    type_filter,
                    graph,
                    rule,
                    labels.as_mut_ptr(),
                    n,
                    counts.as_mut_ptr(),
                )
            },
            self.handle(),
        )?;
        Ok((labels, counts))
    }

    pub fn chill_plus(
        &mut self,
        cutoff: f64,
        graph: i32,
        rule: i32,
    ) -> Result<(Vec<i32>, [i32; 9]), Error> {
        let n = self.count() as usize;
        let mut labels = vec![0; n];
        let mut counts = [0; 9];
        check(
            unsafe {
                ffi::seams_chill_plus(
                    self.ptr,
                    cutoff,
                    0,
                    graph,
                    rule,
                    labels.as_mut_ptr(),
                    counts.as_mut_ptr(),
                )
            },
            self.handle(),
        )?;
        Ok((labels, counts))
    }

    pub fn steinhardt(
        &mut self,
        cutoff: f64,
        order: i32,
        average: bool,
    ) -> Result<(Vec<f64>, Option<Vec<f64>>), Error> {
        let n = self.count() as usize;
        let mut ql = vec![0.0; n];
        let mut bar = if average { vec![0.0; n] } else { Vec::new() };
        let (bar_ptr, bar_len) = if average {
            (bar.as_mut_ptr(), n)
        } else {
            (std::ptr::null_mut(), 0)
        };
        check(
            unsafe {
                ffi::seams_steinhardt(
                    self.ptr,
                    cutoff,
                    0,
                    order,
                    ql.as_mut_ptr(),
                    n,
                    bar_ptr,
                    bar_len,
                )
            },
            self.handle(),
        )?;
        Ok((ql, if average { Some(bar) } else { None }))
    }

    pub fn q4q6(&mut self, cutoff: f64) -> Result<(Vec<f64>, Vec<f64>), Error> {
        let n = self.count() as usize;
        let mut q4 = vec![0.0; n];
        let mut q6 = vec![0.0; n];
        check(
            unsafe { ffi::seams_q4q6(self.ptr, cutoff, 0, q4.as_mut_ptr(), q6.as_mut_ptr()) },
            self.handle(),
        )?;
        Ok((q4, q6))
    }

    pub fn neighbour_degree(&mut self, cutoff: f64) -> Result<Vec<i32>, Error> {
        let mut nb = std::ptr::null_mut();
        check(
            unsafe {
                ffi::seams_neighbours_build(
                    self.ptr,
                    cutoff,
                    0,
                    ffi::SEAMS_GRAPH_CUTOFF,
                    &mut nb,
                )
            },
            self.handle(),
        )?;
        let mut offsets = std::ptr::null();
        let mut cols = std::ptr::null();
        let mut n_offsets = 0;
        let mut nnz = 0;
        let code = unsafe {
            ffi::seams_neighbours_csr(nb, &mut offsets, &mut n_offsets, &mut cols, &mut nnz)
        };
        if code != ffi::SEAMS_OK {
            unsafe { ffi::seams_neighbours_destroy(nb) };
            return check(code, nb as *const c_void).map(|()| Vec::new());
        }
        let mut degree = Vec::with_capacity((n_offsets - 1) as usize);
        unsafe {
            for i in 0..(n_offsets - 1) {
                degree.push(*offsets.add(i as usize + 1) - *offsets.add(i as usize));
            }
            ffi::seams_neighbours_destroy(nb);
        }
        Ok(degree)
    }

    pub fn knn(&mut self, k: i32, cutoff: f64, mutual: bool) -> Result<Neighbours, Error> {
        let mut out = std::ptr::null_mut();
        check(
            unsafe {
                ffi::seams_neighbours_knn(
                    self.ptr,
                    k,
                    cutoff,
                    0,
                    if mutual { 1 } else { 0 },
                    &mut out,
                )
            },
            self.handle(),
        )?;
        Ok(Neighbours { ptr: out })
    }

    pub fn rings(&mut self, cutoff: f64, max_depth: i32) -> Result<Rings, Error> {
        let mut out = std::ptr::null_mut();
        check(
            unsafe {
                ffi::seams_rings_find(
                    self.ptr,
                    cutoff,
                    0,
                    ffi::SEAMS_GRAPH_CUTOFF,
                    max_depth,
                    &mut out,
                )
            },
            self.handle(),
        )?;
        Ok(Rings { ptr: out })
    }

    pub fn six_ring_fraction(&mut self, cutoff: f64) -> Result<(i32, f64), Error> {
        let rings = self.rings(cutoff, 6)?;
        let mut hist = [0; 7];
        let mut fraction = 0.0;
        let handle = rings.ptr as *const c_void;
        let code = unsafe { ffi::seams_rings_histogram(rings.ptr, hist.as_mut_ptr(), 7) };
        let code2 = unsafe { ffi::seams_rings_fraction(rings.ptr, 6, &mut fraction) };
        check(code, handle)?;
        check(code2, handle)?;
        Ok((hist[6], fraction))
    }

    pub fn cages(&mut self, cutoff: f64, graph: i32, k: i32, complete: bool) -> Result<Cages, Error> {
        let mut out = std::ptr::null_mut();
        check(
            unsafe {
                ffi::seams_cages(
                    self.ptr,
                    cutoff,
                    0,
                    graph,
                    k,
                    if complete { 1 } else { 0 },
                    &mut out,
                )
            },
            self.handle(),
        )?;
        Ok(Cages { ptr: out })
    }

    pub fn signature_cages(
        &mut self,
        cutoff: f64,
        graph: i32,
        spec: &str,
    ) -> Result<(usize, usize), Error> {
        let text = c_text(spec)?;
        let mut n_cages = 0usize;
        let mut n_atoms = 0usize;
        check(
            unsafe {
                ffi::seams_signature_cages(
                    self.ptr,
                    cutoff,
                    0,
                    graph,
                    text.as_ptr(),
                    &mut n_cages,
                    &mut n_atoms,
                )
            },
            self.handle(),
        )?;
        Ok((n_cages, n_atoms))
    }

    pub fn rdf(&mut self, rmax: f64, bins: i32) -> Result<(Vec<f64>, Vec<f64>, f64), Error> {
        let mut n_out = 0usize;
        let mut volume = 0.0;
        check(
            unsafe {
                ffi::seams_rdf(
                    self.ptr,
                    rmax,
                    bins,
                    0,
                    0,
                    std::ptr::null_mut(),
                    std::ptr::null_mut(),
                    std::ptr::null_mut(),
                    0,
                    &mut n_out,
                    &mut volume,
                )
            },
            self.handle(),
        )?;
        let mut r = vec![0.0; n_out];
        let mut g = vec![0.0; n_out];
        check(
            unsafe {
                ffi::seams_rdf(
                    self.ptr,
                    rmax,
                    bins,
                    0,
                    0,
                    r.as_mut_ptr(),
                    g.as_mut_ptr(),
                    std::ptr::null_mut(),
                    n_out,
                    &mut n_out,
                    std::ptr::null_mut(),
                )
            },
            self.handle(),
        )?;
        r.truncate(n_out);
        g.truncate(n_out);
        Ok((r, g, volume))
    }

    pub fn running_cn(&mut self, rmax: f64, bins: i32) -> Result<(Vec<f64>, Vec<f64>), Error> {
        let mut n_out = 0usize;
        check(
            unsafe {
                ffi::seams_running_cn(
                    self.ptr,
                    rmax,
                    bins,
                    0,
                    0,
                    std::ptr::null_mut(),
                    std::ptr::null_mut(),
                    0,
                    &mut n_out,
                )
            },
            self.handle(),
        )?;
        let mut r = vec![0.0; n_out];
        let mut cn = vec![0.0; n_out];
        check(
            unsafe {
                ffi::seams_running_cn(
                    self.ptr,
                    rmax,
                    bins,
                    0,
                    0,
                    r.as_mut_ptr(),
                    cn.as_mut_ptr(),
                    n_out,
                    &mut n_out,
                )
            },
            self.handle(),
        )?;
        r.truncate(n_out);
        cn.truncate(n_out);
        Ok((r, cn))
    }

    pub fn coordination(&mut self, rmax: f64, bins: i32) -> Result<f64, Error> {
        let mut cn = 0.0;
        check(
            unsafe { ffi::seams_coordination(self.ptr, rmax, bins, 0, 0, &mut cn) },
            self.handle(),
        )?;
        Ok(cn)
    }

    pub fn hbonds(
        &mut self,
        cutoff: f64,
        oxygen_type: i32,
        hydrogen_type: i32,
        dist: f64,
        angle: f64,
    ) -> Result<i32, Error> {
        let mut n_bonds = 0;
        check(
            unsafe {
                ffi::seams_hbonds(
                    self.ptr,
                    cutoff,
                    oxygen_type,
                    hydrogen_type,
                    dist,
                    angle,
                    &mut n_bonds,
                )
            },
            self.handle(),
        )?;
        Ok(n_bonds)
    }

    pub fn density_axis(
        &mut self,
        type_filter: i32,
        bins: i32,
        axis: i32,
    ) -> Result<(Vec<f64>, Vec<f64>), Error> {
        let mut n_out = 0usize;
        check(
            unsafe {
                ffi::seams_density_axis(
                    self.ptr,
                    type_filter,
                    bins,
                    axis,
                    std::ptr::null_mut(),
                    std::ptr::null_mut(),
                    0,
                    &mut n_out,
                )
            },
            self.handle(),
        )?;
        let mut z = vec![0.0; n_out];
        let mut rho = vec![0.0; n_out];
        check(
            unsafe {
                ffi::seams_density_axis(
                    self.ptr,
                    type_filter,
                    bins,
                    axis,
                    z.as_mut_ptr(),
                    rho.as_mut_ptr(),
                    n_out,
                    &mut n_out,
                )
            },
            self.handle(),
        )?;
        z.truncate(n_out);
        rho.truncate(n_out);
        Ok((z, rho))
    }

    pub fn fingerprint(&mut self, cutoff: f64, hops: i32) -> Result<Fingerprint, Error> {
        let mut out = std::ptr::null_mut();
        check(
            unsafe {
                ffi::seams_fingerprint(
                    self.ptr,
                    cutoff,
                    0,
                    4,
                    ffi::SEAMS_GRAPH_CUTOFF,
                    hops,
                    0,
                    &mut out,
                )
            },
            self.handle(),
        )?;
        Ok(Fingerprint { ptr: out })
    }

    pub fn contact_pairs(&mut self, spec: &str) -> Result<(usize, i32, i32), Error> {
        let text = c_text(spec)?;
        let mut n_pairs = 0usize;
        let mut n_cation = 0;
        let mut n_anion = 0;
        check(
            unsafe {
                ffi::seams_contact_pairs(
                    self.ptr,
                    text.as_ptr(),
                    &mut n_pairs,
                    &mut n_cation,
                    &mut n_anion,
                )
            },
            self.handle(),
        )?;
        Ok((n_pairs, n_cation, n_anion))
    }

    pub fn domains(&mut self, cutoff: f64, spec: &str, subset: i32) -> Result<(i32, i32, f64), Error> {
        let text = c_text(spec)?;
        let mut n_subset = 0;
        let mut largest = 0;
        let mut percolation = 0.0;
        check(
            unsafe {
                ffi::seams_domains(
                    self.ptr,
                    cutoff,
                    text.as_ptr(),
                    subset,
                    &mut n_subset,
                    &mut largest,
                    &mut percolation,
                )
            },
            self.handle(),
        )?;
        Ok((n_subset, largest, percolation))
    }

    pub fn prisms(&mut self, cutoff: f64) -> Result<([i32; 7], [i32; 7]), Error> {
        let mut perfect = [0; 7];
        let mut imperfect = [0; 7];
        check(
            unsafe {
                ffi::seams_prisms(
                    self.ptr,
                    cutoff,
                    0,
                    perfect.as_mut_ptr(),
                    imperfect.as_mut_ptr(),
                    7,
                )
            },
            self.handle(),
        )?;
        Ok((perfect, imperfect))
    }

    pub fn monolayer(&mut self, cutoff: f64, max_depth: i32, sheet_area: f64) -> Result<Vec<i32>, Error> {
        let len = (max_depth + 1) as usize;
        let mut counts = vec![0; len];
        check(
            unsafe {
                ffi::seams_monolayer(
                    self.ptr,
                    cutoff,
                    0,
                    max_depth,
                    sheet_area,
                    counts.as_mut_ptr(),
                    len,
                )
            },
            self.handle(),
        )?;
        Ok(counts)
    }

    pub fn site_table(&mut self, spec: &str) -> Result<Vec<i32>, Error> {
        let text = c_text(spec)?;
        let n = self.count() as usize;
        let mut kinds = vec![0; n];
        check(
            unsafe { ffi::seams_site_table(self.ptr, text.as_ptr(), kinds.as_mut_ptr(), n) },
            self.handle(),
        )?;
        Ok(kinds)
    }
}
