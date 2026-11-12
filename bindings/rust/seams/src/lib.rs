/* Copyright (c) d-SEAMS developers */
//! Safe calls into the d-SEAMS C API.

use seams_sys as ffi;
use std::ffi::CStr;
use std::os::raw::c_int;

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

fn check(code: c_int) -> Result<(), Error> {
    if code == ffi::SEAMS_OK {
        return Ok(());
    }
    let message = unsafe { CStr::from_ptr(ffi::seams_last_error()) }
        .to_string_lossy()
        .into_owned();
    Err(Error { code, message })
}

pub fn version() -> String {
    unsafe { CStr::from_ptr(ffi::seams_version_string()) }
        .to_string_lossy()
        .into_owned()
}

pub fn set_num_threads(n: i32) -> Result<(), Error> {
    check(unsafe { ffi::seams_set_num_threads(n) })
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

impl System {
    pub fn new() -> Result<Self, Error> {
        let ptr = unsafe { ffi::seams_system_create() };
        if ptr.is_null() {
            return check(ffi::SEAMS_ERR_NOMEM).map(|_| unreachable!());
        }
        Ok(Self { ptr })
    }

    pub fn set_ortho(&mut self, xyz: &[f64], lx: f64, ly: f64, lz: f64) -> Result<(), Error> {
        if xyz.len() % 3 != 0 {
            return Err(Error {
                code: ffi::SEAMS_ERR_ARG,
                message: "coordinates are not a multiple of 3".into(),
            });
        }
        let n = (xyz.len() / 3) as c_int;
        check(unsafe {
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
        })
    }

    pub fn set_frame(&mut self, xyz: &[f64], box9: &[f64; 9]) -> Result<(), Error> {
        if xyz.len() % 3 != 0 {
            return Err(Error {
                code: ffi::SEAMS_ERR_ARG,
                message: "coordinates are not a multiple of 3".into(),
            });
        }
        let n = (xyz.len() / 3) as c_int;
        check(unsafe {
            ffi::seams_system_set_frame(
                self.ptr,
                xyz.as_ptr(),
                n,
                box9.as_ptr(),
                std::ptr::null(),
                std::ptr::null(),
            )
        })
    }

    pub fn count(&self) -> i32 {
        unsafe { ffi::seams_system_count(self.ptr) }
    }

    pub fn chill_plus(&mut self, cutoff: f64, graph: i32, rule: i32) -> Result<(Vec<i32>, [i32; 9]), Error> {
        let n = self.count() as usize;
        let mut labels = vec![0; n];
        let mut counts = [0; 9];
        check(unsafe {
            ffi::seams_chill_plus(
                self.ptr,
                cutoff,
                0,
                graph,
                rule,
                labels.as_mut_ptr(),
                counts.as_mut_ptr(),
            )
        })?;
        Ok((labels, counts))
    }

    pub fn q4q6(&mut self, cutoff: f64) -> Result<(Vec<f64>, Vec<f64>), Error> {
        let n = self.count() as usize;
        let mut q4 = vec![0.0; n];
        let mut q6 = vec![0.0; n];
        check(unsafe { ffi::seams_q4q6(self.ptr, cutoff, 0, q4.as_mut_ptr(), q6.as_mut_ptr()) })?;
        Ok((q4, q6))
    }

    pub fn neighbour_degree(&mut self, cutoff: f64) -> Result<Vec<i32>, Error> {
        let mut nb = std::ptr::null_mut();
        check(unsafe {
            ffi::seams_neighbours_build(self.ptr, cutoff, 0, ffi::SEAMS_GRAPH_CUTOFF, &mut nb)
        })?;
        let mut offsets = std::ptr::null();
        let mut cols = std::ptr::null();
        let mut n_offsets = 0;
        let mut nnz = 0;
        let code = unsafe { ffi::seams_neighbours_csr(nb, &mut offsets, &mut n_offsets, &mut cols, &mut nnz) };
        if code != ffi::SEAMS_OK {
            unsafe { ffi::seams_neighbours_destroy(nb) };
            return check(code).map(|()| Vec::new());
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

    pub fn six_ring_fraction(&mut self, cutoff: f64) -> Result<(i32, f64), Error> {
        let mut rings = std::ptr::null_mut();
        check(unsafe {
            ffi::seams_rings_find(self.ptr, cutoff, 0, ffi::SEAMS_GRAPH_CUTOFF, 6, &mut rings)
        })?;
        let mut hist = [0; 7];
        let mut fraction = 0.0;
        let code = unsafe { ffi::seams_rings_histogram(rings, hist.as_mut_ptr(), 7) };
        let code2 = unsafe { ffi::seams_rings_fraction(rings, 6, &mut fraction) };
        unsafe { ffi::seams_rings_destroy(rings) };
        check(code)?;
        check(code2)?;
        Ok((hist[6], fraction))
    }
}
