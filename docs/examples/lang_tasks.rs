//! Neighbours, rings, and CHILL+ from Rust.
//! Each library call is a `Result`. The histograms are iterators.

use std::env;
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_void};
use std::ptr;

#[link(name = "seams_lang")]
extern "C" {
    fn seams_frame_read(path: *const c_char, frame: c_int, type_filter: c_int) -> *mut c_void;
    fn seams_frame_free(frame: *mut c_void);
    fn seams_frame_nop(frame: *mut c_void) -> c_int;
    fn seams_frame_neighbours(
        frame: *mut c_void,
        cutoff: f64,
        type_filter: c_int,
        graph: *const c_char,
    ) -> c_int;
    fn seams_frame_neighbour_pairs(frame: *const c_void) -> c_int;
    fn seams_frame_rings(
        frame: *mut c_void,
        cutoff: f64,
        type_filter: c_int,
        graph: *const c_char,
        max_depth: c_int,
        counts: *mut c_int,
        counts_len: c_int,
    ) -> c_int;
    fn seams_frame_chill_plus(
        frame: *mut c_void,
        cutoff: f64,
        type_filter: c_int,
        graph: *const c_char,
        labels: *mut c_int,
        counts9: *mut c_int,
    ) -> c_int;
    fn seams_last_error() -> *const c_char;
}

struct Frame {
    ptr: *mut c_void,
}

impl Drop for Frame {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { seams_frame_free(self.ptr) }
            self.ptr = ptr::null_mut();
        }
    }
}

fn last_error() -> String {
    unsafe {
        let text = seams_last_error();
        if text.is_null() {
            "seams call failed".to_string()
        } else {
            CStr::from_ptr(text).to_string_lossy().into_owned()
        }
    }
}

fn check(code: c_int) -> Result<(), String> {
    if code == 0 { Ok(()) } else { Err(last_error()) }
}

impl Frame {
    fn read(path: &str) -> Result<Self, String> {
        let path = CString::new(path).map_err(|err| err.to_string())?;
        let ptr = unsafe { seams_frame_read(path.as_ptr(), 1, -1) };
        if ptr.is_null() {
            Err(last_error())
        } else {
            Ok(Frame { ptr })
        }
    }

    fn pairs(&self) -> Result<i32, String> {
        let graph = CString::new("cutoff").unwrap();
        unsafe {
            check(seams_frame_neighbours(self.ptr, 3.5, -1, graph.as_ptr()))?;
            Ok(seams_frame_neighbour_pairs(self.ptr))
        }
    }

    fn rings(&self) -> Result<Vec<i32>, String> {
        let graph = CString::new("cutoff").unwrap();
        let mut counts = vec![0; 7];
        unsafe {
            check(seams_frame_rings(
                self.ptr,
                3.5,
                -1,
                graph.as_ptr(),
                6,
                counts.as_mut_ptr(),
                counts.len() as c_int,
            ))?;
        }
        Ok(counts)
    }

    fn chill(&self) -> Result<Vec<i32>, String> {
        let graph = CString::new("cutoff").unwrap();
        let n = unsafe { seams_frame_nop(self.ptr) };
        let mut labels = vec![0; n as usize];
        let mut counts = vec![0; 9];
        unsafe {
            check(seams_frame_chill_plus(
                self.ptr,
                3.5,
                -1,
                graph.as_ptr(),
                labels.as_mut_ptr(),
                counts.as_mut_ptr(),
            ))?;
        }
        Ok(counts)
    }
}

fn joined(values: &[i32]) -> String {
    values.iter().map(|value| value.to_string()).collect::<Vec<_>>().join(" ")
}

fn main() -> Result<(), String> {
    let root = env::var("DSEAMS_ROOT").map_err(|_| "DSEAMS_ROOT is unset".to_string())?;
    let frame = Frame::read(&format!("{root}/input/traj/ice1h.lammpstrj"))?;
    let pairs = frame.pairs()?;
    let rings = frame.rings()?;
    let chill = frame.chill()?;
    println!("pairs {pairs} rings {} chill {}", joined(&rings), joined(&chill));
    Ok(())
}
