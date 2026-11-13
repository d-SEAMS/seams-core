/* Copyright (c) d-SEAMS developers */
//! Safe Rust for the d-SEAMS C API.
//!
//! `seams-sys` stays a raw FFI crate. This crate owns every handle, frees
//! it in [`Drop`], and returns [`Result<T, Error>`](Result). The error text
//! is the string from `seams_handle_error`. Coordinates and cutoffs are
//! angstroms. A cell is nine LAMMPS numbers: origin, edge lengths, then
//! the tilts `xy`, `xz`, `yz`.
//!
//! Analyses with several knobs take an options struct. Neighbour lists and
//! rings yield iterators. A handle may move to another thread ([`Send`]).
//! It is not [`Sync`]: the C library requires exclusive use of a handle,
//! and the Steinhardt scratch is process-global.
//!
//! # Quick example
//!
//! ```
//! let (xyz, edge) = seams::diamond8();
//! seams::set_num_threads(1)?;
//! let mut system = seams::System::new()?;
//! system.set_ortho(&xyz, edge, edge, edge)?;
//! assert_eq!(system.count(), 8);
//! let ice = system.chill_plus(&seams::ChillOptions::default())?;
//! assert_eq!(ice.counts[seams::Label::Cubic as usize], 8);
//! let orders = system.q4q6(3.5)?;
//! assert!(orders.q6.iter().all(|value| value.is_finite()));
//! # Ok::<(), seams::Error>(())
//! ```
//!
//! Build with `SEAMS_LIB_DIR` pointing at the directory that contains
//! `libseams.so`.

#![deny(missing_docs)]

use seams_sys as ffi;
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_void};

/// Status codes. The numeric values match `seams_status`.
///
/// ```
/// assert_eq!(seams::Status::Ok as i32, 0);
/// assert_eq!(seams::Status::State as i32, 5);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum Status {
    /// The call succeeded.
    Ok = 0,
    /// A required pointer is null.
    Null = 1,
    /// A number is outside the accepted range.
    Range = 2,
    /// The cell is missing or has no volume.
    Cell = 3,
    /// A caller argument is the wrong shape.
    Argument = 4,
    /// The handle has no frame yet.
    State = 5,
    /// An allocation failed.
    Memory = 6,
    /// The library hit an internal failure.
    Internal = 7,
}

impl Status {
    fn from_code(code: i32) -> Option<Self> {
        Some(match code {
            0 => Self::Ok,
            1 => Self::Null,
            2 => Self::Range,
            3 => Self::Cell,
            4 => Self::Argument,
            5 => Self::State,
            6 => Self::Memory,
            7 => Self::Internal,
            _ => return None,
        })
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Ok => "ok",
            Self::Null => "null pointer",
            Self::Range => "out of range",
            Self::Cell => "bad cell",
            Self::Argument => "bad argument",
            Self::State => "no frame",
            Self::Memory => "out of memory",
            Self::Internal => "internal error",
        }
    }
}

/// Neighbour graph.
///
/// ```
/// assert_eq!(seams::Graph::Cutoff as i32, 0);
/// assert_eq!(seams::Graph::Mutual as i32, 1);
/// assert_eq!(seams::Graph::Seeded as i32, 2);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum Graph {
    /// Every pair inside the cutoff.
    Cutoff = ffi::SEAMS_GRAPH_CUTOFF,
    /// Mutual four-nearest neighbours.
    Mutual = ffi::SEAMS_GRAPH_MUTUAL,
    /// Cage graph: mutual seeds and union completions.
    Seeded = ffi::SEAMS_GRAPH_SEEDED,
}

/// CHILL+ threshold table.
///
/// Published edges are inclusive. OVITO edges are open.
///
/// ```
/// assert_eq!(seams::ChillRule::Published as i32, 0);
/// assert_eq!(seams::ChillRule::Ovito as i32, 1);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum ChillRule {
    /// Nguyen and Molinero, 2015.
    Published = ffi::SEAMS_CHILL_PUBLISHED,
    /// OVITO Chill+ modifier.
    Ovito = ffi::SEAMS_CHILL_OVITO,
}

/// CHILL and CHILL+ atom label. `Classification::counts[label as usize]` is the population.
///
/// ```
/// assert_eq!(seams::Label::Hexagonal as i32, 1);
/// assert_eq!(seams::label_name(seams::Label::Hexagonal), "hexagonal");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum Label {
    /// Cubic ice.
    Cubic = ffi::SEAMS_LABEL_CUBIC,
    /// Hexagonal ice.
    Hexagonal = ffi::SEAMS_LABEL_HEXAGONAL,
    /// Liquid water.
    Water = ffi::SEAMS_LABEL_WATER,
    /// Interfacial ice.
    Interfacial = ffi::SEAMS_LABEL_INTERFACIAL,
    /// Clathrate.
    Clathrate = ffi::SEAMS_LABEL_CLATHRATE,
    /// Interfacial clathrate.
    Interclathrate = ffi::SEAMS_LABEL_INTERCLATHRATE,
    /// No class matched.
    Unclassified = ffi::SEAMS_LABEL_UNCLASSIFIED,
    /// Reclassified cubic.
    Recubic = ffi::SEAMS_LABEL_RECUBIC,
    /// Reclassified hexagonal.
    Rehex = ffi::SEAMS_LABEL_REHEX,
}

/// Primitive-ring symmetry class.
///
/// ```
/// assert_eq!(seams::ring_symmetry_name(seams::RingSymmetry::Hexagon), "hexagon");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum RingSymmetry {
    /// No class yet.
    Unclassified = ffi::SEAMS_RING_UNCLASSIFIED,
    /// Double-diamond six-ring.
    Ddc = ffi::SEAMS_RING_DDC,
    /// Basal face of a hexagonal cage.
    HcBasal = ffi::SEAMS_RING_HC_BASAL,
    /// Prismatic face of a hexagonal cage.
    HcPrismatic = ffi::SEAMS_RING_HC_PRISMATIC,
    /// Basal face of both cage types.
    BothBasal = ffi::SEAMS_RING_BOTH_BASAL,
    /// Prismatic face of both cage types.
    BothPrismatic = ffi::SEAMS_RING_BOTH_PRISMATIC,
    /// Prismatic six-ring.
    Prism = ffi::SEAMS_RING_PRISM,
    /// Deformed prismatic six-ring.
    DeformedPrism = ffi::SEAMS_RING_DEFORMED_PRISM,
    /// Prism mixing the two cage types.
    MixedPrism = ffi::SEAMS_RING_MIXED_PRISM,
    /// Five-membered ring.
    Pentagon = ffi::SEAMS_RING_PENTAGON,
    /// Six-ring that is not a cage face.
    Hexagon = ffi::SEAMS_RING_HEXAGON,
    /// Neither a pentagon nor a hexagon face.
    Other = ffi::SEAMS_RING_OTHER,
}

/// Per-atom cage membership. The four codes are disjoint.
///
/// ```
/// assert_eq!(seams::CageLabel::Hex as i32, 1);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum CageLabel {
    /// In neither cage.
    Water = ffi::SEAMS_CAGE_WATER,
    /// Hexagonal cage only.
    Hex = ffi::SEAMS_CAGE_HEX,
    /// Double-diamond cage only.
    Cubic = ffi::SEAMS_CAGE_CUBIC,
    /// Both cages.
    Both = ffi::SEAMS_CAGE_BOTH,
}

/// Site-table kind.
///
/// ```
/// assert_eq!(seams::SiteKind::WaterO as i32, 8);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum SiteKind {
    /// The table did not name this atom.
    Unspecified = ffi::SEAMS_SITE_UNSPECIFIED,
    /// Cation head.
    Cation = ffi::SEAMS_SITE_CATION,
    /// Anion.
    Anion = ffi::SEAMS_SITE_ANION,
    /// Tail atom.
    Tail = ffi::SEAMS_SITE_TAIL,
    /// Donor hydrogen.
    DonorH = ffi::SEAMS_SITE_DONOR_H,
    /// Acceptor.
    Acceptor = ffi::SEAMS_SITE_ACCEPTOR,
    /// Polar atom.
    Polar = ffi::SEAMS_SITE_POLAR,
    /// Apolar atom.
    Apolar = ffi::SEAMS_SITE_APOLAR,
    /// Water oxygen.
    WaterO = ffi::SEAMS_SITE_WATER_O,
    /// Water hydrogen.
    WaterH = ffi::SEAMS_SITE_WATER_H,
    /// Solvent atom.
    Solvent = ffi::SEAMS_SITE_SOLVENT,
}

/// Cell axis.
///
/// ```
/// assert_eq!(seams::Axis::Z as i32, 2);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum Axis {
    /// First cell edge.
    X = 0,
    /// Second cell edge.
    Y = 1,
    /// Third cell edge.
    Z = 2,
}

/// Domain subset. The values match the C API.
///
/// ```
/// assert_eq!(seams::DomainSubset::Polar as i32, 1);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum DomainSubset {
    /// Apolar site subset.
    Apolar = 0,
    /// Polar site subset.
    Polar = 1,
}

/// A failed call. The message is `seams_handle_error` when the call had a handle.
///
/// ```
/// let mut system = seams::System::new()?;
/// let err = system
///     .chill_plus(&seams::ChillOptions::default())
///     .expect_err("an empty frame has no coordinates");
/// assert_eq!(err.status(), Some(seams::Status::State));
/// assert!(!err.message().is_empty());
/// # Ok::<(), seams::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The library rejected the call with a known status.
    Library {
        /// Status code.
        status: Status,
        /// Detail text. May be empty.
        message: String,
    },
    /// The library returned a status this crate does not name.
    Other {
        /// Raw status integer.
        code: i32,
        /// Detail text.
        message: String,
    },
    /// The argument cannot be passed to the C API.
    Argument {
        /// What was wrong.
        message: String,
    },
}

impl Error {
    fn argument(message: impl Into<String>) -> Self {
        Self::Argument {
            message: message.into(),
        }
    }

    fn from_status(code: i32, message: String) -> Self {
        match Status::from_code(code) {
            Some(status) => Self::Library { status, message },
            None => Self::Other { code, message },
        }
    }

    /// Status when the failure came from the library and the code is known.
    ///
    /// ```
    /// let err = seams::nframes("\0").expect_err("a nul cannot be a path");
    /// assert_eq!(err.status(), None);
    /// assert_eq!(err.code(), seams_sys::SEAMS_ERR_ARG);
    /// ```
    #[must_use]
    pub fn status(&self) -> Option<Status> {
        match self {
            Self::Library { status, .. } => Some(*status),
            Self::Other { .. } | Self::Argument { .. } => None,
        }
    }

    /// Status integer. An argument error uses `SEAMS_ERR_ARG`.
    ///
    /// ```
    /// let err = seams::nframes("\0").unwrap_err();
    /// assert_eq!(err.code(), 4);
    /// ```
    #[must_use]
    pub fn code(&self) -> i32 {
        match self {
            Self::Library { status, .. } => *status as i32,
            Self::Other { code, .. } => *code,
            Self::Argument { .. } => ffi::SEAMS_ERR_ARG,
        }
    }

    /// Detail text.
    ///
    /// ```
    /// let err = seams::nframes("\0").unwrap_err();
    /// assert_eq!(err.message(), "text contains a nul");
    /// ```
    #[must_use]
    pub fn message(&self) -> &str {
        match self {
            Self::Library { message, .. }
            | Self::Other { message, .. }
            | Self::Argument { message } => message,
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Library { status, message } => {
                if message.is_empty() {
                    write!(f, "{}", status.as_str())
                } else {
                    write!(f, "{message} ({})", status.as_str())
                }
            }
            Self::Other { code, message } => {
                if message.is_empty() {
                    write!(f, "status {code}")
                } else {
                    write!(f, "{message} (status {code})")
                }
            }
            Self::Argument { message } => write!(f, "{message}"),
        }
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
    Err(Error::from_status(code, message_of(handle)))
}

fn c_text(text: &str) -> Result<CString, Error> {
    CString::new(text).map_err(|_| Error::argument("text contains a nul"))
}

fn static_text(ptr: *const c_char) -> String {
    if ptr.is_null() {
        return String::new();
    }
    unsafe { CStr::from_ptr(ptr).to_string_lossy().into_owned() }
}

fn atom_count(n_coords: usize) -> Result<i32, Error> {
    if n_coords % 3 != 0 {
        return Err(Error::argument("coordinates are not a multiple of 3"));
    }
    i32::try_from(n_coords / 3).map_err(|_| Error::argument("too many points"))
}

/// Eight-atom cubic diamond used by the crate examples.
///
/// The lattice constant is `2.75 * 4 / sqrt(3)` angstroms. The returned
/// edge is one cubic cell of that crystal.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// assert_eq!(xyz.len(), 24);
/// assert!(edge > 3.0);
/// ```
#[must_use]
pub fn diamond8() -> (Vec<f64>, f64) {
    let a = 2.75 * 4.0 / 3.0_f64.sqrt();
    let fcc = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let basis = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]];
    let mut xyz = Vec::with_capacity(24);
    for shift in fcc {
        for site in basis {
            xyz.push((shift[0] + site[0]) * a);
            xyz.push((shift[1] + site[1]) * a);
            xyz.push((shift[2] + site[2]) * a);
        }
    }
    (xyz, a)
}

/// Library version text, for example `2.14.0`.
///
/// ```
/// assert!(!seams::version().is_empty());
/// ```
#[must_use]
pub fn version() -> String {
    static_text(unsafe { ffi::seams_version_string() })
}

/// ABI number compiled into this library. The value is 1.
///
/// ```
/// assert_eq!(seams::abi(), 1);
/// ```
#[must_use]
pub fn abi() -> i32 {
    unsafe { ffi::seams_abi() }
}

/// Set the OpenMP width for the next call. Zero restores the library default.
///
/// The width is process-global. Do not call this while another thread is
/// inside the library.
///
/// ```
/// seams::set_num_threads(1)?;
/// assert_eq!(seams::num_threads(), 1);
/// # Ok::<(), seams::Error>(())
/// ```
pub fn set_num_threads(n: i32) -> Result<(), Error> {
    check(unsafe { ffi::seams_set_num_threads(n) }, std::ptr::null())
}

/// OpenMP width the next call will use.
///
/// ```
/// seams::set_num_threads(1)?;
/// assert_eq!(seams::num_threads(), 1);
/// # Ok::<(), seams::Error>(())
/// ```
#[must_use]
pub fn num_threads() -> i32 {
    unsafe { ffi::seams_get_num_threads() }
}

/// Count frames in a trajectory path.
///
/// ```
/// let err = seams::nframes("\0").expect_err("a nul cannot be a path");
/// assert_eq!(err.code(), seams_sys::SEAMS_ERR_ARG);
/// ```
pub fn nframes(path: &str) -> Result<usize, Error> {
    let text = c_text(path)?;
    let mut count = 0usize;
    check(
        unsafe { ffi::seams_path_nframes(text.as_ptr(), &mut count) },
        std::ptr::null(),
    )?;
    Ok(count)
}

/// Short name of a CHILL+ label.
///
/// ```
/// assert_eq!(seams::label_name(seams::Label::Water), "water");
/// ```
#[must_use]
pub fn label_name(label: Label) -> String {
    static_text(unsafe { ffi::seams_label_name(label as c_int) })
}

/// Short name of a ring-symmetry code.
///
/// ```
/// assert_eq!(seams::ring_symmetry_name(seams::RingSymmetry::Pentagon), "pentagon");
/// ```
#[must_use]
pub fn ring_symmetry_name(symmetry: RingSymmetry) -> String {
    static_text(unsafe { ffi::seams_ring_symmetry_name(symmetry as c_int) })
}

/// Knobs for CHILL and CHILL+.
///
/// ```
/// let options = seams::ChillOptions::default();
/// assert_eq!(options.cutoff, 3.5);
/// assert_eq!(options.graph, seams::Graph::Cutoff);
/// assert_eq!(options.rule, seams::ChillRule::Published);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChillOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Species filter. Zero keeps every atom.
    pub species: i32,
    /// Neighbour graph.
    pub graph: Graph,
    /// Threshold table.
    pub rule: ChillRule,
}

impl Default for ChillOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            species: 0,
            graph: Graph::Cutoff,
            rule: ChillRule::Published,
        }
    }
}

/// Knobs for a local Steinhardt parameter on the cutoff shell.
///
/// ```
/// let options = seams::SteinhardtOptions { order: 4, average: true, ..Default::default() };
/// assert_eq!(options.order, 4);
/// assert!(options.average);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SteinhardtOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Species filter. Zero keeps every atom.
    pub species: i32,
    /// Degree. 3, 4, 6, 8, or 12.
    pub order: i32,
    /// Also return the Lechner-Dellago mean.
    pub average: bool,
}

impl Default for SteinhardtOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            species: 0,
            order: 6,
            average: false,
        }
    }
}

/// Knobs for q_l on a named graph.
///
/// ```
/// let options = seams::OrderOptions::default();
/// assert_eq!(options.graph, seams::Graph::Cutoff);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrderOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Species filter.
    pub species: i32,
    /// Degree. 3, 4, 6, 8, or 12.
    pub order: i32,
    /// Cutoff or mutual.
    pub graph: Graph,
    /// Also return the neighbour mean.
    pub average: bool,
}

impl Default for OrderOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            species: 0,
            order: 6,
            graph: Graph::Cutoff,
            average: false,
        }
    }
}

/// Knobs for the q6-bar ice summary.
///
/// ```
/// assert_eq!(seams::Q6BarOptions::default().ice_cut, 0.36);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Q6BarOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Species filter.
    pub species: i32,
    /// Finite threshold. The nucleation census uses 0.36.
    pub ice_cut: f64,
}

impl Default for Q6BarOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            species: 0,
            ice_cut: 0.36,
        }
    }
}

/// Knobs for a neighbour list.
///
/// ```
/// assert_eq!(seams::NeighbourOptions::default().graph, seams::Graph::Cutoff);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NeighbourOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Species filter.
    pub species: i32,
    /// Cutoff or mutual.
    pub graph: Graph,
}

impl Default for NeighbourOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            species: 0,
            graph: Graph::Cutoff,
        }
    }
}

/// Knobs for a k-nearest neighbour list.
///
/// ```
/// let options = seams::KnnOptions { k: 4, mutual: true, ..Default::default() };
/// assert_eq!(options.k, 4);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KnnOptions {
    /// Positive neighbour count.
    pub k: i32,
    /// Nominal cutoff, angstroms. Candidates lie inside cutoff + 1.5.
    pub cutoff: f64,
    /// Species filter.
    pub species: i32,
    /// Keep only pairs that name each other.
    pub mutual: bool,
}

impl Default for KnnOptions {
    fn default() -> Self {
        Self {
            k: 4,
            cutoff: 3.5,
            species: 0,
            mutual: false,
        }
    }
}

/// Knobs for primitive rings.
///
/// ```
/// assert_eq!(seams::RingOptions::default().max_size, 6);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RingOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Species filter.
    pub species: i32,
    /// Cutoff or mutual.
    pub graph: Graph,
    /// Largest ring size, from 3 through 12.
    pub max_size: i32,
}

impl Default for RingOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            species: 0,
            graph: Graph::Cutoff,
            max_size: 6,
        }
    }
}

/// Knobs for hexagonal and double-diamond cages.
///
/// ```
/// assert_eq!(seams::CageOptions::default().graph, seams::Graph::Seeded);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CageOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Species filter.
    pub species: i32,
    /// Cage graph. Seeded is the usual choice.
    pub graph: Graph,
    /// Neighbour count used to seed the graph.
    pub k: i32,
    /// Complete partial cages.
    pub complete: bool,
}

impl Default for CageOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            species: 0,
            graph: Graph::Seeded,
            k: 4,
            complete: false,
        }
    }
}

/// Knobs for a signature cage census.
///
/// ```
/// let options = seams::SignatureOptions { spec: "512".into(), ..Default::default() };
/// assert_eq!(options.spec, "512");
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct SignatureOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Species filter.
    pub species: i32,
    /// Neighbour graph.
    pub graph: Graph,
    /// Ring census text, for example `512`.
    pub spec: String,
}

impl Default for SignatureOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            species: 0,
            graph: Graph::Cutoff,
            spec: "512".to_string(),
        }
    }
}

/// Knobs for a partial RDF, its running coordination, or the integral.
///
/// ```
/// let options = seams::DistributionOptions { rmax: 6.0, bins: 8, ..Default::default() };
/// assert_eq!(options.bins, 8);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistributionOptions {
    /// Outer limit, angstroms.
    pub rmax: f64,
    /// Bin count. Zero selects about `rmax / 0.1`.
    pub bins: i32,
    /// First species. Zero keeps every atom.
    pub type_i: i32,
    /// Second species. Zero keeps every atom.
    pub type_j: i32,
}

impl Default for DistributionOptions {
    fn default() -> Self {
        Self {
            rmax: 6.0,
            bins: 0,
            type_i: 0,
            type_j: 0,
        }
    }
}

/// Knobs for a hydrogen-bond count.
///
/// ```
/// let options = seams::BondOptions::default();
/// assert_eq!(options.distance, 2.42);
/// assert_eq!(options.angle, 30.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondOptions {
    /// Heavy-atom cutoff, angstroms.
    pub cutoff: f64,
    /// Heavy-atom species. Zero keeps every non-hydrogen.
    pub oxygen: i32,
    /// Hydrogen species. Must be positive.
    pub hydrogen: i32,
    /// Acceptor-H distance, angstroms.
    pub distance: f64,
    /// O-O-H angle, degrees.
    pub angle: f64,
}

impl Default for BondOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            oxygen: 1,
            hydrogen: 2,
            distance: 2.42,
            angle: 30.0,
        }
    }
}

/// Knobs for the axial number density.
///
/// ```
/// assert_eq!(seams::DensityOptions::default().axis, seams::Axis::Z);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DensityOptions {
    /// Species filter. Zero keeps every atom.
    pub species: i32,
    /// Bin count. Zero selects about one bin per 0.1 angstrom.
    pub bins: i32,
    /// Cell axis.
    pub axis: Axis,
}

impl Default for DensityOptions {
    fn default() -> Self {
        Self {
            species: 0,
            bins: 0,
            axis: Axis::Z,
        }
    }
}

/// Knobs for a topological fingerprint.
///
/// ```
/// assert_eq!(seams::FingerprintOptions::default().hops, 2);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FingerprintOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Species filter.
    pub species: i32,
    /// Neighbour count on the bonded graph.
    pub k: i32,
    /// Neighbour graph.
    pub graph: Graph,
    /// Neighbourhood depth. Zero selects 2.
    pub hops: i32,
    /// Colour vertices by LAMMPS type.
    pub colour_types: bool,
}

impl Default for FingerprintOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            species: 0,
            k: 4,
            graph: Graph::Cutoff,
            hops: 2,
            colour_types: false,
        }
    }
}

/// Knobs for a domain search.
///
/// ```
/// let options = seams::DomainOptions { spec: "1=waterO".into(), ..Default::default() };
/// assert_eq!(options.subset, seams::DomainSubset::Apolar);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct DomainOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Site-table text, for example `1=waterO`.
    pub spec: String,
    /// Which subset to measure.
    pub subset: DomainSubset,
}

impl Default for DomainOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            spec: String::new(),
            subset: DomainSubset::Apolar,
        }
    }
}

/// Knobs for prism stacks.
///
/// ```
/// assert_eq!(seams::PrismOptions::default().cutoff, 3.5);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrismOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Species filter.
    pub species: i32,
}

impl Default for PrismOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            species: 0,
        }
    }
}

/// Knobs for a monolayer ring histogram.
///
/// ```
/// let options = seams::SheetOptions { sheet_area: 100.0, ..Default::default() };
/// assert_eq!(options.max_size, 6);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SheetOptions {
    /// Positive cutoff, angstroms.
    pub cutoff: f64,
    /// Species filter.
    pub species: i32,
    /// Largest ring size.
    pub max_size: i32,
    /// Sheet area, square angstroms. Must be positive.
    pub sheet_area: f64,
}

impl Default for SheetOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            species: 0,
            max_size: 6,
            sheet_area: 0.0,
        }
    }
}

/// Knobs for a planar ice front.
///
/// ```
/// assert_eq!(seams::BrineOptions::default().axis, seams::Axis::Z);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BrineOptions {
    /// Species kept as water, or 0.
    pub water_type: i32,
    /// Species kept as ions, or 0 for none.
    pub ion_type: i32,
    /// Bins along the axis. Zero selects 20.
    pub bins: i32,
    /// Front axis.
    pub axis: Axis,
}

impl Default for BrineOptions {
    fn default() -> Self {
        Self {
            water_type: 0,
            ion_type: 0,
            bins: 0,
            axis: Axis::Z,
        }
    }
}

/// CHILL or CHILL+ labels and the nine population counts.
///
/// `counts[k]` is the number of atoms whose label value is `k`.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let ice = system.chill_plus(&seams::ChillOptions::default())?;
/// assert_eq!(ice.counts[0], 8);
/// assert_eq!(ice.labels.len(), 8);
/// # Ok::<(), seams::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Classification {
    /// One label code per atom.
    pub labels: Vec<i32>,
    /// Population in label order.
    pub counts: [i32; 9],
}

/// Local Steinhardt parameter, and the average when it was requested.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let q6 = system.steinhardt(&seams::SteinhardtOptions { average: true, ..Default::default() })?;
/// assert_eq!(q6.local.len(), 8);
/// assert_eq!(q6.average.as_ref().map(Vec::len), Some(8));
/// # Ok::<(), seams::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Steinhardt {
    /// One local q_l per atom. Dimensionless.
    pub local: Vec<f64>,
    /// Neighbour mean of q_lm. `None` when the average was not requested.
    pub average: Option<Vec<f64>>,
}

/// Local q4 and q6 on the cutoff shell.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let orders = system.q4q6(3.5)?;
/// assert_eq!(orders.q4.len(), orders.q6.len());
/// # Ok::<(), seams::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Q4Q6 {
    /// One local q4 per atom.
    pub q4: Vec<f64>,
    /// One local q6 per atom.
    pub q6: Vec<f64>,
}

/// Lechner-Dellago q6-bar and the ice count above a cut.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let summary = system.q6bar(&seams::Q6BarOptions::default())?;
/// assert_eq!(summary.values.len(), 8);
/// assert!(summary.mean.is_finite());
/// # Ok::<(), seams::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Q6Bar {
    /// One q6-bar per atom.
    pub values: Vec<f64>,
    /// Mean of the finite values.
    pub mean: f64,
    /// Median of the finite values.
    pub median: f64,
    /// Count strictly above the cut.
    pub n_ice: i32,
}

/// Partial radial distribution.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let g = system.rdf(&seams::DistributionOptions { rmax: 6.0, bins: 4, ..Default::default() })?;
/// assert!(!g.radius.is_empty());
/// assert!(g.volume > 0.0);
/// # Ok::<(), seams::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Distribution {
    /// Bin centres, angstroms.
    pub radius: Vec<f64>,
    /// g(r), dimensionless.
    pub g: Vec<f64>,
    /// Pair count in each bin.
    pub count: Vec<i32>,
    /// Cell volume, cubic angstroms.
    pub volume: f64,
}

/// Running coordination number along a partial RDF.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let series = system.running_cn(&seams::DistributionOptions { rmax: 6.0, bins: 4, ..Default::default() })?;
/// assert_eq!(series.radius.len(), series.coordination.len());
/// # Ok::<(), seams::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct RunningCoordination {
    /// Bin centres, angstroms.
    pub radius: Vec<f64>,
    /// Running coordination, dimensionless.
    pub coordination: Vec<f64>,
}

/// Disjoint cage populations.
///
/// ```
/// let counts = seams::CageCounts { hexagonal: 1, cubic: 0, water: 2, both: 0 };
/// assert_eq!(counts.hexagonal, 1);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CageCounts {
    /// Hexagonal cage only.
    pub hexagonal: i32,
    /// Double-diamond cage only.
    pub cubic: i32,
    /// Neither cage.
    pub water: i32,
    /// Both cages.
    pub both: i32,
}

/// Closed cages whose ring census matches a signature.
///
/// ```
/// let found = seams::SignatureMatch { cages: 2, atoms: 40 };
/// assert_eq!(found.cages, 2);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SignatureMatch {
    /// Number of closed cages.
    pub cages: usize,
    /// Distinct atoms in those cages.
    pub atoms: usize,
}

/// Mutual nearest unlike contacts of a site table.
///
/// ```
/// let contacts = seams::Contacts { pairs: 0, cations: 1, anions: 1 };
/// assert_eq!(contacts.pairs, 0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Contacts {
    /// Contact count.
    pub pairs: usize,
    /// Cation molecules in the ion cloud.
    pub cations: i32,
    /// Anion molecules in the ion cloud.
    pub anions: i32,
}

/// Largest polar or apolar domain.
///
/// ```
/// let domain = seams::DomainReport { subset: 8, largest: 8, percolation: 1.0 };
/// assert_eq!(domain.largest, 8);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DomainReport {
    /// Atoms in the requested subset.
    pub subset: i32,
    /// Largest connected component.
    pub largest: i32,
    /// Percolation measure, dimensionless.
    pub percolation: f64,
}

/// Perfect and imperfect prism counts by ring size.
///
/// Index is the ring size. Sizes 3 through 6 are filled.
///
/// ```
/// let prisms = seams::PrismCensus::default();
/// assert_eq!(prisms.perfect.len(), 7);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PrismCensus {
    /// Perfect stacks.
    pub perfect: [i32; 7],
    /// Imperfect stacks.
    pub imperfect: [i32; 7],
}

/// Number density along one cell axis.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let rho = system.density(&seams::DensityOptions { bins: 2, ..Default::default() })?;
/// assert_eq!(rho.coordinate.len(), rho.density.len());
/// # Ok::<(), seams::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct DensityProfile {
    /// Bin centres, angstroms.
    pub coordinate: Vec<f64>,
    /// Number density.
    pub density: Vec<f64>,
}

/// The two planar crossings of an ice slab.
///
/// ```
/// let fronts = seams::BrineFronts { enter: 1.0, exit: 2.0, thickness: 1.0, found: true };
/// assert!(fronts.found);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrineFronts {
    /// Crossing on entering the ice, angstroms.
    pub enter: f64,
    /// Crossing on leaving the ice, angstroms.
    pub exit: f64,
    /// Ice run along the axis, angstroms.
    pub thickness: f64,
    /// Both crossings exist.
    pub found: bool,
}

/// Molal rejection coefficient.
///
/// ```
/// let rejection = seams::Rejection { ice: 0.0, liquid: 0.0, k: 0.0, defined: false };
/// assert!(!rejection.defined);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rejection {
    /// Bulk ice molality, mol/kg.
    pub ice: f64,
    /// Bulk liquid molality, mol/kg.
    pub liquid: f64,
    /// Ratio of the two molalities.
    pub k: f64,
    /// The ratio is defined.
    pub defined: bool,
}

/// Concentration profile along the brine axis.
///
/// ```
/// let profile = seams::BrineProfile { coordinate: vec![0.0], molality: vec![0.0], ice_fraction: vec![1.0] };
/// assert_eq!(profile.coordinate.len(), 1);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct BrineProfile {
    /// Bin centres, angstroms.
    pub coordinate: Vec<f64>,
    /// mol/kg.
    pub molality: Vec<f64>,
    /// Dimensionless ice fraction.
    pub ice_fraction: Vec<f64>,
}

/// Brine pockets trapped in the ice.
///
/// ```
/// let pockets = seams::Pockets { pockets: 0, reservoir_waters: 8 };
/// assert_eq!(pockets.reservoir_waters, 8);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Pockets {
    /// Pocket count, excluding the reservoir.
    pub pockets: i32,
    /// Waters in the liquid reservoir.
    pub reservoir_waters: i32,
}

/// One directed neighbour pair. The self entry is omitted.
///
/// ```
/// let pair = seams::NeighbourPair { atom: 0, neighbour: 1 };
/// assert_eq!(pair.neighbour, 1);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NeighbourPair {
    /// Atom index in the frame.
    pub atom: i32,
    /// Neighbour index in the frame.
    pub neighbour: i32,
}

/// Iterator of directed neighbour pairs.
///
/// Yields [`NeighbourPair`]. Created by [`Neighbours::pairs`].
pub struct NeighbourPairs<'a> {
    offsets: &'a [i32],
    cols: &'a [i32],
    atom: usize,
    cursor: usize,
}

impl Iterator for NeighbourPairs<'_> {
    type Item = NeighbourPair;

    fn next(&mut self) -> Option<Self::Item> {
        while self.atom + 1 < self.offsets.len() {
            let end = self.offsets[self.atom + 1];
            let end = usize::try_from(end).unwrap_or(0);
            if self.cursor < end {
                let index = self.cursor;
                self.cursor += 1;
                let neighbour = self.cols.get(index).copied().unwrap_or(-1);
                if neighbour != i32::try_from(self.atom).unwrap_or(i32::MAX) {
                    return Some(NeighbourPair {
                        atom: i32::try_from(self.atom).unwrap_or(0),
                        neighbour,
                    });
                }
                continue;
            }
            self.atom += 1;
            if self.atom < self.offsets.len() {
                self.cursor = usize::try_from(self.offsets[self.atom]).unwrap_or(0);
            }
        }
        None
    }
}

/// Owned neighbour list.
///
/// The handle is freed on drop. It is [`Send`] and not [`Sync`].
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let list = system.neighbours(&seams::NeighbourOptions::default())?;
/// assert!(list.pairs().all(|pair| pair.atom != pair.neighbour));
/// assert!(list.degrees()?.iter().all(|degree| *degree == 4));
/// # Ok::<(), seams::Error>(())
/// ```
pub struct Neighbours {
    ptr: *mut ffi::seams_neighbours,
    offsets: Vec<i32>,
    cols: Vec<i32>,
}

// Safety: the C++ object is a heap allocation of vectors and a string.
// It has no thread-local members. Moving the owner to another thread is
// sound. Sharing it is not: the C API requires exclusive use.
unsafe impl Send for Neighbours {}

impl Drop for Neighbours {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { ffi::seams_neighbours_destroy(self.ptr) };
        }
    }
}

impl Neighbours {
    fn adopt(ptr: *mut ffi::seams_neighbours) -> Result<Self, Error> {
        let handle = ptr as *const c_void;
        let mut n_offsets = 0usize;
        let mut nnz = 0usize;
        check(
            unsafe {
                ffi::seams_neighbours_copy(
                    ptr,
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
        let mut offsets = vec![0; n_offsets];
        let mut cols = vec![0; nnz];
        check(
            unsafe {
                ffi::seams_neighbours_copy(
                    ptr,
                    offsets.as_mut_ptr(),
                    offsets.len(),
                    &mut n_offsets,
                    cols.as_mut_ptr(),
                    cols.len(),
                    &mut nnz,
                )
            },
            handle,
        )?;
        offsets.truncate(n_offsets);
        cols.truncate(nnz);
        Ok(Self { ptr, offsets, cols })
    }

    /// Number of neighbour pairs, omitting self.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let list = system.knn(&seams::KnnOptions { k: 4, mutual: true, ..Default::default() })?;
    /// assert!(list.nnz()? > 0);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn nnz(&self) -> Result<usize, Error> {
        Ok(self.pairs().count())
    }

    /// Degree of each atom. The self entry is not counted.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let list = system.neighbours(&seams::NeighbourOptions::default())?;
    /// assert_eq!(list.degrees()?.len(), 8);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn degrees(&self) -> Result<Vec<i32>, Error> {
        let mut degree = Vec::new();
        if self.offsets.len() >= 2 {
            for window in self.offsets.windows(2) {
                degree.push(window[1] - window[0]);
            }
        }
        Ok(degree)
    }

    /// Directed pairs. Self is omitted.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let list = system.neighbours(&seams::NeighbourOptions::default())?;
    /// assert_eq!(list.pairs().count(), list.nnz()?);
    /// # Ok::<(), seams::Error>(())
    /// ```
    #[must_use]
    pub fn pairs(&self) -> NeighbourPairs<'_> {
        let cursor = self
            .offsets
            .first()
            .copied()
            .and_then(|v| usize::try_from(v).ok())
            .unwrap_or(0);
        NeighbourPairs {
            offsets: &self.offsets,
            cols: &self.cols,
            atom: 0,
            cursor,
        }
    }
}

/// Owned primitive-ring set.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let rings = system.rings(&seams::RingOptions::default())?;
/// assert!(rings.size() > 0);
/// assert!(!rings.vertices(0)?.is_empty());
/// # Ok::<(), seams::Error>(())
/// ```
pub struct Rings {
    ptr: *mut ffi::seams_rings,
}

unsafe impl Send for Rings {}

impl Drop for Rings {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { ffi::seams_rings_destroy(self.ptr) };
        }
    }
}

impl Rings {
    /// Number of primitive rings.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let rings = system.rings(&seams::RingOptions::default())?;
    /// assert_eq!(rings.iter().count(), rings.size() as usize);
    /// # Ok::<(), seams::Error>(())
    /// ```
    #[must_use]
    pub fn size(&self) -> i32 {
        unsafe { ffi::seams_rings_count(self.ptr) }
    }

    /// Atom indices of one ring.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let rings = system.rings(&seams::RingOptions::default())?;
    /// assert!(!rings.vertices(0)?.is_empty());
    /// # Ok::<(), seams::Error>(())
    /// ```
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
                ffi::seams_rings_vertices(self.ptr, index, atoms.as_mut_ptr(), n_out, &mut n_out)
            },
            handle,
        )?;
        atoms.truncate(n_out);
        Ok(atoms)
    }

    /// Histogram of ring size. `counts[k]` is the number of rings of size k.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let rings = system.rings(&seams::RingOptions::default())?;
    /// let hist = rings.histogram(6)?;
    /// assert_eq!(hist.len(), 7);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn histogram(&self, up_to: i32) -> Result<Vec<i32>, Error> {
        let len = usize::try_from(up_to).map_err(|_| Error::argument("ring size is negative"))? + 1;
        let mut counts = vec![0; len];
        check(
            unsafe { ffi::seams_rings_histogram(self.ptr, counts.as_mut_ptr(), up_to + 1) },
            self.ptr as *const c_void,
        )?;
        Ok(counts)
    }

    /// Fraction of atoms that sit on a ring of one size. Zero counts every ring.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let rings = system.rings(&seams::RingOptions::default())?;
    /// assert!(rings.fraction(6)?.is_finite());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn fraction(&self, ring_size: i32) -> Result<f64, Error> {
        let mut value = 0.0;
        check(
            unsafe { ffi::seams_rings_fraction(self.ptr, ring_size, &mut value) },
            self.ptr as *const c_void,
        )?;
        Ok(value)
    }

    /// Walk the rings. Each item is the vertex list of one ring.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let rings = system.rings(&seams::RingOptions::default())?;
    /// assert!(rings.iter().all(|ring| ring.is_ok()));
    /// # Ok::<(), seams::Error>(())
    /// ```
    #[must_use]
    pub fn iter(&self) -> RingIter<'_> {
        RingIter {
            rings: self,
            index: 0,
        }
    }
}

/// Iterator of primitive rings. Each item is a [`Result`] of vertex indices.
pub struct RingIter<'a> {
    rings: &'a Rings,
    index: i32,
}

impl Iterator for RingIter<'_> {
    type Item = Result<Vec<i32>, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.rings.size() {
            return None;
        }
        let index = self.index;
        self.index += 1;
        Some(self.rings.vertices(index))
    }
}

/// Owned hexagonal and double-diamond membership.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let cages = system.cages(&seams::CageOptions::default())?;
/// let counts = cages.counts()?;
/// assert_eq!(counts.hexagonal + counts.cubic + counts.water + counts.both, 8);
/// # Ok::<(), seams::Error>(())
/// ```
pub struct Cages {
    ptr: *mut ffi::seams_cage_set,
}

unsafe impl Send for Cages {}

impl Drop for Cages {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { ffi::seams_cages_destroy(self.ptr) };
        }
    }
}

impl Cages {
    /// Disjoint population counts.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let cages = system.cages(&seams::CageOptions::default())?;
    /// assert!(cages.counts()?.water >= 0);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn counts(&self) -> Result<CageCounts, Error> {
        let mut out = CageCounts::default();
        check(
            unsafe {
                ffi::seams_cages_counts(
                    self.ptr,
                    &mut out.hexagonal,
                    &mut out.cubic,
                    &mut out.water,
                    &mut out.both,
                )
            },
            self.ptr as *const c_void,
        )?;
        Ok(out)
    }

    /// One cage code per atom.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let cages = system.cages(&seams::CageOptions::default())?;
    /// assert_eq!(cages.labels(8)?.len(), 8);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn labels(&self, n_atoms: usize) -> Result<Vec<i32>, Error> {
        let mut per_atom = vec![0; n_atoms];
        check(
            unsafe { ffi::seams_cages_labels(self.ptr, per_atom.as_mut_ptr(), n_atoms) },
            self.ptr as *const c_void,
        )?;
        Ok(per_atom)
    }
}

/// Owned topological fingerprint.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let fp = system.fingerprint(&seams::FingerprintOptions::default())?;
/// assert!(!fp.key()?.is_empty());
/// assert_eq!(fp.rings()?.len(), 13);
/// # Ok::<(), seams::Error>(())
/// ```
pub struct Fingerprint {
    ptr: *mut ffi::seams_topo,
}

unsafe impl Send for Fingerprint {}

impl Drop for Fingerprint {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { ffi::seams_fingerprint_destroy(self.ptr) };
        }
    }
}

impl Fingerprint {
    /// Frame fingerprint text.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let fp = system.fingerprint(&seams::FingerprintOptions::default())?;
    /// assert!(fp.key()?.contains('-') || !fp.key()?.is_empty());
    /// # Ok::<(), seams::Error>(())
    /// ```
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
                ffi::seams_fingerprint_key(self.ptr, buf.as_mut_ptr().cast(), buf.len(), &mut n_out)
            },
            handle,
        )?;
        buf.truncate(n_out);
        Ok(String::from_utf8_lossy(&buf).into_owned())
    }

    /// Local key of one atom.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let fp = system.fingerprint(&seams::FingerprintOptions::default())?;
    /// assert!(!fp.atom(0)?.is_empty());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn atom(&self, index: i32) -> Result<String, Error> {
        let handle = self.ptr as *const c_void;
        let mut n_out = 0usize;
        check(
            unsafe {
                ffi::seams_fingerprint_atom(self.ptr, index, std::ptr::null_mut(), 0, &mut n_out)
            },
            handle,
        )?;
        let mut buf = vec![0u8; n_out + 1];
        check(
            unsafe {
                ffi::seams_fingerprint_atom(
                    self.ptr,
                    index,
                    buf.as_mut_ptr().cast(),
                    buf.len(),
                    &mut n_out,
                )
            },
            handle,
        )?;
        buf.truncate(n_out);
        Ok(String::from_utf8_lossy(&buf).into_owned())
    }

    /// Primitive-ring census stored on the fingerprint.
    ///
    /// `counts[k]` is the number of rings of size k.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let fp = system.fingerprint(&seams::FingerprintOptions::default())?;
    /// assert!(fp.rings()?.iter().sum::<i32>() > 0);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn rings(&self) -> Result<Vec<i32>, Error> {
        let mut counts = vec![0; 13];
        check(
            unsafe { ffi::seams_fingerprint_rings(self.ptr, counts.as_mut_ptr(), counts.len()) },
            self.ptr as *const c_void,
        )?;
        Ok(counts)
    }
}

/// Owned planar ice-front measurement.
///
/// ```
/// let (xyz, edge) = seams::diamond8();
/// let mut system = seams::System::new()?;
/// system.set_ortho(&xyz, edge, edge, edge)?;
/// let ice = vec![1; 8];
/// let front = system.brine(&ice, &seams::BrineOptions::default())?;
/// let pockets = front.pockets()?;
/// assert!(pockets.pockets >= 0);
/// # Ok::<(), seams::Error>(())
/// ```
pub struct Brine {
    ptr: *mut ffi::seams_brine,
}

unsafe impl Send for Brine {}

impl Drop for Brine {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { ffi::seams_brine_destroy(self.ptr) };
        }
    }
}

impl Brine {
    /// The two planar crossings. Lengths are angstroms.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let front = system.brine(&[1; 8], &seams::BrineOptions::default())?;
    /// let crossings = front.fronts()?;
    /// assert!(crossings.thickness.is_finite());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn fronts(&self) -> Result<BrineFronts, Error> {
        let mut out = BrineFronts {
            enter: 0.0,
            exit: 0.0,
            thickness: 0.0,
            found: false,
        };
        let mut found = 0;
        check(
            unsafe {
                ffi::seams_brine_fronts(
                    self.ptr,
                    &mut out.enter,
                    &mut out.exit,
                    &mut out.thickness,
                    &mut found,
                )
            },
            self.ptr as *const c_void,
        )?;
        out.found = found != 0;
        Ok(out)
    }

    /// Molal rejection coefficient. Molality is mol/kg.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let front = system.brine(&[1; 8], &seams::BrineOptions::default())?;
    /// let rejection = front.rejection()?;
    /// assert!(rejection.k.is_finite() || !rejection.defined);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn rejection(&self) -> Result<Rejection, Error> {
        let mut out = Rejection {
            ice: 0.0,
            liquid: 0.0,
            k: 0.0,
            defined: false,
        };
        let mut defined = 0;
        check(
            unsafe {
                ffi::seams_brine_rejection(
                    self.ptr,
                    &mut out.ice,
                    &mut out.liquid,
                    &mut out.k,
                    &mut defined,
                )
            },
            self.ptr as *const c_void,
        )?;
        out.defined = defined != 0;
        Ok(out)
    }

    /// Concentration profile along the front axis.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let front = system.brine(&[1; 8], &seams::BrineOptions::default())?;
    /// let profile = front.profile()?;
    /// assert_eq!(profile.coordinate.len(), profile.molality.len());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn profile(&self) -> Result<BrineProfile, Error> {
        let handle = self.ptr as *const c_void;
        let mut n_out = 0usize;
        check(
            unsafe {
                ffi::seams_brine_profile(
                    self.ptr,
                    std::ptr::null_mut(),
                    std::ptr::null_mut(),
                    std::ptr::null_mut(),
                    0,
                    &mut n_out,
                )
            },
            handle,
        )?;
        let mut coordinate = vec![0.0; n_out];
        let mut molality = vec![0.0; n_out];
        let mut ice_fraction = vec![0.0; n_out];
        check(
            unsafe {
                ffi::seams_brine_profile(
                    self.ptr,
                    coordinate.as_mut_ptr(),
                    molality.as_mut_ptr(),
                    ice_fraction.as_mut_ptr(),
                    n_out,
                    &mut n_out,
                )
            },
            handle,
        )?;
        coordinate.truncate(n_out);
        molality.truncate(n_out);
        ice_fraction.truncate(n_out);
        Ok(BrineProfile {
            coordinate,
            molality,
            ice_fraction,
        })
    }

    /// Pockets trapped in the ice, and the liquid reservoir.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let front = system.brine(&[1; 8], &seams::BrineOptions::default())?;
    /// assert!(front.pockets()?.reservoir_waters >= 0);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn pockets(&self) -> Result<Pockets, Error> {
        let mut out = Pockets::default();
        check(
            unsafe {
                ffi::seams_brine_pockets(self.ptr, &mut out.pockets, &mut out.reservoir_waters)
            },
            self.ptr as *const c_void,
        )?;
        Ok(out)
    }
}

/// One frame: coordinates, a cell, and the analyses that read them.
///
/// The handle is freed on drop. It is [`Send`] and not [`Sync`].
///
/// ```
/// let mut system = seams::System::new()?;
/// assert_eq!(system.count(), 0);
/// let err = system.read("missing-trajectory.lammpstrj", 1, -1).expect_err("no such file");
/// assert!(!err.message().is_empty());
/// # Ok::<(), seams::Error>(())
/// ```
pub struct System {
    ptr: *mut ffi::seams_system,
}

unsafe impl Send for System {}

impl Drop for System {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { ffi::seams_system_destroy(self.ptr) };
        }
    }
}

impl System {
    /// Allocate an empty frame.
    ///
    /// ```
    /// let system = seams::System::new()?;
    /// assert_eq!(system.count(), 0);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn new() -> Result<Self, Error> {
        let ptr = unsafe { ffi::seams_system_create() };
        if ptr.is_null() {
            return Err(Error::from_status(
                ffi::SEAMS_ERR_NOMEM,
                message_of(std::ptr::null()),
            ));
        }
        Ok(Self { ptr })
    }

    fn handle(&self) -> *const c_void {
        self.ptr as *const c_void
    }

    /// Store points in an orthorhombic cell whose origin is zero.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// assert_eq!(system.count(), 8);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn set_ortho(&mut self, xyz: &[f64], lx: f64, ly: f64, lz: f64) -> Result<(), Error> {
        let n = atom_count(xyz.len())?;
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

    /// Store points and a general nine-number cell.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut box9 = [0.0; 9];
    /// box9[3] = edge;
    /// box9[4] = edge;
    /// box9[5] = edge;
    /// box9[6] = 0.3;
    /// let mut system = seams::System::new()?;
    /// system.set_frame(&xyz, &box9)?;
    /// assert_eq!(system.count(), 8);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn set_frame(&mut self, xyz: &[f64], box9: &[f64; 9]) -> Result<(), Error> {
        let n = atom_count(xyz.len())?;
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

    /// Store points, a cell, and optional species and ids.
    ///
    /// `box9` is `ox, oy, oz, lx, ly, lz, xy, xz, yz`. Empty `types` and
    /// `ids` leave those columns unset.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let box9 = [0.0, 0.0, 0.0, edge, edge, edge, 0.0, 0.0, 0.0];
    /// let types = [1; 8];
    /// let mut system = seams::System::new()?;
    /// system.set_atoms(&xyz, &box9, &types, &[])?;
    /// assert_eq!(system.count(), 8);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn set_atoms(
        &mut self,
        xyz: &[f64],
        box9: &[f64],
        types: &[i32],
        ids: &[i32],
    ) -> Result<(), Error> {
        if xyz.len() % 3 != 0 {
            return Err(Error::argument("coordinates are not a multiple of 3"));
        }
        let n = xyz.len() / 3;
        if !types.is_empty() && types.len() != n {
            return Err(Error::argument(
                "types length does not match the atom count",
            ));
        }
        if !ids.is_empty() && ids.len() != n {
            return Err(Error::argument("ids length does not match the atom count"));
        }
        check(
            unsafe {
                ffi::seams_system_set_atoms(
                    self.ptr,
                    xyz.as_ptr(),
                    n,
                    box9.as_ptr(),
                    box9.len(),
                    if types.is_empty() {
                        std::ptr::null()
                    } else {
                        types.as_ptr()
                    },
                    if ids.is_empty() {
                        std::ptr::null()
                    } else {
                        ids.as_ptr()
                    },
                )
            },
            self.handle(),
        )
    }

    /// Read one frame from a trajectory path.
    ///
    /// `frame` is 1-based. A negative `type_filter` keeps every atom.
    /// The path may be XYZ, CON, a chemfiles format, or a LAMMPS dump.
    ///
    /// ```
    /// let mut system = seams::System::new()?;
    /// let err = system.read("no-such-frame.xyz", 1, -1).unwrap_err();
    /// assert!(!err.message().is_empty());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn read(&mut self, path: &str, frame: i32, type_filter: i32) -> Result<(), Error> {
        let text = c_text(path)?;
        check(
            unsafe { ffi::seams_system_read(self.ptr, text.as_ptr(), frame, type_filter) },
            self.handle(),
        )
    }

    /// Number of points in the frame.
    ///
    /// ```
    /// let system = seams::System::new()?;
    /// assert_eq!(system.count(), 0);
    /// # Ok::<(), seams::Error>(())
    /// ```
    #[must_use]
    pub fn count(&self) -> i32 {
        unsafe { ffi::seams_system_count(self.ptr) }
    }

    /// Original CHILL labels and nine population counts.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let ice = system.chill(&seams::ChillOptions::default())?;
    /// assert_eq!(ice.labels.len(), 8);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn chill(&mut self, options: &ChillOptions) -> Result<Classification, Error> {
        self.classify(false, options)
    }

    /// CHILL+ labels and nine population counts.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let ice = system.chill_plus(&seams::ChillOptions {
    ///     graph: seams::Graph::Mutual,
    ///     rule: seams::ChillRule::Ovito,
    ///     ..seams::ChillOptions::default()
    /// })?;
    /// assert_eq!(ice.counts.iter().sum::<i32>() as usize, ice.labels.len());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn chill_plus(&mut self, options: &ChillOptions) -> Result<Classification, Error> {
        self.classify(true, options)
    }

    /// Local q_l on the cutoff shell.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let q3 = system.steinhardt(&seams::SteinhardtOptions { order: 3, ..Default::default() })?;
    /// assert!(q3.average.is_none());
    /// assert_eq!(q3.local.len(), 8);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn steinhardt(&mut self, options: &SteinhardtOptions) -> Result<Steinhardt, Error> {
        let n = self.count() as usize;
        let mut local = vec![0.0; n];
        let mut average = vec![0.0; n];
        let bar = if options.average {
            average.as_mut_ptr()
        } else {
            std::ptr::null_mut()
        };
        check(
            unsafe {
                ffi::seams_steinhardt(
                    self.ptr,
                    options.cutoff,
                    options.species,
                    options.order,
                    local.as_mut_ptr(),
                    n,
                    bar,
                    n,
                )
            },
            self.handle(),
        )?;
        Ok(Steinhardt {
            local,
            average: if options.average { Some(average) } else { None },
        })
    }

    /// Local q_l on a named graph.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let q4 = system.ql(&seams::OrderOptions { order: 4, ..Default::default() })?;
    /// assert_eq!(q4.local.len(), 8);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn ql(&mut self, options: &OrderOptions) -> Result<Steinhardt, Error> {
        let n = self.count() as usize;
        let mut local = vec![0.0; n];
        let mut average = vec![0.0; n];
        let bar = if options.average {
            average.as_mut_ptr()
        } else {
            std::ptr::null_mut()
        };
        check(
            unsafe {
                ffi::seams_ql(
                    self.ptr,
                    options.cutoff,
                    options.species,
                    options.order,
                    options.graph as c_int,
                    local.as_mut_ptr(),
                    n,
                    bar,
                    n,
                )
            },
            self.handle(),
        )?;
        Ok(Steinhardt {
            local,
            average: if options.average { Some(average) } else { None },
        })
    }

    /// Local q4 and q6 on the cutoff shell.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let orders = system.q4q6(3.5)?;
    /// let mean = orders.q6.iter().sum::<f64>() / orders.q6.len() as f64;
    /// assert!(mean > 0.3);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn q4q6(&mut self, cutoff: f64) -> Result<Q4Q6, Error> {
        let n = self.count() as usize;
        let mut q4 = vec![0.0; n];
        let mut q6 = vec![0.0; n];
        check(
            unsafe { ffi::seams_q4q6(self.ptr, cutoff, 0, q4.as_mut_ptr(), q6.as_mut_ptr()) },
            self.handle(),
        )?;
        Ok(Q4Q6 { q4, q6 })
    }

    /// q6-bar, its mean and median, and the ice count above the cut.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let summary = system.q6bar(&seams::Q6BarOptions::default())?;
    /// assert!(summary.n_ice >= 0);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn q6bar(&mut self, options: &Q6BarOptions) -> Result<Q6Bar, Error> {
        let n = self.count() as usize;
        let mut values = vec![0.0; n];
        let mut mean = 0.0;
        let mut median = 0.0;
        let mut n_ice = 0;
        check(
            unsafe {
                ffi::seams_q6bar(
                    self.ptr,
                    options.cutoff,
                    options.species,
                    options.ice_cut,
                    values.as_mut_ptr(),
                    n,
                    &mut mean,
                    &mut median,
                    &mut n_ice,
                )
            },
            self.handle(),
        )?;
        Ok(Q6Bar {
            values,
            mean,
            median,
            n_ice,
        })
    }

    /// Cutoff or mutual neighbour list.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let list = system.neighbours(&seams::NeighbourOptions::default())?;
    /// assert_eq!(list.degrees()?.iter().sum::<i32>(), list.nnz()? as i32);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn neighbours(&mut self, options: &NeighbourOptions) -> Result<Neighbours, Error> {
        let mut out = std::ptr::null_mut();
        check(
            unsafe {
                ffi::seams_neighbours_build(
                    self.ptr,
                    options.cutoff,
                    options.species,
                    options.graph as c_int,
                    &mut out,
                )
            },
            self.handle(),
        )?;
        Neighbours::adopt(out)
    }

    /// k-nearest neighbour list.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let list = system.knn(&seams::KnnOptions { mutual: true, ..Default::default() })?;
    /// assert!(list.nnz()? > 0);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn knn(&mut self, options: &KnnOptions) -> Result<Neighbours, Error> {
        let mut out = std::ptr::null_mut();
        check(
            unsafe {
                ffi::seams_neighbours_knn(
                    self.ptr,
                    options.k,
                    options.cutoff,
                    options.species,
                    i32::from(options.mutual),
                    &mut out,
                )
            },
            self.handle(),
        )?;
        Neighbours::adopt(out)
    }

    /// Degree of each atom on the cutoff graph.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let degree = system.neighbour_degree(3.5)?;
    /// assert!(degree.iter().all(|d| *d == 4));
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn neighbour_degree(&mut self, cutoff: f64) -> Result<Vec<i32>, Error> {
        self.neighbours(&NeighbourOptions {
            cutoff,
            ..NeighbourOptions::default()
        })?
        .degrees()
    }

    /// Primitive rings up to a maximum size.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let rings = system.rings(&seams::RingOptions::default())?;
    /// assert!(rings.fraction(0)?.is_finite());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn rings(&mut self, options: &RingOptions) -> Result<Rings, Error> {
        let mut out = std::ptr::null_mut();
        check(
            unsafe {
                ffi::seams_rings_find(
                    self.ptr,
                    options.cutoff,
                    options.species,
                    options.graph as c_int,
                    options.max_size,
                    &mut out,
                )
            },
            self.handle(),
        )?;
        Ok(Rings { ptr: out })
    }

    /// Count of six-rings and the fraction of atoms that sit on one.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let (six, fraction) = system.six_ring_fraction(3.5)?;
    /// assert!(six >= 0);
    /// assert!(fraction.is_finite());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn six_ring_fraction(&mut self, cutoff: f64) -> Result<(i32, f64), Error> {
        let rings = self.rings(&RingOptions {
            cutoff,
            max_size: 6,
            ..RingOptions::default()
        })?;
        let hist = rings.histogram(6)?;
        let fraction = rings.fraction(6)?;
        Ok((hist[6], fraction))
    }

    /// Hexagonal and double-diamond cage membership.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let cages = system.cages(&seams::CageOptions::default())?;
    /// let _ = cages.counts()?;
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn cages(&mut self, options: &CageOptions) -> Result<Cages, Error> {
        let mut out = std::ptr::null_mut();
        check(
            unsafe {
                ffi::seams_cages(
                    self.ptr,
                    options.cutoff,
                    options.species,
                    options.graph as c_int,
                    options.k,
                    i32::from(options.complete),
                    &mut out,
                )
            },
            self.handle(),
        )?;
        Ok(Cages { ptr: out })
    }

    /// Closed cages of one ring census.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let err = system.signature_cages(&seams::SignatureOptions::default()).unwrap_err();
    /// assert_eq!(err.status(), Some(seams::Status::Argument));
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn signature_cages(&mut self, options: &SignatureOptions) -> Result<SignatureMatch, Error> {
        let text = c_text(&options.spec)?;
        let mut out = SignatureMatch::default();
        check(
            unsafe {
                ffi::seams_signature_cages(
                    self.ptr,
                    options.cutoff,
                    options.species,
                    options.graph as c_int,
                    text.as_ptr(),
                    &mut out.cages,
                    &mut out.atoms,
                )
            },
            self.handle(),
        )?;
        Ok(out)
    }

    /// Partial radial distribution and the cell volume.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let g = system.rdf(&seams::DistributionOptions { rmax: 6.0, bins: 4, ..Default::default() })?;
    /// assert_eq!(g.radius.len(), g.g.len());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn rdf(&mut self, options: &DistributionOptions) -> Result<Distribution, Error> {
        let mut n_out = 0usize;
        let mut volume = 0.0;
        check(
            unsafe {
                ffi::seams_rdf(
                    self.ptr,
                    options.rmax,
                    options.bins,
                    options.type_i,
                    options.type_j,
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
        let mut radius = vec![0.0; n_out];
        let mut g = vec![0.0; n_out];
        let mut count = vec![0; n_out];
        check(
            unsafe {
                ffi::seams_rdf(
                    self.ptr,
                    options.rmax,
                    options.bins,
                    options.type_i,
                    options.type_j,
                    radius.as_mut_ptr(),
                    g.as_mut_ptr(),
                    count.as_mut_ptr(),
                    n_out,
                    &mut n_out,
                    std::ptr::null_mut(),
                )
            },
            self.handle(),
        )?;
        radius.truncate(n_out);
        g.truncate(n_out);
        count.truncate(n_out);
        Ok(Distribution {
            radius,
            g,
            count,
            volume,
        })
    }

    /// Running coordination number along the partial RDF.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let series = system.running_cn(&seams::DistributionOptions { rmax: 6.0, bins: 4, ..Default::default() })?;
    /// assert!(!series.coordination.is_empty());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn running_cn(
        &mut self,
        options: &DistributionOptions,
    ) -> Result<RunningCoordination, Error> {
        let mut n_out = 0usize;
        check(
            unsafe {
                ffi::seams_running_cn(
                    self.ptr,
                    options.rmax,
                    options.bins,
                    options.type_i,
                    options.type_j,
                    std::ptr::null_mut(),
                    std::ptr::null_mut(),
                    0,
                    &mut n_out,
                )
            },
            self.handle(),
        )?;
        let mut radius = vec![0.0; n_out];
        let mut coordination = vec![0.0; n_out];
        check(
            unsafe {
                ffi::seams_running_cn(
                    self.ptr,
                    options.rmax,
                    options.bins,
                    options.type_i,
                    options.type_j,
                    radius.as_mut_ptr(),
                    coordination.as_mut_ptr(),
                    n_out,
                    &mut n_out,
                )
            },
            self.handle(),
        )?;
        radius.truncate(n_out);
        coordination.truncate(n_out);
        Ok(RunningCoordination {
            radius,
            coordination,
        })
    }

    /// Coordination number integrated to rmax.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let cn = system.coordination(&seams::DistributionOptions { rmax: 3.5, bins: 8, ..Default::default() })?;
    /// assert!(cn.is_finite());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn coordination(&mut self, options: &DistributionOptions) -> Result<f64, Error> {
        let mut cn = 0.0;
        check(
            unsafe {
                ffi::seams_coordination(
                    self.ptr,
                    options.rmax,
                    options.bins,
                    options.type_i,
                    options.type_j,
                    &mut cn,
                )
            },
            self.handle(),
        )?;
        Ok(cn)
    }

    /// Hydrogen-bond count.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let bonds = system.hbonds(&seams::BondOptions::default())?;
    /// assert_eq!(bonds, 0);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn hbonds(&mut self, options: &BondOptions) -> Result<i32, Error> {
        let mut n_bonds = 0;
        check(
            unsafe {
                ffi::seams_hbonds(
                    self.ptr,
                    options.cutoff,
                    options.oxygen,
                    options.hydrogen,
                    options.distance,
                    options.angle,
                    &mut n_bonds,
                )
            },
            self.handle(),
        )?;
        Ok(n_bonds)
    }

    /// Number density along one cell axis.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let rho = system.density(&seams::DensityOptions { bins: 2, ..Default::default() })?;
    /// assert_eq!(rho.density.len(), 2);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn density(&mut self, options: &DensityOptions) -> Result<DensityProfile, Error> {
        let mut n_out = 0usize;
        check(
            unsafe {
                ffi::seams_density_axis(
                    self.ptr,
                    options.species,
                    options.bins,
                    options.axis as c_int,
                    std::ptr::null_mut(),
                    std::ptr::null_mut(),
                    0,
                    &mut n_out,
                )
            },
            self.handle(),
        )?;
        let mut coordinate = vec![0.0; n_out];
        let mut density = vec![0.0; n_out];
        check(
            unsafe {
                ffi::seams_density_axis(
                    self.ptr,
                    options.species,
                    options.bins,
                    options.axis as c_int,
                    coordinate.as_mut_ptr(),
                    density.as_mut_ptr(),
                    n_out,
                    &mut n_out,
                )
            },
            self.handle(),
        )?;
        coordinate.truncate(n_out);
        density.truncate(n_out);
        Ok(DensityProfile {
            coordinate,
            density,
        })
    }

    /// Topological fingerprint of the bonded graph.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let fp = system.fingerprint(&seams::FingerprintOptions::default())?;
    /// assert!(!fp.key()?.is_empty());
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn fingerprint(&mut self, options: &FingerprintOptions) -> Result<Fingerprint, Error> {
        let mut out = std::ptr::null_mut();
        check(
            unsafe {
                ffi::seams_fingerprint(
                    self.ptr,
                    options.cutoff,
                    options.species,
                    options.k,
                    options.graph as c_int,
                    options.hops,
                    i32::from(options.colour_types),
                    &mut out,
                )
            },
            self.handle(),
        )?;
        Ok(Fingerprint { ptr: out })
    }

    /// Mutual nearest unlike contacts. `spec` is site-table text.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let contacts = system.contact_pairs("1=cationHead")?;
    /// assert_eq!(contacts.pairs, 0);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn contact_pairs(&mut self, spec: &str) -> Result<Contacts, Error> {
        let text = c_text(spec)?;
        let mut out = Contacts::default();
        check(
            unsafe {
                ffi::seams_contact_pairs(
                    self.ptr,
                    text.as_ptr(),
                    &mut out.pairs,
                    &mut out.cations,
                    &mut out.anions,
                )
            },
            self.handle(),
        )?;
        Ok(out)
    }

    /// Largest polar or apolar domain.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let domain = system.domains(&seams::DomainOptions { spec: "1=waterO".into(), ..Default::default() })?;
    /// assert!(domain.largest >= 0);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn domains(&mut self, options: &DomainOptions) -> Result<DomainReport, Error> {
        let text = c_text(&options.spec)?;
        let mut subset = 0;
        let mut largest = 0;
        let mut percolation = 0.0;
        check(
            unsafe {
                ffi::seams_domains(
                    self.ptr,
                    options.cutoff,
                    text.as_ptr(),
                    options.subset as c_int,
                    &mut subset,
                    &mut largest,
                    &mut percolation,
                )
            },
            self.handle(),
        )?;
        Ok(DomainReport {
            subset,
            largest,
            percolation,
        })
    }

    /// Perfect and imperfect prism counts by ring size.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let prisms = system.prisms(&seams::PrismOptions::default())?;
    /// assert_eq!(prisms.perfect.len(), 7);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn prisms(&mut self, options: &PrismOptions) -> Result<PrismCensus, Error> {
        let mut out = PrismCensus::default();
        check(
            unsafe {
                ffi::seams_prisms(
                    self.ptr,
                    options.cutoff,
                    options.species,
                    out.perfect.as_mut_ptr(),
                    out.imperfect.as_mut_ptr(),
                    out.perfect.len(),
                )
            },
            self.handle(),
        )?;
        Ok(out)
    }

    /// Primitive-ring histogram of a monolayer.
    ///
    /// `counts[k]` is the number of rings of size k.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let sheet = system.monolayer(&seams::SheetOptions { sheet_area: 10.0, ..Default::default() })?;
    /// assert_eq!(sheet.len(), 7);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn monolayer(&mut self, options: &SheetOptions) -> Result<Vec<i32>, Error> {
        let len = usize::try_from(options.max_size)
            .map_err(|_| Error::argument("ring size is negative"))?
            + 1;
        let mut counts = vec![0; len];
        check(
            unsafe {
                ffi::seams_monolayer(
                    self.ptr,
                    options.cutoff,
                    options.species,
                    options.max_size,
                    options.sheet_area,
                    counts.as_mut_ptr(),
                    len,
                )
            },
            self.handle(),
        )?;
        Ok(counts)
    }

    /// Per-atom site kind. `spec` is site-table text such as `1=waterO`.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let kinds = system.site_table("1=waterO")?;
    /// assert_eq!(kinds.len(), 8);
    /// # Ok::<(), seams::Error>(())
    /// ```
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

    /// Per-atom ice-shell label.
    ///
    /// A non-positive `shell_rmsd` selects the library ceiling. The label
    /// is 0 for ice Ih, 1 for ice Ic, 2 for stacking disorder, or -1 when
    /// the shell does not match.
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let labels = system.symmetry(0.0)?;
    /// assert_eq!(labels.len(), 8);
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn symmetry(&mut self, shell_rmsd: f64) -> Result<Vec<i32>, Error> {
        let n = self.count() as usize;
        let mut labels = vec![0; n];
        check(
            unsafe { ffi::seams_system_symmetry(self.ptr, shell_rmsd, labels.as_mut_ptr(), n) },
            self.handle(),
        )?;
        Ok(labels)
    }

    /// Planar ice front from an ice assignment.
    ///
    /// `ice[i]` is non-zero when atom i is ice. The slice length must
    /// equal [`count`](Self::count).
    ///
    /// ```
    /// let (xyz, edge) = seams::diamond8();
    /// let mut system = seams::System::new()?;
    /// system.set_ortho(&xyz, edge, edge, edge)?;
    /// let front = system.brine(&[1; 8], &seams::BrineOptions::default())?;
    /// let _ = front.profile()?;
    /// # Ok::<(), seams::Error>(())
    /// ```
    pub fn brine(&mut self, ice: &[i32], options: &BrineOptions) -> Result<Brine, Error> {
        let mut out = std::ptr::null_mut();
        check(
            unsafe {
                ffi::seams_brine_analyse(
                    self.ptr,
                    ice.as_ptr(),
                    ice.len(),
                    options.water_type,
                    options.ion_type,
                    options.bins,
                    options.axis as c_int,
                    &mut out,
                )
            },
            self.handle(),
        )?;
        Ok(Brine { ptr: out })
    }

    fn classify(&mut self, plus: bool, options: &ChillOptions) -> Result<Classification, Error> {
        let n = self.count() as usize;
        let mut labels = vec![0; n];
        let mut counts = [0; 9];
        if plus {
            check(
                unsafe {
                    ffi::seams_chill_plus(
                        self.ptr,
                        options.cutoff,
                        options.species,
                        options.graph as c_int,
                        options.rule as c_int,
                        labels.as_mut_ptr(),
                        counts.as_mut_ptr(),
                    )
                },
                self.handle(),
            )?;
        } else {
            check(
                unsafe {
                    ffi::seams_chill(
                        self.ptr,
                        options.cutoff,
                        options.species,
                        options.graph as c_int,
                        options.rule as c_int,
                        labels.as_mut_ptr(),
                        n,
                        counts.as_mut_ptr(),
                    )
                },
                self.handle(),
            )?;
        }
        Ok(Classification { labels, counts })
    }
}
