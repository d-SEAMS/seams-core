/* Copyright (c) d-SEAMS developers */
#ifndef SEAMS_HPP_
#define SEAMS_HPP_

/**
 * @file seams.hpp
 * @brief Idiomatic C++ for the d-SEAMS C API.
 *
 * @ingroup seams_cpp
 *
 * One header, no library of its own. It includes `seams.h` and links
 * `libseams`. Handles are move-only and free themselves. A status other
 * than success throws `seams::Error`, and the message is
 * `seams_handle_error`. Coordinates and cutoffs are Angstrom. A cell is
 * nine LAMMPS numbers: origin, edge lengths, then the tilts `xy`, `xz`,
 * `yz`.
 *
 * Analyses with several knobs take an options aggregate. Results come
 * back as values, `std::vector`, or `std::string`. The caller does not
 * pass output buffers.
 */

#include "seams.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/**
 * @defgroup seams_cpp C++ API
 * @brief Move-only handles and value results over `seams.h`.
 */

namespace seams {

/**
 * @brief Status codes. Zero is success.
 * @ingroup seams_cpp
 *
 * The numeric values match `seams_status`.
 */
enum class Status : int {
  Ok = SEAMS_OK,             /**< The call succeeded. */
  Null = SEAMS_ERR_NULL,     /**< A required pointer is null. */
  Range = SEAMS_ERR_RANGE,   /**< A number is outside the accepted range. */
  Cell = SEAMS_ERR_BOX,      /**< The cell is missing or has no volume. */
  Argument = SEAMS_ERR_ARG,  /**< A caller argument is the wrong shape. */
  State = SEAMS_ERR_STATE,   /**< The handle has no frame yet. */
  Memory = SEAMS_ERR_NOMEM,  /**< An allocation failed. */
  Internal = SEAMS_ERR_INTERNAL /**< The library hit an internal failure. */
};

/**
 * @brief Neighbour graph.
 * @ingroup seams_cpp
 *
 * Cutoff keeps every pair inside the cutoff. Mutual keeps the four
 * nearest neighbours inside the cutoff plus 1.5 Angstrom, and only the
 * pairs that name each other. Seeded is the cage graph: mutual seeds
 * and union completions.
 */
enum class Graph : int {
  Cutoff = SEAMS_GRAPH_CUTOFF, /**< Every pair inside the cutoff. */
  Mutual = SEAMS_GRAPH_MUTUAL, /**< Mutual four-nearest neighbours. */
  Seeded = SEAMS_GRAPH_SEEDED  /**< Cage graph of seeds and completions. */
};

/**
 * @brief CHILL+ threshold table.
 * @ingroup seams_cpp
 *
 * Published edges are inclusive. OVITO edges are open.
 */
enum class ChillRule : int {
  Published = SEAMS_CHILL_PUBLISHED, /**< Nguyen and Molinero, 2015. */
  Ovito = SEAMS_CHILL_OVITO          /**< OVITO Chill+ modifier. */
};

/**
 * @brief CHILL and CHILL+ atom label.
 * @ingroup seams_cpp
 *
 * The nine codes are the population order in `Classification::counts`.
 */
enum class Label : int {
  Cubic = SEAMS_LABEL_CUBIC,                 /**< Cubic ice. */
  Hexagonal = SEAMS_LABEL_HEXAGONAL,         /**< Hexagonal ice. */
  Water = SEAMS_LABEL_WATER,                 /**< Liquid water. */
  Interfacial = SEAMS_LABEL_INTERFACIAL,     /**< Interfacial ice. */
  Clathrate = SEAMS_LABEL_CLATHRATE,         /**< Clathrate. */
  Interclathrate = SEAMS_LABEL_INTERCLATHRATE, /**< Interfacial clathrate. */
  Unclassified = SEAMS_LABEL_UNCLASSIFIED,   /**< No class matched. */
  Recubic = SEAMS_LABEL_RECUBIC,             /**< Reclassified cubic. */
  Rehex = SEAMS_LABEL_REHEX                  /**< Reclassified hexagonal. */
};

/**
 * @brief Primitive-ring symmetry class.
 * @ingroup seams_cpp
 */
enum class RingSymmetry : int {
  Unclassified = SEAMS_RING_UNCLASSIFIED,       /**< No class yet. */
  Ddc = SEAMS_RING_DDC,                         /**< Double-diamond six-ring. */
  HcBasal = SEAMS_RING_HC_BASAL,                /**< Basal face of a hexagonal cage. */
  HcPrismatic = SEAMS_RING_HC_PRISMATIC,        /**< Prismatic face of a hexagonal cage. */
  BothBasal = SEAMS_RING_BOTH_BASAL,            /**< Basal face of both cage types. */
  BothPrismatic = SEAMS_RING_BOTH_PRISMATIC,    /**< Prismatic face of both cage types. */
  Prism = SEAMS_RING_PRISM,                     /**< Prismatic six-ring. */
  DeformedPrism = SEAMS_RING_DEFORMED_PRISM,    /**< Deformed prismatic six-ring. */
  MixedPrism = SEAMS_RING_MIXED_PRISM,          /**< Prism mixing the two cage types. */
  Pentagon = SEAMS_RING_PENTAGON,               /**< Five-membered ring. */
  Hexagon = SEAMS_RING_HEXAGON,                 /**< Six-ring that is not a cage face. */
  Other = SEAMS_RING_OTHER                      /**< Neither a pentagon nor a hexagon face. */
};

/**
 * @brief Per-atom cage membership. The four codes are disjoint.
 * @ingroup seams_cpp
 */
enum class CageLabel : int {
  Water = SEAMS_CAGE_WATER, /**< In neither cage. */
  Hex = SEAMS_CAGE_HEX,     /**< Hexagonal cage only. */
  Cubic = SEAMS_CAGE_CUBIC, /**< Double-diamond cage only. */
  Both = SEAMS_CAGE_BOTH    /**< Both cages. */
};

/**
 * @brief Site-table kind.
 * @ingroup seams_cpp
 */
enum class SiteKind : int {
  Unspecified = SEAMS_SITE_UNSPECIFIED, /**< The table did not name this atom. */
  Cation = SEAMS_SITE_CATION,           /**< Cation head. */
  Anion = SEAMS_SITE_ANION,             /**< Anion. */
  Tail = SEAMS_SITE_TAIL,               /**< Tail atom. */
  DonorH = SEAMS_SITE_DONOR_H,          /**< Donor hydrogen. */
  Acceptor = SEAMS_SITE_ACCEPTOR,       /**< Acceptor. */
  Polar = SEAMS_SITE_POLAR,             /**< Polar atom. */
  Apolar = SEAMS_SITE_APOLAR,           /**< Apolar atom. */
  WaterO = SEAMS_SITE_WATER_O,          /**< Water oxygen. */
  WaterH = SEAMS_SITE_WATER_H,          /**< Water hydrogen. */
  Solvent = SEAMS_SITE_SOLVENT          /**< Solvent atom. */
};

/**
 * @brief Cell axis. Zero is x, one is y, two is z.
 * @ingroup seams_cpp
 */
enum class Axis : int {
  X = 0, /**< First cell edge. */
  Y = 1, /**< Second cell edge. */
  Z = 2  /**< Third cell edge. */
};

/**
 * @brief Domain subset. Apolar is 0 and polar is 1, matching the C API.
 * @ingroup seams_cpp
 */
enum class DomainSubset : int {
  Apolar = 0, /**< Apolar site subset. */
  Polar = 1   /**< Polar site subset. */
};

/**
 * @brief Failure from the C API, or from a rejected argument.
 * @ingroup seams_cpp
 *
 * `what()` is the text from `seams_handle_error`. `code()` is the
 * status integer. The object does not own a handle.
 */
class Error : public std::runtime_error {
 public:
  /**
   * @brief Build an error from a status and a message.
   * @param status Status code. Zero is not thrown by `check`.
   * @param message Detail text. An empty string is allowed.
   */
  Error(int status, std::string message)
      : std::runtime_error(std::move(message)), code_(status) {}

  /**
   * @brief Status as the `Status` enumeration when the code is known.
   * @return The matching enumerator.
   */
  [[nodiscard]] Status status() const noexcept {
    return static_cast<Status>(code_);
  }

  /**
   * @brief Status integer. Same values as `seams_status`.
   * @return The code stored at construction.
   */
  [[nodiscard]] int code() const noexcept { return code_; }

 private:
  int code_;
};

/**
 * @brief CHILL or CHILL+ labels and the nine population counts.
 * @ingroup seams_cpp
 *
 * `counts[k]` is the number of atoms whose label is the enumerator
 * with value `k`.
 */
struct Classification {
  std::vector<int> labels;     /**< One label code per atom. */
  std::array<int, 9> counts{}; /**< Population in label order. */
};

/**
 * @brief Local Steinhardt parameter, and the average when it was asked for.
 * @ingroup seams_cpp
 *
 * `average` is empty when the call did not request the Lechner-Dellago mean.
 */
struct Steinhardt {
  std::vector<double> local;   /**< One local q_l per atom. Dimensionless. */
  std::vector<double> average; /**< Neighbour mean of q_lm. Empty if skipped. */
};

/**
 * @brief Local q4 and q6 on the cutoff shell.
 * @ingroup seams_cpp
 */
struct Q4Q6 {
  std::vector<double> q4; /**< One local q4 per atom. */
  std::vector<double> q6; /**< One local q6 per atom. */
};

/**
 * @brief Lechner-Dellago q6-bar and the ice count above a cut.
 * @ingroup seams_cpp
 */
struct Q6Bar {
  std::vector<double> values; /**< One q6-bar per atom. */
  double mean = 0;            /**< Mean of the finite values. */
  double median = 0;          /**< Median of the finite values. */
  int n_ice = 0;              /**< Count strictly above the cut. */
};

/**
 * @brief Partial radial distribution.
 * @ingroup seams_cpp
 */
struct Distribution {
  std::vector<double> radius; /**< Bin centers, Angstrom. */
  std::vector<double> g;      /**< g(r), dimensionless. */
  std::vector<int> count;     /**< Pair count in each bin. */
  double volume = 0;          /**< Cell volume, cubic Angstrom. */
};

/**
 * @brief Running coordination number along a partial RDF.
 * @ingroup seams_cpp
 */
struct RunningCoordination {
  std::vector<double> radius;       /**< Bin centers, Angstrom. */
  std::vector<double> coordination; /**< Running coordination, dimensionless. */
};

/**
 * @brief Disjoint cage populations.
 * @ingroup seams_cpp
 */
struct CageCounts {
  int hexagonal = 0; /**< Hexagonal cage only. */
  int cubic = 0;     /**< Double-diamond cage only. */
  int water = 0;     /**< Neither cage. */
  int both = 0;      /**< Both cages. */
};

/**
 * @brief Closed cages whose ring census matches a signature.
 * @ingroup seams_cpp
 */
struct SignatureMatch {
  std::size_t cages = 0; /**< Number of closed cages. */
  std::size_t atoms = 0; /**< Distinct atoms in those cages. */
};

/**
 * @brief Mutual nearest unlike contacts of a site table.
 * @ingroup seams_cpp
 */
struct Contacts {
  std::size_t pairs = 0; /**< Contact count. */
  int cations = 0;       /**< Cation molecules in the ion cloud. */
  int anions = 0;        /**< Anion molecules in the ion cloud. */
};

/**
 * @brief Largest polar or apolar domain.
 * @ingroup seams_cpp
 */
struct DomainReport {
  int subset = 0;          /**< Atoms in the requested subset. */
  int largest = 0;         /**< Largest connected component. */
  double percolation = 0;  /**< Percolation measure, dimensionless. */
};

/**
 * @brief Perfect and imperfect prism counts by ring size.
 * @ingroup seams_cpp
 *
 * Index is the ring size. Sizes 3 through 6 are filled.
 */
struct PrismCensus {
  std::array<int, 7> perfect{};   /**< Perfect stacks. */
  std::array<int, 7> imperfect{}; /**< Imperfect stacks. */
};

/**
 * @brief Number density along one cell axis.
 * @ingroup seams_cpp
 */
struct DensityProfile {
  std::vector<double> coordinate; /**< Bin centers, Angstrom. */
  std::vector<double> density;    /**< Number density. */
};

/**
 * @brief Number and mass density of one histogram.
 * @ingroup seams_cpp
 */
struct FractionalProfile {
  std::vector<double> coordinate;     /**< Bin centres. Fractional, or radians for azimuth. */
  std::vector<double> number_density; /**< Atoms per cubic Angstrom. */
  std::vector<double> mass_density;   /**< Caller mass units per cubic Angstrom. Mass 1 when masses is empty.
                                         An amu mass is amu per cubic Angstrom; times 1.66054 it is g/cm^3. */
  std::vector<double> volume;         /**< Bin volume, cubic Angstrom. */
};

/**
 * @brief Knobs for a fractional or directional density.
 * @ingroup seams_cpp
 */
struct FractionalDensityOptions {
  int species = 0;                         /**< Species filter. Zero keeps every atom. */
  int bins = 1;                            /**< Positive bin count. */
  Axis axis = Axis::Z;                     /**< Cell axis. Ignored when use_direction is set. */
  bool use_direction = false;              /**< Bin an arbitrary Cartesian direction. */
  std::array<double, 3> direction{{0, 0, 1}}; /**< Cartesian direction. */
  double origin = 0;                       /**< Fractional origin, wrapped into [0, 1). */
  std::span<const double> masses{};        /**< One mass per stored atom. Empty means mass 1. */
};

/**
 * @brief Pore axis and centre.
 * @ingroup seams_cpp
 */
struct PoreFrame {
  std::array<double, 3> axis{};               /**< Unit Cartesian axis. */
  std::array<double, 3> centre{};             /**< Cartesian centre, Angstrom. */
  std::array<double, 3> centre_fractional{};  /**< Fractional centre. */
  int axis_index = 0;                         /**< 0, 1, or 2. */
};

/**
 * @brief Knobs for detecting a pore axis and centre.
 * @ingroup seams_cpp
 */
struct PoreFrameOptions {
  int framework = 0;                              /**< Species for an automatic centre. Zero keeps every atom. */
  int axis = -1;                                  /**< 0, 1, 2, or -1 for the longest edge. */
  bool use_direction = false;                     /**< Direction replaces axis and must be a cell edge. */
  std::array<double, 3> direction{{0, 0, 1}};     /**< Cartesian direction. */
  bool use_centre = false;                        /**< Use centre_fractional instead of the framework centroid. */
  std::array<double, 3> centre_fractional{{0, 0, 0}}; /**< Fractional centre. */
};

/**
 * @brief Radial, azimuthal, and axial pore profiles.
 * @ingroup seams_cpp
 */
struct PoreProfiles {
  FractionalProfile radial;     /**< Cylindrical shells. */
  FractionalProfile azimuthal;  /**< Azimuthal sectors. */
  FractionalProfile axial;      /**< Slabs along the pore axis. */
};

/**
 * @brief Knobs for cylindrical pore profiles.
 * @ingroup seams_cpp
 */
struct PoreProfileOptions {
  int species = 0;                                /**< Species filter. Zero keeps every atom. */
  int radial_bins = 1;                            /**< Positive radial bin count. */
  int azimuthal_bins = 1;                         /**< Positive azimuthal bin count. */
  int axial_bins = 1;                             /**< Positive axial bin count. */
  double r_max = 0;                               /**< Outer radius, Angstrom. Zero covers the face. */
  int axis = -1;                                  /**< 0, 1, 2, or -1 for the longest edge. */
  bool use_direction = false;                     /**< Direction replaces axis and must be a cell edge. */
  std::array<double, 3> direction{{0, 0, 1}};     /**< Cartesian direction. */
  bool use_centre = false;                        /**< Use centre_fractional instead of the framework centroid. */
  std::array<double, 3> centre_fractional{{0, 0, 0}}; /**< Fractional centre. */
  int framework = 0;                              /**< Species for an automatic centre. */
  double origin = 0;                              /**< Fractional axial origin. */
  std::span<const double> masses{};               /**< One mass per stored atom. Empty means mass 1. */
};

/**
 * @brief Knobs for host-guest contacts.
 * @ingroup seams_cpp
 */
struct HostGuestOptions {
  int guest = 1;          /**< Positive guest species. */
  int framework = 0;      /**< Framework species. Zero keeps every other type. Equal to guest is rejected. */
  double cutoff = 0;      /**< Contact distance, Angstrom. Zero counts a coincident site. */
};

/**
 * @brief Nearest framework site and contact counts.
 * @ingroup seams_cpp
 */
struct HostGuestContacts {
  std::vector<int> guest_index;         /**< Guest atom indices. */
  std::vector<int> nearest_index;       /**< Nearest framework index, or -1. */
  std::vector<int> nearest_type;        /**< Type of that site. */
  std::vector<double> nearest_distance; /**< Distance, Angstrom. Negative when no site exists. */
  std::vector<int> framework_types;     /**< Framework types, sorted. */
  std::vector<int> contacts;            /**< Row-major counts, one row per guest. */
};

/**
 * @brief The two planar crossings of an ice slab.
 * @ingroup seams_cpp
 */
struct BrineFronts {
  double enter = 0;    /**< Crossing on entering the ice, Angstrom. */
  double exit = 0;     /**< Crossing on leaving the ice, Angstrom. */
  double thickness = 0; /**< Ice run along the axis, Angstrom. */
  bool found = false;  /**< Both crossings exist. */
};

/**
 * @brief Molal rejection coefficient.
 * @ingroup seams_cpp
 */
struct Rejection {
  double ice = 0;     /**< Bulk ice molality, mol/kg. */
  double liquid = 0;  /**< Bulk liquid molality, mol/kg. */
  double k = 0;       /**< Ratio of the two molalities. */
  bool defined = false; /**< The ratio is defined. */
};

/**
 * @brief Concentration profile along the brine axis.
 * @ingroup seams_cpp
 */
struct BrineProfile {
  std::vector<double> coordinate;    /**< Bin centers, Angstrom. */
  std::vector<double> molality;      /**< mol/kg. */
  std::vector<double> ice_fraction;  /**< Dimensionless. */
};

/**
 * @brief Brine pockets trapped in the ice.
 * @ingroup seams_cpp
 */
struct Pockets {
  int pockets = 0;           /**< Pocket count, excluding the reservoir. */
  int reservoir_waters = 0;  /**< Waters in the liquid reservoir. */
};

/**
 * @brief Knobs for CHILL and CHILL+.
 * @ingroup seams_cpp
 */
struct ChillOptions {
  double cutoff = 3.5;                         /**< Positive cutoff, Angstrom. */
  int species = 0;                             /**< Species filter. Zero keeps every atom. */
  Graph graph = Graph::Cutoff;                 /**< Neighbour graph. */
  ChillRule rule = ChillRule::Published;       /**< Threshold table. */
};

/**
 * @brief Knobs for a local Steinhardt parameter on the cutoff shell.
 * @ingroup seams_cpp
 */
struct SteinhardtOptions {
  double cutoff = 3.5;  /**< Positive cutoff, Angstrom. */
  int species = 0;      /**< Species filter. Zero keeps every atom. */
  int order = 6;        /**< Degree. 3, 4, 6, 8, or 12. */
  bool average = false; /**< Also return the Lechner-Dellago mean. */
};

/**
 * @brief Knobs for q_l on a named graph.
 * @ingroup seams_cpp
 */
struct OrderOptions {
  double cutoff = 3.5;          /**< Positive cutoff, Angstrom. */
  int species = 0;              /**< Species filter. */
  int order = 6;                /**< Degree. 3, 4, 6, 8, or 12. */
  Graph graph = Graph::Cutoff;  /**< Cutoff or mutual. */
  bool average = false;         /**< Also return the neighbour mean. */
};

/**
 * @brief Knobs for the q6-bar ice summary.
 * @ingroup seams_cpp
 */
struct Q6BarOptions {
  double cutoff = 3.5;   /**< Positive cutoff, Angstrom. */
  int species = 0;       /**< Species filter. */
  double ice_cut = 0.36; /**< Finite threshold. The nucleation census uses 0.36. */
};

/**
 * @brief Knobs for a neighbour list.
 * @ingroup seams_cpp
 */
struct NeighbourOptions {
  double cutoff = 3.5;          /**< Positive cutoff, Angstrom. */
  int species = 0;              /**< Species filter. */
  Graph graph = Graph::Cutoff;  /**< Cutoff or mutual. */
};

/**
 * @brief Knobs for a k-nearest neighbour list.
 * @ingroup seams_cpp
 */
struct KnnOptions {
  int k = 4;            /**< Positive neighbour count. */
  double cutoff = 3.5;  /**< Nominal cutoff, Angstrom. Candidates lie inside cutoff + 1.5. */
  int species = 0;      /**< Species filter. */
  bool mutual = false;  /**< Keep only pairs that name each other. */
};

/**
 * @brief Knobs for primitive rings.
 * @ingroup seams_cpp
 */
struct RingOptions {
  double cutoff = 3.5;          /**< Positive cutoff, Angstrom. */
  int species = 0;              /**< Species filter. */
  Graph graph = Graph::Cutoff;  /**< Cutoff or mutual. */
  int max_size = 6;             /**< Largest ring size, from 3 through 12. */
};

/**
 * @brief Knobs for hexagonal and double-diamond cages.
 * @ingroup seams_cpp
 */
struct CageOptions {
  double cutoff = 3.5;          /**< Positive cutoff, Angstrom. */
  int species = 0;              /**< Species filter. */
  Graph graph = Graph::Seeded;  /**< Cage graph. Seeded is the usual choice. */
  int k = 4;                    /**< Neighbour count used to seed the graph. */
  bool complete = false;        /**< Complete partial cages. */
};

/**
 * @brief Knobs for a signature cage census.
 * @ingroup seams_cpp
 */
struct SignatureOptions {
  double cutoff = 3.5;          /**< Positive cutoff, Angstrom. */
  int species = 0;              /**< Species filter. */
  Graph graph = Graph::Cutoff;  /**< Neighbour graph. */
  std::string spec = "512";     /**< Ring census text, for example `512`. */
};

/**
 * @brief Knobs for a partial RDF, its running coordination, or the integral.
 * @ingroup seams_cpp
 */
struct DistributionOptions {
  double rmax = 6;   /**< Outer limit, Angstrom. */
  int bins = 0;      /**< Bin count. Zero selects about rmax / 0.1. */
  int type_i = 0;    /**< First species. Zero keeps every atom. */
  int type_j = 0;    /**< Second species. Zero keeps every atom. */
};

/**
 * @brief Knobs for a hydrogen-bond count.
 * @ingroup seams_cpp
 */
struct BondOptions {
  double cutoff = 3.5;       /**< Heavy-atom cutoff, Angstrom. */
  int oxygen = 1;            /**< Heavy-atom species. Zero keeps every non-hydrogen. */
  int hydrogen = 2;          /**< Hydrogen species. Must be positive. */
  double distance = 2.42;    /**< Acceptor-H distance, Angstrom. */
  double angle = 30;         /**< O-O-H angle, degrees. */
};

/**
 * @brief One directed hydrogen bond.
 * @ingroup seams_cpp
 *
 * Indices are frame indices. Heavies from set_hydrogen_frame occupy
 * `[0, n_heavy)` and hydrogens follow. `angle` is acceptor-first degrees:
 * 0 is a linear donor-hydrogen-acceptor geometry.
 */
struct DirectedBond {
  int donor = -1;       /**< Donor heavy-atom index. */
  int hydrogen = -1;    /**< Hydrogen index. */
  int acceptor = -1;    /**< Acceptor heavy-atom index. */
  double d_AH = 0;      /**< Acceptor-hydrogen distance, Angstrom. */
  double d_DA = 0;      /**< Donor-acceptor distance, Angstrom. */
  double angle = 0;     /**< Acceptor-first angle, degrees. */
};

/**
 * @brief Knobs for directed hydrogen bonds.
 * @ingroup seams_cpp
 */
struct DirectedBondOptions {
  double cutoff = 3.5;     /**< Donor-acceptor cutoff, Angstrom. Strictly inside. */
  double distance = 2.42;  /**< Acceptor-hydrogen cutoff, Angstrom. Strictly inside. */
  double angle = 30;       /**< Acceptor-first angle, degrees. Kept at or below this value. */
  int hydrogen = 2;        /**< Hydrogen species. Must be positive and distinct from the type sets. */
  bool exclude_same_mol = true; /**< Drop a donor and acceptor that share a molID. */
};

/**
 * @brief Knobs for the axial number density.
 * @ingroup seams_cpp
 */
struct DensityOptions {
  int species = 0;   /**< Species filter. Zero keeps every atom. */
  int bins = 0;      /**< Bin count. Zero selects about one bin per 0.1 Angstrom. */
  Axis axis = Axis::Z; /**< Cell axis. */
};

/**
 * @brief Knobs for a topological fingerprint.
 * @ingroup seams_cpp
 */
struct FingerprintOptions {
  double cutoff = 3.5;          /**< Positive cutoff, Angstrom. */
  int species = 0;              /**< Species filter. */
  int k = 4;                    /**< Neighbour count on the bonded graph. */
  Graph graph = Graph::Cutoff;  /**< Neighbour graph. */
  int hops = 2;                 /**< Neighbourhood depth. Zero selects 2. */
  bool color_types = false;    /**< Color vertices by LAMMPS type. */
};

/**
 * @brief Knobs for a domain search.
 * @ingroup seams_cpp
 */
struct DomainOptions {
  double cutoff = 3.5;                 /**< Positive cutoff, Angstrom. */
  std::string spec;                    /**< Site-table text, for example `1=waterO`. */
  DomainSubset subset = DomainSubset::Apolar; /**< Which subset to measure. */
};

/**
 * @brief Knobs for prism stacks.
 * @ingroup seams_cpp
 */
struct PrismOptions {
  double cutoff = 3.5; /**< Positive cutoff, Angstrom. */
  int species = 0;     /**< Species filter. */
};

/**
 * @brief Knobs for a monolayer ring histogram.
 * @ingroup seams_cpp
 */
struct SheetOptions {
  double cutoff = 3.5;     /**< Positive cutoff, Angstrom. */
  int species = 0;         /**< Species filter. */
  int max_size = 6;        /**< Largest ring size. */
  double sheet_area = 0;   /**< Sheet area, square Angstrom. Must be positive. */
};

/**
 * @brief Knobs for a planar ice front.
 * @ingroup seams_cpp
 */
struct BrineOptions {
  int water_type = 0; /**< Species kept as water, or 0. */
  int ion_type = 0;   /**< Species kept as ions, or 0 for none. */
  int bins = 0;       /**< Bins along the axis. Zero selects 20. */
  Axis axis = Axis::Z; /**< Front axis. */
};

/**
 * @brief Per-molecule intact periodic centres.
 * @ingroup seams_cpp
 */
struct MoleculeCentres {
  std::vector<int> mol_ids;    /**< Unique molecule ids in appearance order. */
  std::vector<double> centres; /**< Interleaved coordinates, 3 * n_mols doubles. */
};

/**
 * @brief Per-water orientation vectors under minimum image.
 * @ingroup seams_cpp
 */
struct WaterVectors {
  std::vector<int> mol_ids;       /**< Water molecule ids. */
  std::vector<double> oh_mid;     /**< Raw O-to-H-midpoint vectors, Angstrom. */
  std::vector<double> bisector;   /**< Unit bisector vectors. */
  std::vector<double> hh;         /**< Raw H1-to-H2 vectors, Angstrom. */
  std::vector<double> hh_unit;    /**< Unit H1-to-H2 vectors. */
};

/**
 * @brief Per-molecule dipole or geometric orientation vectors.
 * @ingroup seams_cpp
 */
struct MoleculeDipoles {
  std::vector<int> mol_ids;          /**< Unique molecule ids. */
  std::vector<double> dipoles;       /**< Raw dipole or bond vectors, 3 * n_mols. */
  std::vector<double> unit_vectors;  /**< Normalized orientation vectors, 3 * n_mols. */
};

/**
 * @brief Axis-projected orientational order parameters P1 and P2.
 * @ingroup seams_cpp
 */
struct OrientationOrder {
  std::vector<double> cos_theta; /**< Per-vector cos(theta). */
  std::vector<double> p2;        /**< Per-vector P2(cos(theta)). */
  double mean_p1 = 0;            /**< Ensemble mean of cos(theta). */
  double mean_p2 = 0;            /**< Ensemble mean of P2(cos(theta)). */
};

/**
 * @brief Accumulated multi-origin MSD curves.
 * @ingroup seams_cpp
 */
struct MsdCurve {
  std::vector<double> lag_time;      /**< Lag times k * dt. */
  std::vector<double> total;         /**< 3D total MSD. */
  std::vector<double> axial;         /**< 1D axial MSD. */
  std::vector<double> lateral;       /**< 2D lateral MSD. */
  std::vector<std::size_t> n_origins; /**< Time-origin count at each lag. */
};

/**
 * @brief Linear Einstein diffusion fit coefficients.
 * @ingroup seams_cpp
 */
struct DiffusionFit {
  double d_total = 0;   /**< 3D diffusion coefficient (slope / 6). */
  double d_axial = 0;   /**< 1D axial diffusion coefficient (slope / 2). */
  double d_lateral = 0; /**< 2D lateral diffusion coefficient (slope / 4). */
  /** Log-log anomalous diffusion exponent, NaN when the fit was skipped. */
  double alpha = std::numeric_limits<double>::quiet_NaN();
  bool alpha_fitted = false; /**< False when the window had < 2 usable log-log lags. */
};

/**
 * @brief Intermittent and continuous H-bond lifetime curves.
 * @ingroup seams_cpp
 */
struct HbondLifetimeCurve {
  std::vector<double> lag_time;       /**< Lag times k * dt. */
  std::vector<double> c_intermittent; /**< Intermittent autocorrelation C_HB. */
  std::vector<double> s_continuous;   /**< Continuous survival S_HB. */
  std::vector<std::size_t> n_origins; /**< Time-origin count at each lag. */
};

/**
 * @brief Integrated and single-exponential H-bond lifetimes.
 * @ingroup seams_cpp
 */
struct HbondLifetimeFit {
  double tau_c_int = 0; /**< Trapezoidal integral of C_HB. */
  double tau_c_exp = 0; /**< Exponential decay time of C_HB. */
  double tau_s_int = 0; /**< Trapezoidal integral of S_HB. */
  double tau_s_exp = 0; /**< Exponential decay time of S_HB. */
};

/**
 * @brief Accumulated C1 and C2 orientational correlation curves.
 * @ingroup seams_cpp
 */
struct ReorientationCurve {
  std::vector<double> lag_time;       /**< Lag times k * dt. */
  std::vector<double> c1;             /**< First Legendre correlation C1. */
  std::vector<double> c2;             /**< Second Legendre correlation C2. */
  std::vector<std::size_t> n_origins; /**< Time-origin count at each lag. */
};

/**
 * @brief Integrated and single-exponential reorientation times.
 * @ingroup seams_cpp
 */
struct ReorientationFit {
  double tau1_int = 0; /**< Trapezoidal integral of C1. */
  double tau1_exp = 0; /**< Exponential decay time of C1. */
  double tau2_int = 0; /**< Trapezoidal integral of C2. */
  double tau2_exp = 0; /**< Exponential decay time of C2. */
};

namespace detail {

inline void check(int status, const void *handle) {
  if (status == SEAMS_OK) {
    return;
  }
  const char *text = seams_handle_error(handle);
  throw Error(status, text == nullptr ? std::string() : std::string(text));
}

template <class T, void (*Destroy)(T *)>
class Handle {
 public:
  Handle() = default;
  explicit Handle(T *ptr) noexcept : ptr_(ptr) {}
  Handle(const Handle &) = delete;
  Handle &operator=(const Handle &) = delete;
  Handle(Handle &&other) noexcept : ptr_(std::exchange(other.ptr_, nullptr)) {}
  Handle &operator=(Handle &&other) noexcept {
    if (this != &other) {
      reset();
      ptr_ = std::exchange(other.ptr_, nullptr);
    }
    return *this;
  }
  ~Handle() { reset(); }

  [[nodiscard]] explicit operator bool() const noexcept { return ptr_ != nullptr; }

 protected:
  [[nodiscard]] T *get() const noexcept { return ptr_; }

  void reset() noexcept {
    if (ptr_ != nullptr) {
      Destroy(ptr_);
      ptr_ = nullptr;
    }
  }
  T *ptr_ = nullptr;
};

// The C trajectory calls take no xyz length and read nframes * 3 * natoms.
inline std::size_t trajectory_coordinates(std::span<const double> xyz, std::size_t nframes,
                                          std::size_t natoms) {
  constexpr std::size_t most = std::numeric_limits<std::size_t>::max();
  if (natoms > most / 3 || (natoms != 0 && nframes > most / (3 * natoms))) {
    throw Error(SEAMS_ERR_RANGE, "nframes * 3 * natoms overflows size_t");
  }
  const std::size_t want = nframes * 3 * natoms;
  if (xyz.size() != want) {
    throw Error(SEAMS_ERR_RANGE, "xyz length is not nframes * 3 * natoms");
  }
  return want;
}

template <class Call>
inline std::string copy_text(Call call, const void *handle) {
  std::size_t n_out = 0;
  check(call(handle, nullptr, 0, &n_out), handle);
  std::string text(n_out + 1, '\0');
  check(call(handle, text.data(), text.size(), &n_out), handle);
  text.resize(n_out);
  return text;
}

}  // namespace detail

/**
 * @brief Chordless cycles of one length.
 * @ingroup seams_cpp
 *
 * `atoms` holds count() rows of `size` atom indices, in the order
 * seams_rings_chordless writes them. `per_atom` is the number of cycles
 * through each atom.
 */
struct ChordlessCycles {
  int size = 0;                /**< Cycle length. */
  std::vector<int> atoms;      /**< count() rows of `size` indices. */
  std::vector<int> per_atom;   /**< Cycles through each atom. */

  /** @brief Number of cycles. */
  [[nodiscard]] std::size_t count() const {
    return size > 0 ? atoms.size() / static_cast<std::size_t>(size) : 0;
  }
  /** @brief Atom indices of one cycle. */
  [[nodiscard]] std::vector<int> cycle(std::size_t index) const {
    const auto first =
        atoms.begin() + static_cast<std::ptrdiff_t>(index * static_cast<std::size_t>(size));
    return std::vector<int>(first, first + size);
  }
};

/**
 * @brief Owned neighbour list.
 * @ingroup seams_cpp
 *
 * Move-only. Destruction frees the C handle. The self entry is omitted.
 * A neighbour list is not safe to share across threads.
 */
class Neighbours : public detail::Handle<seams_neighbours, seams_neighbours_destroy> {
 public:
  using detail::Handle<seams_neighbours, seams_neighbours_destroy>::Handle;

  /**
   * @brief Number of neighbour pairs, omitting self.
   * @return The pair count.
   */
  [[nodiscard]] std::size_t nnz() const {
    std::size_t n_offsets = 0;
    std::size_t count = 0;
    detail::check(seams_neighbours_copy(get(), nullptr, 0, &n_offsets, nullptr, 0, &count), get());
    return count;
  }

  /**
   * @brief Degree of each atom.
   * @return One integer per atom. The self entry is not counted.
   */
  [[nodiscard]] std::vector<int> degrees() const {
    std::size_t n_offsets = 0;
    std::size_t count = 0;
    detail::check(seams_neighbours_copy(get(), nullptr, 0, &n_offsets, nullptr, 0, &count), get());
    std::vector<int> offsets(n_offsets);
    std::vector<int> cols(count);
    detail::check(seams_neighbours_copy(get(), offsets.data(), offsets.size(), &n_offsets,
                                        cols.data(), cols.size(), &count),
                  get());
    std::vector<int> degree;
    if (n_offsets > 0) {
      degree.reserve(n_offsets - 1);
      for (std::size_t i = 0; i + 1 < n_offsets; ++i) {
        degree.push_back(offsets[i + 1] - offsets[i]);
      }
    }
    return degree;
  }

  /**
   * @brief Directed pairs `(atom, neighbour)`.
   * @return One pair per stored column. Self is omitted.
   */
  [[nodiscard]] std::vector<std::pair<int, int>> pairs() const {
    std::size_t n_offsets = 0;
    std::size_t count = 0;
    detail::check(seams_neighbours_copy(get(), nullptr, 0, &n_offsets, nullptr, 0, &count), get());
    std::vector<int> offsets(n_offsets);
    std::vector<int> cols(count);
    detail::check(seams_neighbours_copy(get(), offsets.data(), offsets.size(), &n_offsets,
                                        cols.data(), cols.size(), &count),
                  get());
    std::vector<std::pair<int, int>> out;
    out.reserve(count);
    for (std::size_t atom = 0; atom + 1 < n_offsets; ++atom) {
      for (int cursor = offsets[atom]; cursor < offsets[atom + 1]; ++cursor) {
        const int neighbour = cols[static_cast<std::size_t>(cursor)];
        if (neighbour != static_cast<int>(atom)) {
          out.emplace_back(static_cast<int>(atom), neighbour);
        }
      }
    }
    return out;
  }

  /**
   * @brief Chordless cycles of one length on this graph.
   * @param size Cycle length, from 3 through 12.
   * @return The rows and the per-atom counts of seams_rings_chordless.
   */
  [[nodiscard]] ChordlessCycles chordless(int size) const {
    ChordlessCycles out;
    out.size = size;
    out.per_atom.resize(static_cast<std::size_t>(seams_neighbours_count(get())));
    std::size_t count = 0;
    detail::check(seams_rings_chordless(get(), size, nullptr, 0, &count, out.per_atom.data(),
                                        out.per_atom.size()),
                  get());
    out.atoms.resize(count * static_cast<std::size_t>(size));
    detail::check(seams_rings_chordless(get(), size, out.atoms.data(), out.atoms.size(), &count,
                                        nullptr, 0),
                  get());
    return out;
  }
};

/**
 * @brief Owned primitive-ring set.
 * @ingroup seams_cpp
 *
 * Move-only. `begin` and `end` walk the rings and yield vertex indices.
 */
class Rings : public detail::Handle<seams_rings, seams_rings_destroy> {
 public:
  using detail::Handle<seams_rings, seams_rings_destroy>::Handle;

  /**
   * @brief Number of primitive rings.
   * @return The ring count. A moved-from handle reports 0.
   */
  [[nodiscard]] int size() const { return seams_rings_count(get()); }

  /**
   * @brief Atom indices of one ring.
   * @param index Ring index, from 0 to `size() - 1`.
   * @return Frame indices, in ring order.
   */
  [[nodiscard]] std::vector<int> vertices(int index) const {
    std::size_t n_out = 0;
    detail::check(seams_rings_vertices(get(), index, nullptr, 0, &n_out), get());
    std::vector<int> atoms(n_out);
    detail::check(seams_rings_vertices(get(), index, atoms.data(), n_out, &n_out), get());
    atoms.resize(n_out);
    return atoms;
  }

  /**
   * @brief Histogram of ring size.
   * @param up_to Highest size stored. `counts[k]` is the number of rings of size k.
   * @return The histogram. Length is `up_to + 1`.
   */
  [[nodiscard]] std::vector<int> histogram(int up_to) const {
    std::vector<int> counts(static_cast<std::size_t>(up_to + 1));
    detail::check(seams_rings_histogram(get(), counts.data(), up_to + 1), get());
    return counts;
  }

  /**
   * @brief Fraction of atoms that sit on a ring of one size.
   * @param ring_size Size to count. Zero counts every ring.
   * @return Dimensionless fraction.
   */
  [[nodiscard]] double fraction(int ring_size) const {
    double value = 0;
    detail::check(seams_rings_fraction(get(), ring_size, &value), get());
    return value;
  }

  /** @brief Input iterator. `operator*` copies one ring. */
  class iterator {
   public:
    using value_type = std::vector<int>;
    using difference_type = std::ptrdiff_t;
    iterator(const Rings *owner, int index) : owner_(owner), index_(index) {}
    [[nodiscard]] value_type operator*() const { return owner_->vertices(index_); }
    iterator &operator++() {
      ++index_;
      return *this;
    }
    [[nodiscard]] bool operator==(const iterator &other) const { return index_ == other.index_; }
    [[nodiscard]] bool operator!=(const iterator &other) const { return !(*this == other); }

   private:
    const Rings *owner_ = nullptr;
    int index_ = 0;
  };

  /** @brief First ring. */
  [[nodiscard]] iterator begin() const { return iterator(this, 0); }
  /** @brief Past the last ring. */
  [[nodiscard]] iterator end() const { return iterator(this, size()); }
};

/**
 * @brief Owned hexagonal and double-diamond membership.
 * @ingroup seams_cpp
 */
class Cages : public detail::Handle<seams_cage_set, seams_cages_destroy> {
 public:
  using detail::Handle<seams_cage_set, seams_cages_destroy>::Handle;

  /**
   * @brief Disjoint population counts.
   * @return Hexagonal, cubic, water, and both.
   */
  [[nodiscard]] CageCounts counts() const {
    CageCounts out;
    detail::check(seams_cages_counts(get(), &out.hexagonal, &out.cubic, &out.water, &out.both),
                  get());
    return out;
  }

  /**
   * @brief One cage code per atom.
   * @param n_atoms Length of the frame the cages were built from.
   * @return Codes matching `CageLabel`.
   */
  [[nodiscard]] std::vector<int> labels(std::size_t n_atoms) const {
    std::vector<int> per_atom(n_atoms);
    detail::check(seams_cages_labels(get(), per_atom.data(), per_atom.size()), get());
    return per_atom;
  }
};

/**
 * @brief Owned topological fingerprint.
 * @ingroup seams_cpp
 */
class Fingerprint : public detail::Handle<seams_topo, seams_fingerprint_destroy> {
 public:
  using detail::Handle<seams_topo, seams_fingerprint_destroy>::Handle;

  /**
   * @brief Frame fingerprint text.
   * @return The key, without a trailing NUL.
   */
  [[nodiscard]] std::string key() const {
    return detail::copy_text(
        [](const void *handle, char *buf, std::size_t cap, std::size_t *n_out) {
          return seams_fingerprint_key(static_cast<const seams_topo *>(handle), buf, cap, n_out);
        },
        get());
  }

  /**
   * @brief Local key of one atom.
   * @param index Atom index in the frame.
   * @return The local key.
   */
  [[nodiscard]] std::string atom(int index) const {
    std::size_t n_out = 0;
    detail::check(seams_fingerprint_atom(get(), index, nullptr, 0, &n_out), get());
    std::string text(n_out + 1, '\0');
    detail::check(seams_fingerprint_atom(get(), index, text.data(), text.size(), &n_out), get());
    text.resize(n_out);
    return text;
  }

  /**
   * @brief Primitive-ring census stored on the fingerprint.
   * @return `counts[k]` is the number of rings of size k. Length 13.
   */
  [[nodiscard]] std::vector<int> rings() const {
    std::vector<int> counts(13);
    detail::check(seams_fingerprint_rings(get(), counts.data(), counts.size()), get());
    return counts;
  }
};

/**
 * @brief Owned planar ice-front measurement.
 * @ingroup seams_cpp
 */
class Brine : public detail::Handle<seams_brine, seams_brine_destroy> {
 public:
  using detail::Handle<seams_brine, seams_brine_destroy>::Handle;

  /** @brief The two planar crossings. Lengths are Angstrom. */
  [[nodiscard]] BrineFronts fronts() const {
    BrineFronts out;
    int found = 0;
    detail::check(seams_brine_fronts(get(), &out.enter, &out.exit, &out.thickness, &found), get());
    out.found = found != 0;
    return out;
  }

  /** @brief Molal rejection coefficient. Molality is mol/kg. */
  [[nodiscard]] Rejection rejection() const {
    Rejection out;
    int defined = 0;
    detail::check(seams_brine_rejection(get(), &out.ice, &out.liquid, &out.k, &defined), get());
    out.defined = defined != 0;
    return out;
  }

  /** @brief Concentration profile along the front axis. */
  [[nodiscard]] BrineProfile profile() const {
    std::size_t n_out = 0;
    detail::check(seams_brine_profile(get(), nullptr, nullptr, nullptr, 0, &n_out), get());
    BrineProfile out;
    out.coordinate.resize(n_out);
    out.molality.resize(n_out);
    out.ice_fraction.resize(n_out);
    detail::check(seams_brine_profile(get(), out.coordinate.data(), out.molality.data(),
                                      out.ice_fraction.data(), n_out, &n_out),
                  get());
    out.coordinate.resize(n_out);
    out.molality.resize(n_out);
    out.ice_fraction.resize(n_out);
    return out;
  }

  /** @brief Pockets trapped in the ice, and the liquid reservoir. */
  [[nodiscard]] Pockets pockets() const {
    Pockets out;
    detail::check(seams_brine_pockets(get(), &out.pockets, &out.reservoir_waters), get());
    return out;
  }
};

/**
 * @brief Owned streaming MSD accumulator.
 * @ingroup seams_cpp
 */
class Msd : public detail::Handle<seams_msd, seams_msd_destroy> {
 public:
  /**
   * @brief Allocate a streaming MSD accumulator.
   * @param natoms Particle count per frame.
   * @param max_lag Maximum lag in frames, or 0 for all lags.
   * @param origin_stride Frame stride between origins.
   * @param dt Frame spacing in physical time units.
   * @param axis Axial direction 3-vector.
   */
  Msd(std::size_t natoms, std::size_t max_lag, std::size_t origin_stride, double dt,
      std::span<const double, 3> axis) {
    seams_msd *raw = nullptr;
    detail::check(seams_msd_create(natoms, max_lag, origin_stride, dt, axis.data(), &raw), nullptr);
    *this = Msd(raw);
  }

  using detail::Handle<seams_msd, seams_msd_destroy>::Handle;

  /** @brief Push one unwrapped coordinate frame. */
  void push(std::span<const double> xyz) {
    if (xyz.size() % 3 != 0) {
      throw Error(SEAMS_ERR_ARG, "coordinates are not a multiple of 3");
    }
    detail::check(seams_msd_push(get(), xyz.data(), xyz.size() / 3, nullptr, 0), get());
  }

  /** @brief Push one periodic coordinate frame and unwrap incrementally. */
  void push(std::span<const double> xyz, std::span<const double, 9> box9) {
    if (xyz.size() % 3 != 0) {
      throw Error(SEAMS_ERR_ARG, "coordinates are not a multiple of 3");
    }
    detail::check(seams_msd_push(get(), xyz.data(), xyz.size() / 3, box9.data(), 9), get());
  }

  /** @brief Copy out accumulated total, axial, and lateral MSD curves. */
  [[nodiscard]] MsdCurve curve() const {
    std::size_t nl = 0;
    detail::check(seams_msd_curve(get(), nullptr, nullptr, nullptr, nullptr, nullptr, 0, &nl),
                  get());
    MsdCurve out;
    out.lag_time.resize(nl);
    out.total.resize(nl);
    out.axial.resize(nl);
    out.lateral.resize(nl);
    out.n_origins.resize(nl);
    detail::check(seams_msd_curve(get(), out.lag_time.data(), out.total.data(), out.axial.data(),
                                  out.lateral.data(), out.n_origins.data(), nl, &nl),
                  get());
    return out;
  }

  /** @brief Fit Einstein diffusion coefficients over `[t_start, t_end]`. */
  [[nodiscard]] DiffusionFit fit(double t_start = 0.0, double t_end = 0.0) const {
    DiffusionFit out;
    detail::check(
        seams_msd_fit(get(), t_start, t_end, &out.d_total, &out.d_axial, &out.d_lateral, &out.alpha),
        get());
    out.alpha_fitted = !std::isnan(out.alpha);
    return out;
  }
};

/**
 * @brief Owned streaming hydrogen-bond lifetime accumulator.
 * @ingroup seams_cpp
 */
class HbondTcf : public detail::Handle<seams_hbond_tcf, seams_hbond_tcf_destroy> {
 public:
  /**
   * @brief Allocate an H-bond lifetime accumulator.
   * @param max_lag Maximum lag in frames, or 0 for all lags.
   * @param origin_stride Frame stride between origins.
   * @param dt Frame spacing.
   * @param directed True for directed `(donor, acceptor)` edges.
   */
  HbondTcf(std::size_t max_lag, std::size_t origin_stride, double dt, bool directed = true) {
    seams_hbond_tcf *raw = nullptr;
    detail::check(seams_hbond_tcf_create(max_lag, origin_stride, dt, directed ? 1 : 0, &raw),
                  nullptr);
    *this = HbondTcf(raw);
  }

  using detail::Handle<seams_hbond_tcf, seams_hbond_tcf_destroy>::Handle;

  /** @brief Push one frame's `(donor, acceptor)` arrays. */
  void push(std::span<const int> donors, std::span<const int> acceptors) {
    if (donors.size() != acceptors.size()) {
      throw Error(SEAMS_ERR_ARG, "donors and acceptors differ in length");
    }
    detail::check(
        seams_hbond_tcf_push(get(), donors.empty() ? nullptr : donors.data(),
                             acceptors.empty() ? nullptr : acceptors.data(), donors.size()),
        get());
  }

  /** @brief Read intermittent C_HB and continuous S_HB curves. */
  [[nodiscard]] HbondLifetimeCurve curve() const {
    std::size_t nl = 0;
    detail::check(seams_hbond_tcf_curve(get(), nullptr, nullptr, nullptr, nullptr, 0, &nl), get());
    HbondLifetimeCurve out;
    out.lag_time.resize(nl);
    out.c_intermittent.resize(nl);
    out.s_continuous.resize(nl);
    out.n_origins.resize(nl);
    detail::check(seams_hbond_tcf_curve(get(), out.lag_time.data(), out.c_intermittent.data(),
                                        out.s_continuous.data(), out.n_origins.data(), nl, &nl),
                  get());
    return out;
  }

  /** @brief Fit integrated and single-exponential lifetimes. */
  [[nodiscard]] HbondLifetimeFit fit(double t_start = 0.0, double t_end = 0.0) const {
    HbondLifetimeFit out;
    detail::check(seams_hbond_tcf_fit(get(), t_start, t_end, &out.tau_c_int, &out.tau_c_exp,
                                      &out.tau_s_int, &out.tau_s_exp),
                  get());
    return out;
  }
};

/**
 * @brief Owned streaming orientational correlation accumulator.
 * @ingroup seams_cpp
 */
class Reorientation : public detail::Handle<seams_reorient, seams_reorient_destroy> {
 public:
  /**
   * @brief Allocate a reorientation accumulator.
   * @param n_vectors Molecular vector count per frame.
   * @param max_lag Maximum lag in frames, or 0 for all lags.
   * @param origin_stride Frame stride between origins.
   * @param dt Frame spacing.
   */
  Reorientation(std::size_t n_vectors, std::size_t max_lag, std::size_t origin_stride, double dt) {
    seams_reorient *raw = nullptr;
    detail::check(seams_reorient_create(n_vectors, max_lag, origin_stride, dt, &raw), nullptr);
    *this = Reorientation(raw);
  }

  using detail::Handle<seams_reorient, seams_reorient_destroy>::Handle;

  /** @brief Push one frame of molecular 3-vectors. */
  void push(std::span<const double> vectors) {
    if (vectors.size() % 3 != 0) {
      throw Error(SEAMS_ERR_ARG, "vectors length is not a multiple of 3");
    }
    detail::check(seams_reorient_push(get(), vectors.data(), vectors.size() / 3), get());
  }

  /** @brief Read C1 and C2 correlation curves. */
  [[nodiscard]] ReorientationCurve curve() const {
    std::size_t nl = 0;
    detail::check(seams_reorient_curve(get(), nullptr, nullptr, nullptr, nullptr, 0, &nl), get());
    ReorientationCurve out;
    out.lag_time.resize(nl);
    out.c1.resize(nl);
    out.c2.resize(nl);
    out.n_origins.resize(nl);
    detail::check(seams_reorient_curve(get(), out.lag_time.data(), out.c1.data(), out.c2.data(),
                                       out.n_origins.data(), nl, &nl),
                  get());
    return out;
  }

  /** @brief Fit integrated and single-exponential reorientation times. */
  [[nodiscard]] ReorientationFit fit(double t_start = 0.0, double t_end = 0.0) const {
    ReorientationFit out;
    detail::check(seams_reorient_fit(get(), t_start, t_end, &out.tau1_int, &out.tau1_exp,
                                     &out.tau2_int, &out.tau2_exp),
                  get());
    return out;
  }
};

/**
 * @brief One frame: coordinates, a cell, and the analyses that read them.
 * @ingroup seams_cpp
 *
 * Move-only. The destructor frees the C handle. Methods that read the
 * frame require exclusive use of this object. The Steinhardt scratch
 * inside the library is process-global, so two frames must not run
 * `steinhardt`, `ql`, `q4q6`, or `q6bar` at the same time.
 * `set_num_threads` is also process-global.
 */
class System : public detail::Handle<seams_system, seams_system_destroy> {
 public:
  /**
   * @brief Allocate an empty frame.
   * @throws Error with status `Memory` when the allocation fails.
   */
  System() : detail::Handle<seams_system, seams_system_destroy>(seams_system_create()) {
    if (get() == nullptr) {
      detail::check(SEAMS_ERR_NOMEM, nullptr);
    }
  }

  /**
   * @brief Store points in an orthorhombic cell whose origin is zero.
   * @param xyz Interleaved x, y, z. Length is a multiple of 3. Angstrom.
   * @param lx Edge length along x, Angstrom.
   * @param ly Edge length along y, Angstrom.
   * @param lz Edge length along z, Angstrom.
   */
  void set_ortho(std::span<const double> xyz, double lx, double ly, double lz) {
    const double origin[3] = {0, 0, 0};
    const double lengths[3] = {lx, ly, lz};
    set_ortho(xyz, origin, lengths);
  }

  /**
   * @brief Store points in an orthorhombic cell.
   * @param xyz Interleaved x, y, z. Length is a multiple of 3. Angstrom.
   * @param origin Cell origin, Angstrom.
   * @param lengths Edge lengths, Angstrom. Each must be positive.
   * @param types One species per atom, or an empty span.
   * @param ids One id per atom, or an empty span.
   */
  void set_ortho(std::span<const double> xyz, std::span<const double, 3> origin,
                 std::span<const double, 3> lengths, std::span<const int> types = {},
                 std::span<const int> ids = {}) {
    if (xyz.size() % 3 != 0) {
      throw Error(SEAMS_ERR_ARG, "coordinates are not a multiple of 3");
    }
    const auto n = static_cast<int>(xyz.size() / 3);
    detail::check(seams_system_set_ortho(get(), xyz.data(), n, origin[0], origin[1], origin[2],
                                         lengths[0], lengths[1], lengths[2],
                                         types.empty() ? nullptr : types.data(),
                                         ids.empty() ? nullptr : ids.data()),
                  get());
  }

  /**
   * @brief Store points and a general nine-number cell.
   * @param xyz Interleaved coordinates, Angstrom.
   * @param box9 `ox, oy, oz, lx, ly, lz, xy, xz, yz`.
   * @param types One species per atom, or empty.
   * @param ids One id per atom, or empty.
   */
  void set_frame(std::span<const double> xyz, std::span<const double, 9> box9,
                 std::span<const int> types = {}, std::span<const int> ids = {}) {
    if (xyz.size() % 3 != 0) {
      throw Error(SEAMS_ERR_ARG, "coordinates are not a multiple of 3");
    }
    const auto n = static_cast<int>(xyz.size() / 3);
    detail::check(seams_system_set_frame(get(), xyz.data(), n, box9.data(),
                                         types.empty() ? nullptr : types.data(),
                                         ids.empty() ? nullptr : ids.data()),
                  get());
  }

  std::pair<std::vector<int>, std::vector<int>> chill(double cutoff, int type, int graph,
                                                      int rule) {
    const int n = count();
    std::vector<int> labels(static_cast<std::size_t>(std::max(n, 0)));
    int counts[9] = {};
    detail::check(seams_chill(get(), cutoff, type, graph, rule, labels.data(),
                              static_cast<std::size_t>(std::max(n, 0)), counts),
                  get());
    return {labels, std::vector<int>(counts, counts + 9)};
  }

  std::pair<std::vector<double>, std::vector<double>> steinhardt(double cutoff, int order,
                                                                 bool average) {
    const int n = count();
    std::vector<double> ql(static_cast<std::size_t>(std::max(n, 0)));
    std::vector<double> bar;
    double *bar_ptr = nullptr;
    if (average) {
      bar.resize(static_cast<std::size_t>(std::max(n, 0)));
      bar_ptr = bar.data();
    }
    detail::check(seams_steinhardt(get(), cutoff, 0, order, ql.data(),
                                   static_cast<std::size_t>(std::max(n, 0)), bar_ptr,
                                   average ? static_cast<std::size_t>(std::max(n, 0)) : 0),
                  get());
    return {ql, bar};
  }

  /**
   * @brief Read one frame from a trajectory path.
   * @param path XYZ, CON, chemfiles, or a LAMMPS dump.
   * @param frame 1-based frame index.
   * @param type_filter Species to keep. Negative keeps every atom.
   */
  void read(std::string_view path, int frame, int type_filter) {
    const std::string text(path);
    detail::check(seams_system_read(get(), text.c_str(), frame, type_filter), get());
  }

  /**
   * @brief Number of points in the frame.
   * @return The atom count. An empty frame returns 0.
   */
  [[nodiscard]] int count() const { return seams_system_count(get()); }

  /** @brief Original CHILL labels and nine population counts. */
  [[nodiscard]] Classification chill(const ChillOptions &options) {
    return classify(false, options);
  }

  /** @brief CHILL+ labels and nine population counts. */
  [[nodiscard]] Classification chill_plus(const ChillOptions &options) {
    return classify(true, options);
  }

  /**
   * @brief Local q_l on the cutoff shell.
   * @param options Cutoff, species, degree, and whether to average.
   */
  [[nodiscard]] Steinhardt steinhardt(const SteinhardtOptions &options) {
    const auto n = static_cast<std::size_t>(count());
    Steinhardt out;
    out.local.resize(n);
    if (options.average) {
      out.average.resize(n);
    }
    detail::check(seams_steinhardt(get(), options.cutoff, options.species, options.order,
                                   out.local.data(), n,
                                   options.average ? out.average.data() : nullptr, n),
                  get());
    return out;
  }

  /** @brief Local q_l on a named graph. */
  [[nodiscard]] Steinhardt ql(const OrderOptions &options) {
    const auto n = static_cast<std::size_t>(count());
    Steinhardt out;
    out.local.resize(n);
    if (options.average) {
      out.average.resize(n);
    }
    detail::check(seams_ql(get(), options.cutoff, options.species, options.order,
                           static_cast<int>(options.graph), out.local.data(), n,
                           options.average ? out.average.data() : nullptr, n),
                  get());
    return out;
  }

  /** @brief Local q4 and q6 on the cutoff shell. */
  [[nodiscard]] Q4Q6 q4q6(double cutoff, int species = 0) {
    const auto n = static_cast<std::size_t>(count());
    Q4Q6 out;
    out.q4.resize(n);
    out.q6.resize(n);
    detail::check(seams_q4q6(get(), cutoff, species, out.q4.data(), out.q6.data()), get());
    return out;
  }

  /** @brief q6-bar, its mean and median, and the ice count above the cut. */
  [[nodiscard]] Q6Bar q6bar(const Q6BarOptions &options) {
    const auto n = static_cast<std::size_t>(count());
    Q6Bar out;
    out.values.resize(n);
    detail::check(seams_q6bar(get(), options.cutoff, options.species, options.ice_cut,
                              out.values.data(), n, &out.mean, &out.median, &out.n_ice),
                  get());
    return out;
  }

  /** @brief Cutoff or mutual neighbour list. */
  [[nodiscard]] Neighbours neighbours(const NeighbourOptions &options) {
    seams_neighbours *out = nullptr;
    detail::check(seams_neighbours_build(get(), options.cutoff, options.species,
                                         static_cast<int>(options.graph), &out),
                  get());
    return Neighbours(out);
  }

  /**
   * @brief Chordless cycles of one length on a neighbour graph.
   * @param options Graph of the cycles.
   * @param size Cycle length, from 3 through 12.
   * @return Neighbours::chordless on that graph.
   */
  [[nodiscard]] ChordlessCycles rings_chordless(const NeighbourOptions &options, int size) {
    return neighbours(options).chordless(size);
  }

  /** @brief k-nearest neighbour list. */
  [[nodiscard]] Neighbours knn(const KnnOptions &options) {
    seams_neighbours *out = nullptr;
    detail::check(seams_neighbours_knn(get(), options.k, options.cutoff, options.species,
                                       options.mutual ? 1 : 0, &out),
                  get());
    return Neighbours(out);
  }

  /**
   * @brief Degree of each atom on the cutoff graph.
   * @param cutoff Positive cutoff, Angstrom.
   * @param species Species filter. Zero keeps every atom.
   */
  [[nodiscard]] std::vector<int> neighbour_degree(double cutoff, int species = 0) {
    return neighbours(NeighbourOptions{cutoff, species, Graph::Cutoff}).degrees();
  }

  /** @brief Primitive rings up to a maximum size. */
  [[nodiscard]] Rings rings(const RingOptions &options) {
    seams_rings *out = nullptr;
    detail::check(seams_rings_find(get(), options.cutoff, options.species,
                                   static_cast<int>(options.graph), options.max_size, &out),
                  get());
    return Rings(out);
  }

  /** @brief Hexagonal and double-diamond cage membership. */
  [[nodiscard]] Cages cages(const CageOptions &options) {
    seams_cage_set *out = nullptr;
    detail::check(seams_cages(get(), options.cutoff, options.species,
                              static_cast<int>(options.graph), options.k, options.complete ? 1 : 0,
                              &out),
                  get());
    return Cages(out);
  }

  /** @brief Closed cages of one ring census. */
  [[nodiscard]] SignatureMatch signature_cages(const SignatureOptions &options) {
    SignatureMatch out;
    detail::check(seams_signature_cages(get(), options.cutoff, options.species,
                                        static_cast<int>(options.graph), options.spec.c_str(),
                                        &out.cages, &out.atoms),
                  get());
    return out;
  }

  /** @brief Partial radial distribution and the cell volume. */
  [[nodiscard]] Distribution rdf(const DistributionOptions &options) {
    std::size_t n_out = 0;
    Distribution out;
    detail::check(seams_rdf(get(), options.rmax, options.bins, options.type_i, options.type_j,
                            nullptr, nullptr, nullptr, 0, &n_out, &out.volume),
                  get());
    out.radius.resize(n_out);
    out.g.resize(n_out);
    out.count.resize(n_out);
    detail::check(seams_rdf(get(), options.rmax, options.bins, options.type_i, options.type_j,
                            out.radius.data(), out.g.data(), out.count.data(), n_out, &n_out,
                            nullptr),
                  get());
    out.radius.resize(n_out);
    out.g.resize(n_out);
    out.count.resize(n_out);
    return out;
  }

  /** @brief Running coordination number along the partial RDF. */
  [[nodiscard]] RunningCoordination running_cn(const DistributionOptions &options) {
    std::size_t n_out = 0;
    detail::check(seams_running_cn(get(), options.rmax, options.bins, options.type_i,
                                   options.type_j, nullptr, nullptr, 0, &n_out),
                  get());
    RunningCoordination out;
    out.radius.resize(n_out);
    out.coordination.resize(n_out);
    detail::check(seams_running_cn(get(), options.rmax, options.bins, options.type_i,
                                   options.type_j, out.radius.data(), out.coordination.data(),
                                   n_out, &n_out),
                  get());
    out.radius.resize(n_out);
    out.coordination.resize(n_out);
    return out;
  }

  /** @brief Coordination number integrated to rmax. */
  [[nodiscard]] double coordination(const DistributionOptions &options) {
    double cn = 0;
    detail::check(seams_coordination(get(), options.rmax, options.bins, options.type_i,
                                     options.type_j, &cn),
                  get());
    return cn;
  }

  /** @brief Hydrogen-bond count. */
  [[nodiscard]] int hbonds(const BondOptions &options) {
    int n_bonds = 0;
    detail::check(seams_hbonds(get(), options.cutoff, options.oxygen, options.hydrogen,
                               options.distance, options.angle, &n_bonds),
                  get());
    return n_bonds;
  }

  /**
   * @brief Store heavies and a separate hydrogen cloud.
   *
   * An empty id or mol span uses the C API default for that cloud.
   */
  void set_hydrogen_frame(std::span<const double> heavy_xyz, std::span<const int> heavy_types,
                          std::span<const int> heavy_mols, std::span<const double> hydrogen_xyz,
                          std::span<const int> hydrogen_types, std::span<const int> hydrogen_mols,
                          std::span<const double, 9> box9, std::span<const int> heavy_ids,
                          std::span<const int> hydrogen_ids) {
    if (heavy_xyz.size() % 3 != 0 || hydrogen_xyz.size() % 3 != 0) {
      throw Error(SEAMS_ERR_ARG, "coordinates are not a multiple of 3");
    }
    const auto n_heavy = heavy_xyz.size() / 3;
    const auto n_hydrogen = hydrogen_xyz.size() / 3;
    if (heavy_types.size() != n_heavy || hydrogen_types.size() != n_hydrogen) {
      throw Error(SEAMS_ERR_ARG, "one type per atom");
    }
    if ((!heavy_mols.empty() && heavy_mols.size() != n_heavy) ||
        (!hydrogen_mols.empty() && hydrogen_mols.size() != n_hydrogen) ||
        (!heavy_ids.empty() && heavy_ids.size() != n_heavy) ||
        (!hydrogen_ids.empty() && hydrogen_ids.size() != n_hydrogen)) {
      throw Error(SEAMS_ERR_ARG, "one id and one molID per atom");
    }
    detail::check(seams_system_set_hydrogen_frame(
                      get(), heavy_xyz.data(), n_heavy, heavy_types.data(),
                      heavy_ids.empty() ? nullptr : heavy_ids.data(),
                      heavy_mols.empty() ? nullptr : heavy_mols.data(),
                      n_hydrogen == 0 ? nullptr : hydrogen_xyz.data(), n_hydrogen,
                      n_hydrogen == 0 ? nullptr : hydrogen_types.data(),
                      hydrogen_ids.empty() ? nullptr : hydrogen_ids.data(),
                      hydrogen_mols.empty() ? nullptr : hydrogen_mols.data(), box9.data(), 9),
                  get());
  }

  /**
   * @brief Directed hydrogen bonds for one donor type set and one acceptor type set.
   */
  [[nodiscard]] std::vector<DirectedBond> directed_hbonds(std::span<const int> donor_types,
                                                         std::span<const int> acceptor_types,
                                                         const DirectedBondOptions &options) {
    if (!(options.cutoff > 0.0) || !(options.distance > 0.0) || !(options.angle > 0.0)) {
      throw Error(SEAMS_ERR_RANGE, "cutoffs must be positive");
    }
    std::size_t n_out = 0;
    detail::check(seams_directed_hbonds(get(), options.cutoff, options.distance, options.angle,
                                        donor_types.empty() ? nullptr : donor_types.data(),
                                        donor_types.size(),
                                        acceptor_types.empty() ? nullptr : acceptor_types.data(),
                                        acceptor_types.size(), options.hydrogen,
                                        options.exclude_same_mol ? 1 : 0, nullptr, 0, &n_out),
                  get());
    std::vector<seams_directed_hbond> rows(n_out);
    for (auto &row : rows) {
      row.size = sizeof(seams_directed_hbond);
    }
    if (n_out > 0) {
      detail::check(seams_directed_hbonds(get(), options.cutoff, options.distance, options.angle,
                                          donor_types.empty() ? nullptr : donor_types.data(),
                                          donor_types.size(),
                                          acceptor_types.empty() ? nullptr : acceptor_types.data(),
                                          acceptor_types.size(), options.hydrogen,
                                          options.exclude_same_mol ? 1 : 0, rows.data(), n_out,
                                          &n_out),
                    get());
    }
    std::vector<DirectedBond> out;
    out.reserve(n_out);
    for (const auto &row : rows) {
      out.push_back(DirectedBond{static_cast<int>(row.donor), static_cast<int>(row.hydrogen),
                                 static_cast<int>(row.acceptor), row.d_AH, row.d_DA, row.angle});
    }
    return out;
  }

  /** @brief Number density along one cell axis. */
  [[nodiscard]] DensityProfile density(const DensityOptions &options) {
    std::size_t n_out = 0;
    detail::check(seams_density_axis(get(), options.species, options.bins,
                                     static_cast<int>(options.axis), nullptr, nullptr, 0, &n_out),
                  get());
    DensityProfile out;
    out.coordinate.resize(n_out);
    out.density.resize(n_out);
    detail::check(seams_density_axis(get(), options.species, options.bins,
                                     static_cast<int>(options.axis), out.coordinate.data(),
                                     out.density.data(), n_out, &n_out),
                  get());
    out.coordinate.resize(n_out);
    out.density.resize(n_out);
    return out;
  }

  /** @brief Number and mass density along a cell axis or a Cartesian direction. */
  [[nodiscard]] FractionalProfile fractional_density(const FractionalDensityOptions &options) {
    const double *direction = options.use_direction ? options.direction.data() : nullptr;
    const double *masses = options.masses.empty() ? nullptr : options.masses.data();
    const std::size_t masses_len = options.masses.empty() ? 0 : options.masses.size();
    std::size_t n_out = 0;
    detail::check(seams_density_fractional(get(), options.species, options.bins,
                                           static_cast<int>(options.axis), direction,
                                           direction == nullptr ? 0 : 3, options.origin, masses,
                                           masses_len, nullptr, nullptr, nullptr, nullptr, 0,
                                           &n_out),
                  get());
    FractionalProfile out;
    out.coordinate.resize(n_out);
    out.number_density.resize(n_out);
    out.mass_density.resize(n_out);
    out.volume.resize(n_out);
    detail::check(seams_density_fractional(get(), options.species, options.bins,
                                           static_cast<int>(options.axis), direction,
                                           direction == nullptr ? 0 : 3, options.origin, masses,
                                           masses_len, out.coordinate.data(),
                                           out.number_density.data(), out.mass_density.data(),
                                           out.volume.data(), n_out, &n_out),
                  get());
    out.coordinate.resize(n_out);
    out.number_density.resize(n_out);
    out.mass_density.resize(n_out);
    out.volume.resize(n_out);
    return out;
  }

  /** @brief Pore axis and centre, supplied or taken from the framework. */
  [[nodiscard]] PoreFrame pore_frame(const PoreFrameOptions &options) {
    const double *direction = options.use_direction ? options.direction.data() : nullptr;
    const double *centre = options.use_centre ? options.centre_fractional.data() : nullptr;
    PoreFrame out;
    int index = 0;
    detail::check(seams_pore_frame(get(), options.framework, options.axis, direction,
                                   direction == nullptr ? 0 : 3, centre, centre == nullptr ? 0 : 3,
                                   out.axis.data(), out.centre.data(), out.centre_fractional.data(),
                                   &index),
                  get());
    out.axis_index = index;
    return out;
  }

  /** @brief Radial, azimuthal, and axial densities about a pore axis. */
  [[nodiscard]] PoreProfiles pore_profiles(const PoreProfileOptions &options) {
    const double *direction = options.use_direction ? options.direction.data() : nullptr;
    const double *centre = options.use_centre ? options.centre_fractional.data() : nullptr;
    const double *masses = options.masses.empty() ? nullptr : options.masses.data();
    const std::size_t masses_len = options.masses.empty() ? 0 : options.masses.size();
    std::size_t n_radial = 0;
    std::size_t n_azimuthal = 0;
    std::size_t n_axial = 0;
    detail::check(seams_pore_profiles(get(), options.species, options.radial_bins,
                                      options.azimuthal_bins, options.axial_bins, options.r_max,
                                      options.axis, direction, direction == nullptr ? 0 : 3, centre,
                                      centre == nullptr ? 0 : 3, options.framework, options.origin,
                                      masses, masses_len, nullptr, nullptr, nullptr, nullptr, 0,
                                      &n_radial, nullptr, nullptr, nullptr, nullptr, 0, &n_azimuthal,
                                      nullptr, nullptr, nullptr, nullptr, 0, &n_axial),
                  get());
    PoreProfiles out;
    out.radial.coordinate.resize(n_radial);
    out.radial.number_density.resize(n_radial);
    out.radial.mass_density.resize(n_radial);
    out.radial.volume.resize(n_radial);
    out.azimuthal.coordinate.resize(n_azimuthal);
    out.azimuthal.number_density.resize(n_azimuthal);
    out.azimuthal.mass_density.resize(n_azimuthal);
    out.azimuthal.volume.resize(n_azimuthal);
    out.axial.coordinate.resize(n_axial);
    out.axial.number_density.resize(n_axial);
    out.axial.mass_density.resize(n_axial);
    out.axial.volume.resize(n_axial);
    detail::check(
        seams_pore_profiles(
            get(), options.species, options.radial_bins, options.azimuthal_bins,
            options.axial_bins, options.r_max, options.axis, direction,
            direction == nullptr ? 0 : 3, centre, centre == nullptr ? 0 : 3, options.framework,
            options.origin, masses, masses_len, out.radial.coordinate.data(),
            out.radial.number_density.data(), out.radial.mass_density.data(),
            out.radial.volume.data(), n_radial, &n_radial, out.azimuthal.coordinate.data(),
            out.azimuthal.number_density.data(), out.azimuthal.mass_density.data(),
            out.azimuthal.volume.data(), n_azimuthal, &n_azimuthal, out.axial.coordinate.data(),
            out.axial.number_density.data(), out.axial.mass_density.data(), out.axial.volume.data(),
            n_axial, &n_axial),
        get());
    return out;
  }

  /** @brief Nearest framework site and contact counts for each guest. */
  [[nodiscard]] HostGuestContacts host_guest(const HostGuestOptions &options) {
    std::size_t n_guest = 0;
    std::size_t n_types = 0;
    detail::check(seams_host_guest(get(), options.guest, options.framework, options.cutoff, nullptr,
                                   nullptr, nullptr, nullptr, 0, &n_guest, nullptr, 0, &n_types,
                                   nullptr, 0),
                  get());
    HostGuestContacts out;
    out.guest_index.resize(n_guest);
    out.nearest_index.resize(n_guest);
    out.nearest_type.resize(n_guest);
    out.nearest_distance.resize(n_guest);
    out.framework_types.resize(n_types);
    out.contacts.resize(n_guest * n_types);
    detail::check(seams_host_guest(get(), options.guest, options.framework, options.cutoff,
                                   out.guest_index.data(), out.nearest_index.data(),
                                   out.nearest_type.data(), out.nearest_distance.data(), n_guest,
                                   &n_guest, out.framework_types.data(), n_types, &n_types,
                                   out.contacts.data(), out.contacts.size()),
                  get());
    out.guest_index.resize(n_guest);
    out.nearest_index.resize(n_guest);
    out.nearest_type.resize(n_guest);
    out.nearest_distance.resize(n_guest);
    out.framework_types.resize(n_types);
    out.contacts.resize(n_guest * n_types);
    return out;
  }

  /** @brief Topological fingerprint of the bonded graph. */
  [[nodiscard]] Fingerprint fingerprint(const FingerprintOptions &options) {
    seams_topo *out = nullptr;
    detail::check(seams_fingerprint(get(), options.cutoff, options.species, options.k,
                                    static_cast<int>(options.graph), options.hops,
                                    options.color_types ? 1 : 0, &out),
                  get());
    return Fingerprint(out);
  }

  /** @brief Mutual nearest unlike contacts. `spec` is site-table text. */
  [[nodiscard]] Contacts contact_pairs(std::string_view spec) {
    const std::string text(spec);
    Contacts out;
    detail::check(seams_contact_pairs(get(), text.c_str(), &out.pairs, &out.cations, &out.anions),
                  get());
    return out;
  }

  /** @brief Largest polar or apolar domain. */
  [[nodiscard]] DomainReport domains(const DomainOptions &options) {
    DomainReport out;
    detail::check(seams_domains(get(), options.cutoff, options.spec.c_str(),
                                static_cast<int>(options.subset), &out.subset, &out.largest,
                                &out.percolation),
                  get());
    return out;
  }

  /** @brief Perfect and imperfect prism counts by ring size. */
  [[nodiscard]] PrismCensus prisms(const PrismOptions &options) {
    PrismCensus out;
    detail::check(seams_prisms(get(), options.cutoff, options.species, out.perfect.data(),
                               out.imperfect.data(), out.perfect.size()),
                  get());
    return out;
  }

  /**
   * @brief Primitive-ring histogram of a monolayer.
   * @return `counts[k]` is the number of rings of size k.
   */
  [[nodiscard]] std::vector<int> monolayer(const SheetOptions &options) {
    std::vector<int> counts(static_cast<std::size_t>(options.max_size + 1));
    detail::check(seams_monolayer(get(), options.cutoff, options.species, options.max_size,
                                  options.sheet_area, counts.data(), counts.size()),
                  get());
    return counts;
  }

  /**
   * @brief Per-atom site kind.
   * @param spec Site-table text, for example `1=waterO`.
   * @return One kind code per atom.
   */
  [[nodiscard]] std::vector<int> site_table(std::string_view spec) {
    const std::string text(spec);
    std::vector<int> kinds(static_cast<std::size_t>(count()));
    detail::check(seams_site_table(get(), text.c_str(), kinds.data(), kinds.size()), get());
    return kinds;
  }

  /**
   * @brief Per-atom ice-shell label.
   * @param shell_rmsd Positive residual ceiling, or a non-positive value for the library ceiling.
   * @return 0 ice Ih, 1 ice Ic, 2 stacking disorder, or -1 when the shell does not match.
   */
  [[nodiscard]] std::vector<int> symmetry(double shell_rmsd) {
    std::vector<int> labels(static_cast<std::size_t>(count()));
    detail::check(seams_system_symmetry(get(), shell_rmsd, labels.data(), labels.size()), get());
    return labels;
  }

  /**
   * @brief Planar ice front from an ice assignment.
   * @param ice Non-zero when the atom is ice. Length must equal `count()`.
   * @param options Water species, ion species, bins, and axis.
   */
  [[nodiscard]] Brine brine(std::span<const int> ice, const BrineOptions &options) {
    seams_brine *out = nullptr;
    detail::check(seams_brine_analyse(get(), ice.data(), ice.size(), options.water_type,
                                      options.ion_type, options.bins,
                                      static_cast<int>(options.axis), &out),
                  get());
    return Brine(out);
  }

  /** @brief Copy the stored molecule ids out. */
  [[nodiscard]] std::vector<int> mols() const {
    const auto n = static_cast<std::size_t>(count());
    std::vector<int> out(n);
    detail::check(seams_mols(get(), out.data(), out.size()), get());
    return out;
  }

  /** @brief Replace the stored molecule ids on a loaded frame. */
  void set_mols(std::span<const int> mol_ids) {
    detail::check(seams_set_mols(get(), mol_ids.data(), mol_ids.size()), get());
  }

  /** @brief Intact periodic centres of molecules in this frame. */
  [[nodiscard]] MoleculeCentres molecule_centres(std::span<const int> mol_ids = {},
                                                 std::span<const double> masses = {}) const {
    std::size_t nm = 0;
    detail::check(seams_molecule_centres(get(), mol_ids.empty() ? nullptr : mol_ids.data(),
                                         mol_ids.size(), masses.empty() ? nullptr : masses.data(),
                                         masses.size(), nullptr, nullptr, 0, &nm),
                  get());
    MoleculeCentres out;
    out.mol_ids.resize(nm);
    out.centres.resize(3 * nm);
    detail::check(seams_molecule_centres(get(), mol_ids.empty() ? nullptr : mol_ids.data(),
                                         mol_ids.size(), masses.empty() ? nullptr : masses.data(),
                                         masses.size(), out.mol_ids.data(), out.centres.data(), nm,
                                         &nm),
                  get());
    return out;
  }

  /** @brief Water O-to-H-midpoint and H-to-H vectors under minimum image. */
  [[nodiscard]] WaterVectors water_vectors(int oxygen_type, int hydrogen_type) const {
    std::size_t nw = 0;
    detail::check(seams_water_vectors(get(), oxygen_type, hydrogen_type, nullptr, nullptr, nullptr,
                                      nullptr, nullptr, 0, &nw),
                  get());
    WaterVectors out;
    out.mol_ids.resize(nw);
    out.oh_mid.resize(3 * nw);
    out.bisector.resize(3 * nw);
    out.hh.resize(3 * nw);
    out.hh_unit.resize(3 * nw);
    detail::check(seams_water_vectors(get(), oxygen_type, hydrogen_type, out.mol_ids.data(),
                                      out.oh_mid.data(), out.bisector.data(), out.hh.data(),
                                      out.hh_unit.data(), nw, &nw),
                  get());
    return out;
  }

  /** @brief Per-molecule dipole or geometric orientation vectors. */
  [[nodiscard]] MoleculeDipoles molecule_dipoles(std::span<const int> mol_ids = {},
                                                 std::span<const double> charges = {}) const {
    std::size_t nm = 0;
    detail::check(seams_molecule_dipoles(get(), mol_ids.empty() ? nullptr : mol_ids.data(),
                                         mol_ids.size(),
                                         charges.empty() ? nullptr : charges.data(),
                                         charges.size(), nullptr, nullptr, nullptr, 0, &nm),
                  get());
    MoleculeDipoles out;
    out.mol_ids.resize(nm);
    out.dipoles.resize(3 * nm);
    out.unit_vectors.resize(3 * nm);
    detail::check(seams_molecule_dipoles(get(), mol_ids.empty() ? nullptr : mol_ids.data(),
                                         mol_ids.size(),
                                         charges.empty() ? nullptr : charges.data(),
                                         charges.size(), out.mol_ids.data(), out.dipoles.data(),
                                         out.unit_vectors.data(), nm, &nm),
                  get());
    return out;
  }

 private:
  [[nodiscard]] Classification classify(bool plus, const ChillOptions &options) {
    const auto n = static_cast<std::size_t>(count());
    Classification out;
    out.labels.resize(n);
    if (plus) {
      detail::check(seams_chill_plus(get(), options.cutoff, options.species,
                                     static_cast<int>(options.graph),
                                     static_cast<int>(options.rule), out.labels.data(),
                                     out.counts.data()),
                    get());
    } else {
      detail::check(seams_chill(get(), options.cutoff, options.species,
                                static_cast<int>(options.graph), static_cast<int>(options.rule),
                                out.labels.data(), n, out.counts.data()),
                    get());
    }
    return out;
  }
};

/**
 * @brief Library version text of the built library.
 * @ingroup seams_cpp
 * @return A copy of the static version string.
 */
[[nodiscard]] inline std::string version() {
  const char *text = seams_version_string();
  return text == nullptr ? std::string() : std::string(text);
}

/**
 * @brief ABI number compiled into this library.
 * @ingroup seams_cpp
 * @return `SEAMS_ABI`, which is 1.
 */
[[nodiscard]] inline int abi() { return seams_abi(); }

/**
 * @brief Set the OpenMP width for the next call.
 * @ingroup seams_cpp
 * @param n Thread count. Zero restores the library default.
 */
inline void set_num_threads(int n) { detail::check(seams_set_num_threads(n), nullptr); }

/**
 * @brief OpenMP width the next call will use.
 * @ingroup seams_cpp
 * @return A positive thread count.
 */
[[nodiscard]] inline int num_threads() { return seams_get_num_threads(); }

/**
 * @brief Count frames in a trajectory path.
 * @ingroup seams_cpp
 * @param path Trajectory path.
 * @return The frame count.
 */
[[nodiscard]] inline std::size_t nframes(std::string_view path) {
  const std::string text(path);
  std::size_t count = 0;
  detail::check(seams_path_nframes(text.c_str(), &count), nullptr);
  return count;
}

/**
 * @brief Short name of a CHILL+ label.
 * @ingroup seams_cpp
 * @param label Label code.
 * @return A static name. The caller does not free it.
 */
[[nodiscard]] inline std::string_view label_name(Label label) {
  const char *text = seams_label_name(static_cast<int>(label));
  return text == nullptr ? std::string_view() : std::string_view(text);
}

/**
 * @brief Short name of a ring-symmetry code.
 * @ingroup seams_cpp
 * @param symmetry Symmetry code.
 * @return A static name.
 */
[[nodiscard]] inline std::string_view ring_symmetry_name(RingSymmetry symmetry) {
  const char *text = seams_ring_symmetry_name(static_cast<int>(symmetry));
  return text == nullptr ? std::string_view() : std::string_view(text);
}

/**
 * @brief Unwrap periodic coordinates across frames.
 * @ingroup seams_cpp
 * @throws Error `SEAMS_ERR_RANGE` when `xyz.size()` is not `nframes * 3 * natoms`,
 *         that product overflows `size_t`, or `boxes.size()` is not a multiple of 9.
 */
[[nodiscard]] inline std::vector<double> unwrap_positions(std::span<const double> xyz,
                                                          std::size_t nframes, std::size_t natoms,
                                                          std::span<const double> boxes) {
  const std::size_t total = detail::trajectory_coordinates(xyz, nframes, natoms);
  if (boxes.size() % 9 != 0) {
    throw Error(SEAMS_ERR_RANGE, "boxes length is not a multiple of 9");
  }
  std::vector<double> out(total);
  detail::check(seams_traj_unwrap_positions(xyz.data(), nframes, natoms, boxes.data(),
                                            boxes.size() / 9, out.data(), out.size()),
                nullptr);
  return out;
}

/**
 * @brief Assemble and unwrap periodic molecule centres across frames.
 * @ingroup seams_cpp
 * @throws Error `SEAMS_ERR_RANGE` when `xyz.size()` is not `nframes * 3 * natoms`,
 *         that product overflows `size_t`, or `boxes.size()` is not a multiple of 9.
 */
[[nodiscard]] inline MoleculeCentres unwrap_centres(std::span<const double> xyz,
                                                    std::size_t nframes, std::size_t natoms,
                                                    std::span<const double> boxes,
                                                    std::span<const int> mols,
                                                    std::span<const double> masses = {}) {
  static_cast<void>(detail::trajectory_coordinates(xyz, nframes, natoms));
  if (boxes.size() % 9 != 0) {
    throw Error(SEAMS_ERR_RANGE, "boxes length is not a multiple of 9");
  }
  const std::size_t nboxes = boxes.size() / 9;
  std::size_t nm = 0;
  detail::check(seams_traj_unwrap_centres(xyz.data(), nframes, natoms, boxes.data(), nboxes,
                                          mols.data(), mols.size(),
                                          masses.empty() ? nullptr : masses.data(), masses.size(),
                                          nullptr, 0, nullptr, 0, &nm),
                nullptr);
  MoleculeCentres out;
  out.mol_ids.resize(nm);
  out.centres.resize(nframes * 3 * nm);
  detail::check(seams_traj_unwrap_centres(
                    xyz.data(), nframes, natoms, boxes.data(), nboxes, mols.data(), mols.size(),
                    masses.empty() ? nullptr : masses.data(), masses.size(), out.mol_ids.data(), nm,
                    out.centres.data(), out.centres.size(), &nm),
                nullptr);
  return out;
}

/**
 * @brief Axis-projected orientational order parameters P1 and P2.
 * @ingroup seams_cpp
 */
[[nodiscard]] inline OrientationOrder orientation_order(std::span<const double> vectors,
                                                        std::span<const double, 3> axis) {
  if (vectors.size() % 3 != 0) {
    throw Error(SEAMS_ERR_ARG, "vectors length is not a multiple of 3");
  }
  const std::size_t nv = vectors.size() / 3;
  OrientationOrder out;
  out.cos_theta.resize(nv);
  out.p2.resize(nv);
  detail::check(seams_orientation_order(nv == 0 ? nullptr : vectors.data(), nv, axis.data(),
                                        nv == 0 ? nullptr : out.cos_theta.data(),
                                        nv == 0 ? nullptr : out.p2.data(), nv, &out.mean_p1,
                                        &out.mean_p2),
                nullptr);
  return out;
}

}  // namespace seams

#endif
