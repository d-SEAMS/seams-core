/* Copyright (c) d-SEAMS developers */
#ifndef SEAMS_HPP_
#define SEAMS_HPP_

/**
 * Header-only C++ wrappers for include/seams.h.
 * A failed status throws seams::Error. The message is seams_handle_error.
 */

#include "seams.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace seams {

class Error : public std::runtime_error {
 public:
  int code;
  Error(int status, std::string message)
      : std::runtime_error(std::move(message)), code(status) {}
};

inline void check(int status, const void *handle) {
  if (status == SEAMS_OK) {
    return;
  }
  const char *text = seams_handle_error(handle);
  throw Error(status, text == nullptr ? std::string() : std::string(text));
}

class Neighbours {
 public:
  explicit Neighbours(seams_neighbours *ptr) : ptr_(ptr) {}
  Neighbours(const Neighbours &) = delete;
  Neighbours &operator=(const Neighbours &) = delete;
  Neighbours(Neighbours &&other) noexcept : ptr_(other.ptr_) { other.ptr_ = nullptr; }
  Neighbours &operator=(Neighbours &&other) noexcept {
    if (this != &other) {
      reset();
      ptr_ = other.ptr_;
      other.ptr_ = nullptr;
    }
    return *this;
  }
  ~Neighbours() { reset(); }

  std::size_t nnz() const {
    std::size_t n_offsets = 0;
    std::size_t count = 0;
    check(seams_neighbours_copy(ptr_, nullptr, 0, &n_offsets, nullptr, 0, &count), ptr_);
    return count;
  }

 private:
  void reset() {
    if (ptr_ != nullptr) {
      seams_neighbours_destroy(ptr_);
      ptr_ = nullptr;
    }
  }
  seams_neighbours *ptr_ = nullptr;
};

class Rings {
 public:
  explicit Rings(seams_rings *ptr) : ptr_(ptr) {}
  Rings(const Rings &) = delete;
  Rings &operator=(const Rings &) = delete;
  Rings(Rings &&other) noexcept : ptr_(other.ptr_) { other.ptr_ = nullptr; }
  Rings &operator=(Rings &&other) noexcept {
    if (this != &other) {
      reset();
      ptr_ = other.ptr_;
      other.ptr_ = nullptr;
    }
    return *this;
  }
  ~Rings() { reset(); }

  std::vector<int> vertices(int index) const {
    std::size_t n_out = 0;
    check(seams_rings_vertices(ptr_, index, nullptr, 0, &n_out), ptr_);
    std::vector<int> atoms(n_out);
    check(seams_rings_vertices(ptr_, index, atoms.data(), n_out, &n_out), ptr_);
    atoms.resize(n_out);
    return atoms;
  }

 private:
  void reset() {
    if (ptr_ != nullptr) {
      seams_rings_destroy(ptr_);
      ptr_ = nullptr;
    }
  }
  seams_rings *ptr_ = nullptr;
};

class Cages {
 public:
  explicit Cages(seams_cage_set *ptr) : ptr_(ptr) {}
  Cages(const Cages &) = delete;
  Cages &operator=(const Cages &) = delete;
  Cages(Cages &&other) noexcept : ptr_(other.ptr_) { other.ptr_ = nullptr; }
  Cages &operator=(Cages &&other) noexcept {
    if (this != &other) {
      reset();
      ptr_ = other.ptr_;
      other.ptr_ = nullptr;
    }
    return *this;
  }
  ~Cages() { reset(); }

  void counts(int &n_hex, int &n_cubic, int &n_water, int &n_both) const {
    check(seams_cages_counts(ptr_, &n_hex, &n_cubic, &n_water, &n_both), ptr_);
  }

 private:
  void reset() {
    if (ptr_ != nullptr) {
      seams_cages_destroy(ptr_);
      ptr_ = nullptr;
    }
  }
  seams_cage_set *ptr_ = nullptr;
};

class Fingerprint {
 public:
  explicit Fingerprint(seams_topo *ptr) : ptr_(ptr) {}
  Fingerprint(const Fingerprint &) = delete;
  Fingerprint &operator=(const Fingerprint &) = delete;
  Fingerprint(Fingerprint &&other) noexcept : ptr_(other.ptr_) { other.ptr_ = nullptr; }
  Fingerprint &operator=(Fingerprint &&other) noexcept {
    if (this != &other) {
      reset();
      ptr_ = other.ptr_;
      other.ptr_ = nullptr;
    }
    return *this;
  }
  ~Fingerprint() { reset(); }

  std::string key() const {
    std::size_t n_out = 0;
    check(seams_fingerprint_key(ptr_, nullptr, 0, &n_out), ptr_);
    std::string text(n_out + 1, '\0');
    check(seams_fingerprint_key(ptr_, text.data(), text.size(), &n_out), ptr_);
    text.resize(n_out);
    return text;
  }

  std::vector<int> rings() const {
    std::vector<int> counts(13);
    check(seams_fingerprint_rings(ptr_, counts.data(), counts.size()), ptr_);
    return counts;
  }

 private:
  void reset() {
    if (ptr_ != nullptr) {
      seams_fingerprint_destroy(ptr_);
      ptr_ = nullptr;
    }
  }
  seams_topo *ptr_ = nullptr;
};

class System {
 public:
  System() : ptr_(seams_system_create()) {
    if (ptr_ == nullptr) {
      check(SEAMS_ERR_NOMEM, nullptr);
    }
  }
  System(const System &) = delete;
  System &operator=(const System &) = delete;
  System(System &&other) noexcept : ptr_(other.ptr_) { other.ptr_ = nullptr; }
  System &operator=(System &&other) noexcept {
    if (this != &other) {
      reset();
      ptr_ = other.ptr_;
      other.ptr_ = nullptr;
    }
    return *this;
  }
  ~System() { reset(); }

  void read(const std::string &path, int frame, int type_filter) {
    check(seams_system_read(ptr_, path.c_str(), frame, type_filter), ptr_);
  }

  int count() const { return seams_system_count(ptr_); }

  std::pair<std::vector<int>, std::vector<int>> chill(double cutoff, int type, int graph, int rule) {
    const int n = count();
    std::vector<int> labels(static_cast<std::size_t>(n));
    int counts[9] = {};
    check(seams_chill(ptr_, cutoff, type, graph, rule, labels.data(),
                      static_cast<std::size_t>(n), counts),
          ptr_);
    return {labels, std::vector<int>(counts, counts + 9)};
  }

  std::pair<std::vector<double>, std::vector<double>> steinhardt(double cutoff, int order, bool average) {
    const int n = count();
    std::vector<double> ql(static_cast<std::size_t>(n));
    std::vector<double> bar;
    double *bar_ptr = nullptr;
    if (average) {
      bar.resize(static_cast<std::size_t>(n));
      bar_ptr = bar.data();
    }
    check(seams_steinhardt(ptr_, cutoff, 0, order, ql.data(), static_cast<std::size_t>(n),
                           bar_ptr, average ? static_cast<std::size_t>(n) : 0),
          ptr_);
    return {ql, bar};
  }

  Cages cages(double cutoff, int graph, int k, bool complete) {
    seams_cage_set *out = nullptr;
    check(seams_cages(ptr_, cutoff, 0, graph, k, complete ? 1 : 0, &out), ptr_);
    return Cages(out);
  }

  std::pair<std::size_t, std::size_t> signature_cages(double cutoff, int graph, const std::string &spec) {
    std::size_t n_cages = 0;
    std::size_t n_atoms = 0;
    check(seams_signature_cages(ptr_, cutoff, 0, graph, spec.c_str(), &n_cages, &n_atoms), ptr_);
    return {n_cages, n_atoms};
  }

  Neighbours knn(int k, double cutoff, bool mutual) {
    seams_neighbours *out = nullptr;
    check(seams_neighbours_knn(ptr_, k, cutoff, 0, mutual ? 1 : 0, &out), ptr_);
    return Neighbours(out);
  }

  std::tuple<std::vector<double>, std::vector<double>, double> rdf(double rmax, int bins) {
    std::size_t n_out = 0;
    double volume = 0.0;
    check(seams_rdf(ptr_, rmax, bins, 0, 0, nullptr, nullptr, nullptr, 0, &n_out, &volume), ptr_);
    std::vector<double> radius(n_out);
    std::vector<double> g(n_out);
    check(seams_rdf(ptr_, rmax, bins, 0, 0, radius.data(), g.data(), nullptr, n_out, &n_out, nullptr),
          ptr_);
    radius.resize(n_out);
    g.resize(n_out);
    return {radius, g, volume};
  }

  std::vector<double> running_cn(double rmax, int bins) {
    std::size_t n_out = 0;
    check(seams_running_cn(ptr_, rmax, bins, 0, 0, nullptr, nullptr, 0, &n_out), ptr_);
    std::vector<double> radius(n_out);
    std::vector<double> cn(n_out);
    check(seams_running_cn(ptr_, rmax, bins, 0, 0, radius.data(), cn.data(), n_out, &n_out), ptr_);
    cn.resize(n_out);
    return cn;
  }

  double coordination(double rmax, int bins) {
    double cn = 0.0;
    check(seams_coordination(ptr_, rmax, bins, 0, 0, &cn), ptr_);
    return cn;
  }

  int hbonds(double cutoff, int oxygen_type, int hydrogen_type, double dist, double angle) {
    int n_bonds = 0;
    check(seams_hbonds(ptr_, cutoff, oxygen_type, hydrogen_type, dist, angle, &n_bonds), ptr_);
    return n_bonds;
  }

  std::pair<std::vector<double>, std::vector<double>> density_axis(int type, int bins, int axis) {
    std::size_t n_out = 0;
    check(seams_density_axis(ptr_, type, bins, axis, nullptr, nullptr, 0, &n_out), ptr_);
    std::vector<double> z(n_out);
    std::vector<double> rho(n_out);
    check(seams_density_axis(ptr_, type, bins, axis, z.data(), rho.data(), n_out, &n_out), ptr_);
    rho.resize(n_out);
    return {z, rho};
  }

  Fingerprint fingerprint(double cutoff, int hops) {
    seams_topo *out = nullptr;
    check(seams_fingerprint(ptr_, cutoff, 0, 4, SEAMS_GRAPH_CUTOFF, hops, 0, &out), ptr_);
    return Fingerprint(out);
  }

  std::tuple<std::size_t, int, int> contact_pairs(const std::string &spec) {
    std::size_t n_pairs = 0;
    int n_cation = 0;
    int n_anion = 0;
    check(seams_contact_pairs(ptr_, spec.c_str(), &n_pairs, &n_cation, &n_anion), ptr_);
    return {n_pairs, n_cation, n_anion};
  }

  std::tuple<int, int, double> domains(double cutoff, const std::string &spec, int subset) {
    int n_subset = 0;
    int largest = 0;
    double percolation = 0.0;
    check(seams_domains(ptr_, cutoff, spec.c_str(), subset, &n_subset, &largest, &percolation), ptr_);
    return {n_subset, largest, percolation};
  }

  std::vector<int> prisms(double cutoff) {
    std::vector<int> perfect(7);
    std::vector<int> imperfect(7);
    check(seams_prisms(ptr_, cutoff, 0, perfect.data(), imperfect.data(), perfect.size()), ptr_);
    return perfect;
  }

  std::vector<int> monolayer(double cutoff, int max_depth, double sheet_area) {
    std::vector<int> counts(static_cast<std::size_t>(max_depth + 1));
    check(seams_monolayer(ptr_, cutoff, 0, max_depth, sheet_area, counts.data(), counts.size()), ptr_);
    return counts;
  }

  Rings rings(double cutoff, int max_depth) {
    seams_rings *out = nullptr;
    check(seams_rings_find(ptr_, cutoff, 0, SEAMS_GRAPH_CUTOFF, max_depth, &out), ptr_);
    return Rings(out);
  }

  std::vector<int> site_table(const std::string &spec) {
    std::vector<int> kinds(static_cast<std::size_t>(count()));
    check(seams_site_table(ptr_, spec.c_str(), kinds.data(), kinds.size()), ptr_);
    return kinds;
  }

 private:
  void reset() {
    if (ptr_ != nullptr) {
      seams_system_destroy(ptr_);
      ptr_ = nullptr;
    }
  }
  seams_system *ptr_ = nullptr;
};

}  // namespace seams

#endif
