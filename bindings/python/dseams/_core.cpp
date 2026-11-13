// Copyright (c) d-SEAMS developers
// Nanobind front. Arrays are the caller's buffers; the C API reads them in place.

#include "seams.h"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <stdexcept>
#include <string>

namespace nb = nanobind;

namespace {

void check(seams_status code) {
  if (code == SEAMS_OK) {
    return;
  }
  throw std::runtime_error(std::string(seams_status_string(code)) + ": " +
                           seams_last_error());
}

void box_of(const nb::ndarray<nb::numpy, double, nb::ndim<1>, nb::c_contig> &box,
            double out[9]) {
  for (int i = 0; i < 9; ++i) {
    out[i] = 0.0;
  }
  if (box.shape(0) == 3) {
    out[3] = box(0);
    out[4] = box(1);
    out[5] = box(2);
    return;
  }
  if (box.shape(0) == 9) {
    for (int i = 0; i < 9; ++i) {
      out[i] = box(i);
    }
    return;
  }
  throw std::invalid_argument("box is 3 lengths or 9 LAMMPS cell numbers");
}

struct Frame {
  seams_system *sys = nullptr;
  Frame() {
    sys = seams_system_create();
    if (sys == nullptr) {
      throw std::runtime_error(seams_last_error());
    }
  }
  ~Frame() { seams_system_destroy(sys); }
  Frame(const Frame &) = delete;
  Frame &operator=(const Frame &) = delete;
};

} // namespace

NB_MODULE(_core, m) {
  m.doc() = "d-SEAMS C API.";
  m.attr("ABI") = static_cast<int>(SEAMS_ABI);
  m.attr("PUBLISHED") = static_cast<int>(SEAMS_CHILL_PUBLISHED);
  m.attr("OVITO") = static_cast<int>(SEAMS_CHILL_OVITO);
  m.attr("CUTOFF") = static_cast<int>(SEAMS_GRAPH_CUTOFF);
  m.attr("MUTUAL") = static_cast<int>(SEAMS_GRAPH_MUTUAL);
  m.def("version", []() { return std::string(seams_version_string()); });
  m.def("set_num_threads", [](int n) { check(seams_set_num_threads(n)); });
  m.def("num_threads", []() { return seams_get_num_threads(); });

  nb::class_<Frame>(m, "Frame")
      .def(nb::init<>())
      .def("set_frame",
           [](Frame &frame,
              nb::ndarray<nb::numpy, double, nb::ndim<2>, nb::c_contig> pos,
              nb::ndarray<nb::numpy, double, nb::ndim<1>, nb::c_contig> box) {
             if (pos.shape(1) != 3) {
               throw std::invalid_argument("positions are N by 3");
             }
             double box9[9];
             box_of(box, box9);
             check(seams_system_set_frame(frame.sys, pos.data(),
                                          static_cast<int>(pos.shape(0)), box9, nullptr,
                                          nullptr));
           })
      .def("count", [](const Frame &frame) { return seams_system_count(frame.sys); })
      .def("chill_plus",
           [](Frame &frame, double cutoff, int graph, int rule) {
             const int n = seams_system_count(frame.sys);
             int *labels = new int[static_cast<std::size_t>(n)];
             int counts[9];
             try {
               check(seams_chill_plus(frame.sys, cutoff, 0, graph, rule, labels, counts));
             } catch (...) {
               delete[] labels;
               throw;
             }
             nb::capsule owner(labels, [](void *p) noexcept { delete[] static_cast<int *>(p); });
             int *copy = new int[9];
             for (int i = 0; i < 9; ++i) {
               copy[i] = counts[i];
             }
             nb::capsule cowner(copy, [](void *p) noexcept { delete[] static_cast<int *>(p); });
             return nb::make_tuple(
                 nb::ndarray<nb::numpy, int>(labels, {static_cast<std::size_t>(n)}, owner),
                 nb::ndarray<nb::numpy, int>(copy, {9}, cowner));
           },
           nb::arg("cutoff") = 3.5, nb::arg("graph") = 0, nb::arg("rule") = 0,
           "CHILL+ labels. rule 0 (PUBLISHED) averages q_lm over the four "
           "closest neighbours. rule 1 (OVITO) averages q_lm over every "
           "neighbour inside the cutoff.")
      .def("q4q6",
           [](Frame &frame, double cutoff) {
             const int n = seams_system_count(frame.sys);
             double *q4 = new double[static_cast<std::size_t>(n)];
             double *q6 = new double[static_cast<std::size_t>(n)];
             try {
               check(seams_q4q6(frame.sys, cutoff, 0, q4, q6));
             } catch (...) {
               delete[] q4;
               delete[] q6;
               throw;
             }
             nb::capsule o4(q4, [](void *p) noexcept { delete[] static_cast<double *>(p); });
             nb::capsule o6(q6, [](void *p) noexcept { delete[] static_cast<double *>(p); });
             return nb::make_tuple(
                 nb::ndarray<nb::numpy, double>(q4, {static_cast<std::size_t>(n)}, o4),
                 nb::ndarray<nb::numpy, double>(q6, {static_cast<std::size_t>(n)}, o6));
           },
           nb::arg("cutoff") = 3.5)
      .def("neighbours",
           [](Frame &frame, double cutoff, int graph) {
             seams_neighbours *nb = nullptr;
             check(seams_neighbours_build(frame.sys, cutoff, 0, graph, &nb));
             const int *offsets = nullptr;
             const int *cols = nullptr;
             int n_offsets = 0;
             int nnz = 0;
             check(seams_neighbours_csr(nb, &offsets, &n_offsets, &cols, &nnz));
             auto *off = new int[static_cast<std::size_t>(n_offsets)];
             auto *col = new int[static_cast<std::size_t>(nnz > 0 ? nnz : 1)];
             for (int i = 0; i < n_offsets; ++i) {
               off[i] = offsets[i];
             }
             for (int i = 0; i < nnz; ++i) {
               col[i] = cols[i];
             }
             seams_neighbours_destroy(nb);
             nb::capsule oo(off, [](void *p) noexcept { delete[] static_cast<int *>(p); });
             nb::capsule oc(col, [](void *p) noexcept { delete[] static_cast<int *>(p); });
             return nb::make_tuple(
                 nb::ndarray<nb::numpy, int>(off, {static_cast<std::size_t>(n_offsets)}, oo),
                 nb::ndarray<nb::numpy, int>(col, {static_cast<std::size_t>(nnz)}, oc));
           },
           nb::arg("cutoff") = 3.5, nb::arg("graph") = 0)
      .def("rings", [](Frame &frame, double cutoff, int max_depth) {
        seams_rings *rings = nullptr;
        check(seams_rings_find(frame.sys, cutoff, 0, SEAMS_GRAPH_CUTOFF, max_depth, &rings));
        const int nr = seams_rings_count(rings);
        int hist[13] = {};
        check(seams_rings_histogram(rings, hist, 13));
        auto *sizes = new int[static_cast<std::size_t>(nr > 0 ? nr : 1)];
        auto *sym = new int[static_cast<std::size_t>(nr > 0 ? nr : 1)];
        if (nr > 0) {
          check(seams_rings_labels(rings, sizes, sym));
        }
        double fraction = 0.0;
        check(seams_rings_fraction(rings, 6, &fraction));
        seams_rings_destroy(rings);
        auto *counts = new int[13];
        for (int i = 0; i < 13; ++i) {
          counts[i] = hist[i];
        }
        nb::capsule os(sizes, [](void *p) noexcept { delete[] static_cast<int *>(p); });
        nb::capsule oy(sym, [](void *p) noexcept { delete[] static_cast<int *>(p); });
        nb::capsule oc(counts, [](void *p) noexcept { delete[] static_cast<int *>(p); });
        return nb::make_tuple(
            nb::ndarray<nb::numpy, int>(sizes, {static_cast<std::size_t>(nr)}, os),
            nb::ndarray<nb::numpy, int>(sym, {static_cast<std::size_t>(nr)}, oy),
            nb::ndarray<nb::numpy, int>(counts, {13}, oc), fraction);
      },
           nb::arg("cutoff") = 3.5, nb::arg("max_depth") = 6);
}
