#include <vector>

#include <cpp11.hpp>
#include <cpp11/data_frame.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/strings.hpp>

#include "seams_lang.h"

using namespace cpp11::literals;

using frame_ptr = cpp11::external_pointer<seams_frame, seams_frame_free>;

namespace {

const char *ice_name(int label) {
  static const char *names[9] = {
      "cubic",     "hexagonal", "water",         "interfacial", "clathrate",
      "interClathrate", "unclassified", "reCubic", "reHex"};
  if (label < 0 || label > 8) {
    return "unclassified";
  }
  return names[label];
}

seams_frame *need(SEXP ptr) {
  frame_ptr xp(ptr);
  seams_frame *frame = xp.get();
  if (frame == nullptr) {
    cpp11::stop("empty dseams frame");
  }
  return frame;
}

} // namespace

[[cpp11::register]]
cpp11::sexp read_frame_(std::string path, int frame, int type_filter) {
  seams_frame *held = seams_frame_read(path.c_str(), frame, type_filter);
  if (held == nullptr) {
    cpp11::stop("%s", seams_last_error());
  }
  return cpp11::sexp(static_cast<SEXP>(frame_ptr(held)));
}

[[cpp11::register]]
int nframes_(std::string path) {
  const int n = seams_path_nframes(path.c_str());
  if (n < 0) {
    cpp11::stop("%s", seams_last_error());
  }
  return n;
}

[[cpp11::register]]
int nop_(SEXP ptr) { return seams_frame_nop(need(ptr)); }

[[cpp11::register]]
cpp11::doubles positions_(SEXP ptr) {
  seams_frame *frame = need(ptr);
  const int n = seams_frame_nop(frame);
  const double *xyz = seams_frame_positions(frame);
  cpp11::writable::doubles out(static_cast<R_xlen_t>(n) * 3);
  for (R_xlen_t i = 0; i < out.size(); ++i) {
    out[i] = xyz[static_cast<std::size_t>(i)];
  }
  return out;
}

[[cpp11::register]]
cpp11::data_frame chill_plus_(SEXP ptr, double cutoff, int type,
                              std::string graph, std::string rule) {
  seams_frame *frame = need(ptr);
  const int n = seams_frame_nop(frame);
  std::vector<int> labels(static_cast<std::size_t>(n));
  int counts[9] = {};
  if (seams_frame_chill_plus_rule(frame, cutoff, type, graph.c_str(),
                                  rule.c_str(), labels.data(), counts) != 0) {
    cpp11::stop("%s", seams_last_error());
  }
  const int *types = seams_frame_types(frame);
  const int *ids = seams_frame_ids(frame);
  cpp11::writable::integers atom(n);
  cpp11::writable::integers type_col(n);
  cpp11::writable::integers label(n);
  cpp11::writable::strings name(n);
  for (int i = 0; i < n; ++i) {
    atom[i] = ids != nullptr ? ids[i] : i + 1;
    type_col[i] = types != nullptr ? types[i] : NA_INTEGER;
    label[i] = labels[static_cast<std::size_t>(i)];
    name[i] = ice_name(label[i]);
  }
  cpp11::writable::data_frame out(
      {"atom"_nm = atom, "type"_nm = type_col, "label"_nm = label,
       "name"_nm = name});
  out.attr("counts") = cpp11::writable::integers(counts, counts + 9);
  return out;
}

[[cpp11::register]]
cpp11::data_frame rings_(SEXP ptr, double cutoff, int type, std::string graph,
                         int max_depth) {
  seams_frame *frame = need(ptr);
  std::vector<int> counts(static_cast<std::size_t>(max_depth + 1), 0);
  if (seams_frame_rings(frame, cutoff, type, graph.c_str(), max_depth,
                        counts.data(), static_cast<int>(counts.size())) != 0) {
    cpp11::stop("%s", seams_last_error());
  }
  int rows = 0;
  for (int size = 3; size <= max_depth; ++size) {
    if (counts[static_cast<std::size_t>(size)] > 0) {
      ++rows;
    }
  }
  cpp11::writable::integers size_col(rows);
  cpp11::writable::integers n_col(rows);
  int row = 0;
  for (int size = 3; size <= max_depth; ++size) {
    const int n = counts[static_cast<std::size_t>(size)];
    if (n == 0) {
      continue;
    }
    size_col[row] = size;
    n_col[row] = n;
    ++row;
  }
  return cpp11::writable::data_frame({"size"_nm = size_col, "n"_nm = n_col});
}

[[cpp11::register]]
cpp11::data_frame steinhardt_(SEXP ptr, double cutoff, int type, int order_l) {
  seams_frame *frame = need(ptr);
  const int n = seams_frame_nop(frame);
  std::vector<double> ql(static_cast<std::size_t>(n));
  std::vector<double> qlbar(static_cast<std::size_t>(n));
  if (seams_frame_steinhardt(frame, cutoff, type, order_l, ql.data(),
                             qlbar.data()) != 0) {
    cpp11::stop("%s", seams_last_error());
  }
  const int *ids = seams_frame_ids(frame);
  cpp11::writable::integers atom(n);
  cpp11::writable::doubles ql_col(n);
  cpp11::writable::doubles bar_col(n);
  for (int i = 0; i < n; ++i) {
    atom[i] = ids != nullptr ? ids[i] : i + 1;
    ql_col[i] = ql[static_cast<std::size_t>(i)];
    bar_col[i] = qlbar[static_cast<std::size_t>(i)];
  }
  cpp11::writable::data_frame out(
      {"atom"_nm = atom, "ql"_nm = ql_col, "qlbar"_nm = bar_col});
  out.attr("order_l") = order_l;
  return out;
}

[[cpp11::register]]
cpp11::data_frame cages_(SEXP ptr, double cutoff, int type, int k,
                         std::string graph) {
  seams_frame *frame = need(ptr);
  const int n = seams_frame_nop(frame);
  std::vector<int> per(static_cast<std::size_t>(n));
  int n_hex = 0, n_cubic = 0, n_water = 0, n_both = 0;
  if (seams_frame_cages(frame, cutoff, type, k, graph.c_str(), 0, per.data(),
                        &n_hex, &n_cubic, &n_water, &n_both) != 0) {
    cpp11::stop("%s", seams_last_error());
  }
  const int *ids = seams_frame_ids(frame);
  cpp11::writable::integers atom(n);
  cpp11::writable::integers cage(n);
  static const char *cage_name[] = {"water", "hexagonal", "cubic", "both"};
  cpp11::writable::strings name(n);
  for (int i = 0; i < n; ++i) {
    atom[i] = ids != nullptr ? ids[i] : i + 1;
    const int code = per[static_cast<std::size_t>(i)];
    cage[i] = code;
    name[i] = cage_name[code >= 0 && code <= 3 ? code : 0];
  }
  cpp11::writable::data_frame out(
      {"atom"_nm = atom, "cage"_nm = cage, "name"_nm = name});
  out.attr("n_ih") = n_hex;
  out.attr("n_ic") = n_cubic;
  out.attr("n_water") = n_water;
  out.attr("n_both") = n_both;
  return out;
}
