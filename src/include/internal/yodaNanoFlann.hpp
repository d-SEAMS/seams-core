///////////////////////////////////////////////////////////////////////////////////////////
//
//
// Structure Factor code
// MIT License

// Copyright (c) 2018 Rohit Goswami, Amrita Goswami
//
// r95g10[at]gmail.com, amrita16thaug646[at]gmail.com
// Extensions to nanoflann for molecular simulations
//
///////////////////////////////////////////////////////////////////////////////////////////

/** Squared Euclidean (L2) distance functor for periodic boxes (suitable for low-dimensionality
 * datasets, like 2D or 3D point clouds) Corresponding distance traits:
 * nanoflann::metric_L2_Simple \tparam T Type of the elements (e.g. double,
 * float, uint8_t) \tparam _DistanceType Type of distance variables (must be
 * signed) (e.g. float, double, int64_t)
 */
#ifndef __YODANANOFLANN_H_
#define __YODANANOFLANN_H_

#include <nanoflann.hpp>

// Use the nanoflann namespace
namespace nanoflann {

template <class T, class DataSource, typename _DistanceType = T>
struct L2_Simple_Adaptor_MD {
  typedef T ElementType;
  typedef _DistanceType DistanceType;

  const DataSource &data_source;

  L2_Simple_Adaptor_MD(const DataSource &_data_source)
      : data_source(_data_source) {}

  inline DistanceType evalMetric(const T *a, const size_t b_idx,
                                 size_t size) const {
    DistanceType result = DistanceType();
    for (size_t i = 0; i < size; ++i) {
      const DistanceType diff = a[i] - data_source.kdtree_get_pt(b_idx, i);
      const DistanceType del =
          diff - data_source.box[i] * round(diff / data_source.box[i]);
      result += del * del;
    }
    return result;
  }

  template <typename U, typename V>
  inline DistanceType accum_dist(const U a, const V b, const size_t) const {
    return (a - b) * (a - b);
  }
};
/** Metaprogramming helper traits class for the L2_simple (Euclidean) molecular simluation metric */
struct metric_L2_Simple_MD : public Metric {
  template <class T, class DataSource> struct traits {
    typedef L2_Simple_Adaptor_MD<T, DataSource> distance_t;
  };
};
} // namespace nanoflann

#endif // __YODANANOFLANN_H_
