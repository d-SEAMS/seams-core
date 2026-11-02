#ifndef SEAMS_OZAKI_H_
#define SEAMS_OZAKI_H_

// Emulated products for the short dots and the one tiny GEMM in d-SEAMS.
//
// Ozaki scheme I splits a floating-point vector into slices whose products
// are exact in FP64 (Ozaki, Ogita and Oishi, Numer. Algorithms 59, 95, 2012).
// Ozaki scheme II replaces that slicing by an integer matrix product modulo
// several coprime moduli and a Chinese-remainder reconstruction (Ozaki,
// Uchino and Imamura, arXiv:2504.08009, 2025). The INT8 products use
// AVX-512 VNNI or AMX-INT8 when the shape fills a tile; shorter products
// use an INT32 accumulator, which is exact for the lengths that occur here.
//
// TwoProd / TwoSum and pairwise summation are the error-free transformations
// of Dekker, Knuth, and Ogita, Rump and Oishi (SIAM J. Sci. Comput. 26,
// 2005). Pairwise summation is the reduction whose error Higham bounds by
// gamma_{ceil(log2 n)} (Accuracy and Stability of Numerical Algorithms,
// SIAM, 2002). It does not depend on a left-to-right fold.

#include <complex>
#include <cstddef>
#include <vector>

namespace seams::product {

enum class Mode {
  Fp64 = 0,
  Dot2 = 1,
  OzakiI = 2,
  OzakiII = 3,
};

enum class Reduction {
  Left = 0,
  Pairwise = 1,
  Neumaier = 2,
};

enum class Int8Backend {
  Auto = 0,
  Scalar = 1,
  Vnni = 2,
  Amx = 3,
};

struct RealDot {
  double value;
  double bound;
  bool certified;
  int work;
  const char *backend;
};

struct ComplexDot {
  double re;
  double im;
  double normA;
  double normB;
  double boundRe;
  double boundIm;
  double boundNormA;
  double boundNormB;
  bool certified;
  const char *backend;
};

Mode mode();
void setMode(Mode m);
const char *modeName(Mode m);

Reduction reduction();
void setReduction(Reduction r);
const char *reductionName(Reduction r);

bool amxAvailable();
bool vnniAvailable();

void twoSum(double a, double b, double &s, double &e);
void twoProd(double a, double b, double &hi, double &lo);

double pairwiseSum(const double *x, int n);
double neumaierSum(const double *x, int n);
double leftSum(const double *x, int n);
double reduce(const double *x, int n, Reduction r);

// Absolute error bound for pairwise summation of n finite terms whose
// absolute values sum to absSum (Higham's gamma_ceil(log2 n) bound,
// rounded upward).
double pairwiseBound(int n, long double absSum);

RealDot realDot(const double *x, const double *y, int n, Mode m);

ComplexDot complexDot(const std::complex<double> *a,
                      const std::complex<double> *b, int n, Mode m);

// c = Re(a · conj(b)) / (|a| |b|), with an absolute certificate.
struct Cij {
  double value;
  double bound;
  bool certified;
  const char *backend;
};

Cij bondCorrelation(const std::complex<double> *a, const std::complex<double> *b,
                    int n, Mode m);

// Mean of the finite entries, in input order. Left matches a plain
// accumulator; pairwise and Neumaier do not depend on a chain of +=.
double meanOfFinite(const double *x, int n, Reduction r);

// Row-major C = A(m x k) * B(k x n). bounds may be null. On Ozaki-II each
// bounds[i*n+j] is an absolute certificate. Fp64 leaves bounds at 0 when
// the pointer is non-null (no certificate for the rounded GEMM).
struct GemmInfo {
  bool ok;
  const char *backend;
  int moduli;
  int sliceWidth;
};

GemmInfo gemm(int m, int n, int k, const double *A, int lda, const double *B,
              int ldb, double *C, int ldc, Mode mode, double *bounds,
              Int8Backend backend = Int8Backend::Auto);

// INT8 C += A * B, row-major, int32 accumulator. Used to check the tiles.
// Returns false when the backend cannot run this shape.
bool int8Gemm(int m, int n, int k, const signed char *A, int lda,
              const signed char *B, int ldb, int *C, int ldc,
              Int8Backend backend);

} // namespace seams::product

#endif
