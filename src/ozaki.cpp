#include <ozaki.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>

#if defined(__x86_64__)
#include <immintrin.h>
#endif
#if defined(__x86_64__) && defined(__linux__)
#include <sys/syscall.h>
#include <unistd.h>
#endif

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace seams::product {
namespace {

// __int128 is a GCC/Clang extension; the typedef keeps -pedantic quiet.
__extension__ typedef __int128 i128;

std::atomic<int> gMode{-1};
std::atomic<int> gReduction{-1};

int fromEnvMode() {
  const char *env = std::getenv("SEAMS_PRODUCT");
  if (env == nullptr) {
    return static_cast<int>(Mode::Fp64);
  }
  if (std::strcmp(env, "dot2") == 0) {
    return static_cast<int>(Mode::Dot2);
  }
  if (std::strcmp(env, "ozaki1") == 0 || std::strcmp(env, "ozaki-i") == 0) {
    return static_cast<int>(Mode::OzakiI);
  }
  if (std::strcmp(env, "ozaki2") == 0 || std::strcmp(env, "ozaki-ii") == 0) {
    return static_cast<int>(Mode::OzakiII);
  }
  return static_cast<int>(Mode::Fp64);
}

int fromEnvReduction() {
  const char *env = std::getenv("SEAMS_SUM");
  if (env == nullptr) {
    return static_cast<int>(Reduction::Left);
  }
  if (std::strcmp(env, "pairwise") == 0) {
    return static_cast<int>(Reduction::Pairwise);
  }
  if (std::strcmp(env, "neumaier") == 0) {
    return static_cast<int>(Reduction::Neumaier);
  }
  return static_cast<int>(Reduction::Left);
}

int ceilLog2(int n) {
  if (n <= 1) {
    return 0;
  }
  int l = 0;
  int v = 1;
  while (v < n && l < 30) {
    v <<= 1;
    ++l;
  }
  return l;
}

double upward(long double x) {
  if (!(x > 0)) {
    return 0.0;
  }
  if (x > static_cast<long double>(std::numeric_limits<double>::max())) {
    return std::numeric_limits<double>::infinity();
  }
  const double d = static_cast<double>(x);
  return std::nextafter(d, std::numeric_limits<double>::infinity());
}

bool finitePair(double a, double b) {
  return std::isfinite(a) && std::isfinite(b);
}

#if defined(__x86_64__) && defined(__linux__)
#ifndef ARCH_REQ_XCOMP_PERM
#define ARCH_REQ_XCOMP_PERM 0x1023
#endif
#ifndef XFEATURE_XTILEDATA
#define XFEATURE_XTILEDATA 18
#endif

bool requestAmx() {
  const long rc = syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA);
  return rc == 0;
}
#endif

struct CrtTable {
  int s = 0;
  uint32_t m[16] = {};
  __uint128_t modulus = 0;
  __uint128_t weight[16] = {};
  int sumBits = 0;
};

CrtTable makeCrt() {
  // Largest pairwise-coprime moduli in 2..256, the INT8 set of Ozaki,
  // Uchino and Imamura, arXiv:2504.08009, equation (19). The product is
  // kept inside a signed 128-bit CRT accumulator.
  static const uint32_t cand[] = {256, 255, 253, 251, 247, 239, 233, 229,
                                  227, 223, 217, 211, 199, 197, 193, 191};
  CrtTable t;
  t.modulus = 1;
  for (uint32_t m : cand) {
    if (t.modulus > (~__uint128_t{0}) / m) {
      break;
    }
    const __uint128_t next = t.modulus * m;
    // c * weight must fit in signed int128, |c| <= m/2, weight < modulus.
    if (next > (__uint128_t{1} << 110)) {
      break;
    }
    t.m[t.s] = m;
    t.modulus = next;
    ++t.s;
  }
  auto inv = [](int64_t a, int64_t mod) {
    int64_t t0 = 0;
    int64_t t1 = 1;
    int64_t r0 = mod;
    int64_t r1 = a % mod;
    while (r1 != 0) {
      const int64_t q = r0 / r1;
      const int64_t r = r0 - q * r1;
      const int64_t nt = t0 - q * t1;
      r0 = r1;
      r1 = r;
      t0 = t1;
      t1 = nt;
    }
    if (t0 < 0) {
      t0 += mod;
    }
    return t0;
  };
  for (int i = 0; i < t.s; i++) {
    const __uint128_t Mi = t.modulus / t.m[i];
    const int64_t yi = inv(static_cast<int64_t>(Mi % t.m[i]), t.m[i]);
    t.weight[i] = Mi * static_cast<__uint128_t>(yi);
  }
  __uint128_t lim = t.modulus / 2 - 1;
  int bits = 0;
  while (bits < 120 && (__uint128_t{1} << bits) <= lim) {
    ++bits;
  }
  t.sumBits = bits - 1;
  return t;
}

const CrtTable &crt() {
  static const CrtTable table = makeCrt();
  return table;
}

int symMod(int32_t a, int32_t m) {
  int32_t r = a % m;
  if (r < 0) {
    r += m;
  }
  if (r >= (m + 1) / 2) {
    r -= m;
  }
  return r;
}

i128 remCentered(i128 s, __uint128_t M) {
  const __uint128_t half = M / 2;
  __uint128_t u;
  if (s >= 0) {
    u = static_cast<__uint128_t>(s) % M;
  } else {
    const __uint128_t neg = static_cast<__uint128_t>(-s) % M;
    u = neg == 0 ? 0 : M - neg;
  }
  if (u > half) {
    return static_cast<i128>(u - M);
  }
  return static_cast<i128>(u);
}

double roundI128(i128 v, int extraExp) {
  if (v == 0) {
    return 0.0;
  }
  const bool neg = v < 0;
  __uint128_t u = neg ? static_cast<__uint128_t>(-v) : static_cast<__uint128_t>(v);
  int bits = 0;
  __uint128_t probe = u;
  while (probe > 0) {
    probe >>= 1;
    ++bits;
  }
  int shift = bits - 53;
  __uint128_t mant;
  if (shift <= 0) {
    mant = u;
  } else {
    const __uint128_t half = __uint128_t{1} << (shift - 1);
    const __uint128_t rem = u & ((__uint128_t{1} << shift) - 1);
    mant = u >> shift;
    if (rem > half || (rem == half && (mant & 1))) {
      ++mant;
    }
    if (mant == (__uint128_t{1} << 53)) {
      mant >>= 1;
      ++shift;
    }
    extraExp += shift;
  }
  const double d = std::ldexp(static_cast<double>(static_cast<uint64_t>(mant)), extraExp);
  return neg ? -d : d;
}

int ufpExpOfRecip(double a) {
  const int ea = std::ilogb(a);
  const double p2 = std::ldexp(1.0, ea);
  if (std::fabs(a) == p2) {
    return -ea;
  }
  return -ea - 1;
}

int sliceWidthFor(int n) {
  const int lg = ceilLog2(std::max(n, 1));
  int w = (53 - lg) / 2;
  if (w > 26) {
    w = 26;
  }
  if (w < 8) {
    w = 8;
  }
  return w;
}

__attribute__((optimize("fp-contract=off"))) double exactIntegerDot(
    const double *cx, const double *cy, int n) {
  double acc = 0.0;
  for (int i = 0; i < n; i++) {
    const double p = cx[i] * cy[i];
    acc = acc + p;
  }
  return acc;
}

struct Slice {
  std::vector<double> coeff;
  int unitExp = 0;
};

bool takeSlice(const std::vector<double> &x, int w, Slice &s, std::vector<double> &rem) {
  double amax = 0.0;
  for (double v : x) {
    amax = std::max(amax, std::fabs(v));
  }
  const int n = static_cast<int>(x.size());
  s.coeff.assign(static_cast<std::size_t>(n), 0.0);
  rem.assign(static_cast<std::size_t>(n), 0.0);
  if (!(amax > 0.0) || !std::isfinite(amax)) {
    return false;
  }
  const int e = std::ilogb(amax);
  s.unitExp = e - (w - 1);
  bool any = false;
  for (int i = 0; i < n; i++) {
    const double v = x[static_cast<std::size_t>(i)];
    if (v == 0.0 || !std::isfinite(v)) {
      rem[static_cast<std::size_t>(i)] = v;
      continue;
    }
    if (std::ilogb(std::fabs(v)) < s.unitExp) {
      rem[static_cast<std::size_t>(i)] = v;
      continue;
    }
    const double q = std::ldexp(v, -s.unitExp);
    const double qi = std::trunc(q);
    s.coeff[static_cast<std::size_t>(i)] = qi;
    const double sliced = std::ldexp(qi, s.unitExp);
    rem[static_cast<std::size_t>(i)] = v - sliced;
    if (qi != 0.0) {
      any = true;
    }
  }
  return any;
}

RealDot dotFp64(const double *x, const double *y, int n) {
  double acc = 0.0;
  long double absSum = 0.0;
  for (int i = 0; i < n; i++) {
    const double p = x[i] * y[i];
    acc += p;
    absSum += std::fabs(static_cast<long double>(p));
  }
  RealDot out;
  out.value = acc;
  out.bound = pairwiseBound(std::max(n, 1), absSum);
  out.certified = false;
  out.work = n;
  out.backend = "fp64";
  return out;
}

RealDot dotPairwiseProducts(const double *x, const double *y, int n) {
  std::vector<double> terms(static_cast<std::size_t>(std::max(n, 0)) * 2);
  long double absSum = 0.0;
  bool exact = true;
  for (int i = 0; i < n; i++) {
    double hi = 0.0;
    double lo = 0.0;
    twoProd(x[i], y[i], hi, lo);
    terms[static_cast<std::size_t>(2 * i)] = hi;
    terms[static_cast<std::size_t>(2 * i + 1)] = lo;
    absSum += std::fabs(static_cast<long double>(hi));
    absSum += std::fabs(static_cast<long double>(lo));
    if (!finitePair(x[i], y[i]) || !std::isfinite(hi) || !std::isfinite(lo)) {
      exact = false;
    }
    if (hi == 0.0 && lo == 0.0 && x[i] != 0.0 && y[i] != 0.0) {
      exact = false;
    }
  }
  RealDot out;
  out.value = pairwiseSum(terms.data(), static_cast<int>(terms.size()));
  out.bound = pairwiseBound(static_cast<int>(terms.size()), absSum);
  out.certified = exact;
  out.work = static_cast<int>(terms.size());
  out.backend = "dot2";
  return out;
}

RealDot dotOzakiI(const double *x, const double *y, int n) {
  if (n <= 0) {
    return RealDot{0.0, 0.0, true, 0, "ozaki-i"};
  }
  const int w = sliceWidthFor(n);
  const int lg = ceilLog2(n);
  const bool widthOk = (2 * w + lg) <= 53;
  std::vector<double> rx(x, x + n);
  std::vector<double> ry(y, y + n);
  std::vector<Slice> sx;
  std::vector<Slice> sy;
  for (int guard = 0; guard < 128; ++guard) {
    Slice s;
    std::vector<double> rem;
    if (!takeSlice(rx, w, s, rem)) {
      break;
    }
    sx.push_back(std::move(s));
    rx.swap(rem);
  }
  for (int guard = 0; guard < 128; ++guard) {
    Slice s;
    std::vector<double> rem;
    if (!takeSlice(ry, w, s, rem)) {
      break;
    }
    sy.push_back(std::move(s));
    ry.swap(rem);
  }
  bool exactSplit = widthOk;
  for (double v : rx) {
    if (v != 0.0) {
      exactSplit = false;
    }
  }
  for (double v : ry) {
    if (v != 0.0) {
      exactSplit = false;
    }
  }
  std::vector<double> terms;
  terms.reserve(sx.size() * sy.size());
  long double absSum = 0.0;
  for (const Slice &a : sx) {
    for (const Slice &b : sy) {
      const double idot = exactIntegerDot(a.coeff.data(), b.coeff.data(), n);
      const double piece = std::ldexp(idot, a.unitExp + b.unitExp);
      terms.push_back(piece);
      absSum += std::fabs(static_cast<long double>(piece));
      if (!std::isfinite(piece)) {
        exactSplit = false;
      }
    }
  }
  RealDot out;
  out.value = pairwiseSum(terms.data(), static_cast<int>(terms.size()));
  out.bound = pairwiseBound(static_cast<int>(terms.size()), absSum);
  out.certified = exactSplit && std::isfinite(out.value);
  out.work = static_cast<int>(terms.size());
  out.backend = "ozaki-i";
  if (!out.certified) {
    const RealDot fallback = dotPairwiseProducts(x, y, n);
    out = fallback;
    out.backend = "ozaki-i-dot2";
  }
  return out;
}

bool truncScaled(const double *x, int n, int kBits, std::vector<int64_t> &out,
                 int &exp) {
  out.assign(static_cast<std::size_t>(n), 0);
  double amax = 0.0;
  for (int i = 0; i < n; i++) {
    if (!std::isfinite(x[i])) {
      return false;
    }
    amax = std::max(amax, std::fabs(x[i]));
  }
  if (amax == 0.0) {
    exp = 0;
    return true;
  }
  const int ue = ufpExpOfRecip(amax);
  exp = (kBits - 1) + ue;
  const double limit = std::ldexp(1.0, kBits);
  for (int i = 0; i < n; i++) {
    const double scaled = std::ldexp(x[i], exp);
    if (!std::isfinite(scaled) || std::fabs(scaled) >= limit) {
      return false;
    }
    out[static_cast<std::size_t>(i)] = static_cast<int64_t>(std::trunc(scaled));
  }
  return true;
}

RealDot dotOzakiII(const double *x, const double *y, int n) {
  RealDot out{0.0, std::numeric_limits<double>::infinity(), false, 0, "ozaki-ii"};
  if (n <= 0) {
    out.value = 0.0;
    out.bound = 0.0;
    out.certified = true;
    return out;
  }
  const CrtTable &table = crt();
  if (table.s < 2 || table.sumBits <= ceilLog2(std::max(n, 1))) {
    return out;
  }
  const int kSum = table.sumBits - ceilLog2(std::max(n, 1));
  const int kA = kSum / 2;
  const int kB = kSum - kA;
  std::vector<int64_t> xp;
  std::vector<int64_t> yp;
  int ex = 0;
  int ey = 0;
  if (!truncScaled(x, n, kA, xp, ex) || !truncScaled(y, n, kB, yp, ey)) {
    return out;
  }
  const int modCount = table.s;
  std::vector<int> residues(static_cast<std::size_t>(modCount));
  for (int t = 0; t < modCount; t++) {
    const int m = static_cast<int>(table.m[t]);
    int32_t acc = 0;
    for (int i = 0; i < n; i++) {
      const int a = symMod(static_cast<int32_t>(xp[static_cast<std::size_t>(i)] % m), m);
      const int b = symMod(static_cast<int32_t>(yp[static_cast<std::size_t>(i)] % m), m);
      acc += static_cast<int32_t>(a) * static_cast<int32_t>(b);
    }
    residues[static_cast<std::size_t>(t)] = symMod(acc, m);
  }
  i128 acc = 0;
  for (int t = 0; t < modCount; t++) {
    const i128 term =
        static_cast<i128>(residues[static_cast<std::size_t>(t)]) *
        static_cast<i128>(table.weight[t]);
    acc = remCentered(acc + term, table.modulus);
  }
  out.value = roundI128(acc, -(ex + ey));
  long double trunc = 0.0;
  const long double invDx = std::ldexp(1.0L, -ex);
  const long double invDy = std::ldexp(1.0L, -ey);
  for (int i = 0; i < n; i++) {
    trunc += invDx * std::fabs(static_cast<long double>(y[i]));
    trunc += std::fabs(static_cast<long double>(xp[static_cast<std::size_t>(i)])) * invDx *
             invDy;
  }
  trunc *= 1.0L + 16.0L * static_cast<long double>(n) *
                       std::numeric_limits<long double>::epsilon();
  const double roundErr =
      std::isfinite(out.value) ? std::ldexp(std::fabs(out.value), -53) : 0.0;
  out.bound = upward(trunc + static_cast<long double>(roundErr));
  out.certified = std::isfinite(out.value);
  out.work = modCount;
  out.backend = "ozaki-ii-scalar";
  return out;
}

void packRealImag(const std::complex<double> *z, int n, std::vector<double> &re,
                  std::vector<double> &im) {
  re.resize(static_cast<std::size_t>(n));
  im.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    re[static_cast<std::size_t>(i)] = z[i].real();
    im[static_cast<std::size_t>(i)] = z[i].imag();
  }
}

#if defined(__x86_64__)
__attribute__((target("avx512f,avx512bw,avx512vnni"))) bool int8GemmVnni(
    int m, int n, int k, const signed char *A, int lda, const signed char *B,
    int ldb, int *C, int ldc) {
  if (!vnniAvailable()) {
    return false;
  }
  const int kp = (k + 3) & ~3;
  const int np = (n + 15) & ~15;
  std::vector<uint8_t> Au(static_cast<std::size_t>(m) * kp);
  std::vector<int8_t> Bp(static_cast<std::size_t>(kp) * np);
  for (int i = 0; i < m; i++) {
    for (int t = 0; t < k; t++) {
      const int s = static_cast<int>(A[i * lda + t]);
      Au[static_cast<std::size_t>(i) * kp + t] = static_cast<uint8_t>(s + 128);
    }
  }
  for (int t = 0; t < k; t++) {
    for (int j = 0; j < n; j++) {
      Bp[static_cast<std::size_t>(t) * np + j] = B[t * ldb + j];
    }
  }
  std::vector<int32_t> colSum(static_cast<std::size_t>(n), 0);
  for (int t = 0; t < k; t++) {
    for (int j = 0; j < n; j++) {
      colSum[static_cast<std::size_t>(j)] += Bp[static_cast<std::size_t>(t) * np + j];
    }
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (m >= 32)
#endif
  for (int i = 0; i < m; i++) {
    std::vector<int32_t> row(static_cast<std::size_t>(np), 0);
    for (int j0 = 0; j0 < n; j0 += 16) {
      __m512i acc = _mm512_setzero_si512();
      for (int t = 0; t < kp; t += 4) {
        uint32_t a4 = 0;
        std::memcpy(&a4, Au.data() + static_cast<std::size_t>(i) * kp + t, 4);
        const __m512i va = _mm512_set1_epi32(static_cast<int>(a4));
        alignas(64) int8_t packed[64];
        for (int lane = 0; lane < 4; lane++) {
          for (int c = 0; c < 16; c++) {
            const int j = j0 + c;
            int8_t b = 0;
            if (j < np && t + lane < kp) {
              b = Bp[static_cast<std::size_t>(t + lane) * np + j];
            }
            packed[c * 4 + lane] = b;
          }
        }
        const __m512i vb = _mm512_load_si512(packed);
        acc = _mm512_dpbusd_epi32(acc, va, vb);
      }
      alignas(64) int32_t tmp[16];
      _mm512_store_si512(tmp, acc);
      for (int c = 0; c < 16 && j0 + c < n; c++) {
        const int j = j0 + c;
        const int32_t corrected =
            tmp[c] - 128 * colSum[static_cast<std::size_t>(j)];
        row[static_cast<std::size_t>(j)] = corrected;
      }
    }
    for (int j = 0; j < n; j++) {
      C[i * ldc + j] += row[static_cast<std::size_t>(j)];
    }
  }
  return true;
}

#if defined(__linux__)
__attribute__((target("amx-tile,amx-int8"))) bool int8GemmAmx(
    int m, int n, int k, const signed char *A, int lda, const signed char *B,
    int ldb, int *C, int ldc) {
  if (!amxAvailable()) {
    return false;
  }
  const int mp = (m + 15) & ~15;
  const int np = (n + 15) & ~15;
  const int kp = (k + 63) & ~63;
  std::vector<int8_t> Ap(static_cast<std::size_t>(mp) * kp, 0);
  std::vector<int8_t> Bp(static_cast<std::size_t>(kp) * np, 0);
  for (int i = 0; i < m; i++) {
    for (int t = 0; t < k; t++) {
      Ap[static_cast<std::size_t>(i) * kp + t] = A[i * lda + t];
    }
  }
  for (int t = 0; t < k; t++) {
    for (int j = 0; j < n; j++) {
      Bp[static_cast<std::size_t>(t) * np + j] = B[t * ldb + j];
    }
  }
  struct TileCfg {
    uint8_t palette_id;
    uint8_t start_row;
    uint8_t reserved[14];
    uint16_t colsb[16];
    uint8_t rows[16];
  };
  // One thread. On this KVM guest, several threads issuing AMX tiles
  // at once occasionally returned a wrong 16x16 tile (3 of 20 runs at
  // 8 threads, n = 512). A single tile stream matched the scalar product.
  void *raw = nullptr;
  if (posix_memalign(&raw, 64, 16 * 64 + 16 * 64 + 16 * 16 * sizeof(int32_t)) != 0) {
    return false;
  }
  if (!requestAmx()) {
    std::free(raw);
    return false;
  }
  TileCfg cfg;
  std::memset(&cfg, 0, sizeof cfg);
  cfg.palette_id = 1;
  cfg.rows[0] = 16;
  cfg.colsb[0] = 64;
  cfg.rows[1] = 16;
  cfg.colsb[1] = 64;
  cfg.rows[2] = 16;
  cfg.colsb[2] = 64;
  _tile_loadconfig(&cfg);
  auto *base = static_cast<unsigned char *>(raw);
  int8_t *at = reinterpret_cast<int8_t *>(base);
  int8_t *bt = at + 16 * 64;
  int32_t *ct = reinterpret_cast<int32_t *>(bt + 16 * 64);
  for (int i0 = 0; i0 < mp; i0 += 16) {
        for (int j0 = 0; j0 < np; j0 += 16) {
          _tile_zero(2);
          for (int t0 = 0; t0 < kp; t0 += 64) {
            for (int r = 0; r < 16; r++) {
              std::memcpy(at + r * 64,
                          Ap.data() + static_cast<std::size_t>(i0 + r) * kp + t0, 64);
            }
            for (int g = 0; g < 16; g++) {
              for (int c = 0; c < 16; c++) {
                for (int lane = 0; lane < 4; lane++) {
                  const int kk = t0 + g * 4 + lane;
                  bt[(g * 16 + c) * 4 + lane] =
                      Bp[static_cast<std::size_t>(kk) * np + (j0 + c)];
                }
              }
            }
            _tile_loadd(0, at, 64);
            _tile_loadd(1, bt, 64);
            _tile_dpbssd(2, 0, 1);
          }
          _tile_stored(2, ct, 64);
          for (int r = 0; r < 16; r++) {
            const int i = i0 + r;
            if (i >= m) {
              break;
            }
            for (int c = 0; c < 16; c++) {
              const int j = j0 + c;
              if (j >= n) {
                break;
              }
              C[i * ldc + j] += ct[r * 16 + c];
            }
          }
        }
      }
  _tile_release();
  std::free(raw);
  return true;
}
#else
// AMX tiles need the Linux arch_prctl permission request.
bool int8GemmAmx(int, int, int, const signed char *, int, const signed char *, int, int *,
                 int) {
  return false;
}
#endif

__attribute__((target("avx512f,fma"))) void gemmAvx512Rows(
    int m, int n, int k, const double *A, int lda, const double *B, int ldb,
    double *C, int ldc) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (m >= 32)
#endif
  for (int i = 0; i < m; i++) {
    int j = 0;
    for (; j + 8 <= n; j += 8) {
      __m512d acc = _mm512_setzero_pd();
      for (int t = 0; t < k; t++) {
        const __m512d b = _mm512_loadu_pd(B + static_cast<std::ptrdiff_t>(t) * ldb + j);
        acc = _mm512_fmadd_pd(_mm512_set1_pd(A[static_cast<std::ptrdiff_t>(i) * lda + t]), b, acc);
      }
      _mm512_storeu_pd(C + static_cast<std::ptrdiff_t>(i) * ldc + j, acc);
    }
    for (; j < n; j++) {
      double acc = 0.0;
      for (int t = 0; t < k; t++) {
        acc += A[static_cast<std::ptrdiff_t>(i) * lda + t] *
               B[static_cast<std::ptrdiff_t>(t) * ldb + j];
      }
      C[static_cast<std::ptrdiff_t>(i) * ldc + j] = acc;
    }
  }
}
#endif

void gemmScalar(int m, int n, int k, const double *A, int lda, const double *B,
                int ldb, double *C, int ldc) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (m >= 64)
#endif
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double acc = 0.0;
      for (int t = 0; t < k; t++) {
        acc += A[i * lda + t] * B[t * ldb + j];
      }
      C[i * ldc + j] = acc;
    }
  }
}

void int8GemmScalar(int m, int n, int k, const signed char *A, int lda,
                    const signed char *B, int ldb, int *C, int ldc) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (m >= 32)
#endif
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      int32_t acc = 0;
      for (int t = 0; t < k; t++) {
        acc += static_cast<int32_t>(A[i * lda + t]) *
               static_cast<int32_t>(B[t * ldb + j]);
      }
      C[i * ldc + j] += acc;
    }
  }
}

} // namespace

Mode mode() {
  int v = gMode.load(std::memory_order_relaxed);
  if (v < 0) {
    const int got = fromEnvMode();
    int expected = -1;
    if (!gMode.compare_exchange_strong(expected, got)) {
      v = expected;
    } else {
      v = got;
    }
  }
  return static_cast<Mode>(v);
}

void setMode(Mode m) { gMode.store(static_cast<int>(m), std::memory_order_relaxed); }

const char *modeName(Mode m) {
  switch (m) {
  case Mode::Fp64:
    return "fp64";
  case Mode::Dot2:
    return "dot2";
  case Mode::OzakiI:
    return "ozaki-i";
  case Mode::OzakiII:
    return "ozaki-ii";
  }
  return "fp64";
}

Reduction reduction() {
  int v = gReduction.load(std::memory_order_relaxed);
  if (v < 0) {
    const int got = fromEnvReduction();
    int expected = -1;
    if (!gReduction.compare_exchange_strong(expected, got)) {
      v = expected;
    } else {
      v = got;
    }
  }
  return static_cast<Reduction>(v);
}

void setReduction(Reduction r) {
  gReduction.store(static_cast<int>(r), std::memory_order_relaxed);
}

const char *reductionName(Reduction r) {
  switch (r) {
  case Reduction::Left:
    return "left";
  case Reduction::Pairwise:
    return "pairwise";
  case Reduction::Neumaier:
    return "neumaier";
  }
  return "left";
}

bool amxAvailable() {
#if defined(__x86_64__) && defined(__linux__)
  static const bool ok = requestAmx();
  return ok;
#else
  return false;
#endif
}

bool vnniAvailable() {
#if defined(__x86_64__)
  return __builtin_cpu_supports("avx512vnni");
#else
  return false;
#endif
}

void twoSum(double a, double b, double &s, double &e) {
  s = a + b;
  const double a1 = s - b;
  const double b1 = s - a1;
  const double da = a - a1;
  const double db = b - b1;
  e = da + db;
}

void twoProd(double a, double b, double &hi, double &lo) {
  hi = a * b;
  lo = std::fma(a, b, -hi);
}

double leftSum(const double *x, int n) {
  double acc = 0.0;
  for (int i = 0; i < n; i++) {
    acc += x[i];
  }
  return acc;
}

double pairwiseSum(const double *x, int n) {
  if (n <= 0) {
    return 0.0;
  }
  std::vector<double> buf(x, x + n);
  int len = n;
  while (len > 1) {
    const int half = len / 2;
    for (int i = 0; i < half; i++) {
      buf[static_cast<std::size_t>(i)] =
          buf[static_cast<std::size_t>(2 * i)] + buf[static_cast<std::size_t>(2 * i + 1)];
    }
    if (len % 2 == 1) {
      buf[static_cast<std::size_t>(half)] = buf[static_cast<std::size_t>(len - 1)];
      len = half + 1;
    } else {
      len = half;
    }
  }
  return buf[0];
}

double neumaierSum(const double *x, int n) {
  if (n <= 0) {
    return 0.0;
  }
  double s = x[0];
  double c = 0.0;
  for (int i = 1; i < n; i++) {
    const double v = x[i];
    const double t = s + v;
    if (std::fabs(s) >= std::fabs(v)) {
      c += (s - t) + v;
    } else {
      c += (v - t) + s;
    }
    s = t;
  }
  return s + c;
}

double reduce(const double *x, int n, Reduction r) {
  switch (r) {
  case Reduction::Pairwise:
    return pairwiseSum(x, n);
  case Reduction::Neumaier:
    return neumaierSum(x, n);
  case Reduction::Left:
    return leftSum(x, n);
  }
  return leftSum(x, n);
}

double pairwiseBound(int n, long double absSum) {
  if (!(absSum > 0) || n <= 1) {
    return 0.0;
  }
  const int levels = std::max(ceilLog2(n), 1);
  constexpr long double u = 0x1p-53L;
  const long double ku = static_cast<long double>(levels) * u;
  if (!(ku < 0.5L)) {
    return std::numeric_limits<double>::infinity();
  }
  const long double gamma = ku / (1.0L - ku);
  long double b = gamma * absSum;
  b *= 1.0L + 8.0L * static_cast<long double>(levels) *
                   std::numeric_limits<long double>::epsilon();
  return upward(b);
}

RealDot realDot(const double *x, const double *y, int n, Mode m) {
  if (n <= 0) {
    return RealDot{0.0, 0.0, true, 0, modeName(m)};
  }
  switch (m) {
  case Mode::Dot2:
    return dotPairwiseProducts(x, y, n);
  case Mode::OzakiI:
    return dotOzakiI(x, y, n);
  case Mode::OzakiII:
    return dotOzakiII(x, y, n);
  case Mode::Fp64:
    return dotFp64(x, y, n);
  }
  return dotFp64(x, y, n);
}

ComplexDot complexDot(const std::complex<double> *a, const std::complex<double> *b,
                      int n, Mode m) {
  ComplexDot out{};
  out.backend = modeName(m);
  out.certified = true;
  if (n <= 0) {
    out.certified = true;
    return out;
  }
  std::vector<double> ar, ai, br, bi;
  packRealImag(a, n, ar, ai);
  packRealImag(b, n, br, bi);
  const RealDot re1 = realDot(ar.data(), br.data(), n, m);
  const RealDot re2 = realDot(ai.data(), bi.data(), n, m);
  const RealDot im1 = realDot(ai.data(), br.data(), n, m);
  const RealDot im2 = realDot(ar.data(), bi.data(), n, m);
  const RealDot nA1 = realDot(ar.data(), ar.data(), n, m);
  const RealDot nA2 = realDot(ai.data(), ai.data(), n, m);
  const RealDot nB1 = realDot(br.data(), br.data(), n, m);
  const RealDot nB2 = realDot(bi.data(), bi.data(), n, m);
  out.re = re1.value + re2.value;
  out.im = im1.value - im2.value;
  out.normA = nA1.value + nA2.value;
  out.normB = nB1.value + nB2.value;
  out.boundRe = upward(static_cast<long double>(re1.bound) + re2.bound +
                       0x1p-53L * std::fabs(static_cast<long double>(out.re)));
  out.boundIm = upward(static_cast<long double>(im1.bound) + im2.bound +
                       0x1p-53L * std::fabs(static_cast<long double>(out.im)));
  out.boundNormA = upward(static_cast<long double>(nA1.bound) + nA2.bound +
                          0x1p-53L * std::fabs(static_cast<long double>(out.normA)));
  out.boundNormB = upward(static_cast<long double>(nB1.bound) + nB2.bound +
                          0x1p-53L * std::fabs(static_cast<long double>(out.normB)));
  out.certified = re1.certified && re2.certified && im1.certified && im2.certified &&
                  nA1.certified && nA2.certified && nB1.certified && nB2.certified;
  out.backend = re1.backend;
  return out;
}

Cij bondCorrelation(const std::complex<double> *a, const std::complex<double> *b,
                    int n, Mode m) {
  const ComplexDot d = complexDot(a, b, n, m);
  Cij out;
  out.backend = d.backend;
  out.certified = false;
  out.bound = std::numeric_limits<double>::infinity();
  const double prod = d.normA * d.normB;
  const double denom = std::sqrt(prod);
  out.value = d.re / denom;
  if (!(denom > 0.0) || !std::isfinite(out.value) || !d.certified) {
    return out;
  }
  const long double nA = d.normA;
  const long double nB = d.normB;
  const long double dProd =
      static_cast<long double>(d.boundNormA) * std::fabs(nB) +
      std::fabs(nA) * static_cast<long double>(d.boundNormB) +
      static_cast<long double>(d.boundNormA) * d.boundNormB;
  const long double dDen = dProd / (2.0L * static_cast<long double>(denom));
  const long double dRe = d.boundRe;
  const long double absC = std::fabs(static_cast<long double>(out.value));
  long double dc = dRe / static_cast<long double>(denom) +
                   absC * dDen / static_cast<long double>(denom);
  dc *= 1.0L + 16.0L * std::numeric_limits<long double>::epsilon();
  out.bound = upward(dc);
  out.certified = std::isfinite(out.bound);
  return out;
}

double meanOfFinite(const double *x, int n, Reduction r) {
  std::vector<double> kept;
  kept.reserve(static_cast<std::size_t>(std::max(n, 0)));
  for (int i = 0; i < n; i++) {
    if (std::isfinite(x[i])) {
      kept.push_back(x[i]);
    }
  }
  if (kept.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return reduce(kept.data(), static_cast<int>(kept.size()), r) /
         static_cast<double>(kept.size());
}

bool int8Gemm(int m, int n, int k, const signed char *A, int lda,
              const signed char *B, int ldb, int *C, int ldc, Int8Backend backend) {
  if (m <= 0 || n <= 0 || k <= 0) {
    return true;
  }
  if (backend == Int8Backend::Auto) {
    if (k >= 64 && m >= 16 && n >= 16 && amxAvailable()) {
      backend = Int8Backend::Amx;
    } else if (k >= 16 && vnniAvailable()) {
      backend = Int8Backend::Vnni;
    } else {
      backend = Int8Backend::Scalar;
    }
  }
#if defined(__x86_64__)
  if (backend == Int8Backend::Amx) {
    return int8GemmAmx(m, n, k, A, lda, B, ldb, C, ldc);
  }
  if (backend == Int8Backend::Vnni) {
    return int8GemmVnni(m, n, k, A, lda, B, ldb, C, ldc);
  }
#else
  (void)backend;
#endif
  int8GemmScalar(m, n, k, A, lda, B, ldb, C, ldc);
  return true;
}

namespace {

void fp64Gemm(int m, int n, int k, const double *A, int lda, const double *B,
              int ldb, double *C, int ldc) {
#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx512f")) {
    gemmAvx512Rows(m, n, k, A, lda, B, ldb, C, ldc);
    return;
  }
#endif
  gemmScalar(m, n, k, A, lda, B, ldb, C, ldc);
}

GemmInfo gemmOzakiII(int m, int n, int k, const double *A, int lda, const double *B,
                     int ldb, double *C, int ldc, double *bounds, Int8Backend backend) {
  GemmInfo info{false, "ozaki-ii", 0, 0};
  const CrtTable &table = crt();
  info.moduli = table.s;
  if (m <= 0 || n <= 0 || k <= 0) {
    info.ok = true;
    return info;
  }
  const int kSum = table.sumBits - ceilLog2(std::max(k, 1));
  if (kSum < 4) {
    return info;
  }
  const int kA = kSum / 2;
  const int kB = kSum - kA;
  std::vector<int64_t> Ap(static_cast<std::size_t>(m) * k);
  std::vector<int64_t> Bp(static_cast<std::size_t>(k) * n);
  std::vector<int> rowExp(static_cast<std::size_t>(m));
  std::vector<int> colExp(static_cast<std::size_t>(n));
  for (int i = 0; i < m; i++) {
    std::vector<int64_t> row;
    int exp = 0;
    if (!truncScaled(A + static_cast<std::ptrdiff_t>(i) * lda, k, kA, row, exp)) {
      return info;
    }
    rowExp[static_cast<std::size_t>(i)] = exp;
    for (int t = 0; t < k; t++) {
      Ap[static_cast<std::size_t>(i) * k + t] = row[static_cast<std::size_t>(t)];
    }
  }
  for (int j = 0; j < n; j++) {
    std::vector<double> col(static_cast<std::size_t>(k));
    for (int t = 0; t < k; t++) {
      col[static_cast<std::size_t>(t)] = B[static_cast<std::ptrdiff_t>(t) * ldb + j];
    }
    std::vector<int64_t> scaled;
    int exp = 0;
    if (!truncScaled(col.data(), k, kB, scaled, exp)) {
      return info;
    }
    colExp[static_cast<std::size_t>(j)] = exp;
    for (int t = 0; t < k; t++) {
      Bp[static_cast<std::size_t>(t) * n + j] = scaled[static_cast<std::size_t>(t)];
    }
  }
  std::vector<i128> acc(static_cast<std::size_t>(m) * n, 0);
  std::vector<signed char> Ai(static_cast<std::size_t>(m) * k);
  std::vector<signed char> Bi(static_cast<std::size_t>(k) * n);
  std::vector<int> partial(static_cast<std::size_t>(m) * n, 0);
  for (int mod = 0; mod < table.s; mod++) {
    const int modulus = static_cast<int>(table.m[mod]);
    for (int i = 0; i < m * k; i++) {
      Ai[static_cast<std::size_t>(i)] = static_cast<signed char>(
          symMod(static_cast<int32_t>(Ap[static_cast<std::size_t>(i)] % modulus), modulus));
    }
    for (int i = 0; i < k * n; i++) {
      Bi[static_cast<std::size_t>(i)] = static_cast<signed char>(
          symMod(static_cast<int32_t>(Bp[static_cast<std::size_t>(i)] % modulus), modulus));
    }
    std::fill(partial.begin(), partial.end(), 0);
    Int8Backend use = backend;
    if (use == Int8Backend::Auto) {
      if (k >= 64 && m >= 16 && n >= 16 && amxAvailable()) {
        use = Int8Backend::Amx;
      } else if (vnniAvailable() && k >= 16 && n >= 16) {
        use = Int8Backend::Vnni;
      } else {
        use = Int8Backend::Scalar;
      }
    }
    if (!int8Gemm(m, n, k, Ai.data(), k, Bi.data(), n, partial.data(), n, use)) {
      return info;
    }
    info.backend = use == Int8Backend::Amx
                       ? "ozaki-ii-amx"
                       : (use == Int8Backend::Vnni ? "ozaki-ii-vnni" : "ozaki-ii-scalar");
    for (int i = 0; i < m * n; i++) {
      const int r = symMod(partial[static_cast<std::size_t>(i)], modulus);
      const i128 term =
          static_cast<i128>(r) * static_cast<i128>(table.weight[mod]);
      acc[static_cast<std::size_t>(i)] =
          remCentered(acc[static_cast<std::size_t>(i)] + term, table.modulus);
    }
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      const int idx = i * n + j;
      const double value =
          roundI128(acc[static_cast<std::size_t>(idx)],
                    -(rowExp[static_cast<std::size_t>(i)] + colExp[static_cast<std::size_t>(j)]));
      C[static_cast<std::ptrdiff_t>(i) * ldc + j] = value;
      if (bounds != nullptr) {
        long double trunc = 0.0;
        const long double invD = std::ldexp(1.0L, -rowExp[static_cast<std::size_t>(i)]);
        const long double invE = std::ldexp(1.0L, -colExp[static_cast<std::size_t>(j)]);
        for (int t = 0; t < k; t++) {
          const double b = B[static_cast<std::ptrdiff_t>(t) * ldb + j];
          trunc += invD * std::fabs(static_cast<long double>(b));
          trunc += std::fabs(static_cast<long double>(Ap[static_cast<std::size_t>(i) * k + t])) *
                   invD * invE;
        }
        trunc *= 1.0L + 16.0L * static_cast<long double>(k) *
                            std::numeric_limits<long double>::epsilon();
        const long double roundErr =
            std::ldexp(std::fabs(static_cast<long double>(value)), -53);
        bounds[idx] = upward(trunc + roundErr);
      }
    }
  }
  info.ok = true;
  info.sliceWidth = kA;
  return info;
}

} // namespace

GemmInfo gemm(int m, int n, int k, const double *A, int lda, const double *B, int ldb,
              double *C, int ldc, Mode which, double *bounds, Int8Backend backend) {
  GemmInfo info{true, "fp64-avx512", 0, 0};
  if (which == Mode::OzakiII) {
    return gemmOzakiII(m, n, k, A, lda, B, ldb, C, ldc, bounds, backend);
  }
  if (which == Mode::Fp64) {
    fp64Gemm(m, n, k, A, lda, B, ldb, C, ldc);
    if (bounds != nullptr) {
      std::fill(bounds, bounds + static_cast<std::ptrdiff_t>(m) * n, 0.0);
    }
    info.backend = "fp64-avx512";
    return info;
  }
  info.backend = modeName(which);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (m * n >= 256)
#endif
  for (int i = 0; i < m; i++) {
    std::vector<double> row(static_cast<std::size_t>(k));
    for (int t = 0; t < k; t++) {
      row[static_cast<std::size_t>(t)] = A[static_cast<std::ptrdiff_t>(i) * lda + t];
    }
    for (int j = 0; j < n; j++) {
      std::vector<double> col(static_cast<std::size_t>(k));
      for (int t = 0; t < k; t++) {
        col[static_cast<std::size_t>(t)] = B[static_cast<std::ptrdiff_t>(t) * ldb + j];
      }
      const RealDot d = realDot(row.data(), col.data(), k, which);
      C[static_cast<std::ptrdiff_t>(i) * ldc + j] = d.value;
      if (bounds != nullptr) {
        bounds[static_cast<std::size_t>(i) * n + j] = d.bound;
      }
    }
  }
  info.ok = true;
  return info;
}

} // namespace seams::product
