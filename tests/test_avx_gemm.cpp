#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <ozaki.hpp>

#include <vector>

// The AVX-512 row kernel is selected for 32 or more rows. Its product has to
// match the scalar product, including when the OpenMP region runs serially.

TEST_CASE("wide fp64 gemm matches the scalar product", "[ozaki]") {
  constexpr int m = 40;
  constexpr int n = 8;
  constexpr int k = 6;
  std::vector<double> a(static_cast<std::size_t>(m * k));
  std::vector<double> b(static_cast<std::size_t>(k * n));
  std::vector<double> c(static_cast<std::size_t>(m * n), 0.0);
  for (int i = 0; i < m * k; ++i) {
    a[static_cast<std::size_t>(i)] = static_cast<double>((i % 7) - 3);
  }
  for (int i = 0; i < k * n; ++i) {
    b[static_cast<std::size_t>(i)] = 0.5 * static_cast<double>(i % 5);
  }

  const auto info = seams::product::gemm(
      m, n, k, a.data(), k, b.data(), n, c.data(), n, seams::product::Mode::Fp64,
      nullptr);
  REQUIRE(info.ok);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (int t = 0; t < k; ++t) {
        sum += a[static_cast<std::size_t>(i * k + t)] *
               b[static_cast<std::size_t>(t * n + j)];
      }
      REQUIRE(c[static_cast<std::size_t>(i * n + j)] == Catch::Approx(sum).margin(1e-9));
    }
  }
}

TEST_CASE("vnni int8 tile matches the scalar product", "[ozaki]") {
  if (!seams::product::vnniAvailable()) {
    return;
  }
  constexpr int m = 32;
  constexpr int n = 32;
  constexpr int k = 16;
  std::vector<signed char> a(static_cast<std::size_t>(m * k));
  std::vector<signed char> b(static_cast<std::size_t>(k * n));
  std::vector<int> fast(static_cast<std::size_t>(m * n), 0);
  std::vector<int> slow(static_cast<std::size_t>(m * n), 0);
  for (int i = 0; i < m * k; ++i) {
    a[static_cast<std::size_t>(i)] = static_cast<signed char>((i % 11) - 5);
  }
  for (int i = 0; i < k * n; ++i) {
    b[static_cast<std::size_t>(i)] = static_cast<signed char>((i % 9) - 4);
  }
  REQUIRE(seams::product::int8Gemm(m, n, k, a.data(), k, b.data(), n, fast.data(), n,
                                   seams::product::Int8Backend::Vnni));
  REQUIRE(seams::product::int8Gemm(m, n, k, a.data(), k, b.data(), n, slow.data(), n,
                                   seams::product::Int8Backend::Scalar));
  REQUIRE(fast == slow);
}
