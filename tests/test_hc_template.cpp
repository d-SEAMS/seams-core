#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <pntCorrespondence.hpp>
#include <ring.hpp>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

// The HCbasal cage reads templates/hc.xyz from the working directory and
// returns those coordinates. Discarding the read left a matrix with no rows.

TEST_CASE("HCbasal cage keeps the template coordinates", "[pntCorrespondence]") {
  namespace fs = std::filesystem;
  const fs::path file = "templates/hc.xyz";
  fs::create_directories(file.parent_path());
  std::string previous;
  const bool had = fs::exists(file);
  if (had) {
    std::ifstream in(file);
    previous.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
  }
  struct Guard {
    fs::path file;
    bool had;
    std::string previous;
    ~Guard() {
      std::error_code ignore;
      if (!had) {
        fs::remove(file, ignore);
        return;
      }
      std::ofstream out(file);
      out << previous;
    }
  } guard{file, had, std::move(previous)};
  {
    std::ofstream out(file);
    REQUIRE(out.good());
    out << "2\ncomment\n"
        << "O 1.25 2.5 3.75\n"
        << "O 4.5 5.5 6.5\n";
  }

  const Eigen::MatrixXd pts = pntToPnt::getPointSetCage(ring::strucType::HCbasal);
  REQUIRE(pts.rows() == 2);
  REQUIRE(pts.cols() == 3);
  REQUIRE(pts(0, 0) == Catch::Approx(1.25));
  REQUIRE(pts(0, 1) == Catch::Approx(2.5));
  REQUIRE(pts(0, 2) == Catch::Approx(3.75));
  REQUIRE(pts(1, 0) == Catch::Approx(4.5));
  REQUIRE(pts(1, 1) == Catch::Approx(5.5));
  REQUIRE(pts(1, 2) == Catch::Approx(6.5));
}
