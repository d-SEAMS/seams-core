#include <chrono>
#include <iostream>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

int main() {
  std::string traj = "traj/exampleTraj.lammpstrj";
  auto cloud = sinp::readLammpsTrjreduced(traj, 1, nullptr, 2, false, {0,0,0}, {0,0,0});

  auto start = std::chrono::high_resolution_clock::now();
  for (int rep = 0; rep < 100; rep++) {
    auto nList = nneigh::neighListO(3.5, &cloud, 2);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "100 iterations of neighListO: " << ms << " ms" << std::endl;
  std::cout << "Per iteration: " << ms / 100.0 << " ms" << std::endl;
  return 0;
}
