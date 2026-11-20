#include <chemfiles.hpp>

#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "usage: derive_confined_water_fixture INPUT.gro INPUT.xtc "
                 "OUTPUT.gro OUTPUT.xtc\n";
    return 2;
  }
  try {
    chemfiles::Trajectory input(argv[2]);
    input.set_topology(argv[1]);
    chemfiles::Trajectory gro(argv[3], 'w', "GRO");
    chemfiles::Trajectory xtc(argv[4], 'w', "XTC");
    for (std::size_t step = 0; step < 3; ++step) {
      auto frame = input.read();
      frame.set("xtc_precision", 1000.0);
      gro.write(frame);
      xtc.write(frame);
    }
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << error.what() << "\n";
    return EXIT_FAILURE;
  }
}
