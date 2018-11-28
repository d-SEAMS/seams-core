#include "opt_parser.h"

// This creates options
cxxopts::ParseResult parse(int argc, char *argv[]) {
  try {
    cxxopts::Options options(
        argv[0], "Structure calculations for molecular simulations");
    options.positional_help("[optional args]").show_positional_help();
    options.allow_unrecognised_options().add_options()(
        "f,file", "File name",
        cxxopts::value<std::string>()->default_value("input/parameters.txt"));
    auto result = options.parse(argc, argv);
    if (result.arguments().size() == 0) {
      std::cout << "DO error handling" << std::endl;
    }
    return result;
  } catch (const cxxopts::OptionException &e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}
