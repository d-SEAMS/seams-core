//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <opt_parser.h>

#include <iostream>

cxxopts::ParseResult parse(int argc, char *argv[]) {
  try {
    cxxopts::Options options(
        argv[0], "Structure calculations for molecular simulations");
    options.positional_help("[optional args]").show_positional_help();
    options.allow_unrecognised_options().add_options()(
        "c,config", "Yaml Config",
        cxxopts::value<std::string>()->default_value("conf.yml"))("h,help",
                                                                 "Print help");
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << options.help({"", "Group"}) << "\n";
      std::exit(0);
    }
    return result;
  } catch (const cxxopts::OptionException &e) {
    std::cout << "error parsing options: " << e.what() << "\n";
    std::exit(1);
  }
}
