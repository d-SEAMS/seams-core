//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include "opt_parser.h"

// This creates options
cxxopts::ParseResult parse(int argc, char *argv[]) {
  try {
    cxxopts::Options options(
        argv[0], "Structure calculations for molecular simulations");
    options.positional_help("[optional args]").show_positional_help();
    options.allow_unrecognised_options().add_options()(
        "c,config", "Yaml Config",
        cxxopts::value<std::string>()->default_value("conf.yml"))

        ("h,help", "Print help");
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << options.help({"", "Group"}) << std::endl;
      exit(0);
    }
    // No options
    if (result.arguments().size() == 0) {
      std::cout << "DO error handling" << std::endl;
    }
    return result;
  } catch (const cxxopts::OptionException &e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}
