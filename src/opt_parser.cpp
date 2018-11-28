#include "opt_parser.h"

// This creates options
cxxopts::ParseResult parse(int argc, char *argv[]) {
  try {
    cxxopts::Options options(
        argv[0], "Structure calculations for molecular simulations");
    options.positional_help("[optional args]").show_positional_help();
    options.allow_unrecognised_options().add_options()(
        "f,file", "File name",
        cxxopts::value<std::string>()->default_value("input/parameters.txt"))(
        "s,script", "Lua Script",
        cxxopts::value<std::string>()->default_value("main.lua"))(
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
