# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*YAML%20Cpp][YAML Cpp:1]]
# A standard cmake-based build
{ clangStdenv, fetchJSON, cmake }:
clangStdenv.mkDerivation {
  name = "yaml-cpp-master";
  src = fetchJSON ./yaml-cpp.src.json;
  nativeBuildInputs = [ cmake ];
}
# YAML Cpp:1 ends here
