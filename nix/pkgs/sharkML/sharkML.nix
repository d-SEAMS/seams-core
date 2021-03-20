# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*SharkML][SharkML:1]]
# A standard cmake-based build
{ clangStdenv, fetchJSON, cmake, boost, openblas, liblapack }:
clangStdenv.mkDerivation {
  name = "sharkML-master";
  src = fetchJSON ./sharkML.src.json;
  nativeBuildInputs = [ cmake boost openblas liblapack ];
}
# SharkML:1 ends here
