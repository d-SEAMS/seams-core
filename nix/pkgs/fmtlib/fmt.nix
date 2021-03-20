# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*FMT][FMT:1]]
{ clangStdenv, fetchJSON, cmake }:
clangStdenv.mkDerivation rec {
  name = "fmtlib-master";
  src = fetchJSON ./fmt.src.json;
  nativeBuildInputs = [ cmake ];
  meta = {
    description = "{fmt} is an open-source formatting library for C++. It can be used as a safe and fast alternative to (s)printf and IOStreams.";
    homepage = http://fmtlib.net;
  };
}
# FMT:1 ends here
