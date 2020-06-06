# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*FMT][FMT:1]]
{ clangStdenv, fetchJSON, cmake }:
clangStdenv.mkDerivation rec {
  name = "rang-master";
  src = fetchJSON ./rang.src.json;
  nativeBuildInputs = [ cmake ];
  meta = {
    description = "A Minimal, Header only Modern c++ library for terminal goodies";
    homepage = https://agauniyal.github.io/rang/;
  };
}
# FMT:1 ends here
