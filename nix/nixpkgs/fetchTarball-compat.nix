# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*Using%20the%20JSON][Using the JSON:2]]
# fetchTarball version that is compatible between all the versions of Nix
{ url, sha256 }@attrs:
let
  inherit (builtins) lessThan nixVersion fetchTarball;
in
if lessThan nixVersion "1.12" then
  fetchTarball { inherit url; }
else
  fetchTarball attrs
# Using the JSON:2 ends here
