#-----------------------------------------------------------------------------------
# d-SEAMS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# A copy of the GNU General Public License is available at
# http://www.gnu.org/licenses/
#-----------------------------------------------------------------------------------

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
