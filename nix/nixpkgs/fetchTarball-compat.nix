#-----------------------------------------------------------------------------------
# d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
#
# Copyright (c) 2018--present d-SEAMS core team
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the MIT License as published by
# the Open Source Initiative.
#
# A copy of the MIT License is included in the LICENSE file of this repository.
# You should have received a copy of the MIT License along with this program.
# If not, see <https://opensource.org/licenses/MIT>.
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
