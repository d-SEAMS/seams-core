#-----------------------------------------------------------------------------------
# d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
#
# Copyright (c) 2019--present d-SEAMS core team
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the MIT License as published by
# the Open Source Initiative.
#
# A copy of the MIT License is included in the LICENSE file of this repository.
# You should have received a copy of the MIT License along with this program.
# If not, see <https://opensource.org/licenses/MIT>.
#-----------------------------------------------------------------------------------

# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*Project%20Root][Project Root:1]]
  # Define
  let
  # Import
    buildpkgs = import ./nix {};
  in
  # Pass to
  buildpkgs.yodaStruct
# Project Root:1 ends here
