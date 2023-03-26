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

# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*Helper%20Functions][Helper Functions:1]]
{ fetchFromGitHub }:

# fetches the source described in the <package>.src.json

# path to the .src.json
path:
let
  data = builtins.fromJSON (builtins.readFile path);
in
  fetchFromGitHub { inherit (data) owner repo rev sha256; }
# Helper Functions:1 ends here
