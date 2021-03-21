#-----------------------------------------------------------------------------------
# d-SEAMS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# A copy of the GNU General Public License is available at
# http://www.gnu.org/licenses/
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
