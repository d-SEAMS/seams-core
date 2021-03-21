#-----------------------------------------------------------------------------------
# d-SEAMS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# A copy of the GNU General Public License is available at
# http://www.gnu.org/licenses/
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
