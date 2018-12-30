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
