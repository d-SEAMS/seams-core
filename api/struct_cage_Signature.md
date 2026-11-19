<!-- https://docs.dseams.info/api/struct_cage_Signature.md -->

# struct cage::Signature {#struct-cage-signature}

## Overview {#overview}

Ring-size census of a polyhedral cage: size -> number of faces. [More…](#details-structcage-1-1-signature)

## Detailed Documentation {#detailed-documentation}

Ring-size census of a polyhedral cage: size -> number of faces.

Parse a comma list (`4:6,6:8`) or a named table entry (`sodalite`, `alpha`, `512`, `51262`, `51264`, `51268`, `sh`, `sH`, `hc`, `ddc`). Counts of a repeated size add. Every size and count must be positive. The names `hc` and `ddc` keep the TUM ice finders (five hexagons and seven hexagons); a raw list `4:6,6:2` is the hexagonal prism.
