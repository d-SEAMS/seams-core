``cage::Signature`` names a polyhedron by its ring-size census
(``sodalite`` is ``{4:6, 6:8}``). ``cage::findBySignature`` grows
face-sharing primitive rings that close (every edge shared by exactly
two faces). ``seams cages --signature 4:6,6:8`` or a named table
entry (``sodalite|alpha|512|51262|hc|ddc``) prints the cage and atom
counts. HC ``{4:6, 6:2}`` and DDC ``{6:7}`` match ``findHC`` /
``findDDC`` vertex sets on the in-tree fixtures.
