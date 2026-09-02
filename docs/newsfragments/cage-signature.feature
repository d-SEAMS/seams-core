``cage::Signature`` names a polyhedron by its ring-size census
(``sodalite`` is ``{4:6, 6:8}``). ``cage::findBySignature`` grows
face-sharing primitive rings that close (every edge shared by exactly
two faces). ``seams cages --signature 4:6,6:8`` or a named table
entry (``sodalite|alpha|512|51262|hc|ddc``) prints the cage and atom
counts. Named ``hc`` and ``ddc`` call ``findHC`` / ``findDDC`` so the
vertex sets match those finders; a raw list ``4:6,6:2`` is the
hexagonal prism.
