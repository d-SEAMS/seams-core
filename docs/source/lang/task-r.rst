R
=

``read_frame``, ``chill_plus``, and ``primitive_rings`` are the R front.
The package returns positions rather than the neighbour list, so the
example counts directed cutoff pairs of those positions under the
minimum image of the dump's box. The ring table omits empty sizes; the
printed histogram still has one integer per size from 0 through 6.

.. literalinclude:: ../../examples/lang_tasks.R
   :language: r

``docs_lang_r`` runs this file. It exits 77 when R or the installed
package is absent.
