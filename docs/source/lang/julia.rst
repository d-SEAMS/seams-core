Julia
=====

Documenter.jl builds ``bindings/julia/Dseams`` and
``bindings/julia/Seams``. ``checkdocs = :exports`` refuses an export
without a docstring, and every exported function has a ``jldoctest``.

Ice Ih is read with ``type = -1`` and the cutoff graph. The Dseams
doctests assert 128 hexagonal waters and 256 six-rings:

.. code-block:: julia

   f = read_frame(Dseams.ice_ih_path(); type=-1)
   chill_plus(f; graph="cutoff").counts["hexagonal"] == 128
   primitive_rings(f; graph="cutoff", max_depth=6)[6] == 256

The docs job runs ``julia bindings/julia/Dseams/docs/make.jl`` and
``julia bindings/julia/Seams/docs/make.jl`` with ``SEAMS_LANG_LIB`` and
``SEAMS_LIB`` pointed at the build.
