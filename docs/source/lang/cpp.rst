C++
===

Doxygen reads ``bindings/cxx/seams_lang.h`` with ``docs/Doxyfile.lang``.
``WARN_AS_ERROR`` is on, so a warning fails the docs job. The public
declarations sit in the groups ``lang_frame``, ``lang_order``,
``lang_rings``, ``lang_cages``, ``lang_more``, and ``lang_errors``.
HTML is written to ``docs/lang-cpp/html``.

``docs/examples/lang_ice.c`` is the Ice Ih census. It prints
``hexagonal 128`` and ``rings6 256`` and returns 0 only when both match.
``meson test docs_lang_ice`` compiles and runs it.
