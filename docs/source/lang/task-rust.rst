Rust
====

Each library call returns a ``Result``. The printed histograms are
iterators over the count buffers, and ``Drop`` frees the frame.

.. literalinclude:: ../../examples/lang_tasks.rs
   :language: rust

``docs_lang_rust`` compiles and runs this file. It exits 77 when
``rustc`` is absent.
