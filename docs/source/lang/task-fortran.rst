Fortran
=======

The program binds the language front with ``iso_c_binding``. Every call
returns a status code, and the frame pointer is released on the success
path and on the error path.

.. literalinclude:: ../../examples/lang_tasks.f90
   :language: fortran

``docs_lang_fortran`` compiles and runs this file. It exits 77 when
``gfortran`` is absent.
