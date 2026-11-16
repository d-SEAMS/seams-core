C
=

The C front returns a status code from every call. The program owns the
frame, the neighbour list, the ring handle, and the label buffer, and it
frees each of them.

.. literalinclude:: ../../examples/lang_tasks_c.c
   :language: c

``docs_lang_c`` compiles and runs this file.
