C++
===

``Frame`` owns the language-front pointer and frees it when the object
leaves scope, including when a call throws.

.. literalinclude:: ../../examples/lang_tasks_cpp.cpp
   :language: cpp

Compile this file as C++17 and link ``libseams``.
