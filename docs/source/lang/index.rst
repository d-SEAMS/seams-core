Language APIs
=============

seams-core documents the C API and the C++ front here. Every other
language binding is its own package with its own documentation.
Each one links ``libseams`` or ``libseams_lang``.

.. list-table::
   :header-rows: 1

   * - Language
     - Repository
     - Package
   * - Python
     - `d-SEAMS/PydSEAMSlib <https://github.com/d-SEAMS/PydSEAMSlib>`_
     - ``pydseamslib`` (``import dseams``)
   * - R
     - `d-SEAMS/dseams-r <https://github.com/d-SEAMS/dseams-r>`_
     - ``dseams``
   * - Julia
     - `d-SEAMS/Dseams.jl <https://github.com/d-SEAMS/Dseams.jl>`_
     - ``Dseams``
   * - Rust
     - `d-SEAMS/dseams-rs <https://github.com/d-SEAMS/dseams-rs>`_
     - ``dseams-sys``, ``dseams``
   * - Fortran
     - `d-SEAMS/dseams-fortran <https://github.com/d-SEAMS/dseams-fortran>`_
     - fpm ``dseams``
   * - Lua, Fennel
     - `d-SEAMS/yodaStruct <https://github.com/d-SEAMS/yodaStruct>`_
     - LuaRocks package ``dseams``

.. toctree::
   :maxdepth: 1

   tasks
   cpp
