PLUMED
======

Get the d-SEAMS PLUMED plugin v0.4.1 from the
:doc:`PLUMED engine page <../engines/plumed>`.

``engines/plumed/SeamsColvar.cpp`` is the ``SEAMS`` action. The
standalone PLUMED action (``dseams-plumed`` v0.4.1) is indexed in
`PLUMED-NEST egg 26.018 <https://www.plumed-nest.org/eggs/26/018/>`_.
``engines/plumed/plumed.dat`` loads it on ``engines/plumed/ice.xyz``,
restrains ``seams.q6``, and prints ``seams.ice`` and ``seams.ring``.
Coordinates in the xyz file are nanometers. ``calculate`` multiplies
by 10 before the library call. ``ice`` is the cubic fraction and
``ring`` is the fraction of atoms on a primitive six-ring. Only ``q6``
carries derivatives, which is why the deck can restrain it.

.. code-block:: text

   plumed ice 1.000000 ring 1.000000 hexagonal 0 six 128

When ``plumed`` and ``libSeamsPlumed.so`` are both present, ``plumed
driver`` has to print the same ice and ring fractions.

.. figure:: figures/plumed.svg
   :alt: PLUMED SEAMS action: ice fraction 1 and six-ring fraction 1.

   The action's ice and ring components.

What this tells you
-------------------

Every atom is cubic and every atom sits on a six-ring, so both
fractions are 1 and the hexagonal count stays 0. The restraint is
legal only because ``q6`` is the component with derivatives.
