===============================
Secondary Structure Constraints
===============================

Secondary structure constraints provide an easy control of which structures
the prediction algorithms actually include into their solution space and how
these structures are evaluated.

.. toctree::
   :maxdepth: 1

   constraints/intro
   constraints/hard
   constraints/soft


High Level Constraints Interfaces
---------------------------------

Simplified interfaces to the soft constraints framework can be obtained
by the implementations in the submodules

  * :ref:`SHAPE_reactivities` and
  * :ref:`constraints_ligand`.

An implementation that generates soft constraints for unpaired nucleotides
by minimizing the discrepancy between their predicted and expected pairing
probability is available in submodule :ref:`perturbation`.

General API symbols
-------------------

.. doxygengroup:: constraints

