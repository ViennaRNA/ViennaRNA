===============================
Secondary Structure Constraints
===============================

Secondary structure constraints provide an easy control of which structures
the prediction algorithms actually include into their solution space and how
these structures are evaluated.

.. toctree::
   :maxdepth: 1
   :caption: Specialized Modules:

   constraints/hard
   constraints/soft

Introduction
------------

Secondary Structure constraints can be subdivided into two groups:

  * :doc:`/grammar/constraints/hard`, and
  * :doc:`/grammar/constraints/soft`.

While *hard constraints* directly influence the production rules used in the folding
recursions by allowing, disallowing, or enforcing certain decomposition steps,
*soft constraints* are used to change position specific contributions
in the recursions by adding bonuses/penalties in form of pseudo free energies
to certain loop configurations.


Secondary structure constraints are always applied at decomposition level, i.e.
in each step of the recursive structure decomposition, for instance during MFE
prediction. Below is a visualization of the decomposition scheme

.. image:: /gfx/recursions.svg

For :doc:`/grammar/constraints/hard` the following option flags may be used to constrain
the pairing behavior of single, or pairs of nucleotides:

  * :c:macro:`VRNA_CONSTRAINT_CONTEXT_EXT_LOOP`
  * :c:macro:`VRNA_CONSTRAINT_CONTEXT_HP_LOOP`
  * :c:macro:`VRNA_CONSTRAINT_CONTEXT_INT_LOOP`
  * :c:macro:`VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC`
  * :c:macro:`VRNA_CONSTRAINT_CONTEXT_MB_LOOP`
  * :c:macro:`VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC`
  * :c:macro:`VRNA_CONSTRAINT_CONTEXT_ENFORCE`
  * :c:macro:`VRNA_CONSTRAINT_CONTEXT_NO_REMOVE`
  * :c:macro:`VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS`

However, for :doc:`/grammar/constraints/soft` we do not allow for simple loop type
dependent constraining. But soft constraints are equipped with generic constraint
support. This enables the user to pass arbitrary callback functions that
return auxiliary energy contributions for evaluation the evaluation of any
decomposition.

The callback will then always be notified about the type of decomposition
that is happening, and the corresponding delimiting sequence positions. The
following decomposition steps are distinguished, and should be captured by
the user's implementation of the callback:

  * :c:macro:`VRNA_DECOMP_PAIR_HP`
  * :c:macro:`VRNA_DECOMP_PAIR_IL`
  * :c:macro:`VRNA_DECOMP_PAIR_ML`
  * :c:macro:`VRNA_DECOMP_ML_ML_ML`
  * :c:macro:`VRNA_DECOMP_ML_STEM`
  * :c:macro:`VRNA_DECOMP_ML_ML`
  * :c:macro:`VRNA_DECOMP_ML_UP`
  * :c:macro:`VRNA_DECOMP_ML_ML_STEM`
  * :c:macro:`VRNA_DECOMP_ML_COAXIAL`
  * :c:macro:`VRNA_DECOMP_EXT_EXT`
  * :c:macro:`VRNA_DECOMP_EXT_UP`
  * :c:macro:`VRNA_DECOMP_EXT_STEM`
  * :c:macro:`VRNA_DECOMP_EXT_EXT_EXT`
  * :c:macro:`VRNA_DECOMP_EXT_STEM_EXT`
  * :c:macro:`VRNA_DECOMP_EXT_STEM_OUTSIDE`
  * :c:macro:`VRNA_DECOMP_EXT_EXT_STEM`
  * :c:macro:`VRNA_DECOMP_EXT_EXT_STEM1`


General API symbols
-------------------

.. doxygengroup:: constraints
    :no-title:

High Level Constraints Interfaces
---------------------------------

High-level interfaces that build upon the soft constraints framework
can be obtained by the implementations in the submodules:

* :doc:`/modified_bases`
* :doc:`/probing/SHAPE`
* :doc:`/ligands/constraints`

An implementation that generates soft constraints for unpaired nucleotides
by minimizing the discrepancy between their predicted and expected pairing
probability is available in submodule :doc:`/probing/perturbation`.

