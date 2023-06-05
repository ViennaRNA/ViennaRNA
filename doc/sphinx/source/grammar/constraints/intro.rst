============
Introduction
============

.. contents:: Table of Contents
   :local:
   :depth: 1

Secondary Structure constraints can be subdivided into two groups:

  * :ref:`hard_constraints`, and
  * :ref:`soft_constraints`.

While Hard-Constraints directly influence the production rules used in the folding
recursions by allowing, disallowing, or enforcing certain decomposition steps,
Soft-constraints on the other hand are used to change position specific contributions
in the recursions by adding bonuses/penalties in form of pseudo free energies
to certain loop configurations.


Secondary structure constraints are always applied at decomposition level, i.e.
in each step of the recursive structure decomposition, for instance during MFE
prediction. Below is a visualization of the decomposition scheme

.. image:: /gfx/recursions.svg

For :ref:`hard_constraints` the following option flags may be used to constrain
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

However, for :ref:`soft_constraints` we do not allow for simple loop type
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


.. doxygengroup:: constraints

