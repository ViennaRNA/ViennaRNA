===============================
Secondary Structure Constraints
===============================

Introduction
------------

Secondary Structure constraints can be subdivided into two groups:

* :ref:`advanced-constraints-hard`, and
* :ref:`advanced-constraints-soft`.

While **hard constraints** directly influence the production rules used in the folding
recursions by allowing, disallowing, or enforcing certain decomposition steps,
**soft constraints** are used to change position specific contributions
in the recursions by adding bonuses/penalties in form of pseudo free energies
to certain loop configurations.


.. note::

  Secondary structure constraints are always applied at decomposition level, i.e.
  in each step of the recursive structure decomposition, for instance during MFE
  prediction. Below is a visualization of the decomposition scheme

  .. image:: /gfx/recursions.png


.. seealso::

   :ref:`api-grammar-constraints` lists all related constraints related API symbols.


.. _advanced-constraints-hard:

Hard Constraints
----------------

Hard constraints as implemented in our library can be specified for individual loop
types, i.e. the atomic derivations of the RNA folding grammar rules. Hence, the pairing
behavior of both, single nucleotides and pairs of bases, can be constrained in every
loop context separately. Additionally, an abstract implementation using a callback
mechanism allows for full control of more complex hard constraints.

The following option flags may be used to constrain the pairing behavior of single,
or pairs of nucleotides:

* :c:macro:`VRNA_CONSTRAINT_CONTEXT_EXT_LOOP`
* :c:macro:`VRNA_CONSTRAINT_CONTEXT_HP_LOOP`
* :c:macro:`VRNA_CONSTRAINT_CONTEXT_INT_LOOP`
* :c:macro:`VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC`
* :c:macro:`VRNA_CONSTRAINT_CONTEXT_MB_LOOP`
* :c:macro:`VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC`
* :c:macro:`VRNA_CONSTRAINT_CONTEXT_ENFORCE`
* :c:macro:`VRNA_CONSTRAINT_CONTEXT_NO_REMOVE`
* :c:macro:`VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS`

.. seealso::

    A full list of available functions and data structures for hard constraints can
    be found under :ref:`api-grammar-constraints-hard`.


.. _advanced-constraints-soft:

Soft Constraints
----------------

Soft-constraints are used to change position specific contributions in the recursions
by adding bonuses/penalties in form of pseudo free energies to certain loop configurations.

For the sake of memory efficiency, we do not implement a loop context aware version of
soft constraints. The *static* soft constraints as implemented only distinguish unpaired
from paired nucleotides. This is usually sufficient for most use-case scenarios.
However, similar to hard constraints, an abstract soft constraints implementation using
a callback mechanism exists, that allows for any soft constraint that is compatible with
the RNA folding grammar. Thus, loop contexts and even individual decomposition rules can
be addressed separately for maximum flexibility in soft-constraints application.

Soft constraint callback functions will always be notified about the type of decomposition
that is happening, and the corresponding delimiting sequence positions. The following
decomposition steps are distinguished, and should be captured by the user's implementation
of the callback:

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


.. seealso::

    A full list of available functions and data structures for soft constraints can
    be found under :ref:`api-grammar-constraints-soft`.



High Level Constraints Interfaces
---------------------------------

High-level interfaces that build upon the soft constraints framework
can be obtained by the implementations in the submodules:

* :doc:`/modified_bases`
* :doc:`/probing_data`
* :doc:`/probing/SHAPE`
* :doc:`/ligands/constraints`

An implementation that generates soft constraints for unpaired nucleotides
by minimizing the discrepancy between their predicted and expected pairing
probability is available in submodule :doc:`/probing/perturbation`.

