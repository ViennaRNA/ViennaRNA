Soft Constraints
================

Functions and data structures for secondary structure soft constraints.

.. contents:: Table of Contents
    :local:
    :depth: 2


Introduction
------------

Soft-constraints are used to change position specific contributions
in the recursions by adding bonuses/penalties in form of pseudo free energies
to certain loop configurations.

.. note::
  For the sake of memory efficiency, we do not implement a loop context aware version of
  soft constraints. The *static* soft constraints as implemented only distinguish unpaired
  from paired nucleotides. This is usually sufficient for most use-case scenarios.
  However, similar to hard constraints, an abstract soft constraints implementation using
  a callback mechanism exists, that allows for any soft constraint that is compatible with
  the RNA folding grammar. Thus, loop contexts and even individual derivation rules can
  be addressed separately for maximum flexibility in soft-constraints application.


Energy Precision and Accumulation
--------------------------------

For single-sequence folding, :c:func:`vrna_sc_add_up` and
:c:func:`vrna_sc_add_bp` accumulate pseudo free energy contributions for each
unpaired nucleotide or base pair without first rounding them to 0.01 kcal/mol.
The corresponding :c:func:`vrna_sc_set_up` and :c:func:`vrna_sc_set_bp`
functions replace the stored contributions and also retain fractional values.

Partition function calculations use the stored values without rounding to a
fixed energy increment. This allows ensemble energies and base-pair
probabilities to respond continuously to these soft constraints. MFE
calculations retain the existing integer energy scale: each accumulated
nucleotide or base-pair contribution is rounded to the nearest 0.01 kcal/mol.
For an unpaired stretch, the rounded contributions of its nucleotides are
summed.

For example, adding 0.003 kcal/mol twice to the same nucleotide gives a stored
contribution of 0.006 kcal/mol and an MFE contribution of 0.01 kcal/mol:

.. code-block:: python

    import RNA

    fc = RNA.fold_compound("AAAAAA")
    fc.sc_add_up(3, 0.003)
    fc.sc_add_up(3, 0.003)
    print(f"{fc.mfe()[1]:.2f}")  # 0.01 kcal/mol
    print(f"{fc.pf()[1]:.3f}")   # 0.006 kcal/mol

These rules also apply to negative contributions. The floating-point
precision of PF tables still depends on whether ViennaRNA was configured
with ``--enable-floatpf``. Comparative folding APIs and stack soft constraints
retain their existing rounding behavior.

MFE rounding now uses the stored double value directly. Values very close to
a rounding boundary may therefore produce a different result from older
versions, which first converted the contribution to float. For example,
0.0049999999 kcal/mol rounds to 0.00 kcal/mol. This also applies to soft
constraints derived from probing data.


Common API symbols
------------------

.. doxygengroup:: soft_constraints
    :no-title:

Constraints for Unpaired Positions
----------------------------------

.. doxygengroup:: soft_constraints_up
    :no-title:


Constraints for Base Pairs
--------------------------

.. doxygengroup:: soft_constraints_bp
    :no-title:


Constraints for Stacked Base Pairs
----------------------------------

.. doxygengroup:: soft_constraints_st
    :no-title:

Generic implementation
----------------------

.. doxygengroup:: soft_constraints_generic
    :no-title:

