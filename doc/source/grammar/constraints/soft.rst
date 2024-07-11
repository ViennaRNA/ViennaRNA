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

