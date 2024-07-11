Hard Constraints
================

This module covers all functionality for hard constraints in secondary
structure prediction.

.. contents:: Table of Contents
    :local:
    :depth: 2


Introduction
------------

Hard constraints as implemented in our library can be specified for individual loop
types, i.e. the atomic derivations of the RNA folding grammar rules. Hence, the pairing
behavior of both, single nucleotides and pairs of bases, can be constrained in every
loop context separately. Additionally, an abstract implementation using a callback
mechanism allows for full control of more complex hard constraints.


Hard Constraints API
--------------------

.. doxygengroup:: hard_constraints
    :no-title:

