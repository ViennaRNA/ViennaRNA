.. _api-grammar-constraints:

Secondary Structure Constraints
===============================

Secondary structure constraints provide an easy control of which structures
the prediction algorithms actually include into their solution space and how
these structures are evaluated.


General API symbols
-------------------

.. doxygengroup:: constraints
    :no-title:


.. _api-grammar-constraints-hard:

Hard Constraints API
--------------------

Functions and data structures for secondary structure hard constraints

.. doxygengroup:: hard_constraints
    :no-title:


.. _api-grammar-constraints-soft:

Soft Constraints API
--------------------

Functions and data structures for secondary structure soft constraints.

Common API symbols
^^^^^^^^^^^^^^^^^^

.. doxygengroup:: soft_constraints
    :no-title:

Constraints for Unpaired Positions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: soft_constraints_up
    :no-title:


Constraints for Base Pairs
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: soft_constraints_bp
    :no-title:


Constraints for Stacked Base Pairs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: soft_constraints_st
    :no-title:

Generic implementation
^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: soft_constraints_generic
    :no-title:

