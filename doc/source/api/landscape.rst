The RNA Secondary Structure Landscape
=====================================

Basic Landscape API
-------------------

.. doxygengroup:: landscape
    :no-title:

Neighborhood and Move Sets
--------------------------

Functions to generate structural neighbors of a secondary
structure according to a particular Move Set.

.. doxygengroup:: neighbors
    :no-title:


Folding Paths
-------------

API for various RNA (re-)folding path algorithms that also allow for determination
of saddle points, energy barriers, and local minima.

.. doxygengroup:: paths
    :no-title:


Direct Refolding Paths
----------------------

Heuristics to explore direct, optimal (re-)folding paths between
two secondary structures.

.. doxygengroup:: paths_direct
    :no-title:


Other Folding Paths
-------------------

Implementation of gradient- and random walks starting from
a single secondary structure.

.. doxygengroup:: paths_walk
    :no-title:
