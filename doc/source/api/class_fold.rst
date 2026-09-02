Classified Dynamic Programming Variants
=======================================

This part describes algorithms that partition the structure space *a priori* and
compute MFE or partition function for each of the resulting partitions. This
approach directly leads to *Classified Dynamic Programming*.

Basic Interface
---------------

.. doxygengroup:: class_fold


Distance Based Partitioning
---------------------------

The secondary structure space is divided into partitions according to the base
pair distance to two given reference structures and all relevant properties are
calculated for each of the resulting partitions.

.. admonition:: See also...

  For further details, we refer to :cite:t:`lorenz:2009`

General
^^^^^^^

.. doxygengroup:: kl_neighborhood
    :no-title:


MFE Variants
^^^^^^^^^^^^

Compute the minimum free energy (MFE) and secondary structures for a partitioning
of the secondary structure space according to the base pair distance to two fixed
reference structures basepair distance to two fixed reference structures.

.. doxygengroup:: kl_neighborhood_mfe
    :no-title:


Partition Function Variants
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute the partition function and stochastically sample secondary structures for
a partitioning of the secondary structure space according to the base pair distance
to two fixed reference structures.

.. doxygengroup:: kl_neighborhood_pf
    :no-title:


Stochastic Backtracking
^^^^^^^^^^^^^^^^^^^^^^^

Functions related to stochastic backtracking from a specified distance class.

.. doxygengroup:: kl_neighborhood_stochbt
    :no-title:


Density of States
-----------------

.. doxygengroup:: dos
    :no-title:

