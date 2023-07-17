Distance Based Partitioning of the Secondary Structure Space
============================================================

The secondary structure space is divided into partitions according to the base
pair distance to two given reference structures and all relevant properties are
calculated for each of the resulting partitions.

.. admonition:: See also...

  For further details, we refer to :cite:t:`lorenz:2009`

.. contents:: Table of Contents
    :local:

General
-------

.. doxygengroup:: kl_neighborhood
   :content-only:


MFE Variants
------------

Compute the minimum free energy (MFE) and secondary structures for a partitioning
of the secondary structure space according to the base pair distance to two fixed
reference structures basepair distance to two fixed reference structures.

.. doxygengroup:: kl_neighborhood_mfe


Partition Function Variants
---------------------------


Compute the partition function and stochastically sample secondary structures for
a partitioning of the secondary structure space according to the base pair distance
to two fixed reference structures.

.. doxygengroup:: kl_neighborhood_pf


Stochastic Backtracking
-----------------------

Functions related to stochastic backtracking from a specified distance class.

.. doxygengroup:: kl_neighborhood_stochbt
