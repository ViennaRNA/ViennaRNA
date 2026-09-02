Classified Dynamic Programming Variants
=======================================

Usually, thermodynamic properties using the basic recursions for :doc:`/advanced/mfe`,
:doc:`/advanced/partfunc`, and so forth, are computed over the entire structure space.
However, sometimes it is desired to partition the structure space *a priori* and compute
the above properties for each of the resulting partitions. This approach directly leads
to **Classified Dynamic Programming**.


Distance Based Partitioning
---------------------------

The secondary structure space is divided into partitions according to the base
pair distance to two given reference structures and all relevant properties are
calculated for each of the resulting partitions.

.. admonition:: See also...

  For further details, we refer to :cite:t:`lorenz:2009`

MFE Variants
^^^^^^^^^^^^

Compute the minimum free energy (MFE) and secondary structures for a partitioning
of the secondary structure space according to the base pair distance to two fixed
reference structures basepair distance to two fixed reference structures.



Partition Function Variants
^^^^^^^^^^^^^^^^^^^^^^^^^^^


Compute the partition function and stochastically sample secondary structures for
a partitioning of the secondary structure space according to the base pair distance
to two fixed reference structures.


Stochastic Backtracking
^^^^^^^^^^^^^^^^^^^^^^^

Functions related to stochastic backtracking from a specified distance class.


Density of States
-----------------

