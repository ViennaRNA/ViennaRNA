RNA-RNA interaction as a stepwise process
=========================================

In this approach to cofolding the interaction between two RNA molecules is
seen as a stepwise process. In a first step, the target molecule has to
adopt a structure in which a binding site is accessible. In a second step,
the ligand molecule will hybridize with a region accessible to an
interaction. Consequently the algorithm is designed as a two step process:
The first step is the calculation of the probability
that a region within the target is unpaired, or equivalently, the
calculation of the free energy needed to expose a region. In the second step
we compute the free energy of an interaction for every possible binding site.

.. doxygengroup:: up_cofold
