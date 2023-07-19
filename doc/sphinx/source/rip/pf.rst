Partition Function for Two Hybridized Sequences
===============================================

To simplify the implementation the partition function computation is done
internally in a null model that does not include the duplex initiation
energy, i.e. the entropic penalty for producing a dimer from two
monomers). The resulting free energies and pair probabilities are initially
relative to that null model. In a second step the free energies can be
corrected to include the dimerization penalty, and the pair probabilities
can be divided into the conditional pair probabilities given that a re
dimer is formed or not formed. See :cite:t:`bernhart:2006` for further details.

As for folding one RNA molecule, this computes the partition function
of all possible structures and the base pair probabilities. Uses the
same global #pf_scale variable to avoid overflows.

After computing the partition functions of all possible dimeres one
can compute the probabilities of base pairs, the concentrations out of
start concentrations and sofar and soaway.

Dimer formation is inherently concentration dependent. Given the free
energies of the monomers A and B and dimers AB, AA, and BB one can compute
the equilibrium concentrations, given input concentrations of A and B, see
e.g. Dimitrov & Zuker (2004)

.. doxygengroup:: pf_cofold
    :content-only:
