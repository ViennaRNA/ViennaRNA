RNA-RNA Interaction
===================

The function of an RNA molecule often depends on its interaction with
other RNAs. The following routines therefore allows one to predict structures
formed by two RNA molecules upon hybridization.


Concatenating RNA sequences
---------------------------

One approach to co-folding two RNAs consists of concatenating the two
sequences and keeping track of the concatenation point in all energy
evaluations. Correspondingly, many of the cofold() and
co_pf_fold() routines take one sequence string as argument
and use the the global variable #cut_point to mark the concatenation
point. Note that while the *RNAcofold* program uses the ``&`` character
to mark the chain break in its input.

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

After computing the partition functions of all possible dimers one
can compute the probabilities of base pairs, the concentrations out of
start concentrations and sofar and soaway.

Dimer formation is inherently concentration dependent. Given the free
energies of the monomers A and B and dimers AB, AA, and BB one can compute
the equilibrium concentrations, given input concentrations of A and B, see
e.g. Dimitrov & Zuker (2004)


RNA-RNA interaction as a Stepwise Process
-----------------------------------------

In this approach to cofolding the interaction between two RNA molecules is
seen as a stepwise process. In a first step, the target molecule has to
adopt a structure in which a binding site is accessible. In a second step,
the ligand molecule will hybridize with a region accessible to an
interaction. Consequently the algorithm is designed as a two step process:
The first step is the calculation of the probability
that a region within the target is unpaired, or equivalently, the
calculation of the free energy needed to expose a region. In the second step
we compute the free energy of an interaction for every possible binding site.

This approach is implemented for the interaction between a long
target sequence and a short ligand RNA. Function pf_unstru() calculates
the partition function over all unpaired regions in the input
sequence. Function pf_interact(), which calculates the
partition function over all possible interactions between two
sequences, needs both sequence as separate strings as input.

