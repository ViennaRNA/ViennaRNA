Partition Function and Equilibrium Probability
==============================================

Introduction
------------

The principal idea behind computing the partition function for RNA secondary structures
can be summarized as follows:

  In equilibrium, statistical mechanics and polymer theory tells us that the frequency
  or probability :math:`p(s)` of a particular state :math:`s` depends on its energy
  :math:`E(s)` and follows a Boltzmann distribution.

This can be formally expressed as

.. math::

  p(s) \propto e^{-\beta E(s)} \text{ with } \beta = \frac{1}{kT}

where :math:`k \approx 1.987 \cdot 10^{-3} \frac{kcal}{mol~K}` is the Boltzmann constant,
and :math:`T` the thermodynamic temperature. From that relation, the actual probability
of state :math:`s` can then be obtained using a proper scaling factor, the *canonical
partition function*

.. math::

  Z = \sum_{s \in \Omega} e^{-\beta E(s)}

where :math:`\Omega` is the finite set of all states. Finally, the equilibrium probability
of state :math:`s` can be computed as

.. math::

  p(s) = \frac{e^{-\beta E(s)}}{Z}

The first application of the above principles to an ensemble :math:`\Omega` of RNA
secondary structures :math:`s \in \Omega` has been described by :cite:t:`mccaskill:1990`.
Instead of enumerating all states exhaustively to compute :math:`Z` the basic idea is to
make use of the recursive decompositions described in :doc:`grammar`. This in turn allows
for an efficient computation of :math:`Z` in :math:`\mathcal{O}(n^3)`, i.e. cubic time.
An **outside** variant of the same recursions is then used to compute probabilities for
base pairs, stretches of consecutive unpaired nucleotides, or structural motifs.

Implementation
--------------

Our implementations provide access to the algorithm(s) made available by :cite:t:`mccaskill:1990`
and its variations. This allows us to efficiently compute

* the partition function, and from that
* various equilibrium probabilities, for instance base pair probabilities,
  probabilities of individual structure motifs, and many more.

Similar to our :doc:`/advanced/mfe`, we provide two different flavors for partition
function computations, **global** and **local** (sliding-window) algorithms, as well
as **simplified** and **generic** interfaces.

Comparative predictions for multiple sequence alignment (MSA) input compute an average
of the free energy contributions plus an additional covariance pseudo-energy term, exactly
as we do for the :doc:`/advanced/mfe`.

.. warning::

    Local (consensus) partition function computations for multiple sequence alignment (MSA) input
    is currently not supported!


Boltzmann weights for the free energy contributions of individual loops can be
found in :ref:`api-eval-loops`.

Our library also implements a stochastic (Boltzmann) sampling procedure to draw random
RNA secondary structures from the ensemble according to their equilibrium probabilty.

.. seealso::

    Further details of the Partition function and Base Pair Probability algorithm
    can be obtained from :cite:t:`mccaskill:1990`. The description for the basic
    partition function API can be found at :ref:`api-partfunc`.


Global Partition Function and Equilibrium Probabilities
-------------------------------------------------------

We provide implementations of the global partition function algorithm for

* Single sequences,
* Multiple sequence alignments (MSA), and
* RNA-RNA hybrids



.. seealso::

    A full list of available API symbols is available at :ref:`api-partfunc-global`.


Local Partition Function and Equilibrium Probabilities
------------------------------------------------------

Scanning version using a sliding window approach to compute equilibrium
probabilities are also available. They allow for efficient prediction of
local secondary structure features, e.g. for very long RNA sequences.

The following ``Python`` example demonstrates how to use the local partition
function API to compute **base pair probabilities**:

.. tab-set::

    .. tab-item:: Return value based

        .. code:: python

            import RNA

            # a rather long RNA sequence
            longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU"

            # get a list of base pair probabilities with window
            # size of 200nt, maximum base pair span of 150nt
            # and a probability threshold of 0.01
            for p in RNA.pfl_fold(longseq, 200, 150, 0.01):
                print(p)

    .. tab-item:: Callback function based

        .. code:: python

            import RNA

            # a rather long RNA sequence
            longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU"

            # callback function to store the probabilities
            def pf_window_callback(v, v_size, i, maxsize, what, data=None):
                if what & RNA.PROBS_WINDOW_BPP:
                    data.extend([{'i': i, 'j': j, 'p': p} for j, p in enumerate(v) if (p is not None) and (p >= 0.01)])

            # Prepare fold compound for local structure prediction
            # with window size of 200nt, maximum base pair span of
            # 150nt
            fc = RNA.fold_compound(longseq,
                                   RNA.md(max_bp_span = 150, window_size = 200),
                                   RNA.OPTION_WINDOW)
            # storage container for probabilities 
            data = []
            # run local structure prediction callback API
            fc.probs_window(10, RNA.PROBS_WINDOW_BPP, pf_window_callback, data)

            for p in data:
                print(p)

**Unpaired probabilities** can be predicted similarly:

.. tab-set::

    .. tab-item:: Return value based

        .. code:: python

            import RNA

            # a rather long RNA sequence
            longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU"

            # get a list of unpaired probabilities for subsequences
            # up to the size of 10, computed with window size of
            # 200nt and maximum base pair span of 150nt
            for up in RNA.pfl_fold_up(10, longseq, 200, 150):
                print(p)

    .. tab-item:: Callback function based

        .. code:: python

            import RNA

            # a rather long RNA sequence
            longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU"

            # callback function to store the unpaired probabilities
            def up_callback(v, v_size, i, maxsize, what, data):
                if what & RNA.PROBS_WINDOW_UP:
                    data.append({ 'i': i, 'up': v})

            # Prepare fold compound for local structure prediction
            # with window size of 200nt, maximum base pair span of
            # 150nt
            fc = RNA.fold_compound(longseq,
                                   RNA.md(max_bp_span = 150, window_size = 200),
                                   RNA.OPTION_WINDOW)
            # storage container for probabilities 
            data = []
            # run local structure prediction callback API
            fc.probs_window(10, RNA.PROBS_WINDOW_UP, up_callback, data)

            for up in data:
                print(up)

.. seealso::

    A full list of available API symbols is available at :ref:`api-partfunc-local`.


Predicting various Thermodynamic Properties
-------------------------------------------

As stated above, many thermodynamic properties can be derived from the :math:`Z`.
In particular, within our model of RNA secondary structures, the equilibrium
probabilty :math:`p(F)` of a particular structural feature :math:`F` follows
Boltzmanns law, i.e.:

.. math::

  p(F) \propto \sum_{s \mid F \in s} e^{\frac{-E(s)}{kT}}.

The actual probabilities can then be obtained from the ratio of those
structures containing :math:`F` and *all* structures, i.e.

.. math::

  p(F) = \frac{1}{Z} \sum_{s \mid F \in s} e^{\frac{-E(s)}{kT}}.


Efficient dynamic programming algorithms exist to compute the equilibrium
probabilities

.. math::

  p_{ij} = \frac{1}{Z} \sum_{s \mid (i,j) \in s} e^{\frac{-E(s)}{kT}}

of base pairs :math:`(i,j)` without the need for exhaustive enumeration
of :math:`s`.

.. seealso::

    A full list of available API symbols is available at :ref:`api-partfunc-probabilities`.
