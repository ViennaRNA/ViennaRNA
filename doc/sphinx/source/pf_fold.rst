Partition Function and Equilibrium Properties
=============================================

In contrast to methods that compute the property of a single structure in the ensemble,
e.g. :doc:`/mfe`, the partition function algorithms always consider the *entire*
equilibrium ensemble. For that purpose, the algorithm(s) made available by
:cite:t:`mccaskill:1990` and its variants can be used to efficiently compute

- the partition function, and from that
- various equilibrium probabilities, for instance base pair probabilities,
  probabilities of individual structure motifs, and many more.

The principal idea behind this approach is that in equilibrium, statistical mechanics
and polymer theory tells us that the frequency or probability :math:`p(s)` of a particular
state :math:`s` depends on its energy :math:`E(s)` and follows a Boltzmann distribution, i.e.

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

Instead of enumerating all states exhaustively to compute :math:`Z` one can apply the
:ref:`grammar:secondary structure folding recurrences` again for an efficient computation
in cubic time. An *outside* variant of the same recursions is then used to compute
probabilities for base pairs, stretches of consecutive unpaired nucleotides, or structural
motifs.

.. admonition:: See also...

  Further details of the Partition function and Base Pair Probability algorithm
  can be obtained from :cite:t:`mccaskill:1990`


.. toctree::
   :maxdepth: 1
   :caption: Specialized Modules:

   partfunc/global
   partfunc/window
   partfunc/thermodynamics
   partfunc/global_deprecated
   partfunc/window_deprecated


Partition Function API
----------------------

Similar to our :doc:`/mfe`, we provide two different flavors for partition
function computations:

* :doc:`/partfunc/global` - to compute the partition function for a full length sequence
* :doc:`/partfunc/window` - to compute the partition function of each window using a sliding window approach

While the global partition function approach supports predictions using single sequences as
well as consensus partition functions for multiple sequence alignments (MSA), we currently
do not support MSA input for the local variant.

Comparative prediction computes an average of the free energy contributions plus an additional
covariance pseudo-energy term, exactly as we do for the @ref mfe implementation.

Boltzmann weights for the free energy contributions of individual loops can be
found in :doc:`/eval/eval_loops`.

Our implementations also provide a stochastic backtracking procedure to draw
@ref subopt_stochbt according to their equilibrium probabilty.

General Partition Function API
------------------------------

.. doxygengroup:: pf_fold
