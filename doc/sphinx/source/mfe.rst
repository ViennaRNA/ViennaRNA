Minimum Free Energy (MFE) Algorithms
====================================

Computing the Minimum Free Energy (MFE), i.e. the most stable conformation
in thermodynamic equilibrium.

.. toctree::
   :maxdepth: 1
   :caption: Specialized Modules:

   mfe/global
   mfe/global_deprecated
   mfe/window
   mfe/window_deprecated
   mfe/backtracking


Zuker's Algorithm
-----------------

Our library provides fast dynamic programming Minimum Free Energy (MFE)
folding algorithms derived from the decomposition scheme as described by
:cite:t:`zuker:1981`.

MFE for circular RNAs
---------------------

Folding of *circular* RNA sequences is handled as a post-processing step
of the forward recursions. See :cite:t:`hofacker:2006` for further details.


MFE Algorithm API
-----------------

Predicting the Minimum Free Energy (MFE) and a corresponding (consensus)
secondary structure.

In a nutshell we provide two different flavors for MFE prediction:

* :doc:`/mfe/global` - to compute the MFE for the entire sequence
* :doc:`/mfe/window` - to compute MFEs for each window using a sliding window approach

Each of these flavors, again, provides two implementations to either compute the MFE based on

*  single RNA (DNA) sequence(s), or
*  multiple sequences interacting with each other, or
*  a comparative approach using multiple sequence alignments (MSA).

For the latter, a consensus secondary structure is predicted and our implementations compute
an average of free energies for each sequence in the MSA plus an additional covariance
pseudo-energy term.

The implementations for @ref mfe_backtracking are generally agnostic with respect to whether
local or global structure prediction is in place.
