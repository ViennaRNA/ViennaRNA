Predicting various Thermodynamic Properties
===========================================

Compute various thermodynamic properties using the partition function.

Many thermodynamic properties can be derived from the partition function

.. math::

  Z = \sum_{s \in \omega} e^{\frac{-E(s)}{kT}}.

In particular, for nucleic acids in equilibrium the probabilty :math:`p(F)`
of a particular structural feature :math:`F` follows Boltzmanns law, i.e.:

.. math::

  p(F) \propto \sum_{s \mid F \in s} e^{\frac{-E(s)}{kT}}.

The actual probabilities can then be obtained from the ratio of those
structures containing :math:`F` and *all* structures, i.e.

.. math::

  p(F) = \frac{1}{Z} \sum_{s \mid F \in s} e^{\frac{-E(s)}{kT}}.


Consequently, a particular secondary structure :math:`s` has equilibrium
probability

.. math::

  p(s) = \frac{1}{Z} e^{\frac{-E(s)}{kT}}

which can be easily computed once :math:`Z` and :math:`E(s)` are known.

Efficient dynamic programming algorithms exist to compute the equilibrium
probabilities

.. math::

  p_{ij} = \frac{1}{Z} \sum_{s \mid (i,j) \in s} e^{\frac{-E(s)}{kT}}

of base pairs :math:`(i,j)` without the need for exhaustive enumeration
of :math:`s`.

This interface provides the functions for all thermodynamic property
computations implemented in *RNAlib*.

Thermodynamic Properties API
----------------------------

.. doxygengroup:: thermodynamics
