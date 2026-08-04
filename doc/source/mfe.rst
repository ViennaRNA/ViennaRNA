Minimum Free Energy (MFE) Algorithms
====================================

General Remarks
---------------

According to the laws of Thermodynamics, the most probable
conformation (here, secondary structure :math:`s`) is the
one with least free energy (:math:`E(s)`), since:

.. math::

  p(s) \propto e^{-\beta E(s)} \text{ with } \beta = \frac{1}{RT}

where :math:`R` is the gas constant and :math:`T` the thermodynamic
temperature. Hence, a reasonable RNA secondary structure representative
is the structure with minimum free energy (MFE). The basic idea behind
an efficient dynamic programming MFE prediction algorithm was first
described by :cite:t:`zuker:1981`. Our library implements various
variants of the original algorithmic idea that can be applied to
single sequences, multiple sequence alignments, multiple interacting
sequences and circular sequences.

All algorithms are accessible via a :ref:`simplified-mfe-api` that requires
minimal input and user interaction, but, at the same time is limited
to default setting. For a more fine grained selection of algorithmic
details, e.g. parameter settings or constraints, the :ref:`generic-mfe-api`
should be used. The latter provides functions that operate on
:c:type:`vrna_fold_compound_t` objects, which, in turn, can be fully
customized to the specific needs for the use-case.

The entire API for MFE prediction is basically subdivided into two
different flavors,

* :doc:`/mfe/global` - to compute the MFE for the entire sequence
* :doc:`/mfe/window` - to compute MFEs for each window using a sliding window approach


Each of these flavors, again, provides implementations to either
compute the MFE based on

*  single RNA (DNA) sequence(s), or
*  multiple sequences interacting with each other, or
*  a comparative approach using multiple sequence alignments (MSA).

For the latter, a consensus secondary structure is predicted and
our implementations compute an average of free energies for each
sequence in the MSA plus an additional covariance pseudo-energy
term.

The implementations for :doc:`/mfe/backtracking` are generally
agnostic with respect to whether local or global structure
prediction is in place.


.. _simplified-mfe-api:

Simplified Interface
--------------------

The most simple interface to these implementations might be the
:c:func:`vrna_fold` function. It only requires an input RNA sequence
in uppercase RNA alphabet (AUCG) to return the MFE in :math:`kcal \cdot mol^{-1}`
and a corresponding secondary structure in :ref:`dot-bracket-notation`.

.. admonition:: Example

  .. dropdown:: Single Sequence MFE prediction
      :open:

      .. tab-set::

        .. tab-item:: C

            .. literalinclude:: ../../examples/helloworld_mfe.c
                :encoding: latin-1
                :language: c
                :emphasize-lines: 18

        .. tab-item:: Python

            .. literalinclude:: ../../examples/Python/helloworld_mfe.py
               :encoding: latin-1
               :language: python
               :emphasize-lines: 7

        .. tab-item:: Perl

            .. literalinclude:: ../../examples/Perl5/helloworld_mfe.pl
               :encoding: latin-1
               :language: perl
               :emphasize-lines: 7

      The result of the above code should look like:

      .. code-block:: bash

          GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA
          ..(((((........)))))(((((((...)))))))..... [  -8.80 ]

  .. dropdown:: Consensus Structure Prediction for an MSA

      .. tab-set::

        .. tab-item:: C

            .. literalinclude:: ../../examples/helloworld_mfe_comparative.c
                :encoding: latin-1
                :language: c
                :emphasize-lines: 27

        .. tab-item:: Python

            .. literalinclude:: ../../examples/Python/helloworld_mfe_comparative.py
               :encoding: latin-1
               :language: python
               :emphasize-lines: 14

        .. tab-item:: Perl

            .. literalinclude:: ../../examples/Perl5/helloworld_mfe_comparative.pl
               :encoding: latin-1
               :language: perl
               :emphasize-lines: 14
    

.. note::

  Further simplified functions that compute the MFE for more special cases
  include:

  * :c:func:`vrna_circfold` - Treats the RNA sequence as circular sequence (:cite:t:`hofacker:2006`)
  * :c:func:`vrna_alifold` - Uses a Multiple Sequence Alignment as input to compute consensus MFE
  * :c:func:`vrna_circalifold` - Treats the aligned sequences as circular sequences


.. _generic-mfe-api:

Generic Interface
-----------------

The generic MFE interface requires one to create a :c:type:`vrna_fold_compound_t` object
first before calling the :c:func:`vrna_mfe` function (or method). Adaptations to the
algorithm parameters can be either introduced through the second or third
:c:type:`vrna_fold_compound_t` constructor parameter (see also the
:c:func:`vrna_fold_compound` function), or later via the various
:c:type:`vrna_fold_compound_t` associated methods (or functions in C).

.. admonition:: Example

  This is an example where we change some of the model parameters before
  actually computing the MFE (structure) under the non-default settings.


  .. tab-set::

    .. tab-item:: C

        .. literalinclude:: ../../examples/fold_compound_md.c
            :encoding: latin-1
            :language: c
            :emphasize-lines: 35

    .. tab-item:: Python

        .. literalinclude:: ../../examples/Python/helloworld_nondefault.py
           :encoding: latin-1
           :language: python
           :emphasize-lines: 18

    .. tab-item:: Perl

        .. literalinclude:: ../../examples/Perl5/helloworld_nondefault.pl
           :encoding: latin-1
           :language: perl
           :emphasize-lines: 18

  This should produce the following output:

  .. code-block:: bash
  
    GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA
    ..(((((........)))))(((((((...)))))))..... [  -5.64 ]



Full MFE API Documentation
--------------------------

* :doc:`mfe/global`
* :doc:`mfe/window`
* :doc:`mfe/backtracking`
* :doc:`mfe/global_deprecated`
* :doc:`mfe/window_deprecated`


.. toctree::
   :maxdepth: 1
   :caption: :
   :hidden:

   mfe/global
   mfe/window
   mfe/backtracking
   mfe/global_deprecated
   mfe/window_deprecated

