Minimum Free Energy (MFE)
=========================

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
variants of the original algorithmic idea that can be applied to:

* single sequences,
* multiple sequence alignments (MSA),
* multiple interacting sequences (RNA-RNA interaction), and
* circular sequences

All algorithms are usually accessible in two forms:

* a :ref:`simplified-mfe-api` that requires minimal input and user
  interaction, but, at the same time is limited to default setting, and
* a :ref:`generic-mfe-api` that allows for a fine grained selection of
  algorithmic details, e.g. parameter settings or constraints
  
The latter provides functions that operate on :c:type:`vrna_fold_compound_t`
objects, which, in turn, can be fully customized to the specific needs
for the use-case.

.. seealso::

    The full list of available API symbols for global MFE prediction can be found in
    :ref:`api-mfe-global`.


.. _simplified-mfe-api:

Simplified MFE Interface
------------------------

The most simple interface to these implementations might be the
:c:func:`vrna_fold` function. It only requires an input RNA sequence
in uppercase RNA alphabet (AUCG) to return the MFE in :math:`kcal \cdot mol^{-1}`
and a corresponding secondary structure in :ref:`dot-bracket-notation`.

Single Sequence MFE
^^^^^^^^^^^^^^^^^^^

Below, you'll find an example that illustrates how to predict the MFE
and a corresponding secondary structure for a single (linear) RNA sequence:

.. tab-set::

    .. tab-item:: C

        .. literalinclude:: ../../../examples/helloworld_mfe.c
            :encoding: latin-1
            :language: c
            :emphasize-lines: 18

    .. tab-item:: Python

        .. literalinclude:: ../../../examples/Python/helloworld_mfe.py
            :encoding: latin-1
            :language: python
            :emphasize-lines: 7

    .. tab-item:: Perl

        .. literalinclude:: ../../../examples/Perl5/helloworld_mfe.pl
            :encoding: latin-1
            :language: perl
            :emphasize-lines: 7


The result of the above code should look like:

.. code-block:: bash

    GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA
    ..(((((........)))))(((((((...)))))))..... [  -8.80 ]


Consensus MFE
^^^^^^^^^^^^^

In cases where you want to predict the consensus MFE structure for a set
of aligned RNA sequences (MSA), consider the following example:

.. tab-set::

    .. tab-item:: C

        .. literalinclude:: ../../../examples/helloworld_mfe_comparative.c
            :encoding: latin-1
            :language: c
            :emphasize-lines: 27

    .. tab-item:: Python

        .. literalinclude:: ../../../examples/Python/helloworld_mfe_comparative.py
           :encoding: latin-1
           :language: python
           :emphasize-lines: 14

    .. tab-item:: Perl

        .. literalinclude:: ../../../examples/Perl5/helloworld_mfe_comparative.pl
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

The following example demonstrates how to change some of the model parameters
before actually computing the MFE (structure) under the non-default settings.

.. tab-set::

    .. tab-item:: C

        .. literalinclude:: ../../../examples/fold_compound_md.c
            :encoding: latin-1
            :language: c
            :emphasize-lines: 35

    .. tab-item:: Python

        .. literalinclude:: ../../../examples/Python/helloworld_nondefault.py
            :encoding: latin-1
            :language: python
            :emphasize-lines: 18

    .. tab-item:: Perl

        .. literalinclude:: ../../../examples/Perl5/helloworld_nondefault.pl
            :encoding: latin-1
            :language: perl
            :emphasize-lines: 18


This should produce the following output:

.. code-block:: bash

    GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA
    ..(((((........)))))(((((((...)))))))..... [  -5.64 ]



Backtracking
------------

The implementations for :doc:`mfe/backtracking` are generally
agnostic with respect to whether local or global structure
prediction is in place.

.. seealso::

    The full list of available API symbols for backtracking of MFE structures
    can be found in :ref:`api-mfe-backtracking`.


Local MFE Prediction
--------------------

In addition to the global MFE prediction algorithm(s) above, we also provide implementations
for the local (sliding window) MFE prediction, in particular for

* Single sequences,
* Multiple sequence alignments (MSA), and

Note, that our implementation scans an RNA sequence (or MSA) from the 3' to the 5'
end, and reports back locally optimal (consensus) structures, the corresponding free
energy, and the position of the sliding window in global coordinates.

For any particular RNA sequence (or MSA) multiple locally optimal (consensus)
secondary structures may be predicted. Thus, we tried to implement an interface that
allows for an effortless conversion of the corresponding hits into any target data
structure. As a consequence, we provide two distinct ways to retrieve the corresponding
predictions, either

* through directly writing to an open ``FILE`` stream on-the-fly, or
* through a callback function mechanism.

The latter allows one to store the results in any possible target data structure. Our
implementations then pass the results through the user-implemented callback as soon as
the prediction for a particular window is finished.

Return Valued Based Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is some example code that shows how to use the simple local MFE prediction API:

.. tab-set::

    .. tab-item:: C

        .. code-block:: c

            #include <stdlib.h>
            #include <stdio.h>
            #include <string.h>

            #include <ViennaRNA/mfe/local.h>
            #include <ViennaRNA/utils/basic.h>

            int
            main()
            {
              /* a long RNA sequence input */
              char  *seq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCA" \
                           "GGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUC" \
                           "AUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUC" \
                           "AGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGC" \
                           "CAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCC" \
                           "AGGUAAUUCAAAGGCCUCAAGUCCU";

              /* allocate memory for MFE structure (length + 1) */
              char  *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));

              /* predict local MFE structures with maximum base pair distance of 50 nt */
              float e = vrna_Lfold(seq, 50, NULL);

              /* cleanup memory */
              free(structure);

              return 0;
            }

    .. tab-item:: Python

        .. code-block:: python

            import RNA

            # a long RNA sequence input
            longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCA" \
                      "GGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUC" \
                      "AUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUC" \
                      "AGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGC" \
                      "CAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCC" \
                      "AGGUAAUUCAAAGGCCUCAAGUCCU"

            # predict local MFE structures with maximum base pair distance of 50 nt
            e = RNA.Lfold(longseq, 50)


    .. tab-item:: Perl

        .. code-block:: perl

            use RNA;

            # a long RNA sequence input
            my $longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCA" .
                          "GGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUC" .
                          "AUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUC" .
                          "AGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGC" .
                          "CAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCC" .
                          "AGGUAAUUCAAAGGCCUCAAGUCCU";

            # predict local MFE structures with maximum base pair distance of 50 nt
            my $e = RNA::Lfold($longseq, 50);

In the example above, we didn't provide an output file stream as third argument
for :c:func:`vrna_Lfold`. Hence, the output was automatically written to
``stdout``.

Callback Based Example
^^^^^^^^^^^^^^^^^^^^^^

In contrast to the previous example, here is a variant that uses the callback
based local MFE prediction API:

.. tab-set::

    .. tab-item:: C

        .. code-block:: c
            
            #include <stdlib.h>
            #include <stdio.h>
            #include <string.h>

            #include <ViennaRNA/mfe/local.h>
            #include <ViennaRNA/utils/basic.h>

            /* a callback function that simply prints the resulting
             * local MFE prediction
             */
            void
            callback(unsigned int start,
                     unsigned int end,
                     const char *structure,
                     float en,
                     void *data)
            {
              if (structure) {
                printf("%.*s\n%s (%6.2f) - [%u:%u]\n",
                       end - start + 1,
                       (char *)data) + start - 1,
                       structure,
                       en,
                       start,
                       end);
              }
            }
            
            int
            main()
            {
              /* a long RNA sequence input */
              char  *seq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCA" \
                           "GGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUC" \
                           "AUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUC" \
                           "AGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGC" \
                           "CAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCC" \
                           "AGGUAAUUCAAAGGCCUCAAGUCCU";

              /* predict local MFE structures with maximum base pair distance of 50 nt */
              float e = vrna_Lfold_cb(seq, 50, callback, (void *)seq);

              return 0;
            }

    .. tab-item:: Python

        .. code-block:: python

            import RNA

            # a long RNA sequence input
            longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCA" \
                      "GGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUC" \
                      "AUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUC" \
                      "AGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGC" \
                      "CAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCC" \
                      "AGGUAAUUCAAAGGCCUCAAGUCCU"

            def callback(start, end, structure, en, data = None):
                if structure:
                    print(f"{data[start - 1:end]}\n{structure} ({en:6.2f}) - [{start}-{end}]")

            # predict local MFE structures with maximum base pair distance of 50 nt
            e = RNA.Lfold_cb(longseq, 50, callback, longseq)


    .. tab-item:: Perl

        .. code-block:: perl

            use RNA;

            # a long RNA sequence input
            my $longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCA" .
                          "GGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUC" .
                          "AUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUC" .
                          "AGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGC" .
                          "CAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCC" .
                          "AGGUAAUUCAAAGGCCUCAAGUCCU";

            sub callback {
              my ($start, $end, $structure, $en, $data) = @_;
              if ($structure) {
                printf(substr($data, start - 1, end - start + 1) .
                       $structure .
                       "(%6.2f) - [%u:%u]",
                       en,
                       start,
                       end);
              }
            }

            # predict local MFE structures with maximum base pair distance of 50 nt
            my $e = RNA::Lfold_cb($longseq, 50, \&callback, \&longseq);

.. seealso::

    The full list of available API symbols for local (sliding-window) MFE prediction can
    be found in :ref:`api-mfe-local`.

