Local (sliding window) MFE Prediction
=====================================

Variations of the local (sliding window) Minimum Free Energy (MFE) prediction algorithm.

We provide implementations for the local (sliding window) MFE prediction algorithm for

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

.. tip::

  .. dropdown:: Simple local MFE prediction
      :open:

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

  .. dropdown:: Callback based local MFE prediction

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



Local MFE API Documentation
---------------------------

.. doxygengroup:: mfe_window
    :no-title:
