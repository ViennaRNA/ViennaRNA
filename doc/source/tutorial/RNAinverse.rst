==========================
The Program ``RNAinverse``
==========================

.. contents:: Table of Contents
    :maxdepth: 1
    :local:


Introduction
============

``RNAinverse`` searches for sequences folding into a predefined structure,
thereby inverting the folding algorithm. Input consists of the target
structures (in dot-bracket notation) and a starting sequence, which is
optional.

Lower case characters in the start sequence indicate fixed positions,
i.e. they can be used to add sequence constraints. ``N``'s in the
starting sequence will be replaced by a random nucleotide.
For each search the best sequence found and its Hamming distance to the
start sequence are printed to *stdout*. If the the search was
unsuccessful a structure distance to the target is appended.

By default the program stops as soon as it finds a sequence that has the
target as MFE structure. The option ``-Fp`` switches ``RNAinverse`` to
the partition function mode where the probability of the target structure
:math:`\exp(-E(S)/RT) / Z` is maximized. This tends to produce sequences
with a more well-defined structure.

This probability is written in dot-brackets after the found sequence
and Hamming distance. With the option ``-R`` you can specify how often
the search should be repeated.


Sequence Design
===============

- Prepare an input file ``inv.in`` containing the target structure and
  sequence constraints::

    (((.(((....))).)))
    NNNgNNNNNNNNNNaNNN

- Design sequences using ``RNAinverse``::

    $ RNAinverse < inv.in
      GGUgUUGGAUCCGAaACC    5

  or design even more sequences with::

    $ RNAinverse -R5 -Fp < inv.in
      GGUgUGAACCCUCGaACC    5
      GGCgCCCUUUUGGGaGCC   12  (0.967418)
      CUCgAUCUCACGAUaGGG    6
      GGCgCCCGAAAGGGaGCC   13  (0.967548)
      GUUgAGCCCAUGCUaAGC    6
      GGCgCCCUUAUGGGaGCC   10  (0.967418)
      CGGgUGUUGUGACAaCCG    5
      GCGgGUCGAAAGGCaCGC   12  (0.925482)
      GCCgUAUCCGGGUGaGGC    6
      GGCgCCCUUUUGGGaGCC   13  (0.967418)


The output consists of the calculated sequence and the number of mutations
needed to get the MFE-structure from the start sequence (start sequence not
shown). Additionaly, with the partition function folding (``-Fp``) set, the
second output is another refinement so that the ensemble preferes the MFE
and folds into your given structure with a distinct probability, shown in
brackets.

Other RNA design tools
======================

Another useful program for inverse folding is ``RNA designer``, see
http://www.rnasoft.ca. ``RNA Designer`` takes a secondary structure description
as input and returns an RNA strand that is likely to fold in the given
secondary structure.

The ``sequence design application`` of the ``ViennaRNA Design Webservices``,
see http://nibiru.tbi.univie.ac.at/rnadesign/index.html, uses a different approach,
allowing for more than one secondary structure as input. For more detail read
the online Documentation and the next section of this tutorial.
