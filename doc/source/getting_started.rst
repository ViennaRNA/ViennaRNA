Getting Started
===============

Here you find some more or less elaborate tutorials and manuals on how
to use our software.


Global RNA Secondary Structure Prediction
-----------------------------------------

Several tools for structure prediction of single RNA sequences are
available within the ``ViennaRNA Package``, each with its own special
subset of implemented algorithms.

.. toctree::
   :maxdepth: 1

   tutorial/RNAfold
   tutorial/RNApvmin
   tutorial/RNAsubopt

Consensus Structure Prediction
------------------------------

Sequence co-variations are a direct consequence of RNA base pairing
rules and can be deduced to alignments. RNA helices normally contain 
only 6 out of the 16 possible combinations: the Watson-Crick pairs
``GC``, ``CG``, ``AU``, ``UA``, and the somewhat weaker wobble pairs
``GU`` and ``UG``. Mutations in helical regions therefore have to be
correlated. In particular we often find *compensatory mutations*
where a mutation on one side of the helix is compensated by a second
mutation on the other side, e.g. a ``CG`` pair changes into a
``UA`` pair. Mutations where only one pairing partner changes (such
as ``CG`` to ``UG`` are termed *consistent mutations*.


.. toctree::
   :maxdepth: 1

   tutorial/RNAalifold


RNA-RNA interaction
-------------------

A common problem is the prediction of binding sites between two RNAs, as in
the case of miRNA-mRNA interactions. Following tools of the ``ViennaRNA Package``
can be used to calculate base pairing probabilities.

.. toctree::
   :maxdepth: 1

   tutorial/RNAcofold
   tutorial/RNAduplex


Plotting Structures
-------------------

.. toctree::
   :maxdepth: 1

   tutorial/RNAplot


RNA Design
----------

.. toctree::
   :maxdepth: 1

   tutorial/RNAinverse
   tutorial/switch


RNA folding kinetics
--------------------

RNA folding kinetics describes the dynamical process of how a RNA molecule
approaches to its unique folded biological active conformation (often
referred to as the native state) starting from an initial ensemble of
disordered conformations e.g. the unfolded open chain. The key for
resolving the dynamical behavior of a folding RNA chain lies in the
understanding of the ways in which the molecule explores its astronomically
large free energy landscape, a rugged and complex hyper-surface established
by all the feasible base pairing patterns a RNA sequence can form. The
challenge is to understand how the interplay of formation and break up of
base pairing interactions along the RNA chain can lead to an efficient
search in the energy landscape which reaches the native state of the
molecule on a biologically meaningful time scale.

.. toctree::
   :maxdepth: 1

   tutorial/RNA2Dfold
   tutorial/barriers_treekin


Other Utilities
---------------

.. toctree::
   :maxdepth: 1

   utilities

