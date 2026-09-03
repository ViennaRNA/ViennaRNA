.. _programs:

==================
Command-Line Tools
==================

The ViennaRNA Package comes with a number of executable
programs that provide command line interfaces to the most
important algorithms implemented in *RNAlib*.

Find an overview of these programs and their corresponding
manual pages below.

Main Programs
-------------

.. only:: latex

    ======================== =====================
    Program                  Description
    ======================== =====================
    :doc:`man/RNA2Dfold`     Compute MFE structure, partition function and
                             representative sample structures of k,l neighborhoods
    :doc:`man/RNAaliduplex`  Predict conserved RNA-RNA interactions between two alignments
    :doc:`man/RNAalifold`    Calculate secondary structures for a set of aligned RNA sequences
    :doc:`man/RNAcofold`     Calculate secondary structures of two RNAs with dimerization
    :doc:`man/RNAdistance`   Calculate distances between RNA secondary structures
    :doc:`man/RNAdos`        Calculate the density of states for each energy band of an RNA
    :doc:`man/RNAduplex`     Compute the structure upon hybridization of two RNA strands
    :doc:`man/RNAeval`       Evaluate free energy of RNA sequences with given secondary structure
    :doc:`man/RNAfold`       Calculate minimum free energy secondary structures and partition function of RNAs
    :doc:`man/RNAheat`       Calculate the specific heat (melting curve) of an RNA sequence
    :doc:`man/RNAinverse`    Find RNA sequences with given secondary structure (sequence design)
    :doc:`man/RNALalifold`   Calculate locally stable secondary structures for a set of aligned RNAs
    :doc:`man/RNALfold`      Calculate locally stable secondary structures of long RNAs
    :doc:`man/RNAmultifold`  Compute thermodynamic properties for interaction complexes of multiple RNAs
    :doc:`man/RNApaln`       RNA alignment based on sequence base pairing propensities
    :doc:`man/RNApdist`      Calculate distances between thermodynamic RNA secondary structures ensembles
    :doc:`man/RNAparconv`    Convert energy parameter files from ViennaRNA 1.8 to 2.0 format
    :doc:`man/RNAPKplex`     Predict RNA secondary structures including pseudoknots
    :doc:`man/RNAplex`       Find targets of a query RNA
    :doc:`man/RNAplfold`     Calculate average pair probabilities for locally stable secondary structures
    :doc:`man/RNAplot`       Draw RNA Secondary Structures in PostScript, SVG, or GML
    :doc:`man/RNApvmin`      Calculate a perturbation vector that minimizes discrepancies between predicted and observed pairing probabilities
    :doc:`man/RNAsnoop`      Find targets of a query H/ACA snoRNA
    :doc:`man/RNAsubopt`     Calculate suboptimal secondary structures of RNAs
    :doc:`man/RNAup`         Calculate the thermodynamics of RNA-RNA interactions
    ======================== =====================

.. only:: html

    .. grid:: 1 2 3 3
        :gutter: 1 1 2 3

        .. grid-item-card::

            **RNA2Dfold**
            ^^^
            Compute MFE structure, partition function and
            representative sample structures of k,l neighborhoods
            +++
            .. button-ref:: tutorial/RNA2Dfold
                :ref-type: doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNA2Dfold Tutorial

                RNA2Dfold Tutorial

            .. button-ref:: man/RNA2Dfold
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNA2Dfold Manpage

                Manpage

        .. grid-item-card::

            **RNAaliduplex**
            ^^^
            Predict conserved RNA-RNA interactions between two alignments
            +++
            .. button-ref:: man/RNAaliduplex
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAaliduplex Manpage

                Manpage

        .. grid-item-card::

            **RNAalifold**
            ^^^
            Calculate secondary structures for a set of aligned RNA sequences
            +++
            .. button-ref:: tutorial/RNAalifold
                :ref-type: doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAalifold Tutorial

                RNAalifold Tutorial

            .. button-ref:: man/RNAalifold
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAalifold Manpage

                Manpage

        .. grid-item-card::

            **RNAcofold**
            ^^^
            Calculate secondary structures of two RNAs with dimerization
            +++
            .. button-ref:: tutorial/RNAcofold
                :ref-type: doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAcofold Tutorial

                RNAcofold Tutorial

            .. button-ref:: man/RNAcofold
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAcofold Manpage

                Manpage

        .. grid-item-card::

            **RNAdistance**
            ^^^
            Calculate distances between RNA secondary structures
            +++
            .. button-ref:: man/RNAdistance
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAdistance Manpage

                Manpage

        .. grid-item-card::

            **RNAdos**
            ^^^
            Calculate the density of states for each energy band of an RNA
            +++
            .. button-ref:: man/RNAdos
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAdos Manpage

                Manpage

        .. grid-item-card::

            **RNAduplex**
            ^^^
            Compute the structure upon hybridization of two RNA strands
            +++
            .. button-ref:: tutorial/RNAduplex
                :ref-type: doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAduplex Tutorial

                RNAduplex Tutorial

            .. button-ref:: man/RNAduplex
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAduplex Manpage

                Manpage

        .. grid-item-card::

            **RNAeval**
            ^^^
            Evaluate free energy of RNA sequences with given secondary structure
            +++
            .. button-ref:: man/RNAeval
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAeval Manpage

                Manpage

        .. grid-item-card::
            :text-align: justify

            **RNAfold**
            ^^^
            Compute the minimum free energy structure
            and base pair probabilities for an RNA sequence.
            +++
            .. button-ref:: tutorial/RNAfold
                :ref-type: doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAfold Tutorial

                RNAfold Tutorial

            .. button-ref:: man/RNAfold
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAfold Manpage

                Manpage

        .. grid-item-card::

            **RNAheat**
            ^^^
            Calculate the specific heat (melting curve) of an RNA sequence
            +++
            .. button-ref:: man/RNAheat
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAheat Manpage

                Manpage

        .. grid-item-card::

            **RNAinverse**
            ^^^
            Find RNA sequences with given secondary structure (sequence design)
            +++
            .. button-ref:: tutorial/RNAinverse
                :ref-type: doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAinverse Tutorial

                RNAinverse Tutorial

            .. button-ref:: man/RNAinverse
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAinverse Manpage

                Manpage

        .. grid-item-card::

            **RNALalifold**
            ^^^
            Calculate locally stable secondary structures for a set of aligned RNAs
            +++
            .. button-ref:: man/RNALalifold
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNALalifold Manpage

                Manpage

        .. grid-item-card::

            **RNALfold**
            ^^^
            Calculate locally stable secondary structures of long RNAs
            +++
            .. button-ref:: man/RNALfold
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNALfold Manpage

                Manpage

        .. grid-item-card::

            **RNAmultifold**
            ^^^
            Compute thermodynamic properties for interaction complexes of multiple RNAs
            +++
            .. button-ref:: man/RNAmultifold
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAmultifold Manpage

                Manpage

        .. grid-item-card::

            **RNApaln**
            ^^^
            RNA alignment based on sequence base pairing propensities
            +++
            .. button-ref:: man/RNApaln
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNApaln Manpage

                Manpage

        .. grid-item-card::

            **RNApdist**
            ^^^
            Calculate distances between thermodynamic RNA secondary structures ensembles
            +++
            .. button-ref:: man/RNApdist
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNApdist Manpage

                Manpage

        .. grid-item-card::

            **RNAparconv**
            ^^^
            Convert energy parameter files from ViennaRNA 1.8 to 2.0 format
            +++
            .. button-ref:: man/RNAparconv
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAparconv Manpage

                Manpage

        .. grid-item-card::

            **RNAPKplex**
            ^^^
            Predict RNA secondary structures including pseudoknots
            +++
            .. button-ref:: man/RNAPKplex
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAPKplex Manpage

                Manpage

        .. grid-item-card::

            **RNAplex**
            ^^^
            Find targets of a query RNA
            +++
            .. button-ref:: man/RNAplex
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAplex Manpage

                Manpage

        .. grid-item-card::

            **RNAplfold**
            ^^^
            Calculate average pair probabilities for locally stable secondary structures
            +++
            .. button-ref:: man/RNAplfold
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAplfold Manpage

                Manpage

        .. grid-item-card::

            **RNAplot**
            ^^^
            Draw RNA Secondary Structures in PostScript, SVG, or GML
            +++
            .. button-ref:: tutorial/RNAplot
                :ref-type: doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAplot Tutorial

                RNAplot Tutorial

            .. button-ref:: man/RNAplot
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAplot Manpage

                Manpage

        .. grid-item-card::

            **RNApvmin**
            ^^^
            Calculate a perturbation vector that minimizes discrepancies between predicted and observed pairing probabilities
            +++
            .. button-ref:: tutorial/RNApvmin
                :ref-type: doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNApvmin Tutorial

                RNApvmin Tutorial

            .. button-ref:: man/RNApvmin
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNApvmin Manpage

                Manpage

        .. grid-item-card::

            **RNAsnoop**
            ^^^
            Find targets of a query H/ACA snoRNA
            +++
            .. button-ref:: man/RNAsnoop
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAsnoop Manpage

                Manpage

        .. grid-item-card::

            **RNAsubopt**
            ^^^
            Calculate suboptimal secondary structures of RNAs
            +++
            .. button-ref:: tutorial/RNAsubopt
                :ref-type: doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAsubopt Tutorial

                RNAsubopt Tutorial

            .. button-ref:: man/RNAsubopt
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAsubopt Manpage

                Manpage

        .. grid-item-card::

            **RNAup**
            ^^^
            Calculate the thermodynamics of RNA-RNA interactions
            +++
            .. button-ref:: man/RNAup
                :ref-type: doc
                :color: muted
                :align: center
                :expand:
                :outline:
                :tooltip: To the RNAup Manpage

                Manpage


Additional Programs
-------------------

We include the following additional programs in our distribution of the
ViennaRNA Package. Whether or not they are installed together with the
ViennaRNA Package depends on its :doc:`/configuration`.

.. only:: latex

    ========================= =====================
    Program                   Description
    ========================= =====================
    AnalyseDists              Analyse a distance matrix
    AnalyseSeqs               Analyse a set of sequences of common length
    Kinfold                   Simulate kinetic folding of RNA secondary
                              structures
    kinwalker                 Predict RNA folding trajectories
    RNAforester               Compare RNA secondary structures via forest alignment
    RNAlocmin                 Calculate local minima from structures via gradient walks
    RNAxplorer                Explore the RNA conformation space through
                              sampling and other techniques
    ========================= =====================

.. only:: html

    .. grid:: 1 2 3 4
        :gutter: 1 1 2 3

        .. grid-item-card:: AnalyseDists

            Analyse a distance matrix

        .. grid-item-card:: AnalyseSeqs

            Analyse a set of sequences of common length

        .. grid-item-card:: Kinfold

            Simulate kinetic folding of RNA secondary structures

        .. grid-item-card:: kinwalker

            Predict RNA folding trajectories

        .. grid-item-card:: RNAforester

            Compare RNA secondary structures via forest alignment

        .. grid-item-card:: RNAlocmin

            Calculate local minima from structures via gradient walks

        .. grid-item-card:: RNAxplorer

            Explore the RNA conformation space through sampling and other techniques


Utilities
---------

We also ship a number of small utilities, many of them to manipulate the
PostScript files produced by the structure prediction programs :doc:`man/RNAfold`
and :doc:`man/RNAalifold`.

.. tip::

    Most of the `Perl 5` utilities contain embedded `pod` documentation. Type e.g.

    .. code:: bash

        perldoc relplot.pl

    for detailed instructions.


==============  =======================
Tool Name       Description
==============  =======================
ct2db           | Produce dot bracket notation of an RNA secondary structure
                | given as mfold ``.ct`` file
b2mt.pl         | Produce a mountain representation of a secondary structure
                | from it's dot-bracket notation, as produced by :doc:`man/RNAfold`.
                | Output consists of simple ``x y`` data suitable as input to a
                | plotting program. The mountain representation is a xy plot with
                | sequence position on the x-axis and the number of base pairs
                | enclosing that position on the y-axis.

                Example:

                .. code:: bash

                    RNAfold < my.seq | b2mt.pl | xmgrace -pipe

                .. image:: ../gfx/5smnt.png
                    :alt: mountain plot for 5S RNA
                    :align: center

cmount.pl       | Produce a PostScript mountain plot from a color dot plot as
                | created by :doc:`man/RNAalifold` -p or ``alidot -p``. Each base pair
                | is represented by a trapez whose color encodes the number of
                | consistent and compensatory mutations supporting that pair:
                | Red marks pairs with no sequence variation; ochre, green, turquoise,
                | blue, and violet mark pairs with 2,3,4,5,6 different types of pairs,
                | respectively.

                Example:

                .. code:: bash
                
                    cmount.pl alidot.ps > cmount.ps

                .. image:: ../gfx/TAR_ss.png
                    :alt: color anotated TAR hairpin
                    :align: center

coloraln.pl     | Reads a sequence alignment in ``CLUSTAL`` format and a consensus
                | secondary structure (which it extracts from a secondary structure
                | plot as produced by :doc:`man/RNAalifold`), and produces a postscript figure
                | of the alignment annotated using the consensus structure, coloring
                | base pair using the same color scheme as ``cmount.pl``, :doc:`man/RNAalifold`
                | and ``alidot``.

                Example:

                .. code:: bash

                    coloraln.pl -s alirna.ps file.aln > coloraln.ps

                .. image:: ../gfx/TAR_aln.png
                    :alt: colored TAR alignment
                    :align: center

colorrna.pl     | Reads a consensus secondary structure plot and a color dot plot
                | as produced by :doc:`man/RNAalifold` -p, and writes a new secondary
                | structure plot in which base pairs a colored using the color
                | information from the dot plot.

                Example:

                .. code:: bash
                
                    colorrna.pl alirna.ps alidot.ps > colorRNA.ps

mountain.pl     | Similar to ``b2mt.pl``, but produces a mountain plot from the
                | pair probabilities contained in a PostScript dot plot. It write
                | 3 sets of x y data, suitable as input for a plot program. The
                | first two sets containing the mountain representation from pair
                | probabilities and MFE structure, the third set is the
                | "positional entropy" a measure of structural well-definedness.

                Example:

                .. code:: bash

                    mountain.pl dot.ps | xmgrace -pipe

                .. image:: ../gfx/TAR_cmt.png
                    :alt: color mountain plot for HIV TAR hairpin
                    :align: center

refold.pl       Refold using consensus structure as constraint
relplot.pl      | Reads a postscript secondary structure plot and a dot plot
                | containing pair probabilities as produced by :doc:`man/RNAfold` -p,
                | and produces a new structure plot, color annotated with reliability
                | information in the form of either pair probabilities or positional
                | entropy (default).

                Example:

                .. code:: bash

                    relplot.pl foo_ss.ps foo_dp.ps > foo_rss.ps

                .. image:: ../gfx/5srel.png
                    :alt: structure of a 5s rRNA with reliability annotation
                    :align: center

rotate_ss.pl    | Reads a postscript secondary structure plot as produced by
                | :doc:`man/RNAfold` and produces a new rotated and/or mirrored
                | structure plot.

                Example:

                .. code:: bash

                    rotate_ss.pl -a 30 -m foo_ss.ps > foo_new_ss.ps

switch.pl       | Design sequences that can adopt two different structure, i.e.
                | design RNA switches. The program will sample the set of sequences
                | compatible with two input structures in order to find sequences with
                | desired thermodynamic and kinetic properties. In particular it is
                | possible to specify two different temperatures such that structure
                | 1 is favored at T1 and structure 2 at T2 to design temperature
                | sensitive switches (RNA thermometers). The desired height of the
                | energy barrier separating the two structures, thus determining
                | the refolding time between meta-stable states.
RNAdesign.pl    | Flexible design of multi-stable RNA molecules. An initially random
                | sequence is iteratively mutated and evaluated according to an
                | objective function (see Option: ``--optfun``). Whenever a better
                | scoring sequence has been found, the mutation is accepted, the
                | algorithm terminates once a local minimum is found. This script
                | makes heavy use of the ``RNA::Design`` sub-package that comes
                | with the ViennaRNA Package ``Perl 5`` interface.
==============  =======================
