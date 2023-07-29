Utilities
=========

We also ship a number of small utilities, many of them to manipulate the
PostScript files produced by the structure prediction programs :doc:`/man/RNAfold`
and :doc:`/man/RNAalifold`.

Most of the `Perl 5` utilities contain embedded `pod` documentation. Type e.g.

.. code:: bash

    perldoc relplot.pl

for detailed instructions.

Available Tools
---------------

==============  =======================
Tool Name       Description
==============  =======================
ct2db           | Produce dot bracket notation of an RNA secondary structure
                | given as mfold ``.ct`` file
b2mt.pl         | Produce a mountain representation of a secondary structure
                | from it's dot-bracket notation, as produced by :doc:`/man/RNAfold`.
                | Output consists of simple ``x y`` data suitable as input to a
                | plotting program. The mountain representation is a xy plot with
                | sequence position on the x-axis and the number of base pairs
                | enclosing that position on the y-axis.

                Example:

                .. code:: bash

                    RNAfold < my.seq | b2mt.pl | xmgrace -pipe

                .. image:: gfx/5smnt.png
                    :alt: mountain plot for 5S RNA
                    :align: center

cmount.pl       | Produce a PostScript mountain plot from a color dot plot as
                | created by :doc:`/man/RNAalifold` -p or ``alidot -p``. Each base pair
                | is represented by a trapez whose color encodes the number of
                | consistent and compensatory mutations supporting that pair:
                | Red marks pairs with no sequence variation; ochre, green, turquoise,
                | blue, and violet mark pairs with 2,3,4,5,6 different types of pairs,
                | respectively.

                Example:

                .. code:: bash
                
                    cmount.pl alidot.ps > cmount.ps

                .. image:: gfx/TAR_ss.png
                    :alt: color anotated TAR hairpin
                    :align: center

coloraln.pl     | Reads a sequence alignment in ``CLUSTAL`` format and a consensus
                | secondary structure (which it extracts from a secondary structure
                | plot as produced by :doc:`/man/RNAalifold`), and produces a postscript figure
                | of the alignment annotated using the consensus structure, coloring
                | base pair using the same color scheme as ``cmount.pl``, :doc:`/man/RNAalifold`
                | and ``alidot``.

                Example:

                .. code:: bash

                    coloraln.pl -s alirna.ps file.aln > coloraln.ps

                .. image:: gfx/TAR_aln.png
                    :alt: colored TAR alignment
                    :align: center

colorrna.pl     | Reads a consensus secondary structure plot and a color dot plot
                | as produced by :doc:`/man/RNAalifold` -p, and writes a new secondary
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

                .. image:: gfx/TAR_cmt.png
                    :alt: color mountain plot for HIV TAR hairpin
                    :align: center

refold.pl       Refold using consensus structure as constraint
relplot.pl      | Reads a postscript secondary structure plot and a dot plot
                | containing pair probabilities as produced by :doc:`/man/RNAfold` -p,
                | and produces a new structure plot, color annotated with reliability
                | information in the form of either pair probabilities or positional
                | entropy (default).

                Example:

                .. code:: bash

                    relplot.pl foo_ss.ps foo_dp.ps > foo_rss.ps

                .. image:: gfx/5srel.png
                    :alt: structure of a 5s rRNA with reliability annotation
                    :align: center

rotate_ss.pl    | Reads a postscript secondary structure plot as produced by
                | :doc:`/man/RNAfold` and produces a new rotated and/or mirrored
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
