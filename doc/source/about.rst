About the ViennaRNA Package
===========================

The core of the **ViennaRNA Package** (:cite:t:`lorenz:2011,hofacker:1994`)
is formed by a collection of algorithms for the **prediction and comparison**
of RNA secondary structures.

**Physics based RNA secondary structure prediction** is the main functionality
in the package. We provide variants of three different kinds of dynamic
programming algorithms for structure prediction, (i) the **minimum free energy**
algorithm of :cite:t:`zuker:1981` which yields a single optimal structure,
(ii) the **partition function** algorithm of :cite:t:`mccaskill:1990` that allows
for **base pair probability** computation or stochastic backtracking of structure
samples, and (iii) the **suboptimal folding** algorithm of :cite:t:`wuchty:1999`
which generates all suboptimal structures within a given energy range of the
optimal energy. For secondary structure comparison, the package contains
several measures of distance (dissimilarities) using either **string alignment**
or **tree-editing** :cite:t:`shapiro:1990`. Finally, we provide an algorithm to
design sequences with a predefined structure, i.e. **inverse folding**.


History
-------

Over the past decades since the ViennaRNA Package first saw the light of day as part
of Ivo L. Hofackers PhD project, several authors contributed more and more
algorithm implementations. In 2008, Ronny Lorenz took over the extensive task
to harmonize and simplify the already existing implementations for the sake
of easier feature addition. This eventually lead to version 2.0 of the ViennaRNA
Package. Since then, he (re-)implemented a large portion of the currently existing
library features, such as the new, generalized constraints framework, RNA folding
grammar domain extensions, and the major part of the scripting language interface.

.. grid:: 1 1 2 2

    .. grid-item-card:: Code contributors

        * Daniel Wiegreffe (RNAturtle and RNApuzzler secondary structure layouts)
        * Andreas Gruber (first approach on RNALfold Z-score filtering)
        * Juraj Michalik (non-redundant Boltzmann sampling)
        * Gregor Entzian (neighbor, walk)
        * Mario Koestl (worked on SWIG interface and related unit testing)
        * Dominik Luntzer (pertubation fold)
        * Stefan Badelt (cofold evaluation, RNAdesign.pl, cofold findpath extensions)
        * Stefan Hammer (parts of SWIG interface and corresponding unit tests)
        * Ronny Lorenz (circfold, version 2.0, generic constraints, grammar extensions, and much more)
        * Hakim Tafer (RNAplex, RNAsnoop)
        * Ulrike Mueckstein (RNAup)
        * Stephan Bernhart (RNAcofold, RNAplfold, unpaired probabilities, alifold, and so many more)
        * Stefan Wuchty (RNAsubopt)
        * Ivo Hofacker, Peter Stadler, and Christoph Flamm (almost every implementation up to version 1.8.5)

    .. grid-item-card:: Additional thanks to

        * Sebastian Bonhoeffer's implementation of partition function folding served as a precursor to our part_func.c
        * Manfred Tacker hacked constrained folding into fold.c for the first time
        * Martin Fekete made the first attempts at "alignment folding"
        * Andrea Tanzer and Martin Raden (Mann) for not stopping to report bugs found through comprehensive usage of our applications and RNAlib

        Thanks also to everyone else who helped testing and finding bugs, especially
        Christoph Flamm, Martijn Huynen, Baerbel Krakhofer, and many more.

