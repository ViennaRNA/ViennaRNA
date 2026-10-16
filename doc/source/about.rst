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

Contributors
------------

This is an open source project that has only been made possible with the help of
many researchers and users that contributed code and documentation.

Here we want to thank everyone who spend time and effort to help thrive of the
ViennaRNA Package through discussions, code contributions, testing or bug finding,
in particular:

..  contributors:: ../../
    :provider: local
    :avatars:
    :names:
    :contributions: emoji
    :order: byCommit
