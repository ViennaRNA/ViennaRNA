Manpages
========

The ViennaRNA Package comes with a number of executable
programs that provide command line interfaces to the most
important algorithms implemented in *RNAlib*.

Find an overview of these programs and their corresponding
manual pages below.

.. toctree::
   :hidden:
   :titlesonly:

   man/RNA2Dfold
   man/RNAaliduplex
   man/RNAalifold
   man/RNAcofold
   man/RNAdistance
   man/RNAdos
   man/RNAduplex
   man/RNAeval
   man/RNAfold
   man/RNAheat
   man/RNAinverse
   man/RNALalifold
   man/RNALfold
   man/RNAmultifold
   man/RNApaln
   man/RNAparconv
   man/RNApdist
   man/RNAPKplex
   man/RNAplex
   man/RNAplfold
   man/RNAplot
   man/RNApvmin
   man/RNAsnoop
   man/RNAsubopt
   man/RNAup


.. rubric:: Main Programs

========================= =====================
Program                   Description
========================= =====================
:doc:`/man/RNA2Dfold`     Compute MFE structure, partition function and
                          representative sample structures of k,l neighborhoods
:doc:`/man/RNAaliduplex`  Predict conserved RNA-RNA interactions between two alignments
:doc:`/man/RNAalifold`    Calculate secondary structures for a set of aligned RNA sequences
:doc:`/man/RNAcofold`     Calculate secondary structures of two RNAs with dimerization
:doc:`/man/RNAdistance`   Calculate distances between RNA secondary structures
:doc:`/man/RNAdos`        Calculate the density of states for each energy band of an RNA
:doc:`/man/RNAduplex`     Compute the structure upon hybridization of two RNA strands
:doc:`/man/RNAeval`       Evaluate free energy of RNA sequences with given secondary structure
:doc:`/man/RNAfold`       Calculate minimum free energy secondary structures and partition function of RNAs
:doc:`/man/RNAheat`       Calculate the specific heat (melting curve) of an RNA sequence
:doc:`/man/RNAinverse`    Find RNA sequences with given secondary structure (sequence design)
:doc:`/man/RNALalifold`   Calculate locally stable secondary structures for a set of aligned RNAs
:doc:`/man/RNALfold`      Calculate locally stable secondary structures of long RNAs
:doc:`/man/RNAmultifold`  Compute thermodynamic properties for interaction complexes of multiple RNAs
:doc:`/man/RNApaln`       RNA alignment based on sequence base pairing propensities
:doc:`/man/RNApdist`      Calculate distances between thermodynamic RNA secondary structures ensembles
:doc:`/man/RNAparconv`    Convert energy parameter files from ViennaRNA 1.8 to 2.0 format
:doc:`/man/RNAPKplex`     Predict RNA secondary structures including pseudoknots
:doc:`/man/RNAplex`       Find targets of a query RNA
:doc:`/man/RNAplfold`     Calculate average pair probabilities for locally stable secondary structures
:doc:`/man/RNAplot`       Draw RNA Secondary Structures in PostScript, SVG, or GML
:doc:`/man/RNApvmin`      Calculate a perturbation vector that minimizes discrepancies between predicted and observed pairing probabilities
:doc:`/man/RNAsnoop`      Find targets of a query H/ACA snoRNA
:doc:`/man/RNAsubopt`     Calculate suboptimal secondary structures of RNAs
:doc:`/man/RNAup`         Calculate the thermodynamics of RNA-RNA interactions
========================= =====================


.. rubric:: Additional Programs

We include the following additional programs in our distribution of the
ViennaRNA Package. Whether or not they are installed together with the
ViennaRNA Package depends on its :doc:`/configuration`.

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


