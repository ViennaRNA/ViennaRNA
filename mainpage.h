/**
  \mainpage Vienna RNA Package
\n

\htmlonly <center> \endhtmlonly

<h2>A Library for folding and comparing RNA secondary structures</h2>

\htmlonly </center> \endhtmlonly

\n

\date     1994-2010
\authors   Ivo Hofacker, Peter Stadler, Ronny Lorenz and some more

<h3>Table of Contents</h3>
<hr>

\li \ref mp_intro
\li \ref mp_fold
\li \ref mp_parse
\li \ref mp_utils
\li \ref mp_example
\li \ref mp_ref

<hr>

\section  mp_intro     Introduction

The core of the Vienna RNA Package is formed by a collection of routines
for the prediction and comparison of RNA secondary structures. These
routines can be accessed through stand-alone programs, such as RNAfold,
RNAdistance etc., which should be sufficient for most users. For those who
wish to develop their own programs we provide a library which can be
linked to your own code.

This document describes the library and will be primarily useful to
programmers. However, it also contains details about the implementation
that may be of interest to advanced users. The stand-alone programs are
described in separate man pages. The latest version of the package
including source code and html versions of the documentation can be found
at
\n\n
http://www.tbi.univie.ac.at/~ivo/RNA/


\page mp_fold Folding Routines - Functions for Folding RNA Secondary Structures

<h3>Table of Contents</h3>
<hr>

\li \ref mp_mfe_Fold
\li \ref mp_PF_Fold
\li \ref mp_Inverse_Fold
\li \ref mp_Suboptimal_folding
\li \ref mp_Cofolding
\li \ref mp_Local_Fold
\li \ref mp_Alignment_Fold
\li \ref mp_Fold_Vars
\li \ref mp_Param_Files

<hr>

\section mp_mfe_Fold              Calculating Minimum Free Energy Structures

The library provides a fast dynamic programming minimum free energy
folding algorithm as described by \ref zuker_81 "Zuker & Stiegler (1981)".

Associated functions are:

\verbatim
float fold (char* sequence, char* structure)
\endverbatim

circfold(), energy_of_structure(), energy_of_circ_structure(),
free_arrays(), update_fold_params() and many more.
\see fold.h for a complete list of available functions.

\section mp_PF_Fold               Calculating Partition Functions and Pair Probabilities

Instead of the minimum free energy structure the partition function of
all possible structures and from that the pairing probability for every
possible pair can be calculated, using a dynamic programming algorithm
as described by \ref mccaskill_90 "McCaskill (1990)". The following
functions are provided:

\verbatim float pf_fold (char* sequence, char* structure)\endverbatim

\section mp_Inverse_Fold          Searching for Predefined Structures

\section mp_Suboptimal_folding    Enumerating Suboptimal Structures

\section mp_Cofolding             Folding of 2 RNA molecules

\section mp_Local_Fold            Predicting local structures of large sequences

\section mp_Alignment_Fold        Predicting Consensus Structures from Alignment

\section mp_Fold_Vars             Global Variables for the Folding Routines

\section mp_Param_Files           Reading Energy Parameters from File

\page  mp_parse     Parsing and Comparing - Functions to Manipulate Structures
\page  mp_utils     Utilities - Odds and Ends
\page  mp_example   Example - A Small Example Program

\page  mp_ref       References

-# \anchor mathews_99 D.H. Mathews, J. Sabina, M. Zuker and H. Turner (1999)\n
   Expanded sequence dependence of thermodynamic parameters provides
   robust prediction of RNA secondary structure, JMB, 288: 911-940
-# \anchor zuker_81 Zuker and P. Stiegler (1981)\n
   Optimal  computer  folding  of large RNA sequences using
   thermodynamic and auxiliary information, Nucl Acid Res 9: 133-148
-# \anchor dimitrov_04 D.A. Dimitrov, M.Zuker(2004)\n
   Prediction of hybridization and melting for double stranded nucleic
   acids, Biophysical J. 87: 215-226,
-# \anchor mccaskill_90 J.S. McCaskill (1990)\n
   The equilibrium partition function and base pair binding
   probabilities for RNA secondary structures, Biopolymers 29: 1105-1119
-# \anchor turner_88 D.H. Turner, N. Sugimoto and S.M. Freier (1988)\n
   RNA structure prediction, Ann Rev Biophys Biophys Chem 17: 167-192
-# \anchor jaeger_89 J.A. Jaeger, D.H. Turner and M. Zuker (1989)\n
   Improved predictions of secondary structures for RNA,
   Proc. Natl. Acad. Sci. 86: 7706-7710
-# \anchor he_91 L. He, R. Kierzek, J. SantaLucia, A.E. Walter and D.H. Turner (1991)\n
   Nearest-Neighbor Parameters For GU Mismatches,
   Biochemistry 30: 11124-11132
-# \anchor peritz_91 A.E. Peritz, R. Kierzek, N, Sugimoto, D.H. Turner (1991)\n
   Thermodynamic Study of Internal Loops in Oligoribonucleotides ... ,
   Biochemistry 30: 6428--6435
-# \anchor walter_94 A. Walter, D. Turner, J. Kim, M. Lyttle, P. M&uuml;ller, D. Mathews and M. Zuker (1994)\n
   Coaxial stacking of helices enhances binding of Oligoribonucleotides..,
   Proc. Natl. Acad. Sci. 91: 9218-9222
-# \anchor shapiro_88 B.A. Shapiro, (1988)\n
   An algorithm for comparing multiple  RNA secondary structures,
   CABIOS 4, 381-393
-# \anchor shapiro_90 B.A. Shapiro and K. Zhang (1990)\n
   Comparing multiple RNA secondary structures using tree comparison,
   CABIOS 6, 309-318
-# \anchor bruccoleri_88 R. Bruccoleri and G. Heinrich (1988)\n
   An improved algorithm for nucleic acid secondary structure display,
   CABIOS 4, 167-173
-# \anchor fontana_93a W. Fontana , D.A.M. Konings, P.F. Stadler, P. Schuster (1993) \n
   Statistics of RNA secondary structures, Biopolymers 33, 1389-1404
-# \anchor fontana_93b W. Fontana, P.F. Stadler, E.G. Bornberg-Bauer, T. Griesmacher, I.L.
   Hofacker, M. Tacker, P. Tarazona, E.D. Weinberger, P. Schuster (1993)\n
   RNA folding and combinatory landscapes, Phys. Rev. E 47: 2083-2099
-# \anchor hofacker_94a I.L. Hofacker, W. Fontana, P.F. Stadler, S. Bonhoeffer, M. Tacker, P.
   Schuster (1994) Fast Folding and Comparison of RNA Secondary Structures.
   Monatshefte f. Chemie 125: 167-188
-# \anchor hofacker_94b I.L. Hofacker (1994) The Rules of the Evolutionary Game for RNA:
   A Statistical Characterization of the Sequence to Structure Mapping in RNA.
   PhD Thesis, University of Vienna.
-# \anchor hofacker_02 I.L. Hofacker, M. Fekete, P.F. Stadler (2002).
   Secondary Structure Prediction for Aligned RNA Sequences.
   J. Mol. Biol. 319:1059-1066
-# \anchor adams_79 D. Adams (1979)\n
   The hitchhiker's guide to the galaxy, Pan Books, London

**/

