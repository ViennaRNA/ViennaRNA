/**
  \mainpage Vienna RNA Package
\n

\htmlonly <center> \endhtmlonly

<h2>A Library for folding and comparing RNA secondary structures</h2>

\htmlonly </center> \endhtmlonly

\n

\date     1994-2010
\authors   Ivo Hofacker, Peter Stadler, Ronny Lorenz and many more

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

\anchor toc

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
float fold (char* sequence, char* structure);
\endverbatim
\copybrief fold()

\verbatim
float circfold (char* sequence, char* structure);
\endverbatim
\copybrief circfold()

\verbatim
float energy_of_structure(const char *string, const char *structure, int verbosity_level);
\endverbatim
\copybrief energy_of_structure()

\verbatim
float energy_of_circ_structure(const char *string, const char *structure, int verbosity_level);
\endverbatim
\copybrief energy_of_circ_structure()

\verbatim
void  update_fold_params(void);
\endverbatim
\copybrief update_fold_params()

\verbatim
void  free_arrays(void);
\endverbatim
\copybrief free_arrays()

\see fold.h, cofold.h, 2Dfold.h, Lfold.h, alifold.h and subopt.h for a complete list of available functions.

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section mp_PF_Fold               Calculating Partition Functions and Pair Probabilities

Instead of the minimum free energy structure the partition function of
all possible structures and from that the pairing probability for every
possible pair can be calculated, using a dynamic programming algorithm
as described by \ref mccaskill_90 "McCaskill (1990)". The following
functions are provided:

\verbatim
float pf_fold (char* sequence, char* structure)
\endverbatim
\copybrief pf_fold()

\verbatim
void free_pf_arrays (void)
\endverbatim
\copybrief free_pf_arrays()

\verbatim
void update_pf_params (int length)
\endverbatim
\copybrief update_pf_params()

\verbatim
char *get_centroid_struct_pl(int length, double *dist, plist *pl);
\endverbatim
\copybrief get_centroid_struct_pl()

\verbatim
char *get_centroid_struct_pr(int length, double *dist, double *pr);
\endverbatim
\copybrief get_centroid_struct_pr()

\verbatim
double mean_bp_distance_pr(int length, double *pr);
\endverbatim
\copybrief mean_bp_distance_pr

\see part_func.h, part_func_co.h, part_func_up.h, 2Dpfold.h, LPfold.h, alifold.h and MEA.h for a complete list of available functions.

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly


\section mp_Inverse_Fold          Searching for Predefined Structures

We provide two functions that search for sequences with a given
structure, thereby inverting the folding routines.

\verbatim
float inverse_fold (char *start, char *target)
\endverbatim
\copybrief inverse_fold()

\verbatim
float inverse_pf_fold (char *start, char *target)
\endverbatim
\copybrief inverse_pf_fold()

The following global variables define the behavior or show the results of the
inverse folding routines:

\verbatim
char *symbolset
\endverbatim
\copybrief symbolset

\see inverse.h for more details and a complete list of available functions.

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section mp_Suboptimal_folding    Enumerating Suboptimal Structures

\verbatim
SOLUTION *subopt (char *sequence, char *constraint, int *delta, FILE *fp)
\endverbatim
\copybrief subopt()

\verbatim
SOLUTION *subopt_circ (char *sequence, char *constraint, int *delta, FILE *fp)
\endverbatim
\copybrief subopt_circ()

\verbatim
SOLUTION  *zukersubopt(const char *string);
\endverbatim
\copybrief zukersubopt()

\verbatim
char  *TwoDpfold_pbacktrack (TwoDpfold_vars *vars, unsigned int d1, unsigned int d2)
\endverbatim
\copybrief TwoDpfold_pbacktrack()

\verbatim
char  *alipbacktrack (double *prob)
\endverbatim
\copybrief alipbacktrack()

\verbatim
char    *pbacktrack(char *sequence);
\endverbatim
\copybrief pbacktrack()

\verbatim
char    *pbacktrack_circ(char *sequence);
\endverbatim
\copybrief pbacktrack_circ()

\see subopt.h, part_func.h, alifold.h and 2Dpfold.h for more detailed descriptions

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section mp_Cofolding             Predicting hybridization structures of two molecules
The function of an RNA molecule often depends on its interaction with
other RNAs. The following routines therefore allow to predict structures
formed by two RNA molecules upon hybridization.\n
One approach to co-folding two RNAs consists of concatenating the two
sequences and keeping track of the concatenation point in all energy
evaluations. Correspondingly, many of the cofold() and
co_pf_fold() routines below take one sequence string as argument
and use the the global variable #cut_point to mark the concatenation
point. Note that while the <i>RNAcofold</i> program uses the '&' character
to mark the chain break in its input, you should not use an '&' when using
the library routines (set #cut_point instead).\n
In a second approach to co-folding two RNAs, cofolding is seen as a
stepwise process. In the first step the probability of an unpaired region
is calculated and in a second step this probability of an unpaired region
is  multiplied with the probability of an interaction between the two RNAs.
This approach is implemented for the interaction between a long
target sequence and a short ligand RNA. Function pf_unstru() calculates
the partition function over all unpaired regions in the input
sequence. Function pf_interact(), which calculates the
partition function over all possible interactions between two
sequences, needs both sequence as separate strings as input.

\verbatim
int cut_point
\endverbatim
\copybrief cut_point

\verbatim
float cofold (char *sequence, char *structure)
\endverbatim
\copybrief cofold()

\verbatim
void  free_co_arrays (void)
\endverbatim
\copybrief free_co_arrays()

<b>Partition Function Cofolding</b>

To simplify the implementation the partition function computation is done
internally in a null model that does not include the duplex initiation
energy, i.e. the entropic penalty for producing a dimer from two
monomers). The resulting free energies and pair probabilities are initially
relative to that null model. In a second step the free energies can be
corrected to include the dimerization penalty, and the pair probabilities
can be divided into the conditional pair probabilities given that a re
dimer is formed or not formed.

\verbatim
cofoldF co_pf_fold(char *sequence, char *structure);
\endverbatim
\copybrief co_pf_fold()

\verbatim
void    free_co_pf_arrays(void);
\endverbatim
\copybrief free_co_pf_arrays()

<b>Cofolding all Dimeres, Concentrations</b>

After computing the partition functions of all possible dimeres one
can compute the probabilities of base pairs, the concentrations out of
start concentrations and sofar and soaway.

\verbatim
void  compute_probabilities(
              double FAB,
              double FEA,
              double FEB,
              struct plist  *prAB,
              struct plist  *prA,
              struct plist  *prB,
              int Alength)
\endverbatim
\copybrief compute_probabilities()

\verbatim
ConcEnt *get_concentrations(
              double FEAB, double FEAA, double FEBB,
              double FEA, double FEB, double * startconc)
\endverbatim
\copybrief get_concentrations()

<b>Partition Function Cofolding as a stepwise process</b>

In this approach to cofolding the interaction between two RNA molecules is
seen as a stepwise process. In a first step, the target molecule has to
adopt a structure in which a binding site is accessible. In a second step,
the ligand molecule will hybridize with a region accessible to an
interaction. Consequently the algorithm is designed as a two step process:
The first step is the calculation of the probability
that a region within the target is unpaired, or equivalently, the
calculation of the free energy needed to expose a region. In the second step
we compute the free energy of an interaction for every possible binding
site.
Associated functions are:


\verbatim
pu_contrib *pf_unstru (char *sequence, int max_w)
\endverbatim
\copybrief pf_unstru()

\verbatim
void  free_pu_contrib_struct (pu_contrib *pu)
\endverbatim
\copybrief free_pu_contrib_struct()

\verbatim
interact *pf_interact(
              const char *s1,
              const char *s2,
              pu_contrib *p_c,
              pu_contrib *p_c2,
              int max_w,
              char *cstruc,
              int incr3,
              int incr5)
\endverbatim
\copybrief pf_interact()

\verbatim
void free_interact (interact *pin)
\endverbatim
\copybrief free_interact()

\see cofold.h, part_func_co.h and part_func_up.h for more details

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section mp_Local_Fold            Predicting local structures of large sequences
Local structures can be predicted by a modified version of the
fold() algorithm that restricts the span of all base pairs.

\verbatim
float Lfold(const char *string, char *structure, int maxdist)
\endverbatim
\copybrief Lfold()

\verbatim
float aliLfold(const char **strings, char *structure, int maxdist)
\endverbatim
\copybrief aliLfold()

\verbatim
float Lfoldz(const char *string, char *structure, int maxdist, int zsc, double min_z)
\endverbatim
\copybrief Lfoldz()

\verbatim
plist *pfl_fold(
            char *sequence,
            int winSize,
            int pairSize,
            float cutoffb,
            double **pU,
            struct plist **dpp2,
            FILE *pUfp,
            FILE *spup)
\endverbatim
\copybrief pfl_fold()

\see Lfold.h and LPfold.h for more details

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section mp_Alignment_Fold        Predicting Consensus Structures from Alignment

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section mp_Fold_Vars             Global Variables for the Folding Routines

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section mp_Param_Files           Reading Energy Parameters from File

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\page  mp_parse     Parsing and Comparing - Functions to Manipulate Structures
\page  mp_utils     Utilities - Odds and Ends
\page  mp_example   Example - A Small Example Program

\page  mp_ref       References

-# \anchor mathews_04 D.H. Mathews, M. D. Disney, J.L. Childs, S.J. Schroeder, M. Zuker, D.H. Turner (2004)\n
   Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of
   RNA secondary structure, Proc Natl Acad Sci U S A, 101(19):7287-92
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

