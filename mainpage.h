/**
  \mainpage ViennaRNA Package core - RNAlib
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
float energy_of_structure(const char *string,
                          const char *structure,
                          int verbosity_level);
\endverbatim
\copybrief energy_of_structure()

\verbatim
float energy_of_circ_structure( const char *string,
                                const char *structure,
                                int verbosity_level);
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
float pf_fold ( char* sequence,
                char* structure)
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
char *get_centroid_struct_pl( int length,
                              double *dist,
                              plist *pl);
\endverbatim
\copybrief get_centroid_struct_pl()

\verbatim
char *get_centroid_struct_pr( int length,
                              double *dist,
                              double *pr);
\endverbatim
\copybrief get_centroid_struct_pr()

\verbatim
double mean_bp_distance_pr( int length,
                            double *pr);
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
float inverse_fold (char *start,
                    char *target)
\endverbatim
\copybrief inverse_fold()

\verbatim
float inverse_pf_fold ( char *start,
                        char *target)
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
SOLUTION *subopt (char *sequence,
                  char *constraint,
                  int *delta,
                  FILE *fp)
\endverbatim
\copybrief subopt()

\verbatim
SOLUTION *subopt_circ ( char *sequence,
                        char *constraint,
                        int *delta,
                        FILE *fp)
\endverbatim
\copybrief subopt_circ()

\verbatim
SOLUTION  *zukersubopt(const char *string);
\endverbatim
\copybrief zukersubopt()

\verbatim
char  *TwoDpfold_pbacktrack ( TwoDpfold_vars *vars,
                              unsigned int d1,
                              unsigned int d2)
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
float cofold (char *sequence,
              char *structure)
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
cofoldF co_pf_fold( char *sequence,
                    char *structure);
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
ConcEnt *get_concentrations(double FEAB,
                            double FEAA,
                            double FEBB,
                            double FEA,
                            double FEB,
                            double * startconc)
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
pu_contrib *pf_unstru ( char *sequence,
                        int max_w)
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
float Lfold ( const char *string,
              char *structure,
              int maxdist)
\endverbatim
\copybrief Lfold()

\verbatim
float aliLfold( const char **strings,
                char *structure,
                int maxdist)
\endverbatim
\copybrief aliLfold()

\verbatim
float Lfoldz (const char *string,
              char *structure,
              int maxdist,
              int zsc,
              double min_z)
\endverbatim
\copybrief Lfoldz()

\verbatim
plist *pfl_fold (
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

Consensus structures can be predicted by a modified version of the
fold() algorithm that takes a set of aligned sequences instead
of a single sequence. The energy function consists of the mean energy
averaged over the sequences, plus a covariance term that favors pairs
with consistent and compensatory mutations and penalizes pairs that
cannot be formed by all structures. For details see \ref hofacker_02 "Hofacker (2002)".

\verbatim
float  alifold (const char **strings,
                char *structure)
\endverbatim
\copybrief alifold()

\verbatim
float  circalifold (const char **strings,
                    char *structure)
\endverbatim
\copybrief circalifold()

\verbatim
void    free_alifold_arrays (void)
\endverbatim
\copybrief free_alifold_arrays()

\verbatim
float   energy_of_alistruct (
            const char **sequences,
            const char *structure,
            int n_seq,
            float *energy)
\endverbatim
\copybrief energy_of_alistruct()

\verbatim
struct pair_info
\endverbatim
\copybrief pair_info

\verbatim
double  cv_fact
\endverbatim
\copybrief cv_fact

\verbatim
double  nc_fact
\endverbatim
\copybrief nc_fact

\see alifold.h for more details

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section mp_Fold_Vars             Global Variables for the Folding Routines

The following global variables change the behavior the folding
algorithms or contain additional information after folding.

\verbatim
int noGU
\endverbatim
\copybrief noGU

\verbatim
int no_closingGU
\endverbatim
\copybrief no_closingGU

\verbatim
int noLonelyPairs
\endverbatim
\copybrief noLonelyPairs

\verbatim
int tetra_loop
\endverbatim
\copybrief tetra_loop

\verbatim
int energy_set
\endverbatim
\copybrief energy_set

\verbatim
float temperature
\endverbatim
\copybrief temperature

\verbatim
int dangles
\endverbatim
\copybrief dangles

\verbatim
char *nonstandards
\endverbatim
\copybrief nonstandards

\verbatim
int cut_point
\endverbatim
\copybrief cut_point

\verbatim
float pf_scale
\endverbatim
\copybrief pf_scale

\verbatim
int fold_constrained
\endverbatim
\copybrief fold_constrained

\verbatim
int do_backtrack
\endverbatim
\copybrief do_backtrack

\verbatim
char backtrack_type
\endverbatim
\copybrief backtrack_type

include fold_vars.h if you want to change any of these variables
from their defaults.

\see fold_vars.h for a more complete and detailed description of all global variables and how to use them

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section mp_Param_Files           Reading Energy Parameters from File

A default set of parameters, identical to the one described in \ref mathews_04 "Mathews et.al. (2004)", is
compiled into the library.\n
Alternately, parameters can be read from and written to a file.

\verbatim
void  read_parameter_file (const char fname[])
\endverbatim
\copybrief read_parameter_file()

\verbatim
void  write_parameter_file (const char fname[])
\endverbatim
\copybrief write_parameter_file

To preserve some backward compatibility the RNAlib also provides functions to convert energy parameter
files from the format used in version 1.4-1.8 into the new format used since version 2.0

\verbatim
void convert_parameter_file (
            const char *iname,
            const char *oname,
            unsigned int options)
\endverbatim
\copybrief convert_parameter_file()

\see read_epars.h and convert_epars.h for detailed description of the available functions

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\ref mp_parse "Next Page: Parsing and Comparing"

\page  mp_parse     Parsing and Comparing - Functions to Manipulate Structures

<h2>Representations of Secondary Structures</h2>

The standard representation of a secondary structure is the <i>bracket
notation</i>, where matching brackets symbolize base pairs and unpaired
bases are shown as dots. Alternatively, one may use two types of node
labels, 'P' for paired and 'U' for unpaired; a dot is then replaced by
'(U)', and each closed bracket is assigned an additional identifier 'P'.
We call this the expanded notation. In \ref fontana_93b "Fontana et al. (1993)" a
condensed
representation of the secondary structure is proposed, the so-called
homeomorphically irreducible tree (HIT) representation. Here a stack is
represented as a single pair of matching brackets labeled 'P' and
weighted by the number of base pairs.  Correspondingly, a contiguous
strain of unpaired bases is shown as one pair of matching brackets
labeled 'U' and weighted by its length.  Generally any string consisting
of matching brackets and identifiers is equivalent to a plane tree with
as many different types of nodes as there are identifiers.

\ref shapiro_88 "Bruce Shapiro (1988)" proposed a coarse grained representation, which,
does not retain the full information of the secondary structure. He
represents the different structure elements by single matching brackets
and labels them as 'H' (hairpin loop), 'I' (interior loop), 'B'
(bulge), 'M' (multi-loop), and 'S' (stack). We extend his alphabet by an
extra letter for external elements 'E'. Again these identifiers may be
followed by a weight corresponding to the number of unpaired bases or
base pairs in the structure element.  All tree representations (except
for the dot-bracket form) can be encapsulated into a virtual root
(labeled 'R'), see the example below.

The following example illustrates the different linear tree representations
used by the package. All lines show the same secondary structure.

\verbatim
a) .((((..(((...)))..((..)))).)).
   (U)(((((U)(U)((((U)(U)(U)P)P)P)(U)(U)(((U)(U)P)P)P)P)(U)P)P)(U)
b) (U)(((U2)((U3)P3)(U2)((U2)P2)P2)(U)P2)(U)
c) (((H)(H)M)B)
   ((((((H)S)((H)S)M)S)B)S)
   (((((((H)S)((H)S)M)S)B)S)E)
d) ((((((((H3)S3)((H2)S2)M4)S2)B1)S2)E2)R)
\endverbatim

Above: Tree representations of secondary structures.  a) Full structure:
the first line shows the more convenient condensed notation which is
used by our programs; the second line shows the rather clumsy expanded
notation for completeness, b) HIT structure, c) different versions of
coarse grained structures: the second line is exactly Shapiro's
representation, the first line is obtained by neglecting the stems.
Since each loop is closed by a unique stem, these two lines are
equivalent.  The third line is an extension taking into account also the
external digits.  d) weighted coarse structure, this time including the
virtual root.

For the output of aligned structures from string editing, different
representations are needed, where we put the label on both sides.
The above examples for tree representations would then look like:

\verbatim
a) (UU)(P(P(P(P(UU)(UU)(P(P(P(UU)(UU)(UU)P)P)P)(UU)(UU)(P(P(UU)(U...
b) (UU)(P2(P2(U2U2)(P2(U3U3)P3)(U2U2)(P2(U2U2)P2)P2)(UU)P2)(UU)
c) (B(M(HH)(HH)M)B)
   (S(B(S(M(S(HH)S)(S(HH)S)M)S)B)S)
   (E(S(B(S(M(S(HH)S)(S(HH)S)M)S)B)S)E)
d) (R(E2(S2(B1(S2(M4(S3(H3)S3)((H2)S2)M4)S2)B1)S2)E2)R)
\endverbatim

Aligned structures additionally contain the gap character '_'.

<h2>Parsing and Coarse Graining of Structures</h2>

Several functions are provided for parsing structures and converting to
different representations.

\verbatim
char  *expand_Full(const char *structure)
\endverbatim
\copybrief expand_Full()

\verbatim
char *b2HIT (const char *structure)
\endverbatim
\copybrief b2HIT()

\verbatim
char *b2C (const char *structure)
\endverbatim
\copybrief b2C()

\verbatim
char *b2Shapiro (const char *structure)
\endverbatim
\copybrief b2Shapiro()

\verbatim
char  *expand_Shapiro (const char *coarse);
\endverbatim
\copybrief expand_Shapiro()

\verbatim
char *add_root (const char *structure)
\endverbatim
\copybrief add_root()

\verbatim
char  *unexpand_Full (const char *ffull)
\endverbatim
\copybrief unexpand_Full()

\verbatim
char  *unweight (const char *wcoarse)
\endverbatim
\copybrief unweight()

\verbatim
void   unexpand_aligned_F (char *align[2])
\endverbatim
\copybrief unexpand_aligned_F()

\verbatim
void   parse_structure (const char *structure)
\endverbatim
\copybrief parse_structure()

\see RNAstruct.h for prototypes and more detailed description

<h2>Distance Measures</h2>

A simple measure of dissimilarity between secondary structures of equal
length is the base pair distance, given by the number of pairs present in
only one of the two structures being compared. I.e. the number of base
pairs that have to be opened or closed to transform one structure into the
other. It is therefore particularly useful for comparing structures on the
same sequence. It is implemented by

\verbatim
int bp_distance(const char *str1,
                const char *str2)
\endverbatim
\copybrief bp_distance()

For other cases a distance measure that allows for gaps is preferable.
We can define distances between structures as edit distances between
trees or their string representations. In the case of string distances
this is the same as "sequence alignment". Given a set of edit operations
and edit costs, the edit distance is given by the minimum sum of the
costs along an edit path converting one object into the other. Edit
distances like these always define a metric. The edit operations used by us
are insertion, deletion and replacement of nodes.
String editing does not pay attention to the matching of brackets, while
in tree editing matching brackets represent a single node of the tree.
Tree editing is therefore usually preferable, although somewhat
slower. String edit distances are always smaller or equal to tree edit
distances.

The different level of detail in the structure representations defined
above naturally leads to different measures of distance. For full
structures we use a cost of 1 for deletion or insertion of an unpaired
base and 2 for a base pair. Replacing an unpaired base for a pair incurs
a cost of 1.

Two cost matrices are provided for coarse grained structures:

\verbatim
/*  Null,   H,   B,   I,   M,   S,   E    */
   {   0,   2,   2,   2,   2,   1,   1},   /* Null replaced */
   {   2,   0,   2,   2,   2, INF, INF},   /* H    replaced */
   {   2,   2,   0,   1,   2, INF, INF},   /* B    replaced */
   {   2,   2,   1,   0,   2, INF, INF},   /* I    replaced */
   {   2,   2,   2,   2,   0, INF, INF},   /* M    replaced */
   {   1, INF, INF, INF, INF,   0, INF},   /* S    replaced */
   {   1, INF, INF, INF, INF, INF,   0},   /* E    replaced */


/* Null,   H,   B,   I,   M,   S,   E   */
   {   0, 100,   5,   5,  75,   5,   5},   /* Null replaced */
   { 100,   0,   8,   8,   8, INF, INF},   /* H    replaced */
   {   5,   8,   0,   3,   8, INF, INF},   /* B    replaced */
   {   5,   8,   3,   0,   8, INF, INF},   /* I    replaced */
   {  75,   8,   8,   8,   0, INF, INF},   /* M    replaced */
   {   5, INF, INF, INF, INF,   0, INF},   /* S    replaced */
   {   5, INF, INF, INF, INF, INF,   0},   /* E    replaced */
\endverbatim

The lower matrix uses the costs given in \ref shapiro_90 "Shapiro (1990)".
All distance functions use the following global variables:

\verbatim
int  cost_matrix;
\endverbatim
\copybrief cost_matrix

\verbatim
int   edit_backtrack;
\endverbatim
\copybrief edit_backtrack

\verbatim
char *aligned_line[4];
\endverbatim
\copybrief aligned_line

\see utils.h, dist_vars.h and stringdist.h for more details

<h3>Functions for Tree Edit Distances</h3>

\verbatim
Tree   *make_tree (char *struc)
\endverbatim
\copybrief make_tree()

\verbatim
float   tree_edit_distance (Tree *T1,
                            Tree *T2) 
\endverbatim
\copybrief tree_edit_distance()

\verbatim
void    free_tree(Tree *t)
\endverbatim
\copybrief free_tree()

\see dist_vars.h and treedist.h for prototypes and more detailed descriptions

<h3>Functions for String Alignment</h3>

\verbatim
swString *Make_swString (char *string)
\endverbatim
\copybrief Make_swString()

\verbatim
float     string_edit_distance (swString *T1,
                                swString *T2)
\endverbatim
\copybrief string_edit_distance()

\see dist_vars.h and stringdist.h for prototypes and more detailed descriptions

<h3>Functions for Comparison of Base Pair Probabilities</h3>

For comparison of base pair probability matrices, the matrices are first
condensed into probability profiles which are the compared by alignment.

\verbatim
float *Make_bp_profile_bppm ( double *bppm,
                              int length)
\endverbatim
\copybrief Make_bp_profile_bppm()

\verbatim
float profile_edit_distance ( const float *T1,
                              const float *T2)
\endverbatim
\copybrief profile_edit_distance()

\see ProfileDist.h for prototypes and more details of the above functions

\ref mp_utils "Next Page: Utilities"

\page  mp_utils     Utilities - Odds and Ends

\anchor toc

<h3>Table of Contents</h3>
<hr>

\li \ref  utils_ss
\li \ref  utils_dot
\li \ref  utils_aln
\li \ref  utils_seq
\li \ref  utils_struc
\li \ref  utils_misc

<hr>

\section utils_ss Producing secondary structure graphs

\verbatim
int PS_rna_plot ( char *string,
                  char *structure,
                  char *file)
\endverbatim
\copybrief PS_rna_plot()

\verbatim
int PS_rna_plot_a (
            char *string,
            char *structure,
            char *file,
            char *pre,
            char *post)
\endverbatim
\copybrief PS_rna_plot_a()

\verbatim
int gmlRNA (char *string,
            char *structure,
            char *ssfile,
            char option)
\endverbatim
\copybrief gmlRNA()

\verbatim
int ssv_rna_plot (char *string,
                  char *structure,
                  char *ssfile)
\endverbatim
\copybrief ssv_rna_plot()

\verbatim
int svg_rna_plot (char *string,
                  char *structure,
                  char *ssfile)
\endverbatim
\copybrief svg_rna_plot()

\verbatim
int xrna_plot ( char *string,
                char *structure,
                char *ssfile)
\endverbatim
\copybrief xrna_plot()

\verbatim
int rna_plot_type
\endverbatim
\copybrief rna_plot_type

Two low-level functions provide direct access to the graph lauyouting
algorithms:

\verbatim
int simple_xy_coordinates ( short *pair_table,
                            float *X,
                            float *Y)
\endverbatim
\copybrief simple_xy_coordinates()

\verbatim
int naview_xy_coordinates ( short *pair_table,
                            float *X,
                            float *Y)
\endverbatim
\copybrief naview_xy_coordinates()

\see PS_dot.h and naview.h for more detailed descriptions.

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section utils_dot Producing (colored) dot plots for base pair probabilities

\verbatim
int PS_color_dot_plot ( char *string,
                        cpair *pi,
                        char *filename)
\endverbatim
\copybrief PS_color_dot_plot()

\verbatim
int PS_color_dot_plot_turn (char *seq,
                            cpair *pi,
                            char *filename,
                            int winSize)
\endverbatim
\copybrief PS_color_dot_plot_turn()

\verbatim
int PS_dot_plot_list (char *seq,
                      char *filename,
                      plist *pl,
                      plist *mf,
                      char *comment)
\endverbatim
\copybrief PS_dot_plot_list()

\verbatim
int PS_dot_plot_turn (char *seq,
                      struct plist *pl,
                      char *filename,
                      int winSize)
\endverbatim
\copybrief PS_dot_plot_turn()

\see PS_dot.h for more detailed descriptions.

\section utils_aln Producing (colored) alignments

\verbatim
int PS_color_aln (
            const char *structure,
            const char *filename,
            const char *seqs[],
            const char *names[])
\endverbatim
\copybrief PS_color_aln()

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section  utils_seq   RNA sequence related utilities

Several functions provide useful applications to RNA sequences

\verbatim
char  *random_string (int l,
                      const char symbols[])
\endverbatim
\copybrief random_string()

\verbatim
int   hamming ( const char *s1,
                const char *s2)
\endverbatim
\copybrief hamming()

\verbatim
void str_DNA2RNA(char *sequence);
\endverbatim
\copybrief str_DNA2RNA()

\verbatim
void str_uppercase(char *sequence);
\endverbatim
\copybrief str_uppercase()

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section utils_struc  RNA secondary structure related utilities

\verbatim
char *pack_structure (const char *struc)
\endverbatim
\copybrief pack_structure()

\verbatim
char *unpack_structure (const char *packed)
\endverbatim
\copybrief unpack_structure()

\verbatim
short *make_pair_table (const char *structure)
\endverbatim
\copybrief make_pair_table()

\verbatim
short *copy_pair_table (const short *pt)
\endverbatim
\copybrief copy_pair_table()

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\section  utils_misc  Miscellaneous Utilities

\verbatim
void print_tty_input_seq (void)
\endverbatim
\copybrief print_tty_input_seq()

\verbatim
void print_tty_constraint_full (void)
\endverbatim
\copybrief print_tty_constraint_full()

\verbatim
void print_tty_constraint (unsigned int option)
\endverbatim
\copybrief print_tty_constraint()

\verbatim
int   *get_iindx (unsigned int length)
\endverbatim
\copybrief get_iindx()

\verbatim
int   *get_indx (unsigned int length)
\endverbatim
\copybrief get_indx()

\verbatim
void constrain_ptypes (
                const char *constraint,
                unsigned int length,
                char *ptype,
                int *BP,
                int min_loop_size,
                unsigned int idx_type)
\endverbatim
\copybrief constrain_ptypes()

\verbatim
char  *get_line(FILE *fp);
\endverbatim
\copybrief get_line()

\verbatim
unsigned int read_record(
                char **header,
                char **sequence,
                char ***rest,
                unsigned int options);
\endverbatim
\copybrief read_record()

\verbatim
char  *time_stamp (void)
\endverbatim
\copybrief time_stamp()

\verbatim
void warn_user (const char message[])
\endverbatim
\copybrief warn_user()

\verbatim
void nrerror (const char message[])
\endverbatim
\copybrief nrerror()

\verbatim
void   init_rand (void)
\endverbatim
\copybrief init_rand()

\verbatim
unsigned short xsubi[3];
\endverbatim
\copybrief xsubi

\verbatim
double urn (void)
\endverbatim
\copybrief urn()

\verbatim
int    int_urn (int from, int to)
\endverbatim
\copybrief int_urn()

\verbatim
void  *space (unsigned size)
\endverbatim
\copybrief space()

\verbatim
void  *xrealloc ( void *p,
                  unsigned size)
\endverbatim
\copybrief xrealloc()

\see  utils.h for a complete overview and detailed description of the utility functions

\htmlonly
<hr>
<a href="#toc">Table of Contents</a>
<hr>
\endhtmlonly

\ref mp_example "Next Page: Examples"

\page  mp_example   Example - A Small Example Program

The following program exercises most commonly used functions of the library.
The program folds two sequences using both the mfe and partition function
algorithms and calculates the tree edit and profile distance of the
resulting structures and base pairing probabilities.

\verbatim
#include  <stdio.h>
#include  <math.h>
#include  "utils.h"
#include  "fold_vars.h"
#include  "fold.h"
#include  "part_func.h"
#include  "inverse.h"
#include  "RNAstruct.h"
#include  "treedist.h"
#include  "stringdist.h"
#include  "ProfileDist.h"

void main()
{
   char *seq1="CGCAGGGAUACCCGCG", *seq2="GCGCCCAUAGGGACGC",
        *struct1,* struct2,* xstruc;
   float e1, e2, tree_dist, string_dist, profile_dist, kT;
   Tree *T1, *T2;
   swString *S1, *S2;
   float **pf1, **pf2;
   FLT_OR_DBL *bppm;
   /* fold at 30C instead of the default 37C */
   temperature = 30.;      /* must be set *before* initializing  */

   /* allocate memory for structure and fold */
   struct1 = (char* ) space(sizeof(char)*(strlen(seq1)+1));
   e1 =  fold(seq1, struct1);

   struct2 = (char* ) space(sizeof(char)*(strlen(seq2)+1));
   e2 =  fold(seq2, struct2);

   free_arrays();     /* free arrays used in fold() */

   /* produce tree and string representations for comparison */
   xstruc = expand_Full(struct1);
   T1 = make_tree(xstruc);
   S1 = Make_swString(xstruc);
   free(xstruc);

   xstruc = expand_Full(struct2);
   T2 = make_tree(xstruc);
   S2 = Make_swString(xstruc);
   free(xstruc);

   /* calculate tree edit distance and aligned structures with gaps */
   edit_backtrack = 1;
   tree_dist = tree_edit_distance(T1, T2);
   free_tree(T1); free_tree(T2);
   unexpand_aligned_F(aligned_line);
   printf("%s\n%s  %3.2f\n", aligned_line[0], aligned_line[1], tree_dist);

   /* same thing using string edit (alignment) distance */
   string_dist = string_edit_distance(S1, S2);
   free(S1); free(S2);
   printf("%s  mfe=%5.2f\n%s  mfe=%5.2f  dist=%3.2f\n",
          aligned_line[0], e1, aligned_line[1], e2, string_dist);

   /* for longer sequences one should also set a scaling factor for
      partition function folding, e.g: */
   kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */
   pf_scale = exp(-e1/kT/strlen(seq1));

   /* calculate partition function and base pair probabilities */
   e1 = pf_fold(seq1, struct1);
   /* get the base pair probability matrix for the previous run of pf_fold() */
   bppm = export_bppm();
   pf1 = Make_bp_profile_bppm(bppm, strlen(seq1));

   e2 = pf_fold(seq2, struct2);
   /* get the base pair probability matrix for the previous run of pf_fold() */
   bppm = export_bppm();
   pf2 = Make_bp_profile(strlen(seq2));

   free_pf_arrays();  /* free space allocated for pf_fold() */

   profile_dist = profile_edit_distance(pf1, pf2);
   printf("%s  free energy=%5.2f\n%s  free energy=%5.2f  dist=%3.2f\n",
          aligned_line[0], e1, aligned_line[1], e2, profile_dist);

   free_profile(pf1); free_profile(pf2);
}
\endverbatim

In a typical Unix environment you would compile this program using:
\verbatim
cc ${OPENMP_CFLAGS} -c example.c -I${hpath}
\endverbatim
and link using
\verbatim
cc ${OPENMP_CFLAGS} -o example -L${lpath} -lRNA -lm
\endverbatim
where \e ${hpath} and \e ${lpath} point to the location of the header
files and library, respectively.
\note As default, the RNAlib is compiled with build-in \e OpenMP multithreading
support. Thus, when linking your own object files to the library you have to pass
the compiler specific \e ${OPENMP_CFLAGS} (e.g. '-fopenmp' for \b gcc) even if your code does not
use openmp specific code. However, in that case the \e OpenMP flags may be ommited when compiling
example.c

\ref mp_ref "Next Page: References"

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

