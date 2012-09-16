/*
#############################################################
# The next comment is used to order the modules correctly   #
#############################################################
*/

/**
 *  \defgroup   folding_routines          RNA Secondary Structure Folding
 *  \brief This module contains all functions related to thermodynamic folding of RNAs
 *
 *  \defgroup   mfe_fold                  Calculating Minimum Free Energy (MFE) Structures
 *  \ingroup    folding_routines
 *
 *  \defgroup   pf_fold                   Calculating Partition Functions and Pair Probabilities
 *  \ingroup    folding_routines
 *
 *  \defgroup   mea_fold                  Compute the structure with maximum expected accuracy (MEA)
 *  \ingroup    pf_fold
 *
 *  \defgroup   centroid_fold             Compute the centroid structure
 *  \ingroup    pf_fold
 *
 *  \defgroup   subopt_fold               Enumerating Suboptimal Structures
 *  \ingroup    folding_routines
 *
 *  \defgroup   subopt_zuker              Suboptimal structures according to Zuker et al. 1989
 *  \ingroup    subopt_fold
 *
 *  \defgroup   subopt_wuchty             Suboptimal structures within an energy band arround the MFE
 *  \ingroup    subopt_fold
 *
 *  \defgroup   subopt_stochbt            Stochastic backtracking in the Ensemble
 *  \ingroup    subopt_fold
 *
 *  \defgroup   cofold                    Calculate Secondary Structures of two RNAs upon Dimerization
 *  \ingroup    folding_routines
 *
 *  \defgroup   mfe_cofold                MFE Structures of two hybridized Sequences
 *  \ingroup    cofold mfe_fold
 *
 *  \defgroup   pf_cofold                 Partition Function for two hybridized Sequences
 *  \ingroup    cofold pf_fold
 *
 *  \defgroup   up_cofold                 Partition Function for two hybridized Sequences as a stepwise Process
 *  \ingroup    cofold pf_fold
 *
 *  \defgroup   consensus_fold            Predicting Consensus Structures from Alignment(s)
 *  \ingroup    folding_routines
 *
 *  \defgroup   consensus_mfe_fold        MFE Consensus Structures for Sequence Alignment(s)
 *  \ingroup    consensus_fold mfe_fold
 *
 *  \defgroup   consensus_pf_fold         Partition Function and Base Pair Probabilities for Sequence Alignment(s)
 *  \ingroup    consensus_fold pf_fold
 *
 *  \defgroup   consensus_stochbt         Stochastic Backtracking of Consensus Structures from Sequence Alignment(s)
 *  \ingroup    consensus_fold subopt_stochbt
 *
 *  \defgroup   local_fold                Predicting Locally stable structures of large sequences
 *  \ingroup    folding_routines
 *
 *  \defgroup   local_mfe_fold            Local MFE structure Prediction and Z-scores
 *  \ingroup    local_fold mfe_fold
 *
 *  \defgroup   local_pf_fold             Partition functions for locally stable secondary structures
 *  \ingroup    local_fold pf_fold
 *
 *  \defgroup   local_consensus_fold      Local MFE consensus structures for Sequence Alignments
 *  \ingroup    local_fold consensus_fold
 *
 *  \defgroup   energy_parameters         Change and Precalculate Energy Parameter Sets and Boltzmann Factors
 *  \ingroup    folding_routines
 *
 *  \defgroup   energy_parameters_rw      Reading/Writing energy parameter sets from/to File
 *  \ingroup    energy_parameters
 *
 *  \defgroup   energy_parameters_convert Converting energy parameter files
 *  \ingroup    energy_parameters_rw
 *
 *  \defgroup   eval                      Energy evaluation
 *  \ingroup    folding_routines
 *
 *  \defgroup   inverse_fold              Searching Sequences for Predefined Structures
 *  \ingroup    folding_routines
 *
 *  \defgroup   class_fold                Classified Dynamic Programming
 *  \ingroup    folding_routines
 *
 *  \defgroup   kl_neighborhood           Distance based partitioning of the Secondary Structure Space
 *  \ingroup    class_fold
 *
 *  \defgroup   kl_neighborhood_mfe       Calculating MFE representatives of a Distance Based Partitioning
 *  \ingroup    kl_neighborhood mfe_fold
 *
 *  \defgroup   kl_neighborhood_pf        Calculate Partition Functions of a Distance Based Partitioning
 *  \ingroup    kl_neighborhood pf_fold
 *
 *  \defgroup   kl_neighborhood_stochbt   Stochastic Backtracking of Structures from Distance Based Partitioning
 *  \ingroup    kl_neighborhood subopt_stochbt
 *
 *  \defgroup   dos                       Compute the Density of States
 *  \ingroup    class_fold
 *
 *  \defgroup   parse                     Parsing and Comparing - Functions to Manipulate Structures
 */


/*
#############################################################
# Now the mainpage text is following                        #
#############################################################
*/

/**
  \mainpage ViennaRNA Package core - RNAlib
\n

\htmlonly <center> \endhtmlonly

<h2>A Library for folding and comparing RNA secondary structures</h2>

\htmlonly </center> \endhtmlonly

\n

\date     1994-2012
\authors   Ivo Hofacker, Peter Stadler, Ronny Lorenz and many more

<h3>Table of Contents</h3>
<hr>

\li \ref mp_intro
\li \ref folding_routines
\li \ref mp_parse
\li \ref mp_utils
\li \ref mp_example
\li \ref mp_ref

<hr>

\section  mp_intro     Introduction

The core of the Vienna RNA Package (\cite lorenz:2011, \cite hofacker:1994)
is formed by a collection of routines
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



\page  mp_parse     Parsing and Comparing - Functions to Manipulate Structures

<h2>Representations of Secondary Structures</h2>

The standard representation of a secondary structure is the <i>bracket
notation</i>, where matching brackets symbolize base pairs and unpaired
bases are shown as dots. Alternatively, one may use two types of node
labels, 'P' for paired and 'U' for unpaired; a dot is then replaced by
'(U)', and each closed bracket is assigned an additional identifier 'P'.
We call this the expanded notation. In \cite fontana:1993b a condensed
representation of the secondary structure is proposed, the so-called
homeomorphically irreducible tree (HIT) representation. Here a stack is
represented as a single pair of matching brackets labeled 'P' and
weighted by the number of base pairs.  Correspondingly, a contiguous
strain of unpaired bases is shown as one pair of matching brackets
labeled 'U' and weighted by its length.  Generally any string consisting
of matching brackets and identifiers is equivalent to a plane tree with
as many different types of nodes as there are identifiers.

Bruce Shapiro proposed a coarse grained representation \cite shapiro:1988,
which, does not retain the full information of the secondary structure. He
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

The lower matrix uses the costs given in \cite shapiro:1990.
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

\code{.c}
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  "utils.h"
#include  "fold_vars.h"
#include  "fold.h"
#include  "part_func.h"
#include  "inverse.h"
#include  "RNAstruct.h"
#include  "treedist.h"
#include  "stringdist.h"
#include  "profiledist.h"

void main()
{
   char *seq1="CGCAGGGAUACCCGCG", *seq2="GCGCCCAUAGGGACGC",
        *struct1,* struct2,* xstruc;
   float e1, e2, tree_dist, string_dist, profile_dist, kT;
   Tree *T1, *T2;
   swString *S1, *S2;
   float *pf1, *pf2;
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
   pf2 = Make_bp_profile_bppm(bppm, strlen(seq2));

   free_pf_arrays();  /* free space allocated for pf_fold() */

   profile_dist = profile_edit_distance(pf1, pf2);
   printf("%s  free energy=%5.2f\n%s  free energy=%5.2f  dist=%3.2f\n",
          aligned_line[0], e1, aligned_line[1], e2, profile_dist);

   free_profile(pf1); free_profile(pf2);
}
\endcode

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


**/

