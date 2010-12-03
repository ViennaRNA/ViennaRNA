#ifndef __VIENNA_RNA_PACKAGE_FOLD_VARS_H__
#define __VIENNA_RNA_PACKAGE_FOLD_VARS_H__

#include "data_structures.h"

/** \file fold_vars.h  */

/* to use floats instead of doubles in pf_fold() comment next line */
#define LARGE_PF
#ifdef  LARGE_PF
#define FLT_OR_DBL double
#else
#define FLT_OR_DBL float
#endif

#define PUBLIC
#define PRIVATE static

/**
 *  \brief Global switch to activate/deactivate folding with structure constraints
 */
extern int    fold_constrained;

/**
 *  \brief Global switch to avoid/allow helices of length 1
 * 
 *  Disallow all pairs which can @strong{only} occur as lonely pairs (i.e. as helix
 *  of length 1). This avoids lonely base pairs in the predicted structures in
 *  most cases.
 */
extern int    noLonelyPairs;

/**
 *  \brief Switch the energy model for dangling end contributions (0, 1, 2, 3)
 * 
 *  If set to 0 no stabilizing energies are assigned to bases adjacent to
 *  helices in free ends and multiloops (so called dangling ends). Normally
 *  (dangles = 1) dangling end energies are assigned only to unpaired
 *  bases and a base cannot participate simultaneously in two dangling ends. In
 *  the partition function algorithm pf_fold() these checks are neglected.
 *  If #dangles is set to 2, all folding routines will follow this convention.
 *  This treatment of dangling ends gives more favorable energies to helices
 *  directly adjacent to one another, which can be beneficial since such
 *  helices often do engage in stabilizing interactions through co-axial
 *  stacking.\n
 *  If dangles = 3 co-axial stacking is explicitly included for
 *  adjacent helices in mutli-loops. The option affects only mfe folding
 *  and energy evaluation (fold() and energy_of_structure()), as
 *  well as suboptimal folding (subopt()) via re-evaluation of energies.
 *  Co-axial stacking with one intervening mismatch is not considered so far.
 * 
 *  Default is 2 in most algorithms, partition function algorithms can only handle 0 and 2
 */
extern int  dangles;

/**
 *  \brief Global switch to forbid/allow GU base pairs at all
 */
extern int  noGU;

/**
 *  \brief GU allowed only inside stacks if set to 1
 */
extern int  no_closingGU;

/**
 *  \brief Include special stabilizing energies for some tri-, tetra- and hexa-loops;
 * 
 *  default is 1.
 */
extern int  tetra_loop;

/**
 *  \brief 0 = BP; 1=any mit GC; 2=any mit AU-parameter
 * 
 *  If set to 1 or 2: fold sequences from an artificial alphabet ABCD..., where A
 *  pairs B, C pairs D, etc. using either GC (1) or AU parameters (2);
 *  default is 0, you probably don't want to change it.
 */
extern int  energy_set;

/**
 *  \brief backward compatibility variable.. this does not effect anything
 */
extern int  circ;

/**
 *  \brief generate comma seperated output
 */
extern int  csv;

extern int oldAliEn;        /* use old alifold energies (with gaps) */
extern int ribo;            /* use ribosum matrices */
extern char *RibosumFile;   /* warning this variable will vanish in the future
                               ribosums will be compiled in instead */
/**
 *  \brief contains allowed non standard base pairs
 * 
 *  Lists additional base pairs that will be allowed to form in addition to
 *  GC, CG, AU, UA, GU and UG. Nonstandard base pairs are given a stacking
 *  energy of 0.
 */
extern char *nonstandards;

/**
 *  \brief Rescale energy parameters to a temperature in degC.
 * 
 *  Default is 37C. You have to call the update_..._params() functions after
 *  changing this parameter.
 */
extern double temperature;

extern int  james_rule;     /* interior loops of size 2 get energy 0.8Kcal and
                               no mismatches, default 1 */
extern int  logML;          /* use logarithmic multiloop energy function */

/**
 *  \brief Marks the position (starting from 1) of the first
 *  nucleotide of the second molecule within the concatenated sequence.
 * 
 *  To evaluate the energy of a duplex structure (a structure formed by two
 *  strands), concatenate the to sequences and set it to the
 *  first base of the second strand in the concatenated sequence.
 *  The default value of -1 stands for single molecule folding. The
 *  cut_point variable is also used by PS_rna_plot() and
 *  PS_dot_plot() to mark the chain break in postscript plots.
 */
extern int  cut_point;

/**
 *  \brief Contains a list of base pairs after a call to fold().
 * 
 *  base_pair[0].i contains the total number of pairs.
 *  \deprecated Do not use this variable anymore!
 */
extern bondT  *base_pair;

/**
 *  \brief A pointer to the base pair probability matrix
 * 
 *  \deprecated Do not use this variable anymore!
 */
extern FLT_OR_DBL *pr;

/**
 *  \brief index array to move through pr.
 * 
 *  The probability for base i and j to form a pair is in pr[iindx[i]-j].
 *  \deprecated Do not use this variable anymore!
 */
extern int   *iindx;

/**
 *  \brief A scaling factor used by pf_fold() to avoid overflows.
 * 
 *  Should be set to approximately \f$exp{((-F/kT)/length)}\f$, where \f$F\f$ is an estimate
 *  for the ensemble free energy, for example the minimum free energy. You must
 *  call update_pf_params() after changing this parameter.\n
 *  If pf_scale is -1 (the default) , an estimate will be provided
 *  automatically when computing partition functions, e.g. pf_fold()
 *  The automatic estimate is usually insufficient for sequences more
 *  than a few hundred bases long.
 */
extern double pf_scale;

/**
 *  \brief do backtracking, i.e. compute secondary structures or base pair probabilities
 * 
 *  If 0, do not calculate pair probabilities in pf_fold(); this is about
 *  twice as fast. Default is 1.
 */
extern int    do_backtrack;

/**
 *  \brief A backtrack array marker for inverse_fold()
 * 
 *  If set to 'C': force (1,N) to be paired,
 *  'M' fold as if the sequence were inside a multi-loop. Otherwise ('F') the
 *  usual mfe structure is computed.
 */
extern char backtrack_type;

char * option_string(void);
#endif
