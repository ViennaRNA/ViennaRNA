#ifndef VIENNA_RNA_PACKAGE_FOLD_VARS_H
#define VIENNA_RNA_PACKAGE_FOLD_VARS_H

#include <ViennaRNA/data_structures.h>
/*  For now, we include model.h by default to provide backwards compatibility
    However, this will most likely change, since fold_vars.h is scheduled to
    vanish from the sources at latest in ViennaRNA Package v3
*/
#include <ViennaRNA/model.h>


/**
 *  \file fold_vars.h
 *  \brief Here all all declarations of the global variables used throughout RNAlib
 */

#ifdef VRNA_WARN_DEPRECATED
# ifdef __GNUC__
#  define DEPRECATED(func) func __attribute__ ((deprecated))
# else
#  define DEPRECATED(func) func
# endif
#else
# define DEPRECATED(func) func
#endif

/**
 *  \brief Global switch to activate/deactivate folding with structure constraints
 */
extern int    fold_constrained;

/**
 *  \brief generate comma seperated output
 */
extern int  csv;

/**
 *  warning this variable will vanish in the future
 *  ribosums will be compiled in instead
 */
extern char *RibosumFile;   

/**
 *  interior loops of size 2 get energy 0.8Kcal and
 *  no mismatches, default 1
 */
extern int  james_rule;

/**
 *  use logarithmic multiloop energy function
 */
extern int  logML;

/**
 *  \brief Marks the position (starting from 1) of the first
 *  nucleotide of the second molecule within the concatenated sequence.
 * 
 *  To evaluate the energy of a duplex structure (a structure formed by two
 *  strands), concatenate the to sequences and set it to the
 *  first base of the second strand in the concatenated sequence.
 *  The default value of -1 stands for single molecule folding. The
 *  cut_point variable is also used by vrna_file_PS_rnaplot() and
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




#endif
