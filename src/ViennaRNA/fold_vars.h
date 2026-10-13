#ifndef VIENNA_RNA_PACKAGE_FOLD_VARS_H
#define VIENNA_RNA_PACKAGE_FOLD_VARS_H


#ifdef VRNA_WARN_DEPRECATED
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

#include <ViennaRNA/datastructures/basic.h>
/*  For now, we include model.h by default to provide backwards compatibility
    However, this will most likely change, since fold_vars.h is scheduled to
    vanish from the sources at latest in ViennaRNA Package v3
*/
#include <ViennaRNA/model.h>


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  \file fold_vars.h
 *  \brief Here all all declarations of the global variables used throughout RNAlib
 */

/**
 *  \brief Global switch to activate/deactivate folding with structure constraints
 */
DEPRECATED(extern int    fold_constrained,
           "Use constraints framework instead!");

/**
 *  \brief generate comma seperated output
 */
DEPRECATED(extern int  csv,
           "This flag is obsolete!");

/**
 *  warning this variable will vanish in the future
 *  ribosums will be compiled in instead
 */
DEPRECATED(extern char *RibosumFile,
           "Use ribosum score API directly!");

/**
 *  internal loops of size 2 get energy 0.8Kcal and
 *  no mismatches, default 1
 */
DEPRECATED(extern int  james_rule,
           "This flag is obsolete!");

/**
 *  use logarithmic multiloop energy function
 */
DEPRECATED(extern int  logML,
           "Use logML flag in vrna_md_t instead!");

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
DEPRECATED(extern int  cut_point,
           "Use strand identifiers in vrna_fold_compound_t instead!");

/**
 *  \brief Contains a list of base pairs after a call to fold().
 * 
 *  base_pair[0].i contains the total number of pairs.
 *  \deprecated Do not use this variable anymore!
 */
DEPRECATED(extern bondT  *base_pair,
           "Use backtracking API directly!");

/**
 *  \brief A pointer to the base pair probability matrix
 * 
 *  \deprecated Do not use this variable anymore!
 */
DEPRECATED(extern FLT_OR_DBL *pr,
           "Access base pair probabilities from vrna_fold_compound_t!");

/**
 *  \brief index array to move through pr.
 * 
 *  The probability for base i and j to form a pair is in pr[iindx[i]-j].
 *  \deprecated Do not use this variable anymore!
 */
DEPRECATED(extern int   *iindx,
           "Use iindx field from vrna_fold_compound_t or vrna_idx_row_wise()!");


#endif


#endif
