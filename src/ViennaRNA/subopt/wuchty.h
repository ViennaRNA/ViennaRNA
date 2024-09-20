/* subopt.h */
#ifndef VIENNA_RNA_PACKAGE_SUBOPT_WUCHTY_H
#define VIENNA_RNA_PACKAGE_SUBOPT_WUCHTY_H

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

/**
 *  @file ViennaRNA/subopt/wuchty.h
 *  @ingroup subopt_and_representatives
 *  @brief RNAsubopt and density of states declarations
 */

#include <stdio.h>

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/subopt/basic.h>


/**
 *  @brief Maximum density of states discretization for subopt
 */
#define MAXDOS                1000
#define VRNA_UNSORTED                           0
#define VRNA_SORT_BY_ENERGY_LEXICOGRAPHIC_ASC   1
#define VRNA_SORT_BY_ENERGY_ASC                 2

/**
 *  @brief Returns list of subopt structures or writes to fp
 *
 *  This function produces <b>all</b> suboptimal secondary structures within
 *  'delta' * 0.01 kcal/mol of the optimum, see @rstinline :cite:t:`wuchty:1999` @endrst. The results
 *  are either directly written to a 'fp' (if 'fp' is not NULL), or
 *  (fp==NULL) returned in a #vrna_subopt_solution_t * list terminated
 *  by an entry were the 'structure' member is NULL.
 *
 *  @ingroup subopt_wuchty
 *
 *  @note This function requires all multibranch loop DP matrices for unique
 *        multibranch loop backtracing. Therefore, the supplied #vrna_fold_compound_t
 *        @p fc (argument 1) must be initialized with #vrna_md_t.uniq_ML = 1, for
 *        instance like this:
 *        @code
 * vrna_md_t md;
 * vrna_md_set_default(&md);
 * md.uniq_ML = 1;
 *
 * vrna_fold_compound_t *fc=vrna_fold_compound("GGGGGGAAAAAACCCCCC", &md, VRNA_OPTION_DEFAULT);
 *        @endcode
 *
 *  @see vrna_subopt_cb(), vrna_subopt_zuker()
 *
 *  @param  fc
 *  @param  delta
 *  @param  sorted  Sort results by energy in ascending order
 *  @param  fp
 *  @return
 */
vrna_subopt_solution_t *
vrna_subopt(vrna_fold_compound_t  *fc,
            int                   delta,
            int                   sorted,
            FILE                  *fp);


/**
 *  @brief  Generate suboptimal structures within an energy band arround the MFE
 *
 *  This is the most generic implementation of the suboptimal structure generator
 *  according to @rstinline :cite:t:`wuchty:1999` @endrst. Identical to vrna_subopt(), it computes all
 *  secondary structures within an energy band @p delta arround the MFE. However,
 *  this function does not print the resulting structures and their corresponding
 *  free energies to a file pointer, or returns them as a list. Instead, it calls
 *  a user-provided callback function which it passes the structure in dot-bracket
 *  format, the corresponding free energy in kcal/mol, and a user-provided data
 *  structure each time a structure was backtracked successfully. This function
 *  indicates the final output, i.e. the end of the backtracking procedure by
 *  passing NULL instead of an actual dot-bracket string to the callback.
 *
 *  @ingroup subopt_wuchty
 *
 *  @note This function requires all multibranch loop DP matrices for unique
 *        multibranch loop backtracing. Therefore, the supplied #vrna_fold_compound_t
 *        @p fc (argument 1) must be initialized with #vrna_md_t.uniq_ML = 1, for
 *        instance like this:
 *        @code
 * vrna_md_t md;
 * vrna_md_set_default(&md);
 * md.uniq_ML = 1;
 *
 * vrna_fold_compound_t *fc=vrna_fold_compound("GGGGGGAAAAAACCCCCC", &md, VRNA_OPTION_DEFAULT);
 *        @endcode
 *
 *  @see vrna_subopt_result_f, vrna_subopt(), vrna_subopt_zuker()
 *
 *  @param  fc      fold compount with the sequence data
 *  @param  delta   Energy band arround the MFE in 10cal/mol, i.e. deka-calories
 *  @param  cb      Pointer to a callback function that handles the backtracked structure and its free energy in kcal/mol
 *  @param  data    Pointer to some data structure that is passed along to the callback
 */
void
vrna_subopt_cb(vrna_fold_compound_t *fc,
               int                  delta,
               vrna_subopt_result_f cb,
               void                 *data);


/**
 *  @brief printing threshold for use with logML
 *
 *  @ingroup subopt_wuchty
 *
 */
extern double print_energy;

/**
 *  @brief Sort output by energy
 *
 *  @ingroup subopt_wuchty
 *
 */
extern int subopt_sorted;

/**
 *  @addtogroup dos
 *  @{
 */

/**
 *  @brief The Density of States
 *
 *  This array contains the density of states for an RNA sequences after a call to subopt_par(),
 *  subopt() or subopt_circ().
 *
 *  @pre  Call one of the functions subopt_par(), subopt() or subopt_circ() prior accessing the contents
 *        of this array
 *
 *  @see  subopt_par(), subopt(), subopt_circ()
 *
 */
extern int density_of_states[MAXDOS + 1];

/** @} */ /* End of group dos */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief Returns list of subopt structures or writes to fp
 *
 *  This function produces <b>all</b> suboptimal secondary structures within
 *  'delta' * 0.01 kcal/mol of the optimum. The results are either
 *  directly written to a 'fp' (if 'fp' is not NULL), or
 *  (fp==NULL) returned in a #SOLUTION * list terminated
 *  by an entry were the 'structure' pointer is NULL.
 *
 *  @ingroup subopt_wuchty
 *
 *  @param  seq
 *  @param  structure
 *  @param  delta
 *  @param  fp
 *  @return
 */
DEPRECATED(SOLUTION * subopt(char *seq, char *structure, int delta, FILE * fp),
           "Use vrna_subopt() or vrna_subopt_cb() instead");

/**
 *  @brief Returns list of subopt structures or writes to fp
 *
 *  @ingroup subopt_wuchty
 */
DEPRECATED(SOLUTION *
           subopt_par(char *seq, char *structure, vrna_param_t * parameters, int delta,
                      int is_constrained,
                      int is_circular, FILE * fp),
           "Use vrna_subopt() or vrna_subopt_cb() instead");

/**
 *  @brief Returns list of circular subopt structures or writes to fp
 *
 *  This function is similar to subopt() but calculates secondary structures
 *  assuming the RNA sequence to be circular instead of linear
 *
 *  @ingroup subopt_wuchty
 *
 *  @param  seq
 *  @param  sequence
 *  @param  delta
 *  @param  fp
 *  @return
 */
DEPRECATED(SOLUTION * subopt_circ(char *seq, char *sequence, int delta, FILE * fp),
           "Use vrna_subopt() or vrna_subopt_cb() instead");

/**
 *  @brief Compute Zuker type suboptimal structures
 *
 *  Compute Suboptimal structures according to M. Zuker, i.e. for every
 *  possible base pair the minimum energy structure containing the resp. base pair.
 *  Returns a list of these structures and their energies.
 *
 *  @ingroup subopt_zuker
 *
 *  @deprecated use vrna_zukersubopt() instead
 *
 *  @param  string  RNA sequence
 *  @return         List of zuker suboptimal structures
 */
DEPRECATED(SOLUTION * zukersubopt(const char *string),
           "Use vrna_subopt_zuker() instead");

/**
 *  @brief Compute Zuker type suboptimal structures
 *
 *  @ingroup subopt_zuker
 *
 *  @deprecated use vrna_zukersubopt() instead
 *
 */
DEPRECATED(SOLUTION * zukersubopt_par(const char *string, vrna_param_t * parameters),
           "Use vrna_subopt_zuker() instead");


#endif

#endif
