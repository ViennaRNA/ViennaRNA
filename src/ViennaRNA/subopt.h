/* subopt.h */
#ifndef VIENNA_RNA_PACKAGE_SUBOPT_H
#define VIENNA_RNA_PACKAGE_SUBOPT_H

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
 *  @file subopt.h
 *  @ingroup subopt_and_representatives
 *  @brief RNAsubopt and density of states declarations
 */

#define VRNA_BACKWARD_COMPAT

/**
 *  @brief Typename for the subopt solution list repesenting data structure #vrna_subopt_sol_s
 */
typedef struct vrna_subopt_sol_s   vrna_subopt_solution_t;

/**
 *  @brief  Callback for vrna_subopt_cb()
 *  @ingroup subopt_wuchty
 */
typedef void (vrna_subopt_callback)(const char *stucture, float energy, void *data);

#ifdef VRNA_BACKWARD_COMPAT

/**
 *  @brief  Backward compatibility typedef for #vrna_subopt_sol_s
 *  @deprecated Use #vrna_subopt_solution_t instead!
 */
typedef struct vrna_subopt_sol_s   SOLUTION;

#endif

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>


/**
 *  @brief  Solution element from subopt.c
 */
struct vrna_subopt_sol_s {
  float energy;       /**< @brief Free Energy of structure in kcal/mol */
  char *structure;    /**< @brief Structure in dot-bracket notation */
};

/**
 *  @brief Maximum density of states discretization for subopt
 */
#define MAXDOS                1000

/**
 *  @addtogroup subopt_wuchty
 *  @{
 *
 *  @}
 */

/**
 *  @brief Returns list of subopt structures or writes to fp
 * 
 *  This function produces <b>all</b> suboptimal secondary structures within
 *  'delta' * 0.01 kcal/mol of the optimum, see @cite wuchty:1999. The results
 *  are either directly written to a 'fp' (if 'fp' is not NULL), or
 *  (fp==NULL) returned in a #vrna_subopt_solution_t * list terminated
 *  by an entry were the 'structure' member is NULL.
 *
 *  @ingroup subopt_wuchty
 *
 *  @note This function requires all multibranch loop DP matrices for unique
 *        multibranch loop backtracing. Therefore, the supplied #vrna_fold_compound_t
 *        @p vc (argument 1) must be initialized with #vrna_md_t.uniq_ML = 1, for
 *        instance like this:
 *        @code
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.uniq_ML = 1;

  vrna_fold_compound_t *vc=vrna_fold_compound("GGGGGGAAAAAACCCCCC", &md, VRNA_OPTION_DEFAULT);
 *        @endcode
 *
 *  @see vrna_subopt_cb(), vrna_subopt_zuker()
 *  @param  vc
 *  @param  delta
 *  @param  sorted  Sort results by energy in ascending order
 *  @param  fp
 *  @return
 */
vrna_subopt_solution_t *
vrna_subopt(vrna_fold_compound_t *vc,
            int delta,
            int sorted,
            FILE *fp);

/**
 *  @brief  Generate suboptimal structures within an energy band arround the MFE
 *
 *  This is the most generic implementation of the suboptimal structure generator
 *  according to Wuchty et al. 1999 @cite wuchty:1999. Identical to vrna_subopt(), it computes all
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
 *        @p vc (argument 1) must be initialized with #vrna_md_t.uniq_ML = 1, for
 *        instance like this:
 *        @code
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.uniq_ML = 1;

  vrna_fold_compound_t *vc=vrna_fold_compound("GGGGGGAAAAAACCCCCC", &md, VRNA_OPTION_DEFAULT);
 *        @endcode
 *
 *  @see vrna_subopt_callback, vrna_subopt(), vrna_subopt_zuker()
 *  @param  vc      fold compount with the sequence data
 *  @param  delta   Energy band arround the MFE in 10cal/mol, i.e. deka-calories
 *  @param  cb      Pointer to a callback function that handles the backtracked structure and its free energy in kcal/mol
 *  @param  data    Pointer to some data structure that is passed along to the callback
 */
void
vrna_subopt_cb( vrna_fold_compound_t *vc,
                int delta,
                vrna_subopt_callback *cb,
                void *data);

/**
 *  @brief Compute Zuker type suboptimal structures
 *
 *  Compute Suboptimal structures according to M. Zuker @cite zuker:1989 , i.e. for every
 *  possible base pair the minimum energy structure containing the resp. base pair.
 *  Returns a list of these structures and their energies.
 *
 *  @note This function internally uses the cofold implementation to compute
 *        the suboptimal structures. For that purpose, the function doubles
 *        the sequence and enlarges the DP matrices, which in fact will grow
 *        by a factor of 4 during the computation!
 *        At the end of the structure prediction, everything will be re-set
 *        to its original requriements, i.e. normal sequence, normal (empty)
 *        DP matrices.
 *
 *  @bug  Due to resizing, any pre-existing constraints will be lost!
 *
 *  @ingroup subopt_zuker
 *
 *  @see vrna_subopt(), zukersubopt(), zukersubopt_par()
 *
 *  @param  vc  fold compound
 *  @return     List of zuker suboptimal structures
 */
vrna_subopt_solution_t *
vrna_subopt_zuker(vrna_fold_compound_t *vc);

/**
 *  @brief printing threshold for use with logML
 * 
 *  @ingroup subopt_wuchty
 *
 */
extern  double  print_energy;

/**
 *  @brief Sort output by energy
 * 
 *  @ingroup subopt_wuchty
 *
 */
extern  int     subopt_sorted;

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
 *  @see  subopt_par(), subopt(), subopt_circ()
 *
 */
extern  int     density_of_states[MAXDOS+1];

/** @} */ /* End of group dos */

#ifdef VRNA_BACKWARD_COMPAT

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
DEPRECATED(SOLUTION *subopt (char *seq, char *structure, int delta, FILE *fp));

/**
 *  @brief Returns list of subopt structures or writes to fp
 * 
 *  @ingroup subopt_wuchty
 */
DEPRECATED(SOLUTION *subopt_par(char *seq, char *structure, vrna_param_t *parameters, int delta, int is_constrained, int is_circular, FILE *fp));

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
DEPRECATED(SOLUTION *subopt_circ(char *seq, char *sequence, int delta, FILE *fp));

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
DEPRECATED(SOLUTION  *zukersubopt(const char *string));

/**
 *  @brief Compute Zuker type suboptimal structures
 *
 *  @ingroup subopt_zuker
 *
 *  @deprecated use vrna_zukersubopt() instead
 *
 */
DEPRECATED(SOLUTION  *zukersubopt_par(const char *string, vrna_param_t *parameters));


#endif

#endif
