/* subopt/basic.h */
#ifndef VIENNA_RNA_PACKAGE_SUBOPT_H
#define VIENNA_RNA_PACKAGE_SUBOPT_H

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
 *  @file subopt/basic.h
 *  @ingroup subopt_and_representatives
 *  @brief Basic RNAsubopt and density of states declarations
 */

/**
 *  @brief Typename for the subopt solution list repesenting data structure #vrna_subopt_sol_s
 */
typedef struct vrna_subopt_sol_s vrna_subopt_solution_t;

/**
 *  @brief  Callback for vrna_subopt_cb()
 *  @ingroup subopt_wuchty
 *
 * @callback
 * @parblock
 * This function will be called for each suboptimal secondary structure that is successfully backtraced.
 * @endparblock
 *
 * @see vrna_subopt_cb()
 *
 * @param structure The suboptimal secondary structure in dot-bracket notation
 * @param energy    The free energy of the secondary structure in kcal/mol
 * @param data      Some arbitrary, auxiliary data address as passed to vrna_subopt_cb()
 */
typedef void (*vrna_subopt_result_f)(const char  *stucture,
                                    float       energy,
                                    void        *data);

DEPRECATED(typedef void (vrna_subopt_callback)(const char  *stucture,
                                    float       energy,
                                    void        *data),
          "Use vrna_subopt_result_f instead!");

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief  Backward compatibility typedef for #vrna_subopt_sol_s
 *  @deprecated Use #vrna_subopt_solution_t instead!
 */
typedef struct vrna_subopt_sol_s SOLUTION;

#endif

/**
 *  @brief  Solution element from subopt.c
 */
struct vrna_subopt_sol_s {
  float energy;       /**< @brief Free Energy of structure in kcal/mol */
  char  *structure;   /**< @brief Structure in dot-bracket notation */
};

#endif
