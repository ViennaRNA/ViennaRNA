#ifndef VIENNA_RNA_PACKAGE_PART_FUNC_CO_H
#define VIENNA_RNA_PACKAGE_PART_FUNC_CO_H

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
 *  @file     part_func_co.h
 *  @ingroup  part_func_global_deprecated
 *  @brief    Partition function for two RNA sequences
 */

/**
 *  @addtogroup pf_cofold
 *  @{
 */

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/partfunc/global.h>
#include <ViennaRNA/probabilities/structures.h>
#include <ViennaRNA/probabilities/basepairs.h>
#include <ViennaRNA/concentrations.h>
#include <ViennaRNA/structures/problist.h>

/**
 *  @brief Toggles no intrabp in 2nd mol
 */
extern int    mirnatog;

/**
 *  @brief Free energies of the two monomers
 */
extern double F_monomer[2];

/**
 *  @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 #################################################
 # DEPRECATED FUNCTIONS                          #
 #################################################
 */

/**
 *  @brief Calculate partition function and base pair probabilities
 *
 *  This is the cofold partition function folding. The second molecule starts
 *  at the #cut_point nucleotide.
 *
 *  @note OpenMP: Since this function relies on the global parameters
 *        #do_backtrack, #dangles, #temperature and #pf_scale it is not
 *        threadsafe according to concurrent changes in these variables!
 *        Use co_pf_fold_par() instead to circumvent this issue.
 *
 *  @deprecated{Use vrna_pf_dimer() instead!}
 *
 *  @ingroup part_func_global_deprecated
 *
 *  @param  sequence  Concatenated RNA sequences
 *  @param  structure Will hold the structure or constraints
 *  @return           vrna_dimer_pf_t structure containing a set of energies needed for
 *                    concentration computations.
 */
DEPRECATED(vrna_dimer_pf_t co_pf_fold(char  *sequence,
                                      char  *structure),
"Use vrna_pf_co_fold() or vrna_pf_dimer() instead");

/**
 *  @brief Calculate partition function and base pair probabilities
 *
 *  This is the cofold partition function folding. The second molecule starts
 *  at the #cut_point nucleotide.
 *
 *  @deprecated Use vrna_pf_dimer() instead!
 *
 *  @see get_boltzmann_factors(), co_pf_fold()
 *
 *  @ingroup part_func_global_deprecated
 *
 *  @param sequence       Concatenated RNA sequences
 *  @param structure      Pointer to the structure constraint
 *  @param parameters     Data structure containing the precalculated Boltzmann factors
 *  @param calculate_bppm Switch to turn Base pair probability calculations on/off (0==off)
 *  @param is_constrained Switch to indicate that a structure contraint is passed via the
 *                        structure argument (0==off)
 *  @return               vrna_dimer_pf_t structure containing a set of energies needed for
 *                        concentration computations.
 */
DEPRECATED(vrna_dimer_pf_t co_pf_fold_par(char              *sequence,
                                          char              *structure,
                                          vrna_exp_param_t  *parameters,
                                          int               calculate_bppm,
                                          int               is_constrained),
"Use the new API and vrna_pf_dimer() instead");

/**
 *  DO NOT USE THIS FUNCTION ANYMORE
 *  @deprecated{ This function is deprecated and will be removed soon!}
 *  use assign_plist_from_pr() instead!
 */
DEPRECATED(vrna_ep_t *get_plist(vrna_ep_t *pl,
                                int       length,
                                double    cut_off),
"Use vrna_plist() and vrna_plist_from_probs() instead");

/**
 *  @brief Compute Boltzmann probabilities of dimerization without homodimers
 *
 *  Given the pair probabilities and free energies (in the null model) for a
 *  dimer AB and the two constituent monomers A and B, compute the conditional pair
 *  probabilities given that a dimer AB actually forms.
 *  Null model pair probabilities are given as a list as produced by
 *  assign_plist_from_pr(), the dimer probabilities 'prAB' are modified in place.
 *
 *  @deprecated{ Use vrna_pf_dimer_probs() instead!}
 *
 *  @ingroup part_func_global_deprecated
 *
 *  @param FAB      free energy of dimer AB
 *  @param FEA      free energy of monomer A
 *  @param FEB      free energy of monomer B
 *  @param prAB     pair probabilities for dimer
 *  @param prA      pair probabilities monomer
 *  @param prB      pair probabilities monomer
 *  @param Alength  Length of molecule A
 */
DEPRECATED(void compute_probabilities(double    FAB,
                                      double    FEA,
                                      double    FEB,
                                      vrna_ep_t *prAB,
                                      vrna_ep_t *prA,
                                      vrna_ep_t *prB,
                                      int       Alength),
"Use vrna_pf_dimer_probs() instead");

/**
 *  DO NOT USE THIS FUNCTION ANYMORE
 *  @deprecated{ This function is deprecated and will be removed soon!}
 *  @ingroup part_func_global_deprecated
 */
DEPRECATED(void   init_co_pf_fold(int length),
"This function is obsolete");

/**
 *  @brief Get a pointer to the base pair probability array
 *
 *  Accessing the base pair probabilities for a pair (i,j) is achieved by
 *  @verbatim FLT_OR_DBL *pr = export_bppm(); pr_ij = pr[iindx[i]-j]; @endverbatim
 *
 *  @deprecated This function is deprecated and will be removed soon! The base pair
 *              probability array is available through the #vrna_fold_compound_t data
 *              structure, and its associated #vrna_mx_pf_t member.
 *
 *  @ingroup part_func_global_deprecated
 *
 *  @see vrna_idx_row_wise()
 *
 *  @return A pointer to the base pair probability array
 */
DEPRECATED(FLT_OR_DBL *export_co_bppm(void),
"Use the new API with vrna_fold_compound_t instead");

/**
 *  @brief Free the memory occupied by co_pf_fold()
 *
 *  @deprecated This function will be removed for the new API soon!
 *              See vrna_pf_dimer(), vrna_fold_compound(), and
 *              vrna_fold_compound_free() for an alternative
 *  @ingroup part_func_global_deprecated
 */
DEPRECATED(void free_co_pf_arrays(void),
"This function is obsolete");

/**
 *  @brief Recalculate energy parameters
 *
 *  This function recalculates all energy parameters given
 *  the current model settings.
 *
 *  @deprecated   Use vrna_exp_params_subst() instead!
 *
 *  @ingroup part_func_global_deprecated
 *
 *  @param    length      Length of the current RNA sequence
 */
DEPRECATED(void update_co_pf_params(int length),
"This function is obsolete");

/**
 *  @brief Recalculate energy parameters
 *
 *  This function recalculates all energy parameters given
 *  the current model settings.
 *  It's second argument can either be NULL or a data structure
 *  containing the precomputed Boltzmann factors. In the first
 *  scenario, the necessary data structure will be created automatically
 *  according to the current global model settings, i.e. this
 *  mode might not be threadsafe.
 *  However, if the provided data structure is not NULL, threadsafety
 *  for the model parameters #dangles, #pf_scale and #temperature is regained, since their
 *  values are taken from this data structure during subsequent calculations.
 *
 *  @deprecated   Use vrna_exp_params_subst() instead!
 *
 *  @ingroup part_func_global_deprecated
 *
 *  @param    length      Length of the current RNA sequence
 *  @param    parameters  data structure containing the precomputed Boltzmann factors
 */
DEPRECATED(void update_co_pf_params_par(int               length,
                                        vrna_exp_param_t  *parameters),
"Use the new API with vrna_fold_compound_t instead");

#endif

#endif
