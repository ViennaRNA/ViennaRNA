#ifndef VIENNA_RNA_PACKAGE_PART_FUNC_CO_H
#define VIENNA_RNA_PACKAGE_PART_FUNC_CO_H

#ifdef VRNA_WARN_DEPRECATED
# ifdef __GNUC__
#  define DEPRECATED(func) func __attribute__ ((deprecated))
# else
#  define DEPRECATED(func) func
# endif
#else
# define DEPRECATED(func) func
#endif

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

/**
 *  @file     part_func_co.h
 *  @ingroup  pf_fold cofold pf_cofold
 *  @brief    Partition function for two RNA sequences
 */

/**
 *  @addtogroup pf_cofold
 *  @brief Partition Function Cofolding
 *
 *  To simplify the implementation the partition function computation is done
 *  internally in a null model that does not include the duplex initiation
 *  energy, i.e. the entropic penalty for producing a dimer from two
 *  monomers). The resulting free energies and pair probabilities are initially
 *  relative to that null model. In a second step the free energies can be
 *  corrected to include the dimerization penalty, and the pair probabilities
 *  can be divided into the conditional pair probabilities given that a re
 *  dimer is formed or not formed. See @cite bernhart:2006 for further details.
 *
 *  As for folding one RNA molecule, this computes the partition function
 *  of all possible structures and the base pair probabilities. Uses the
 *  same global #pf_scale variable to avoid overflows.
 *
 *  To simplify the implementation the partition function computation is done
 *  internally in a null model that does not include the duplex initiation
 *  energy, i.e. the entropic penalty for producing a dimer from two
 *  monomers). The resulting free energies and pair probabilities are initially
 *  relative to that null model. In a second step the free energies can be
 *  corrected to include the dimerization penalty, and the pair probabilities
 *  can be divided into the conditional pair probabilities given that a re
 *  dimer is formed or not formed.
 *
 *  After computing the partition functions of all possible dimeres one
 *  can compute the probabilities of base pairs, the concentrations out of
 *  start concentrations and sofar and soaway.
 *
 *  Dimer formation is inherently concentration dependent. Given the free
 *  energies of the monomers A and B and dimers AB, AA, and BB one can compute
 *  the equilibrium concentrations, given input concentrations of A and B, see
 *  e.g. Dimitrov & Zuker (2004)
 *
 *  @{
 *  @ingroup  pf_cofold
 */

/** @brief Typename for the data structure that stores the dimer partition functions, #vrna_dimer_pf_s, as returned by vrna_pf_dimer() */
typedef struct vrna_dimer_pf_s  vrna_dimer_pf_t;

/** @brief Typename for the data structure that stores the dimer concentrations, #vrna_dimer_conc_s, as required by vrna_pf_dimer_concentration() */
typedef struct vrna_dimer_conc_s  vrna_dimer_conc_t;


#ifdef VRNA_BACKWARD_COMPAT

/**
 *  @brief Backward compatibility typedef for #vrna_dimer_pf_s
 */
typedef struct vrna_dimer_pf_s    cofoldF;

/**
 *  @brief Backward compatibility typedef for #vrna_dimer_conc_s
 */
typedef struct vrna_dimer_conc_s  ConcEnt;

#endif

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>

/**
 *  @brief Toggles no intrabp in 2nd mol
 */
extern int    mirnatog;

/**
 *  @brief Free energies of the two monomers
 */
extern double F_monomer[2];

/**
 *  @brief  Data structure returned by vrna_pf_dimer()
 */
struct vrna_dimer_pf_s {
  /* free energies for: */
  double F0AB;  /**< @brief Null model without DuplexInit */
  double FAB;   /**< @brief all states with DuplexInit correction */
  double FcAB;  /**< @brief true hybrid states only */
  double FA;    /**< @brief monomer A */
  double FB;    /**< @brief monomer B */
};

/**
 *  @brief  Data structure for concentration dependency computations
 */
struct vrna_dimer_conc_s {
  double Ac_start;    /**< @brief start concentration A */
  double Bc_start;    /**< @brief start concentration B */
  double ABc;         /**< @brief End concentration AB */
  double AAc;
  double BBc;
  double Ac;
  double Bc;
};

/**
 *  @brief  Calculate partition function and base pair probabilities of
 *          nucleic acid/nucleic acid dimers
 *
 *  This is the cofold partition function folding.
 *
 *  @see    vrna_fold_compound() for how to retrieve the necessary data structure
 *
 *  @param  vc        the fold compound data structure
 *  @param  structure Will hold the structure or constraints
 *  @return           vrna_dimer_pf_t structure containing a set of energies needed for
 *                    concentration computations.
 */
vrna_dimer_pf_t
vrna_pf_dimer(vrna_fold_compound_t *vc,
              char *structure);

/**
 *  @brief Compute Boltzmann probabilities of dimerization without homodimers
 *
 *  Given the pair probabilities and free energies (in the null model) for a
 *  dimer AB and the two constituent monomers A and B, compute the conditional pair
 *  probabilities given that a dimer AB actually forms.
 *  Null model pair probabilities are given as a list as produced by
 *  vrna_plist_from_probs(), the dimer probabilities 'prAB' are modified in place.
 *
 *  @param FAB        free energy of dimer AB
 *  @param FA         free energy of monomer A
 *  @param FB         free energy of monomer B
 *  @param prAB       pair probabilities for dimer
 *  @param prA        pair probabilities monomer
 *  @param prB        pair probabilities monomer
 *  @param Alength    Length of molecule A
 *  @param exp_params The precomputed Boltzmann factors
 */
void  vrna_pf_dimer_probs(double FAB,
                          double FA,
                          double FB,
                          vrna_plist_t *prAB,
                          const vrna_plist_t *prA,
                          const vrna_plist_t *prB,
                          int Alength,
                          const vrna_exp_param_t *exp_params);

/**
 *  @brief Given two start monomer concentrations a and b, compute the
 *  concentrations in thermodynamic equilibrium of all dimers and the monomers.
 *
 *  This function takes an array  'startconc' of input concentrations with alternating
 *  entries for the initial concentrations of molecules A and B (terminated by
 *  two zeroes), then computes the resulting equilibrium concentrations
 *  from the free energies for the dimers. Dimer free energies should be the
 *  dimer-only free energies, i.e. the FcAB entries from the #vrna_dimer_pf_t struct.
 *
 *  @param FcAB       Free energy of AB dimer (FcAB entry)
 *  @param FcAA       Free energy of AA dimer (FcAB entry)
 *  @param FcBB       Free energy of BB dimer (FcAB entry)
 *  @param FEA        Free energy of monomer A
 *  @param FEB        Free energy of monomer B
 *  @param startconc  List of start concentrations [a0],[b0],[a1],[b1],...,[an][bn],[0],[0]
 *  @param exp_params The precomputed Boltzmann factors
 *  @return vrna_dimer_conc_t array containing the equilibrium energies and start concentrations
 */
vrna_dimer_conc_t *vrna_pf_dimer_concentrations(double FcAB,
                                      double FcAA,
                                      double FcBB,
                                      double FEA,
                                      double FEB,
                                      const double *startconc,
                                      const vrna_exp_param_t *exp_params);

/**
 *  @}
 */

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
 *  @param  sequence  Concatenated RNA sequences
 *  @param  structure Will hold the structure or constraints
 *  @return           vrna_dimer_pf_t structure containing a set of energies needed for
 *                    concentration computations.
 */
DEPRECATED(vrna_dimer_pf_t co_pf_fold( char *sequence, char *structure));

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
 *  @param sequence       Concatenated RNA sequences
 *  @param structure      Pointer to the structure constraint
 *  @param parameters     Data structure containing the precalculated Boltzmann factors
 *  @param calculate_bppm Switch to turn Base pair probability calculations on/off (0==off)
 *  @param is_constrained Switch to indicate that a structure contraint is passed via the
 *                        structure argument (0==off)
 *  @return               vrna_dimer_pf_t structure containing a set of energies needed for
 *                        concentration computations.
 */
DEPRECATED(vrna_dimer_pf_t co_pf_fold_par(char *sequence, char *structure, vrna_exp_param_t *parameters, int calculate_bppm, int is_constrained));

/**
 *  DO NOT USE THIS FUNCTION ANYMORE
 *  @deprecated{ This function is deprecated and will be removed soon!}
 *  use assign_plist_from_pr() instead!
 */
DEPRECATED(vrna_plist_t  *get_plist( vrna_plist_t *pl,
                              int length,
                              double cut_off));

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
 *  @param FAB      free energy of dimer AB
 *  @param FEA      free energy of monomer A
 *  @param FEB      free energy of monomer B
 *  @param prAB     pair probabilities for dimer
 *  @param prA      pair probabilities monomer
 *  @param prB      pair probabilities monomer
 *  @param Alength  Length of molecule A
 */
DEPRECATED(void compute_probabilities(double FAB, double FEA, double FEB, vrna_plist_t  *prAB, vrna_plist_t  *prA, vrna_plist_t  *prB, int Alength));

/**
 *  @brief Given two start monomer concentrations a and b, compute the
 *  concentrations in thermodynamic equilibrium of all dimers and the monomers.
 *
 *  This function takes an array  'startconc' of input concentrations with alternating
 *  entries for the initial concentrations of molecules A and B (terminated by
 *  two zeroes), then computes the resulting equilibrium concentrations
 *  from the free energies for the dimers. Dimer free energies should be the
 *  dimer-only free energies, i.e. the FcAB entries from the #vrna_dimer_pf_t struct.
 *
 *  @deprecated{ Use vrna_pf_dimer_concentrations() instead!}
 *
 *  @param FEAB       Free energy of AB dimer (FcAB entry)
 *  @param FEAA       Free energy of AA dimer (FcAB entry)
 *  @param FEBB       Free energy of BB dimer (FcAB entry)
 *  @param FEA        Free energy of monomer A
 *  @param FEB        Free energy of monomer B
 *  @param startconc  List of start concentrations [a0],[b0],[a1],[b1],...,[an][bn],[0],[0]
 *  @return vrna_dimer_conc_t array containing the equilibrium energies and start concentrations
 */
DEPRECATED(vrna_dimer_conc_t *get_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconc));

/**
 *  DO NOT USE THIS FUNCTION ANYMORE
 *  @deprecated{ This function is deprecated and will be removed soon!}
 */
DEPRECATED(void   init_co_pf_fold(int length));

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
 *  @see vrna_idx_row_wise()
 *  @return A pointer to the base pair probability array
 */
DEPRECATED(FLT_OR_DBL *export_co_bppm(void));

/**
 *  @brief Free the memory occupied by co_pf_fold()
 *
 *  @deprecated This function will be removed for the new API soon!
 *              See vrna_pf_dimer(), vrna_fold_compound(), and
 *              vrna_fold_compound_free() for an alternative
 */
DEPRECATED(void free_co_pf_arrays(void));

/**
 *  @brief Recalculate energy parameters
 *
 *  This function recalculates all energy parameters given
 *  the current model settings.
 *
 *  @deprecated   Use vrna_exp_params_subst() instead!
 *
 *  @param    length      Length of the current RNA sequence
 */
DEPRECATED(void update_co_pf_params(int length));

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
 *  @param    length      Length of the current RNA sequence
 *  @param    parameters  data structure containing the precomputed Boltzmann factors
 */
DEPRECATED(void update_co_pf_params_par(int length, vrna_exp_param_t *parameters));

#endif
