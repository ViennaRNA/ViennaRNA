#ifndef VIENNA_RNA_PACKAGE_PART_FUNC_CO_H
#define VIENNA_RNA_PACKAGE_PART_FUNC_CO_H

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  \addtogroup pf_cofold
 *  \brief Partition Function Cofolding
 *
 *  To simplify the implementation the partition function computation is done
 *  internally in a null model that does not include the duplex initiation
 *  energy, i.e. the entropic penalty for producing a dimer from two
 *  monomers). The resulting free energies and pair probabilities are initially
 *  relative to that null model. In a second step the free energies can be
 *  corrected to include the dimerization penalty, and the pair probabilities
 *  can be divided into the conditional pair probabilities given that a re
 *  dimer is formed or not formed. See \cite bernhart:2006 for further details.
 *  @{
 *  \file part_func_co.h
 * 
 *  \brief Partition function for two RNA sequences
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
 */

/**
 *  \brief Toggles no intrabp in 2nd mol
 */
extern int    mirnatog;

/**
 *  \brief Free energies of the two monomers
 */
extern double F_monomer[2];

/**
 *  \brief  Calculate partition function and base pair probabilities of
 *          nucleic acid/nucleic acid dimers
 *
 *  This is the cofold partition function folding.
 *
 *  \see    vrna_fold_compound() for how to retrieve the necessary data structure
 *
 *  \param  vc        the fold compound data structure
 *  \param  structure Will hold the structure or constraints
 *  \return           cofoldF structure containing a set of energies needed for
 *                    concentration computations.
 */
cofoldF vrna_pf_dimer(vrna_fold_compound_t *vc,
                        char *structure);

/**
 *  \brief Compute Boltzmann probabilities of dimerization without homodimers
 * 
 *  Given the pair probabilities and free energies (in the null model) for a
 *  dimer AB and the two constituent monomers A and B, compute the conditional pair
 *  probabilities given that a dimer AB actually forms.
 *  Null model pair probabilities are given as a list as produced by
 *  assign_plist_from_pr(), the dimer probabilities 'prAB' are modified in place.
 * 
 *  \param FAB        free energy of dimer AB
 *  \param FEA        free energy of monomer A
 *  \param FEB        free energy of monomer B
 *  \param prAB       pair probabilities for dimer
 *  \param prA        pair probabilities monomer
 *  \param prB        pair probabilities monomer
 *  \param Alength    Length of molecule A
 *  \param exp_params The precomputed Boltzmann factors
 */
void  vrna_co_pf_dimer_probs( double FAB,
                              double FA,
                              double FB,
                              struct plist *prAB,
                              const plist *prA,
                              const plist *prB,
                              int Alength,
                              const vrna_exp_param_t *exp_params);

/**
 *  \brief Given two start monomer concentrations a and b, compute the
 *  concentrations in thermodynamic equilibrium of all dimers and the monomers.
 * 
 *  This function takes an array  'startconc' of input concentrations with alternating
 *  entries for the initial concentrations of molecules A and B (terminated by
 *  two zeroes), then computes the resulting equilibrium concentrations
 *  from the free energies for the dimers. Dimer free energies should be the
 *  dimer-only free energies, i.e. the FcAB entries from the #cofoldF struct.
 * 
 *  \param FEAB       Free energy of AB dimer (FcAB entry)
 *  \param FEAA       Free energy of AA dimer (FcAB entry)
 *  \param FEBB       Free energy of BB dimer (FcAB entry)
 *  \param FEA        Free energy of monomer A
 *  \param FEB        Free energy of monomer B
 *  \param startconc  List of start concentrations [a0],[b0],[a1],[b1],...,[an][bn],[0],[0]
 *  \param exp_params The precomputed Boltzmann factors
 *  \return ConcEnt array containing the equilibrium energies and start concentrations
 */
ConcEnt *vrna_co_pf_get_concentrations( double FcAB,
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
 *  \brief Calculate partition function and base pair probabilities
 *
 *  This is the cofold partition function folding. The second molecule starts
 *  at the #cut_point nucleotide.
 *
 *  \note OpenMP: Since this function relies on the global parameters
 *        #do_backtrack, #dangles, #temperature and #pf_scale it is not
 *        threadsafe according to concurrent changes in these variables!
 *        Use co_pf_fold_par() instead to circumvent this issue.
 *
 *  \deprecated{Use vrna_pf_dimer() instead!}
 *
 *  \param  sequence  Concatenated RNA sequences
 *  \param  structure Will hold the structure or constraints
 *  \return           cofoldF structure containing a set of energies needed for
 *                    concentration computations.
 */
DEPRECATED(cofoldF co_pf_fold( char *sequence, char *structure));

/**
 *  \brief Calculate partition function and base pair probabilities
 *
 *  This is the cofold partition function folding. The second molecule starts
 *  at the #cut_point nucleotide.
 *
 *  \see get_boltzmann_factors(), co_pf_fold()
 *
 *  \param sequence       Concatenated RNA sequences
 *  \param structure      Pointer to the structure constraint
 *  \param parameters     Data structure containing the precalculated Boltzmann factors
 *  \param calculate_bppm Switch to turn Base pair probability calculations on/off (0==off)
 *  \param is_constrained Switch to indicate that a structure contraint is passed via the
 *                        structure argument (0==off)
 *  \return               cofoldF structure containing a set of energies needed for
 *                        concentration computations.
 */
DEPRECATED(cofoldF co_pf_fold_par(char *sequence, char *structure, vrna_exp_param_t *parameters, int calculate_bppm, int is_constrained));

/**
 *  DO NOT USE THIS FUNCTION ANYMORE
 *  \deprecated{ This function is deprecated and will be removed soon!}
 *  use \ref assign_plist_from_pr() instead!
 */
DEPRECATED(plist  *get_plist( struct plist *pl,
                              int length,
                              double cut_off));

/**
 *  \brief Compute Boltzmann probabilities of dimerization without homodimers
 * 
 *  Given the pair probabilities and free energies (in the null model) for a
 *  dimer AB and the two constituent monomers A and B, compute the conditional pair
 *  probabilities given that a dimer AB actually forms.
 *  Null model pair probabilities are given as a list as produced by
 *  assign_plist_from_pr(), the dimer probabilities 'prAB' are modified in place.
 *
 *  \deprecated{ Use vrna_co_pf_dimer_probs() instead!}
 * 
 *  \param FAB      free energy of dimer AB
 *  \param FEA      free energy of monomer A
 *  \param FEB      free energy of monomer B
 *  \param prAB     pair probabilities for dimer
 *  \param prA      pair probabilities monomer
 *  \param prB      pair probabilities monomer
 *  \param Alength  Length of molecule A
 */
DEPRECATED(void compute_probabilities(double FAB, double FEA, double FEB, struct plist  *prAB, struct plist  *prA, struct plist  *prB, int Alength));

/**
 *  \brief Given two start monomer concentrations a and b, compute the
 *  concentrations in thermodynamic equilibrium of all dimers and the monomers.
 * 
 *  This function takes an array  'startconc' of input concentrations with alternating
 *  entries for the initial concentrations of molecules A and B (terminated by
 *  two zeroes), then computes the resulting equilibrium concentrations
 *  from the free energies for the dimers. Dimer free energies should be the
 *  dimer-only free energies, i.e. the FcAB entries from the #cofoldF struct.
 *
 *  \deprecated{ Use vrna_co_pf_get_concentrations() instead!}
 * 
 *  \param FEAB       Free energy of AB dimer (FcAB entry)
 *  \param FEAA       Free energy of AA dimer (FcAB entry)
 *  \param FEBB       Free energy of BB dimer (FcAB entry)
 *  \param FEA        Free energy of monomer A
 *  \param FEB        Free energy of monomer B
 *  \param startconc  List of start concentrations [a0],[b0],[a1],[b1],...,[an][bn],[0],[0]
 *  \return ConcEnt array containing the equilibrium energies and start concentrations
 */
DEPRECATED(ConcEnt *get_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconc));

/**
 *  DO NOT USE THIS FUNCTION ANYMORE
 *  \deprecated{ This function is deprecated and will be removed soon!}
 */
DEPRECATED(void   init_co_pf_fold(int length));

/**
 *  \brief Get a pointer to the base pair probability array
 * 
 *  Accessing the base pair probabilities for a pair (i,j) is achieved by
 *  \verbatim FLT_OR_DBL *pr = export_bppm(); pr_ij = pr[iindx[i]-j]; \endverbatim
 * 
 *  \see vrna_get_iindx()
 *  \return A pointer to the base pair probability array
 */
DEPRECATED(FLT_OR_DBL *export_co_bppm(void));

/**
 *  \brief Free the memory occupied by co_pf_fold()
 */
DEPRECATED(void free_co_pf_arrays(void));

/**
 *  \brief Recalculate energy parameters
 *
 *  This function recalculates all energy parameters given
 *  the current model settings.
 *
 *  @deprecated   Use vrna_exp_params_update() instead!
 *
 *  \param    length      Length of the current RNA sequence
 */
DEPRECATED(void update_co_pf_params(int length));

/**
 *  \brief Recalculate energy parameters
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
 *  @deprecated   Use vrna_exp_params_update() instead!
 *
 *  \param    length      Length of the current RNA sequence
 *  \param    parameters  data structure containing the precomputed Boltzmann factors
 */
DEPRECATED(void update_co_pf_params_par(int length, vrna_exp_param_t *parameters));

#endif
