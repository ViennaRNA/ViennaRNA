#ifndef VIENNA_RNA_PACKAGE_PART_FUNC_H
#define VIENNA_RNA_PACKAGE_PART_FUNC_H

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/centroid.h>

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif



/**
 *  \addtogroup pf_fold
 *  \brief This section provides information about all functions and variables related to
 *  the calculation of the partition function and base pair probabilities.
 *
 *  Instead of the minimum free energy structure the partition function of all possible structures
 *  and from that the pairing probability for every possible pair can be calculated, using a dynamic
 *  programming algorithm as described in \cite mccaskill:1990. 
 *
 *  @{
 *    \file part_func.h
 *    \brief Partition function of single RNA sequences
 * 
 *    This file includes (almost) all function declarations within the <b>RNAlib</b> that are related to
 *    Partion function folding...
 *  @}
 */

/**
 *  \brief Flag indicating that auxilary arrays are needed throughout the computations. This is essential for stochastic backtracking
 *
 *  Set this variable to 1 prior to a call of pf_fold() to ensure that all matrices needed for stochastic backtracking
 *  are filled in the forward recursions
 *
 *  \ingroup subopt_stochbt
 *
 *  \see pbacktrack(), pbacktrack_circ
 */
extern  int st_back;

/*
#################################################
# PARTITION FUNCTION COMPUTATION                #
#################################################
*/

/**
 *  \brief Compute the partition function \f$Q\f$ for a given RNA sequence
 *
 *  If \a structure is not a NULL pointer on input, it contains on
 *  return a string consisting of the letters " . , | { } ( ) " denoting
 *  bases that are essentially unpaired, weakly paired, strongly paired without
 *  preference, weakly upstream (downstream) paired, or strongly up-
 *  (down-)stream paired bases, respectively.
 *  If #fold_constrained is not 0, the \a structure string is
 *  interpreted on input as a list of constraints for the folding. The
 *  character "x" marks bases that must be unpaired, matching brackets " ( ) "
 *  denote base pairs, all other characters are ignored. Any pairs
 *  conflicting with the constraint will be forbidden. This is usually sufficient
 *  to ensure the constraints are honored.
 *  If the parameter calculate_bppm is set to 0 base pairing probabilities will not
 *  be computed (saving CPU time), otherwise after calculations took place #pr will
 *  contain the probability that bases \a i and \a j pair.
 * 
 *  \ingroup pf_fold
 *
 *  \note           The global array #pr is deprecated and the user who wants the calculated
 *                  base pair probabilities for further computations is advised to use the function
 *                  export_bppm()
 *  \see            vnra_get_fold_compound(), vrna_db_get_from_pr(), export_bppm(), get_boltzmann_factors(), free_pf_arrays()
 *  \param[in,out]  vc              The fold compound data structure
 *  \param[in,out]  structure       A pointer to a char array where a base pair probability information can be stored in a
 *                                  pseudo-dot-bracket notation (may be NULL, too)
 *  \return         The Gibbs free energy of the ensemble (\f$G = -RT \cdot \log(Q) \f$) in kcal/mol
 */
float vrna_pf_fold( vrna_fold_compound *vc,
                    char *structure);

/**
 *  \brief Compute the partition function \f$Q\f$ for a given RNA sequence
 *
 *  If \a structure is not a NULL pointer on input, it contains on
 *  return a string consisting of the letters " . , | { } ( ) " denoting
 *  bases that are essentially unpaired, weakly paired, strongly paired without
 *  preference, weakly upstream (downstream) paired, or strongly up-
 *  (down-)stream paired bases, respectively.
 *  If #fold_constrained is not 0, the \a structure string is
 *  interpreted on input as a list of constraints for the folding. The
 *  character "x" marks bases that must be unpaired, matching brackets " ( ) "
 *  denote base pairs, all other characters are ignored. Any pairs
 *  conflicting with the constraint will be forbidden. This is usually sufficient
 *  to ensure the constraints are honored.
 *  If the parameter calculate_bppm is set to 0 base pairing probabilities will not
 *  be computed (saving CPU time), otherwise after calculations took place #pr will
 *  contain the probability that bases \a i and \a j pair.
 * 
 *  \ingroup pf_fold
 *
 *  \deprecated Use vrna_pf_fold() instead
 *
 *  \note           The global array #pr is deprecated and the user who wants the calculated
 *                  base pair probabilities for further computations is advised to use the function
 *                  export_bppm()
 *  \post           After successful run the hidden folding matrices are filled with the appropriate Boltzmann factors.
 *                  Depending on whether the global variable #do_backtrack was set the base pair probabilities are already
 *                  computed and may be accessed for further usage via the export_bppm() function.
 *                  A call of free_pf_arrays() will free all memory allocated by this function.
 *                  Successive calls will first free previously allocated memory before starting the computation.
 *  \see            vrna_pf_fold(), bppm_to_structure(), export_bppm(), get_boltzmann_factors(), free_pf_arrays()
 *  \param[in]      sequence        The RNA sequence input
 *  \param[in,out]  structure       A pointer to a char array where a base pair probability information can be stored in a
 *                                  pseudo-dot-bracket notation (may be NULL, too)
 *  \param[in]      parameters      Data structure containing the precalculated Boltzmann factors
 *  \param[in]      calculate_bppm  Switch to Base pair probability calculations on/off (0==off)
 *  \param[in]      is_constrained  Switch to indicate that a structure contraint is passed via the structure argument (0==off)
 *  \param[in]      is_circular     Switch to (de-)activate postprocessing steps in case RNA sequence is circular (0==off)
 *  \return         The Gibbs free energy of the ensemble (\f$G = -RT \cdot \log(Q) \f$) in kcal/mol
 */
DEPRECATED(float   pf_fold_par(  const char *sequence,
                      char *structure,
                      vrna_exp_param_t *parameters,
                      int calculate_bppm,
                      int is_constrained,
                      int is_circular));

/**
 *  \brief Compute the partition function \f$Q\f$ of an RNA sequence
 * 
 *  If \a structure is not a NULL pointer on input, it contains on
 *  return a string consisting of the letters " . , | { } ( ) " denoting
 *  bases that are essentially unpaired, weakly paired, strongly paired without
 *  preference, weakly upstream (downstream) paired, or strongly up-
 *  (down-)stream paired bases, respectively.
 *  If #fold_constrained is not 0, the \a structure string is
 *  interpreted on input as a list of constraints for the folding. The
 *  character "x" marks bases that must be unpaired, matching brackets " ( ) "
 *  denote base pairs, all other characters are ignored. Any pairs
 *  conflicting with the constraint will be forbidden. This is usually sufficient
 *  to ensure the constraints are honored.
 *  If #do_backtrack has been set to 0 base pairing probabilities will not
 *  be computed (saving CPU time), otherwise #pr will contain the probability
 *  that bases \a i and \a j pair.
 * 
 *  \ingroup pf_fold
 *
 *  \note   The global array #pr is deprecated and the user who wants the calculated
 *          base pair probabilities for further computations is advised to use the function
 *          export_bppm().
 *  \note   \b OpenMP:
 *          This function is not entirely threadsafe. While the recursions are working on their
 *          own copies of data the model details for the recursions are determined from the global
 *          settings just before entering the recursions. Consider using pf_fold_par() for a
 *          really threadsafe implementation.
 *  \pre    This function takes its model details from the global variables provided in \e RNAlib
 *  \post   After successful run the hidden folding matrices are filled with the appropriate Boltzmann factors.
 *          Depending on whether the global variable #do_backtrack was set the base pair probabilities are already
 *          computed and may be accessed for further usage via the export_bppm() function.
 *          A call of free_pf_arrays() will free all memory allocated by this function.
 *          Successive calls will first free previously allocated memory before starting the computation.
 *  \see    pf_fold_par(), pf_circ_fold(), bppm_to_structure(), export_bppm()
 *  \param sequence   The RNA sequence input
 *  \param structure  A pointer to a char array where a base pair probability information can be stored in a pseudo-dot-bracket notation (may be NULL, too)
 *  \return           The Gibbs free energy of the ensemble (\f$G = -RT \cdot \log(Q) \f$) in kcal/mol
 */
DEPRECATED(float   pf_fold(const char *sequence,
                char *structure));

/**
 *  \brief Compute the partition function of a circular RNA sequence
 * 
 *  \ingroup pf_fold
 *
 *  \note           The global array #pr is deprecated and the user who wants the calculated
 *                  base pair probabilities for further computations is advised to use the function
 *                  export_bppm().
 *  \note           \b OpenMP:
 *                  This function is not entirely threadsafe. While the recursions are working on their
 *                  own copies of data the model details for the recursions are determined from the global
 *                  settings just before entering the recursions. Consider using pf_fold_par() for a
 *                  really threadsafe implementation.
 *  \pre            This function takes its model details from the global variables provided in \e RNAlib
 *  \post           After successful run the hidden folding matrices are filled with the appropriate Boltzmann factors.
 *                  Depending on whether the global variable #do_backtrack was set the base pair probabilities are already
 *                  computed and may be accessed for further usage via the export_bppm() function.
 *                  A call of free_pf_arrays() will free all memory allocated by this function.
 *                  Successive calls will first free previously allocated memory before starting the computation.
 *  \see            vrna_pf_fold()
 *  \deprecated     Use vrna_pf_fold() instead!
 *  \param[in]      sequence   The RNA sequence input
 *  \param[in,out]  structure  A pointer to a char array where a base pair probability information can be
 *                  stored in a pseudo-dot-bracket notation (may be NULL, too)
 *  \return         The Gibbs free energy of the ensemble (\f$G = -RT \cdot \log(Q) \f$) in kcal/mol
 */
DEPRECATED(float   pf_circ_fold( const char *sequence,
                      char *structure));

/**
 *  \brief Sample a secondary structure from the Boltzmann ensemble according its probability\n
 *
 *  \ingroup subopt_stochbt
 *  \pre    #st_back has to be set to 1 before calling pf_fold() or pf_fold_par()
 *  \pre    pf_fold_par() or pf_fold() have to be called first to fill the partition function matrices
 *
 *  \param  sequence  The RNA sequence
 *  \return           A sampled secondary structure in dot-bracket notation
 */
DEPRECATED(char    *pbacktrack(char *sequence));

DEPRECATED(char    *pbacktrack5(char *sequence, int length));

/**
 *  \brief Sample a secondary structure of a subsequence from the Boltzmann ensemble according its probability\n
 *
 *  \ingroup subopt_stochbt
 *  \pre    The fold compound has to be obtained using the #VRNA_OPTION_HYBRID option in vrna_get_fold_compound()
 *  \pre    vrna_pf_fold() hasto be called first to fill the partition function matrices
 *
 *  \param  vc      The fold compound data structure
 *  \param  length  The length of the subsequence to consider (starting with 5' end)
 *  \return         A sampled secondary structure in dot-bracket notation
 */
char    *vrna_pbacktrack5(vrna_fold_compound *vc, int length);

/**
 *  \brief Sample a secondary structure from the Boltzmann ensemble according its probability\n
 *
 *  \ingroup subopt_stochbt
 *  \pre    The fold compound has to be obtained using the #VRNA_OPTION_HYBRID option in vrna_get_fold_compound()
 *  \pre    vrna_pf_fold() hasto be called first to fill the partition function matrices
 *
 *  \note The function will automagically detect cicular RNAs based on the model_details in exp_params as
 *        provided via the #vrna_fold_compound
 *
 *  \param  vc      The fold compound data structure
 *  \param  length  The length of the subsequence to consider (starting with 5' end)
 *  \return         A sampled secondary structure in dot-bracket notation
 */
char    *vrna_pbacktrack(vrna_fold_compound *vc);

/**
 *  \brief Sample a secondary structure of a circular RNA from the Boltzmann ensemble according its probability
 * 
 *  This function does the same as \ref pbacktrack() but assumes the RNA molecule to be circular
 *
 *  \ingroup subopt_stochbt

 *  \pre    #st_back has to be set to 1 before calling pf_fold() or pf_fold_par()
 *  \pre    pf_fold_par() or pf_circ_fold() have to be called first to fill the partition function matrices
 *
 *  \deprecated Use vrna_pbacktrack() instead.
 *
 *  \param  sequence  The RNA sequence
 *  \return           A sampled secondary structure in dot-bracket notation
 */
DEPRECATED(char    *pbacktrack_circ(char *sequence));

/**
 *  \brief Free arrays for the partition function recursions
 *
 *  Call this function if you want to free all allocated memory associated with
 *  the partition function forward recursion.
 *  \note Successive calls of pf_fold(), pf_circ_fold() already check if they should free
 *  any memory from a previous run.
 *  \note <b>OpenMP notice:</b><br>
 *  This function should be called before leaving a thread in order to avoid leaking memory
 *  
 *  \deprecated See vrna_fold_compound and its related functions for how to free memory
 *  occupied by the dynamic programming matrices
 *
 *  \ingroup pf_fold
 *
 *  \post   All memory allocated by pf_fold_par(), pf_fold() or pf_circ_fold() will be free'd
 *  \see    pf_fold_par(), pf_fold(), pf_circ_fold()
 */
DEPRECATED(void  free_pf_arrays(void));

/**
 *  \brief Recalculate energy parameters
 * 
 *  Call this function to recalculate the pair matrix and energy parameters
 *  after a change in folding parameters like #temperature
 *
 *  \deprecated Use vrna_update_pf_params() instead
 *  \ingroup pf_fold
 *
 */
DEPRECATED(void  update_pf_params(int length));

/**
 *  \brief Recalculate energy parameters
 *
 *  \deprecated Use vrna_update_pf_params() instead
 *  \ingroup pf_fold
 *
 */
DEPRECATED(void update_pf_params_par(int length, vrna_exp_param_t *parameters));


/**
 *  \brief Update the energy parameters for subsequent partition function computations
 *
 *  This function can be used to properly assign new energy parameters to
 *  a vrna_fold_compound. For this purpose, the data of the provided pointer `params`
 *  will be copied into `vc` and a recomputation of the partition function scaling
 *  factor is issued, if the `pf_scale` attribute of `params` is below `1.0`.
 *
 *  \see vrna_rescale_pf_params(), vrna_exp_param_t, vrna_md_t, vrna_exp_params_get()
 *
 *  \ingroup pf_fold
 *  \param  vc      The fold compound data structure
 *  \param  params  A pointer to the new energy parameters
 */
void vrna_update_pf_params( vrna_fold_compound *vc,
                            vrna_exp_param_t *params);

/**
 *  \brief Rescale Boltzmann factors for partition function computations
 *
 *  This function may be used to (automatically) rescale the Boltzmann factors used
 *  in partition function computations. Since partition functions over subsequences
 *  can easily become extremely large, the RNAlib internally rescales them to avoid
 *  numerical over- and/or underflow. Therefore, a proper scaling factor \f$s\f$ needs to
 *  be chosen that in turn is then used to normalize the corresponding
 *  partition functions \f$\hat{q}[i,j] = q[i,j] / s^{(j-i+1)}\f$.
 *
 *  This function provides two ways to automatically adjust the scaling
 *  factor.
 *  1. Automatic guess
 *  2. Automatic adjustment according to MFE
 *
 *  Passing `NULL` as second parameter activates the _automatic guess mode_. Here,
 *  the scaling factor is recomputed according to a mean free energy of `184.3*length` cal
 *  for random sequences.
 *  \note This recomputation only takes place if the `pf_scale` attribute of the
 *  `exp_params` datastructure contained in `vc` has a value below `1.0`.
 *
 *  On the other hand, if the MFE for a sequence is known, it can be used to recompute
 *  a more robust scaling factor, since it represents the lowest free energy of the entire
 *  ensemble of structures, i.e. the highest Boltzmann factor. To activate this second
 *  mode of _automatic adjustment according to MFE_, a pointer to the MFE value needs to
 *  be passed as second argument. This value is then taken to compute the scaling factor
 *  as \f$ s = exp((sfact * MFE) / kT / length )\f$, where sfact is an additional
 *  scaling weight located in the vrna_md_t datastructure of `exp_params` in `vc`.
 *
 *  The computed scaling factor \f$s\f$ will be stored as `pf_scale` attribute of the
 *  `exp_params` datastructure in `vc`.
 *
 *  \see vrna_update_pf_params(), vrna_md_t, vrna_exp_param_t, vrna_fold_compound
 *  \ingroup pf_fold
 *
 *  \param  vc  The fold compound data structure
 *  \param  mfe A pointer to the MFE (in kcal/mol) or NULL
 */
void vrna_rescale_pf_params(vrna_fold_compound *vc,
                            double *mfe);

/**
 *  \brief Get a pointer to the base pair probability array
 *  \ingroup  pf_fold
 *
 *  Accessing the base pair probabilities for a pair (i,j) is achieved by
 *  \code
 *  FLT_OR_DBL *pr  = export_bppm();
 *  pr_ij           = pr[iindx[i]-j];
 *  \endcode
 *
 *  \pre      Call pf_fold_par(), pf_fold() or pf_circ_fold() first to fill the base pair probability array
 *
 *  \see pf_fold(), pf_circ_fold(), vrna_get_iindx()
 *
 *  \return A pointer to the base pair probability array
 */
DEPRECATED(FLT_OR_DBL  *export_bppm(void));

/*
#################################################
# OTHER PARTITION FUNCTION RELATED DECLARATIONS #
#################################################
*/


/**
 *  \brief Get the pointers to (almost) all relavant computation arrays used in partition function computation
 *
 *  \ingroup    pf_fold
 *  \pre        In order to assign meaningful pointers, you have to call pf_fold_par() or pf_fold() first!
 *  \see        pf_fold_par(), pf_fold(), pf_circ_fold()
 *  \param[out] S_p       A pointer to the 'S' array (integer representation of nucleotides)
 *  \param[out] S1_p      A pointer to the 'S1' array (2nd integer representation of nucleotides)
 *  \param[out] ptype_p   A pointer to the pair type matrix
 *  \param[out] qb_p      A pointer to the Q<sup>B</sup> matrix
 *  \param[out] qm_p      A pointer to the Q<sup>M</sup> matrix
 *  \param[out] q1k_p     A pointer to the 5' slice of the Q matrix (\f$q1k(k) = Q(1, k)\f$)
 *  \param[out] qln_p     A pointer to the 3' slice of the Q matrix (\f$qln(l) = Q(l, n)\f$)
 *  \return     Non Zero if everything went fine, 0 otherwise
 */
int get_pf_arrays(short **S_p,
                  short **S1_p,
                  char **ptype_p,
                  FLT_OR_DBL **qb_p,
                  FLT_OR_DBL **qm_p,
                  FLT_OR_DBL **q1k_p,
                  FLT_OR_DBL **qln_p);

/**
 *  \brief Get the free energy of a subsequence from the q[] array
 */
double get_subseq_F(int i, int j);


/**
 *  \brief Get the mean base pair distance of the last partition function computation
 * 
 *  \ingroup pf_fold
 *
 *  \deprecated Use vrna_mean_bp_distance() or vrna_mean_bp_distance_pr() instead!
 *  \see vrna_mean_bp_distance(), vrna_mean_bp_distance_pr()
 * 
 *  \param    length
 *  \return  mean base pair distance in thermodynamic ensemble
 */
DEPRECATED(double  mean_bp_distance(int length));

/**
 *  \brief Get the mean base pair distance in the thermodynamic ensemble
 * 
 *  This is a threadsafe implementation of \ref mean_bp_dist() !
 * 
 *  \f$<d> = \sum_{a,b} p_a p_b d(S_a,S_b)\f$\n
 *  this can be computed from the pair probs \f$p_ij\f$ as\n
 *  \f$<d> = \sum_{ij} p_{ij}(1-p_{ij})\f$
 * 
 *  \deprecated Use vrna_mean_bp_distance() or vrna_mean_bp_distance_pr() instead!
 * 
 *  \ingroup pf_fold
 *
 *  \param length The length of the sequence
 *  \param pr     The matrix containing the base pair probabilities
 *  \return       The mean pair distance of the structure ensemble
 */
DEPRECATED(double mean_bp_distance_pr(int length, FLT_OR_DBL *pr));

/**
 *  \brief Get the mean base pair distance in the thermodynamic ensemble from a probability matrix
 * 
 *  \f$<d> = \sum_{a,b} p_a p_b d(S_a,S_b)\f$\n
 *  this can be computed from the pair probs \f$p_ij\f$ as\n
 *  \f$<d> = \sum_{ij} p_{ij}(1-p_{ij})\f$
 * 
 *  \ingroup pf_fold
 *
 *  \param length The length of the sequence
 *  \param pr     The matrix containing the base pair probabilities
 *  \return       The mean pair distance of the structure ensemble
 */
double vrna_mean_bp_distance_pr(int length, FLT_OR_DBL *pr);

/**
 *  \brief Get the mean base pair distance in the thermodynamic ensemble
 * 
 *  \f$<d> = \sum_{a,b} p_a p_b d(S_a,S_b)\f$\n
 *  this can be computed from the pair probs \f$p_ij\f$ as\n
 *  \f$<d> = \sum_{ij} p_{ij}(1-p_{ij})\f$
 * 
 *  \ingroup pf_fold
 *
 *  \param vc     The fold compound data structure
 *  \return       The mean pair distance of the structure ensemble
 */
double vrna_mean_bp_distance(vrna_fold_compound *vc);

plist *stackProb(double cutoff);


/*
#################################################
# DEPRECATED FUNCTIONS                          #
#################################################
*/

/**
 *  \brief Allocate space for pf_fold()
 * 
 *  \deprecated This function is obsolete and will be removed soon!
 */
DEPRECATED(void init_pf_fold(int length));

/**
 *  \deprecated This function is deprecated and should not be used anymore as it is not threadsafe!
 *  \see get_centroid_struct_pl(), get_centroid_struct_pr()
 */
DEPRECATED(char *centroid(int length, double *dist));

/**
 *  \deprecated This function is deprecated and should not be used anymore as it is not threadsafe!
 *  \see vrna_get_centroid_struct(), vrna_get_centroid_struct_pr(), vrna_get_centroid_struct_pl()
 */
DEPRECATED(char *get_centroid_struct_gquad_pr(int length,
                                  double *dist));

/**
 *  get the mean pair distance of ensemble
 * 
 *  \deprecated This function is not threadsafe and should not be used anymore. Use \ref mean_bp_distance() instead!
 */
DEPRECATED(double mean_bp_dist(int length));

/**
 *  \deprecated Use \ref exp_E_IntLoop() from loop_energies.h instead
 */
DEPRECATED(double expLoopEnergy(int u1,
                                int u2,
                                int type,
                                int type2,
                                short si1,
                                short sj1,
                                short sp1,
                                short sq1));

/**
 *  \deprecated Use exp_E_Hairpin() from loop_energies.h instead
 */
DEPRECATED(double expHairpinEnergy( int u,
                                    int type,
                                    short si1,
                                    short sj1,
                                    const char *string));

/* this doesn't work if free_pf_arrays() is called before */
DEPRECATED(void assign_plist_gquad_from_pr(plist **pl,
                                int length,
                                double cut_off));

#endif
