#ifndef __VIENNA_RNA_PACKAGE_PART_FUNC_H__
#define __VIENNA_RNA_PACKAGE_PART_FUNC_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif


/**
 *  \file part_func.h
 * 
 *  \brief Partition function of single RNA sequences
 * 
 *  This file includes (almost) all function declarations within the <b>RNAlib</b> that are related to
 *  Partion function folding...
 *
 *  \note If you plan on using the functions provided from this section of the RNAlib concurrently
 *  via <b>OpenMP</b> you have to place a <i>COPYIN</i> clause right before your <i>PARALLEL</i>
 *  directive! Otherwise, some functions may not behave as expected.
 *  A complete list of variables that have to be passed to the <i>COPYIN</i> clause can be found
 *  in the detailed description of each function below.
 */

/**
 *  a flag indicating that auxilary arrays are needed throughout the computations which are necessary for stochastic backtracking
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
 *  If tha parameter calculate_bppm is set to 0 base pairing probabilities will not
 *  be computed (saving CPU time), otherwise after calculations took place #pr will
 *  contain the probability that bases \a i and \a j pair.
 *  \note The global array #pr is deprecated and the user who wants the calculated
 *  base pair probabilities for further computations is advised to use the function
 *  export_bppm()
 * 
 *  \see pf_circ_fold(), bppm_to_structure(), export_bppm(), get_boltzmann_factors()
 * 
 *  \param sequence       The RNA sequence input
 *  \param structure      A pointer to a char array where a base pair probability information can be stored in a pseudo-dot-bracket notation (may be NULL, too)
 *  \param parameters     Data structure containing the precalculated Boltzmann factors
 *  \param calculate_bppm Switch to Base pair probability calculations on/off (0==off)
 *  \param is_constrained Switch to indicate that a structure contraint is passed via the structure argument (0==off)
 *  \param is_circular    Switch to (de-)activate postprocessing steps in case RNA sequence is circular (0==off)
 *  \return               The Gibbs free energy of the ensemble (\f$G = -RT \cdot \log(Q) \f$) in kcal/mol
 */
float   pf_fold_par(  const char *sequence,
                      char *structure,
                      pf_paramT *parameters,
                      int calculate_bppm,
                      int is_constrained,
                      int is_circular);

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
 *  \note The global array #pr is deprecated and the user who wants the calculated
 *  base pair probabilities for further computations is advised to use the function
 *  export_bppm()
 * 
 *  \see pf_circ_fold(), bppm_to_structure(), export_bppm()
 * 
 *  \param sequence   The RNA sequence input
 *  \param structure  A pointer to a char array where a base pair probability information can be stored in a pseudo-dot-bracket notation (may be NULL, too)
 *  \return           The Gibbs free energy of the ensemble (\f$G = -RT \cdot \log(Q) \f$) in kcal/mol
 */
float   pf_fold(const char *sequence,
                char *structure);

/**
 *  \brief Compute the partition function of a circular RNA sequence
 * 
 *  \see pf_fold(), pf_fold_par()
 * 
 *  \param sequence   The RNA sequence input
 *  \param structure  A pointer to a char array where a base pair probability information can be stored in a pseudo-dot-bracket notation (may be NULL, too)
 *  \return           The Gibbs free energy of the ensemble (\f$G = -RT \cdot \log(Q) \f$) in kcal/mol
 */
float   pf_circ_fold( const char *sequence,
                      char *structure);

/**
 *  \brief Sample a secondary structure from the Boltzmann ensemble according its probability\n
 *
 *  \note You have to call pf_fold() first in order to fill the partition function matrices
 *
 *  \note <b>OpenMP notice:</b><br>This function relies on passing the following
 *  variables to the appropriate <i>COPYIN</i> clause <i>(additionally to the ones needed by pf_fold())</i>:<br>
 *  pstruc, sequence
 *
 *  \param  sequence  The RNA sequence
 *  \return           A sampled secondary structure in dot-bracket notation
 */
char    *pbacktrack(char *sequence);

/**
 *  \brief Sample a secondary structure of a circular RNA from the Boltzmann ensemble according its probability
 * 
 *  This function does the same as \ref pbacktrack() but assumes the RNA molecule to be circular
 *
 *  \note <b>OpenMP notice:</b><br>This function relies on passing the following
 *  variables to the appropriate <i>COPYIN</i> clause <i>(additionally to the ones needed by pf_fold())</i>:<br>
 *  pstruc, sequence
 *
 *  \param  sequence  The RNA sequence
 *  \return           A sampled secondary structure in dot-bracket notation
 */
char    *pbacktrack_circ(char *sequence);

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
 *  \see pf_fold(), pf_circ_fold()
 */
void  free_pf_arrays(void);

/**
 *  \brief Recalculate energy parameters
 * 
 *  Call this function to recalculate the pair matrix and energy parameters
 *  after a change in folding parameters like #temperature
 */
void  update_pf_params(int length);

void update_pf_params_par(int length, pf_paramT *parameters);

/**
 *  \brief Get a pointer to the base pair probability array
 *
 *  Accessing the base pair probabilities for a pair (i,j) is achieved by
 *  \verbatim FLT_OR_DBL *pr = export_bppm(); pr_ij = pr[iindx[i]-j]; \endverbatim
 *
 *  \note Call pf_fold() before using this function!
 *
 *  \see pf_fold(), pf_circ_fold(), get_iindx()
 *
 *  \return A pointer to the base pair probability array
 */
FLT_OR_DBL  *export_bppm(void);

/*
#################################################
# OTHER PARTITION FUNCTION RELATED DECLARATIONS #
#################################################
*/

/**
 *  \brief Create a plist from a probability matrix
 * 
 *  The probability matrix given is parsed and all pair probabilities above
 *  the given threshold are used to create an entry in the plist
 * 
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 * 
 *  \note This function is threadsafe
 * 
 *  \param pl     A pointer to the plist that is to be created
 *  \param probs  The probability matrix used for creting the plist
 *  \param length The length of the RNA sequence
 *  \param cutoff The cutoff value
 */
void  assign_plist_from_pr( plist **pl,
                            FLT_OR_DBL *probs,
                            int length,
                            double cutoff);

/**
 *  \brief Get the pointers to (almost) all relavant computation arrays used in partition function computation
 *
 *  \note In order to assign meaningful pointers, you have to call pf_fold first!
 *
 *  \see pf_fold(), pf_circ_fold()
 *
 *  \param S_p      A pointer to the 'S' array (integer representation of nucleotides)
 *  \param S1_p     A pointer to the 'S1' array (2nd integer representation of nucleotides)
 *  \param ptype_p  A pointer to the pair type matrix
 *  \param qb_p     A pointer to the Q<sup>B</sup> matrix
 *  \param qm_p     A pointer to the Q<sup>M</sup> matrix
 *  \param q1k_p    A pointer to the 5' slice of the Q matrix (\f$q1k(k) = Q(1, k)\f$)
 *  \param qln_p    A pointer to the 3' slice of the Q matrix (\f$qln(l) = Q(l, n)\f$)
 *  \return         Non Zero if everything went fine, 0 otherwise
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
 *  \brief Get the centroid structure of the ensemble
 * 
 *  This function is a threadsafe replacement for \ref centroid() with a 'plist' input
 * 
 *  The centroid is the structure with the minimal average distance to all other structures
 *  \n \f$ <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij} \f$ \n
 *  Thus, the centroid is simply the structure containing all pairs with \f$p_ij>0.5\f$
 *  The distance of the centroid to the ensemble is written to the memory adressed by \a dist.
 * 
 *  \param length The length of the sequence
 *  \param dist   A pointer to the distance variable where the centroid distance will be written to
 *  \param pl     A pair list containing base pair probability information about the ensemble
 *  \return       The centroid structure of the ensemble in dot-bracket notation
 */
char  *get_centroid_struct_pl(int length,
                              double *dist,
                              plist *pl);

/**
 *  \brief Get the centroid structure of the ensemble
 * 
 *  This function is a threadsafe replacement for \ref centroid() with a probability array input
 * 
 *  The centroid is the structure with the minimal average distance to all other structures
 *  \n \f$ <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij} \f$ \n
 *  Thus, the centroid is simply the structure containing all pairs with \f$p_ij>0.5\f$
 *  The distance of the centroid to the ensemble is written to the memory adressed by \a dist.
 * 
 *  \param length The length of the sequence
 *  \param dist   A pointer to the distance variable where the centroid distance will be written to
 *  \param pr     A upper triangular matrix containing base pair probabilities (access via iindx \ref get_iindx() )
 *  \return       The centroid structure of the ensemble in dot-bracket notation
 */
char  *get_centroid_struct_pr(int length,
                              double *dist,
                              FLT_OR_DBL *pr);

/**
 *  \brief Get the mean base pair distance of the last partition function computation
 * 
 *  \note To ensure thread-safety, use the function mean_bp_distance_pr() instead!
 *
 *  \see mean_bp_distance_pr()
 * 
 *  \param    length
 *  \return  mean base pair distance in thermodynamic ensemble
 */
double  mean_bp_distance(int length);

/**
 *  \brief Get the mean base pair distance in the thermodynamic ensemble
 * 
 *  This is a threadsafe implementation of \ref mean_bp_dist() !
 * 
 *  \f$<d> = \sum_{a,b} p_a p_b d(S_a,S_b)\f$\n
 *  this can be computed from the pair probs \f$p_ij\f$ as\n
 *  \f$<d> = \sum_{ij} p_{ij}(1-p_{ij})\f$
 * 
 *  \note This function is threadsafe
 * 
 *  \param length The length of the sequence
 *  \param pr     The matrix containing the base pair probabilities
 *  \return       The mean pair distance of the structure ensemble
 */
double  mean_bp_distance_pr(int length,
                            FLT_OR_DBL *pr);

/**
 *  \brief Create a dot-bracket like structure string from base pair probability matrix
 */
void  bppm_to_structure(char *structure,
                        FLT_OR_DBL *pr,
                        unsigned int length);

plist *stackProb(double cutoff);

/**
 *  \brief Get a pseudo dot bracket notation for a given probability information
 */
char    bppm_symbol(const float *x);


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
DEPRECATED(char *centroid(int length,
                          double *dist));     /* mean pair distance of ensemble */

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

#endif
