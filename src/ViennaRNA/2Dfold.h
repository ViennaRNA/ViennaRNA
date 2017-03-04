#ifndef VIENNA_RNA_PACKAGE_TWO_D_FOLD_H
#define VIENNA_RNA_PACKAGE_TWO_D_FOLD_H

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

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
 *  @file 2Dfold.h
 *  @ingroup kl_neighborhood
 *  @brief MFE structures for base pair distance classes
 *
 */

/**
 *  @addtogroup kl_neighborhood_mfe
 *  @brief Compute the minimum free energy (MFE) and secondary structures for a partitioning of
 *  the secondary structure space according to the base pair distance to two fixed reference structures
 *  basepair distance to two fixed reference structures
 *  @see For further details, we refer to Lorenz et al. 2009 @cite lorenz:2009
 *
 *  @{
 *  @ingroup  kl_neighborhood_mfe
 *
 */

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>

/**
 *  @brief Solution element returned from vrna_mfe_TwoD()
 *
 *  This element contains free energy and structure for the appropriate
 *  kappa (k), lambda (l) neighborhood
 *  The datastructure contains two integer attributes 'k' and 'l'
 *  as well as an attribute 'en' of type float representing the free energy
 *  in kcal/mol and an attribute 's' of type char* containg the secondary
 *  structure representative,
 *
 *  A value of #INF in k denotes the end of a list
 *
 *  @see  vrna_mfe_TwoD()
 */
typedef struct vrna_sol_TwoD_t {
  int   k;        /**<  @brief  Distance to first reference */
  int   l;        /**<  @brief  Distance to second reference */
  float en;       /**<  @brief  Free energy in kcal/mol */
  char  *s;       /**<  @brief  MFE representative structure in dot-bracket notation */
} vrna_sol_TwoD_t;


/**
 * @brief Compute MFE's and representative for distance partitioning
 *
 * This function computes the minimum free energies and a representative
 * secondary structure for each distance class according to the two references
 * specified in the datastructure 'vars'.
 * The maximum basepair distance to each of both references may be set
 * by the arguments 'distance1' and 'distance2', respectively.
 * If both distance arguments are set to '-1', no restriction is assumed and
 * the calculation is performed for each distance class possible.
 *
 * The returned list contains an entry for each distance class. If a maximum
 * basepair distance to either of the references was passed, an entry with
 * k=l=-1 will be appended in the list, denoting the class where all structures
 * exceeding the maximum will be thrown into.
 * The end of the list is denoted by an attribute value of #INF in
 * the k-attribute of the list entry.
 *
 *  @see  vrna_fold_compound_TwoD(), vrna_fold_compound_free(), vrna_pf_TwoD()
 *        vrna_backtrack5_TwoD(), #vrna_sol_TwoD_t, #vrna_fold_compound_t
 *
 *  @param vc         The datastructure containing all precomputed folding attributes
 *  @param distance1  maximum distance to reference1 (-1 means no restriction)
 *  @param distance2  maximum distance to reference2 (-1 means no restriction)
 *  @return           A list of minimum free energies (and corresponding structures)
 *                    for each distance class
 */
vrna_sol_TwoD_t *
vrna_mfe_TwoD(vrna_fold_compound_t  *vc,
              int                   distance1,
              int                   distance2);


/**
 * @brief Backtrack a minimum free energy structure from a 5' section of specified length
 *
 * This function allows one to backtrack a secondary structure beginning at the 5' end, a specified
 * length and residing in a specific distance class.
 * If the argument 'k' gets a value of -1, the structure that is backtracked is assumed to
 * reside in the distance class where all structures exceeding the maximum basepair distance
 * specified in vrna_mfe_TwoD() belong to.
 * @note The argument 'vars' must contain precalculated energy values in the energy matrices,
 * i.e. a call to vrna_mfe_TwoD() preceding this function is mandatory!
 *
 * @see vrna_mfe_TwoD()
 *
 * @param vc    The datastructure containing all precomputed folding attributes
 * @param j     The length in nucleotides beginning from the 5' end
 * @param k     distance to reference1 (may be -1)
 * @param l     distance to reference2
 */
char *
vrna_backtrack5_TwoD(vrna_fold_compound_t *vc,
                     int                  k,
                     int                  l,
                     unsigned int         j);


#ifdef  VRNA_BACKWARD_COMPAT

#define TwoDfold_solution       vrna_sol_TwoD_t         /* restore compatibility of struct rename */

/**
 *  @brief Variables compound for 2Dfold MFE folding
 *
 *  @deprecated This data structure will be removed from the library soon!
 *              Use #vrna_fold_compound_t and the corresponding functions vrna_fold_compound_TwoD(),
 *              vrna_mfe_TwoD(), and vrna_fold_compound_free() instead!
 */
typedef struct TwoDfold_vars {
  vrna_param_t          *P;           /**<  @brief  Precomputed energy parameters and model details */
  int                   do_backtrack; /**<  @brief  Flag whether to do backtracing of the structure(s) or not */
  char                  *ptype;       /**<  @brief  Precomputed array of pair types */
  char                  *sequence;    /**<  @brief  The input sequence  */
  short                 *S, *S1;      /**<  @brief  The input sequences in numeric form */
  unsigned int          maxD1;        /**<  @brief  Maximum allowed base pair distance to first reference */
  unsigned int          maxD2;        /**<  @brief  Maximum allowed base pair distance to second reference */


  unsigned int          *mm1;       /**<  @brief  Maximum matching matrix, reference struct 1 disallowed */
  unsigned int          *mm2;       /**<  @brief  Maximum matching matrix, reference struct 2 disallowed */

  int                   *my_iindx;  /**<  @brief  Index for moving in quadratic distancy dimensions */

  double                temperature;

  unsigned int          *referenceBPs1; /**<  @brief  Matrix containing number of basepairs of reference structure1 in interval [i,j] */
  unsigned int          *referenceBPs2; /**<  @brief  Matrix containing number of basepairs of reference structure2 in interval [i,j] */
  unsigned int          *bpdist;        /**<  @brief  Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j] */

  short                 *reference_pt1;
  short                 *reference_pt2;
  int                   circ;
  int                   dangles;
  unsigned int          seq_length;

  int                   ***E_F5;
  int                   ***E_F3;
  int                   ***E_C;
  int                   ***E_M;
  int                   ***E_M1;
  int                   ***E_M2;

  int                   **E_Fc;
  int                   **E_FcH;
  int                   **E_FcI;
  int                   **E_FcM;

  int                   **l_min_values;
  int                   **l_max_values;
  int                   *k_min_values;
  int                   *k_max_values;

  int                   **l_min_values_m;
  int                   **l_max_values_m;
  int                   *k_min_values_m;
  int                   *k_max_values_m;

  int                   **l_min_values_m1;
  int                   **l_max_values_m1;
  int                   *k_min_values_m1;
  int                   *k_max_values_m1;

  int                   **l_min_values_f;
  int                   **l_max_values_f;
  int                   *k_min_values_f;
  int                   *k_max_values_f;

  int                   **l_min_values_f3;
  int                   **l_max_values_f3;
  int                   *k_min_values_f3;
  int                   *k_max_values_f3;

  int                   **l_min_values_m2;
  int                   **l_max_values_m2;
  int                   *k_min_values_m2;
  int                   *k_max_values_m2;

  int                   *l_min_values_fc;
  int                   *l_max_values_fc;
  int                   k_min_values_fc;
  int                   k_max_values_fc;

  int                   *l_min_values_fcH;
  int                   *l_max_values_fcH;
  int                   k_min_values_fcH;
  int                   k_max_values_fcH;

  int                   *l_min_values_fcI;
  int                   *l_max_values_fcI;
  int                   k_min_values_fcI;
  int                   k_max_values_fcI;

  int                   *l_min_values_fcM;
  int                   *l_max_values_fcM;
  int                   k_min_values_fcM;
  int                   k_max_values_fcM;

  /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
  int                   *E_F5_rem;
  int                   *E_F3_rem;
  int                   *E_C_rem;
  int                   *E_M_rem;
  int                   *E_M1_rem;
  int                   *E_M2_rem;

  int                   E_Fc_rem;
  int                   E_FcH_rem;
  int                   E_FcI_rem;
  int                   E_FcM_rem;

#ifdef COUNT_STATES
  unsigned long         ***N_F5;
  unsigned long         ***N_C;
  unsigned long         ***N_M;
  unsigned long         ***N_M1;
#endif

  vrna_fold_compound_t  *compatibility;
} TwoDfold_vars;

/**
 *  @brief Get a structure of type TwoDfold_vars prefilled with current global settings
 *
 *  This function returns a datastructure of type TwoDfold_vars.
 *  The data fields inside the TwoDfold_vars are prefilled by global settings and all memory
 *  allocations necessary to start a computation are already done for the convenience of the user
 *
 *  @note Make sure that the reference structures are compatible with the sequence according to Watson-Crick- and Wobble-base pairing
 *
 *  @deprecated Use the new API that relies on #vrna_fold_compound_t and the corresponding functions
 *              vrna_fold_compound_TwoD(), vrna_mfe_TwoD(), and vrna_fold_compound_free() instead!
 *
 *  @param seq          The RNA sequence
 *  @param structure1   The first reference structure in dot-bracket notation
 *  @param structure2   The second reference structure in dot-bracket notation
 *  @param circ         A switch to indicate the assumption to fold a circular instead of linear RNA (0=OFF, 1=ON)
 *  @returns            A datastructure prefilled with folding options and allocated memory
 */
DEPRECATED(TwoDfold_vars *
           get_TwoDfold_variables(const char  *seq,
                                  const char  *structure1,
                                  const char  *structure2,
                                  int         circ));

/**
 *  @brief Destroy a TwoDfold_vars datastructure without memory loss
 *
 *  This function free's all allocated memory that depends on the datastructure given.
 *
 *  @deprecated Use the new API that relies on #vrna_fold_compound_t and the corresponding functions
 *              vrna_fold_compound_TwoD(), vrna_mfe_TwoD(), and vrna_fold_compound_free() instead!
 *
 *  @param our_variables  A pointer to the datastructure to be destroyed
 */
DEPRECATED(void
           destroy_TwoDfold_variables(TwoDfold_vars *our_variables));

/**
 * @brief Compute MFE's and representative for distance partitioning
 *
 * This function computes the minimum free energies and a representative
 * secondary structure for each distance class according to the two references
 * specified in the datastructure 'vars'.
 * The maximum basepair distance to each of both references may be set
 * by the arguments 'distance1' and 'distance2', respectively.
 * If both distance arguments are set to '-1', no restriction is assumed and
 * the calculation is performed for each distance class possible.
 *
 * The returned list contains an entry for each distance class. If a maximum
 * basepair distance to either of the references was passed, an entry with
 * k=l=-1 will be appended in the list, denoting the class where all structures
 * exceeding the maximum will be thrown into.
 * The end of the list is denoted by an attribute value of #INF in
 * the k-attribute of the list entry.
 *
 *  @deprecated Use the new API that relies on #vrna_fold_compound_t and the corresponding functions
 *              vrna_fold_compound_TwoD(), vrna_mfe_TwoD(), and vrna_fold_compound_free() instead!
 *
 * @param vars      the datastructure containing all predefined folding attributes
 * @param distance1 maximum distance to reference1 (-1 means no restriction)
 * @param distance2 maximum distance to reference2 (-1 means no restriction)
 */
DEPRECATED(TwoDfold_solution *
           TwoDfoldList(TwoDfold_vars *vars,
                        int           distance1,
                        int           distance2));

/**
 * @brief Backtrack a minimum free energy structure from a 5' section of specified length
 *
 * This function allows one to backtrack a secondary structure beginning at the 5' end, a specified
 * length and residing in a specific distance class.
 * If the argument 'k' gets a value of -1, the structure that is backtracked is assumed to
 * reside in the distance class where all structures exceeding the maximum basepair distance
 * specified in TwoDfold() belong to.
 * @note The argument 'vars' must contain precalculated energy values in the energy matrices,
 * i.e. a call to TwoDfold() preceding this function is mandatory!
 *
 *  @deprecated Use the new API that relies on #vrna_fold_compound_t and the corresponding functions
 *              vrna_fold_compound_TwoD(), vrna_mfe_TwoD(), vrna_backtrack5_TwoD(), and
 *              vrna_fold_compound_free() instead!
 *
 * @param j     The length in nucleotides beginning from the 5' end
 * @param k     distance to reference1 (may be -1)
 * @param l     distance to reference2
 * @param vars  the datastructure containing all predefined folding attributes
 */
DEPRECATED(char *TwoDfold_backtrack_f5(unsigned int   j,
                                       int            k,
                                       int            l,
                                       TwoDfold_vars  *vars));

/**
 *
 */
DEPRECATED(TwoDfold_solution **TwoDfold(TwoDfold_vars *our_variables,
                                        int           distance1,
                                        int           distance2));


#endif

/**
 *  @}
 */

#endif
