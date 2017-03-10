#ifndef VIENNA_RNA_PACKAGE_TWO_D_PF_FOLD_H
#define VIENNA_RNA_PACKAGE_TWO_D_PF_FOLD_H

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
 *  @file 2Dpfold.h
 *  @ingroup kl_neighborhood
 *  @brief Partition function implementations for base pair distance classes
 *
 */

/**
 *  @addtogroup kl_neighborhood_pf
 *  @brief Compute the partition function and stochastically sample secondary structures for a partitioning of
 *  the secondary structure space according to the base pair distance to two fixed reference structures
 *
 *  @{
 *  @ingroup  kl_neighborhood_pf
 */

#include <ViennaRNA/data_structures.h>

/**
 *  @brief Solution element returned from vrna_pf_TwoD()
 *
 *  This element contains the partition function for the appropriate
 *  kappa (k), lambda (l) neighborhood
 *  The datastructure contains two integer attributes 'k' and 'l'
 *  as well as an attribute 'q' of type #FLT_OR_DBL
 *
 *  A value of #INF in k denotes the end of a list
 *
 *  @see  vrna_pf_TwoD()
 */
typedef struct vrna_sol_TwoD_pf_t {
  int         k;  /**<  @brief  Distance to first reference */
  int         l;  /**<  @brief  Distance to second reference */
  FLT_OR_DBL  q;  /**<  @brief  partition function */
} vrna_sol_TwoD_pf_t;

/**
 * @brief Compute the partition function for all distance classes
 *
 * This function computes the partition functions for all distance classes
 * according the two reference structures specified in the datastructure 'vars'.
 * Similar to vrna_mfe_TwoD() the arguments maxDistance1 and maxDistance2 specify
 * the maximum distance to both reference structures. A value of '-1' in either of
 * them makes the appropriate distance restrictionless, i.e. all basepair distancies
 * to the reference are taken into account during computation.
 * In case there is a restriction, the returned solution contains an entry where
 * the attribute k=l=-1 contains the partition function for all structures exceeding
 * the restriction.
 * A value of #INF in the attribute 'k' of the returned list denotes the end of the list
 *
 * @see vrna_fold_compound_TwoD(), vrna_fold_compound_free(), #vrna_fold_compound
 *      #vrna_sol_TwoD_pf_t
 *
 * @param vc            The datastructure containing all necessary folding attributes and matrices
 * @param maxDistance1  The maximum basepair distance to reference1 (may be -1)
 * @param maxDistance2  The maximum basepair distance to reference2 (may be -1)
 * @returns             A list of partition funtions for the corresponding distance classes
 */
vrna_sol_TwoD_pf_t *
vrna_pf_TwoD(vrna_fold_compound_t *vc,
             int                  maxDistance1,
             int                  maxDistance2);


/** @} */ /* End of group kl_neighborhood_pf */

/**
 *  @addtogroup kl_neighborhood_stochbt
 *  @brief Contains functions related to stochastic backtracking from a specified distance class
 *  @{
 */

/**
 *  @brief Sample secondary structure representatives from a set of distance classes according to their
 *  Boltzmann probability
 *
 *  If the argument 'd1' is set to '-1', the structure will be backtracked in the distance class
 *  where all structures exceeding the maximum basepair distance to either of the references reside.
 *
 *  @pre      The argument 'vars' must contain precalculated partition function matrices,
 *            i.e. a call to vrna_pf_TwoD() preceding this function is mandatory!
 *
 *  @see      vrna_pf_TwoD()
 *
 *  @param[inout]  vc The #vrna_fold_compound_t datastructure containing all necessary folding attributes and matrices
 *  @param[in]  d1    The distance to reference1 (may be -1)
 *  @param[in]  d2    The distance to reference2
 *  @returns    A sampled secondary structure in dot-bracket notation
 */
char *
vrna_pbacktrack_TwoD(vrna_fold_compound_t *vc,
                     int                  d1,
                     int                  d2);


/**
 * @brief Sample secondary structure representatives with a specified length from a set of distance classes according to their
 *  Boltzmann probability
 *
 * This function does essentially the same as vrna_pbacktrack_TwoD() with the only difference that partial structures,
 * i.e. structures beginning from the 5' end with a specified length of the sequence, are backtracked
 *
 * @note      This function does not work (since it makes no sense) for circular RNA sequences!
 * @pre       The argument 'vars' must contain precalculated partition function matrices,
 *            i.e. a call to vrna_pf_TwoD() preceding this function is mandatory!
 *
 * @see       vrna_pbacktrack_TwoD(), vrna_pf_TwoD()
 *
 *  @param[inout] vc    The #vrna_fold_compound_t datastructure containing all necessary folding attributes and matrices
 *  @param[in]  d1      The distance to reference1 (may be -1)
 *  @param[in]  d2      The distance to reference2
 *  @param[in]  length  The length of the structure beginning from the 5' end
 *  @returns            A sampled secondary structure in dot-bracket notation
 */
char *
vrna_pbacktrack5_TwoD(vrna_fold_compound_t  *vc,
                      int                   d1,
                      int                   d2,
                      unsigned int          length);


/**
 *  @}
 */ /* End of group kl_neighborhood_stochbt */


#ifdef  VRNA_BACKWARD_COMPAT

#define TwoDpfold_solution       vrna_sol_TwoD_pf_t         /* restore compatibility of struct rename */

/**
 *  @brief  Variables compound for 2Dfold partition function folding
 *
 *  @deprecated This data structure will be removed from the library soon!
 *              Use #vrna_fold_compound_t and the corresponding functions vrna_fold_compound_TwoD(),
 *              vrna_pf_TwoD(), and vrna_fold_compound_free() instead!
 */
typedef struct {
  unsigned int          alloc;
  char                  *ptype;       /**<  @brief  Precomputed array of pair types */
  char                  *sequence;    /**<  @brief  The input sequence  */
  short                 *S, *S1;      /**<  @brief  The input sequences in numeric form */
  unsigned int          maxD1;        /**<  @brief  Maximum allowed base pair distance to first reference */
  unsigned int          maxD2;        /**<  @brief  Maximum allowed base pair distance to second reference */

  double                temperature;  /* temperature in last call to scale_pf_params */
  double                init_temp;    /* temperature in last call to scale_pf_params */
  FLT_OR_DBL            *scale;
  FLT_OR_DBL            pf_scale;
  vrna_exp_param_t      *pf_params; /* holds all [unscaled] pf parameters */

  int                   *my_iindx;  /**<  @brief  Index for moving in quadratic distancy dimensions */
  int                   *jindx;     /**<  @brief  Index for moving in the triangular matrix qm1 */

  short                 *reference_pt1;
  short                 *reference_pt2;

  unsigned int          *referenceBPs1; /**<  @brief  Matrix containing number of basepairs of reference structure1 in interval [i,j] */
  unsigned int          *referenceBPs2; /**<  @brief  Matrix containing number of basepairs of reference structure2 in interval [i,j] */
  unsigned int          *bpdist;        /**<  @brief  Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j] */

  unsigned int          *mm1;           /**<  @brief  Maximum matching matrix, reference struct 1 disallowed */
  unsigned int          *mm2;           /**<  @brief  Maximum matching matrix, reference struct 2 disallowed */

  int                   circ;
  int                   dangles;
  unsigned int          seq_length;

  FLT_OR_DBL            ***Q;
  FLT_OR_DBL            ***Q_B;
  FLT_OR_DBL            ***Q_M;
  FLT_OR_DBL            ***Q_M1;
  FLT_OR_DBL            ***Q_M2;

  FLT_OR_DBL            **Q_c;
  FLT_OR_DBL            **Q_cH;
  FLT_OR_DBL            **Q_cI;
  FLT_OR_DBL            **Q_cM;

  int                   **l_min_values;
  int                   **l_max_values;
  int                   *k_min_values;
  int                   *k_max_values;

  int                   **l_min_values_b;
  int                   **l_max_values_b;
  int                   *k_min_values_b;
  int                   *k_max_values_b;

  int                   **l_min_values_m;
  int                   **l_max_values_m;
  int                   *k_min_values_m;
  int                   *k_max_values_m;

  int                   **l_min_values_m1;
  int                   **l_max_values_m1;
  int                   *k_min_values_m1;
  int                   *k_max_values_m1;

  int                   **l_min_values_m2;
  int                   **l_max_values_m2;
  int                   *k_min_values_m2;
  int                   *k_max_values_m2;

  int                   *l_min_values_qc;
  int                   *l_max_values_qc;
  int                   k_min_values_qc;
  int                   k_max_values_qc;

  int                   *l_min_values_qcH;
  int                   *l_max_values_qcH;
  int                   k_min_values_qcH;
  int                   k_max_values_qcH;

  int                   *l_min_values_qcI;
  int                   *l_max_values_qcI;
  int                   k_min_values_qcI;
  int                   k_max_values_qcI;

  int                   *l_min_values_qcM;
  int                   *l_max_values_qcM;
  int                   k_min_values_qcM;
  int                   k_max_values_qcM;

  /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
  FLT_OR_DBL            *Q_rem;
  FLT_OR_DBL            *Q_B_rem;
  FLT_OR_DBL            *Q_M_rem;
  FLT_OR_DBL            *Q_M1_rem;
  FLT_OR_DBL            *Q_M2_rem;

  FLT_OR_DBL            Q_c_rem;
  FLT_OR_DBL            Q_cH_rem;
  FLT_OR_DBL            Q_cI_rem;
  FLT_OR_DBL            Q_cM_rem;

  vrna_fold_compound_t  *compatibility;
} TwoDpfold_vars;

/**
 * @brief Get a datastructure containing all necessary attributes and global folding switches
 *
 * This function prepares all necessary attributes and matrices etc which are needed for a call
 * of TwoDpfold() .
 * A snapshot of all current global model switches (dangles, temperature and so on) is done and
 * stored in the returned datastructure. Additionally, all matrices that will hold the partition
 * function values are prepared.
 *
 *  @deprecated Use the new API that relies on #vrna_fold_compound_t and the corresponding functions
 *              vrna_fold_compound_TwoD(), vrna_pf_TwoD(), and vrna_fold_compound_free() instead!
 *
 * @param seq         the RNA sequence in uppercase format with letters from the alphabet {AUCG}
 * @param structure1  the first reference structure in dot-bracket notation
 * @param structure2  the second reference structure in dot-bracket notation
 * @param circ        a switch indicating if the sequence is linear (0) or circular (1)
 * @returns           the datastructure containing all necessary partition function attributes
 */
DEPRECATED(TwoDpfold_vars *
           get_TwoDpfold_variables(const char *seq,
                                   const char *structure1,
                                   char       *structure2,
                                   int        circ));

/**
 * @brief Free all memory occupied by a TwoDpfold_vars datastructure
 *
 * This function free's all memory occupied by a datastructure obtained from from
 * get_TwoDpfold_variabless() or get_TwoDpfold_variables_from_MFE()
 *
 *  @deprecated Use the new API that relies on #vrna_fold_compound_t and the corresponding functions
 *              vrna_fold_compound_TwoD(), vrna_pf_TwoD(), and vrna_fold_compound_free() instead!
 *
 * @see get_TwoDpfold_variables(), get_TwoDpfold_variables_from_MFE()
 *
 * @param vars   the datastructure to be free'd
 */
DEPRECATED(void
           destroy_TwoDpfold_variables(TwoDpfold_vars *vars));

/**
 * @brief Compute the partition function for all distance classes
 *
 * This function computes the partition functions for all distance classes
 * according the two reference structures specified in the datastructure 'vars'.
 * Similar to TwoDfold() the arguments maxDistance1 and maxDistance2 specify
 * the maximum distance to both reference structures. A value of '-1' in either of
 * them makes the appropriate distance restrictionless, i.e. all basepair distancies
 * to the reference are taken into account during computation.
 * In case there is a restriction, the returned solution contains an entry where
 * the attribute k=l=-1 contains the partition function for all structures exceeding
 * the restriction.
 * A values of #INF in the attribute 'k' of the returned list denotes the end of the list
 *
 *  @deprecated Use the new API that relies on #vrna_fold_compound_t and the corresponding functions
 *              vrna_fold_compound_TwoD(), vrna_pf_TwoD(), and vrna_fold_compound_free() instead!
 *
 * @see get_TwoDpfold_variables(), destroy_TwoDpfold_variables(), #vrna_sol_TwoD_pf_t
 *
 * @param vars          the datastructure containing all necessary folding attributes and matrices
 * @param maxDistance1  the maximum basepair distance to reference1 (may be -1)
 * @param maxDistance2  the maximum basepair distance to reference2 (may be -1)
 * @returns             a list of partition funtions for the appropriate distance classes
 */
DEPRECATED(TwoDpfold_solution *
           TwoDpfoldList(TwoDpfold_vars *vars,
                         int            maxDistance1,
                         int            maxDistance2));

/**
 *  @brief Sample secondary structure representatives from a set of distance classes according to their
 *  Boltzmann probability
 *
 *  If the argument 'd1' is set to '-1', the structure will be backtracked in the distance class
 *  where all structures exceeding the maximum basepair distance to either of the references reside.
 *
 *  @pre      The argument 'vars' must contain precalculated partition function matrices,
 *            i.e. a call to TwoDpfold() preceding this function is mandatory!
 *
 *  @deprecated Use the new API that relies on #vrna_fold_compound_t and the corresponding functions
 *              vrna_fold_compound_TwoD(), vrna_pf_TwoD(), vrna_pbacktrack_TwoD(), and
 *              vrna_fold_compound_free() instead!
 *
 *  @see      TwoDpfold()
 *
 *  @param[in]  vars  the datastructure containing all necessary folding attributes and matrices
 *  @param[in]  d1    the distance to reference1 (may be -1)
 *  @param[in]  d2    the distance to reference2
 *  @returns    A sampled secondary structure in dot-bracket notation
 */
DEPRECATED(char *
           TwoDpfold_pbacktrack(TwoDpfold_vars  *vars,
                                int             d1,
                                int             d2));

/**
 * @brief Sample secondary structure representatives with a specified length from a set of distance classes according to their
 *  Boltzmann probability
 *
 * This function does essentially the same as TwoDpfold_pbacktrack() with the only difference that partial structures,
 * i.e. structures beginning from the 5' end with a specified length of the sequence, are backtracked
 *
 * @note      This function does not work (since it makes no sense) for circular RNA sequences!
 * @pre       The argument 'vars' must contain precalculated partition function matrices,
 *            i.e. a call to TwoDpfold() preceding this function is mandatory!
 *
 *  @deprecated Use the new API that relies on #vrna_fold_compound_t and the corresponding functions
 *              vrna_fold_compound_TwoD(), vrna_pf_TwoD(), vrna_pbacktrack5_TwoD(), and
 *              vrna_fold_compound_free() instead!
 *
 * @see       TwoDpfold_pbacktrack(), TwoDpfold()
 *
 *  @param[in]  vars    the datastructure containing all necessary folding attributes and matrices
 *  @param[in]  d1      the distance to reference1 (may be -1)
 *  @param[in]  d2      the distance to reference2
 *  @param[in]  length  the length of the structure beginning from the 5' end
 *  @returns    A sampled secondary structure in dot-bracket notation
 */
DEPRECATED(char *
           TwoDpfold_pbacktrack5(TwoDpfold_vars *vars,
                                 int            d1,
                                 int            d2,
                                 unsigned int   length));

/**
 * @brief
 *
 *
 */
DEPRECATED(FLT_OR_DBL **TwoDpfold(TwoDpfold_vars  *our_variables,
                                  int             maxDistance1,
                                  int             maxDistance2));

/**
 * @brief
 *
 *
 */
DEPRECATED(FLT_OR_DBL **TwoDpfold_circ(TwoDpfold_vars *our_variables,
                                       int            maxDistance1,
                                       int            maxDistance2));

#endif

#endif
