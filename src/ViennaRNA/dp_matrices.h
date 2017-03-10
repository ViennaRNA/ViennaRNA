#ifndef VIENNA_RNA_PACKAGE_DP_MATRICES_H
#define VIENNA_RNA_PACKAGE_DP_MATRICES_H

/**
 *  @file     dp_matrices.h
 *  @ingroup  data_structures
 *  @brief    Functions to deal with standard dynamic programming (DP) matrices
 */

/**
 *  @addtogroup dp_matrices   The Dynamic Programming Matrices
 *  @brief  This module provides interfaces that deal with creation and destruction of
 *          dynamic programming matrices used within the RNAlib.
 *
 *  @{
 *  @ingroup  dp_matrices
 */

/** @brief Typename for the Minimum Free Energy (MFE) DP matrices data structure #vrna_mx_mfe_s */
typedef struct  vrna_mx_mfe_s vrna_mx_mfe_t;
/** @brief Typename for the Partition Function (PF) DP matrices data structure #vrna_mx_pf_s */
typedef struct  vrna_mx_pf_s vrna_mx_pf_t;

#include <ViennaRNA/data_structures.h>

/**
 *  @brief  An enumerator that is used to specify the type of a polymorphic Dynamic Programming (DP)
 *  matrix data structure
 *  @see #vrna_mx_mfe_t, #vrna_mx_pf_t
 */
typedef enum {
  VRNA_MX_DEFAULT,  /**<  @brief  Default DP matrices */
  VRNA_MX_WINDOW,   /**<  @brief  DP matrices suitable for local structure prediction using
                     *    window approach.
                     *    @see    vrna_mfe_window(), vrna_mfe_window_zscore(), pfl_fold()
                     */
  VRNA_MX_2DFOLD    /**<  @brief  DP matrices suitable for distance class partitioned structure prediction
                     *    @see  vrna_mfe_TwoD(), vrna_pf_TwoD()
                     */
} vrna_mx_type_e;

/**
 *  @brief  Minimum Free Energy (MFE) Dynamic Programming (DP) matrices data structure required within the #vrna_fold_compound_t
 */
struct vrna_mx_mfe_s {
  /** @name Common fields for MFE matrices
   *  @{
   */
  vrna_mx_type_e  type;
  unsigned int    length;  /**<  @brief  Length of the sequence, therefore an indicator of the size of the DP matrices */
  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif
  /** @name Default DP matrices
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_DEFAULT @endcode
   * @{
   */
  int *c;           /**<  @brief  Energy array, given that i-j pair */
  int *f5;          /**<  @brief  Energy of 5' end */
  int *f3;          /**<  @brief  Energy of 3' end */
  int *fc;          /**<  @brief  Energy from i to cutpoint (and vice versa if i>cut) */
  int *fML;         /**<  @brief  Multi-loop auxiliary energy array */
  int *fM1;         /**<  @brief  Second ML array, only for unique multibrnach loop decomposition */
  int *fM2;         /**<  @brief  Energy for a multibranch loop region with exactly two stems, extending to 3' end */
  int *ggg;         /**<  @brief  Energies of g-quadruplexes */
  int Fc;           /**<  @brief  Minimum Free Energy of entire circular RNA */
  int FcH;
  int FcI;
  int FcM;
  /**
   * @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif
  /** @name Local Folding DP matrices using window approach
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_WINDOW @endcode
   * @{
   */
  int **c_local;            /**<  @brief  Energy array, given that i-j pair */
  int *f3_local;            /**<  @brief  Energy of 5' end */
  int **fML_local;          /**<  @brief  Multi-loop auxiliary energy array */
  int **ggg_local;          /**<  @brief  Energies of g-quadruplexes */
  /**
   * @}
   */
#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /** @name Distance Class DP matrices
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_2DFOLD @endcode
   * @{
   */
  int           ***E_F5;
  int           **l_min_F5;
  int           **l_max_F5;
  int           *k_min_F5;
  int           *k_max_F5;

  int           ***E_F3;
  int           **l_min_F3;
  int           **l_max_F3;
  int           *k_min_F3;
  int           *k_max_F3;

  int           ***E_C;
  int           **l_min_C;
  int           **l_max_C;
  int           *k_min_C;
  int           *k_max_C;

  int           ***E_M;
  int           **l_min_M;
  int           **l_max_M;
  int           *k_min_M;
  int           *k_max_M;

  int           ***E_M1;
  int           **l_min_M1;
  int           **l_max_M1;
  int           *k_min_M1;
  int           *k_max_M1;

  int           ***E_M2;
  int           **l_min_M2;
  int           **l_max_M2;
  int           *k_min_M2;
  int           *k_max_M2;

  int           **E_Fc;
  int           *l_min_Fc;
  int           *l_max_Fc;
  int           k_min_Fc;
  int           k_max_Fc;

  int           **E_FcH;
  int           *l_min_FcH;
  int           *l_max_FcH;
  int           k_min_FcH;
  int           k_max_FcH;

  int           **E_FcI;
  int           *l_min_FcI;
  int           *l_max_FcI;
  int           k_min_FcI;
  int           k_max_FcI;

  int           **E_FcM;
  int           *l_min_FcM;
  int           *l_max_FcM;
  int           k_min_FcM;
  int           k_max_FcM;

  /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
  int           *E_F5_rem;
  int           *E_F3_rem;
  int           *E_C_rem;
  int           *E_M_rem;
  int           *E_M1_rem;
  int           *E_M2_rem;

  int           E_Fc_rem;
  int           E_FcH_rem;
  int           E_FcI_rem;
  int           E_FcM_rem;

#ifdef COUNT_STATES
  unsigned long ***N_F5;
  unsigned long ***N_C;
  unsigned long ***N_M;
  unsigned long ***N_M1;
#endif

  /**
   * @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
};
#endif
};

/**
 *  @brief  Partition function (PF) Dynamic Programming (DP) matrices data structure required within the #vrna_fold_compound_t
 */
struct vrna_mx_pf_s {
  /** @name Common fields for DP matrices
   *  @{
   */
  vrna_mx_type_e type;
  unsigned int length;
  FLT_OR_DBL *scale;
  FLT_OR_DBL *expMLbase;


  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif

  /** @name Default PF matrices
   *  @note These data fields are available if
   *        @code vrna_mx_pf_t.type == VRNA_MX_DEFAULT @endcode
   *  @{
   */
  FLT_OR_DBL *q;
  FLT_OR_DBL *qb;
  FLT_OR_DBL *qm;
  FLT_OR_DBL *qm1;
  FLT_OR_DBL *probs;
  FLT_OR_DBL *q1k;
  FLT_OR_DBL *qln;
  FLT_OR_DBL *G;

  FLT_OR_DBL qo;
  FLT_OR_DBL *qm2;
  FLT_OR_DBL qho;
  FLT_OR_DBL qio;
  FLT_OR_DBL qmo;

  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /** @name Local Folding DP matrices using window approach
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_WINDOW @endcode
   * @{
   */
  FLT_OR_DBL **q_local;
  FLT_OR_DBL **qb_local;
  FLT_OR_DBL **qm_local;
  FLT_OR_DBL **pR;
  FLT_OR_DBL **qm2_local;
  FLT_OR_DBL **QI5;
  FLT_OR_DBL **q2l;
  FLT_OR_DBL **qmb;
  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /** @name Distance Class DP matrices
   *  @note These data fields are available if
   *        @code vrna_mx_pf_t.type == VRNA_MX_2DFOLD @endcode
   *  @{
   */
  FLT_OR_DBL ***Q;
  int **l_min_Q;
  int **l_max_Q;
  int *k_min_Q;
  int *k_max_Q;


  FLT_OR_DBL ***Q_B;
  int **l_min_Q_B;
  int **l_max_Q_B;
  int *k_min_Q_B;
  int *k_max_Q_B;

  FLT_OR_DBL ***Q_M;
  int **l_min_Q_M;
  int **l_max_Q_M;
  int *k_min_Q_M;
  int *k_max_Q_M;

  FLT_OR_DBL ***Q_M1;
  int **l_min_Q_M1;
  int **l_max_Q_M1;
  int *k_min_Q_M1;
  int *k_max_Q_M1;

  FLT_OR_DBL ***Q_M2;
  int **l_min_Q_M2;
  int **l_max_Q_M2;
  int *k_min_Q_M2;
  int *k_max_Q_M2;

  FLT_OR_DBL **Q_c;
  int *l_min_Q_c;
  int *l_max_Q_c;
  int k_min_Q_c;
  int k_max_Q_c;

  FLT_OR_DBL **Q_cH;
  int *l_min_Q_cH;
  int *l_max_Q_cH;
  int k_min_Q_cH;
  int k_max_Q_cH;

  FLT_OR_DBL **Q_cI;
  int *l_min_Q_cI;
  int *l_max_Q_cI;
  int k_min_Q_cI;
  int k_max_Q_cI;

  FLT_OR_DBL **Q_cM;
  int *l_min_Q_cM;
  int *l_max_Q_cM;
  int k_min_Q_cM;
  int k_max_Q_cM;

  /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
  FLT_OR_DBL *Q_rem;
  FLT_OR_DBL *Q_B_rem;
  FLT_OR_DBL *Q_M_rem;
  FLT_OR_DBL *Q_M1_rem;
  FLT_OR_DBL *Q_M2_rem;

  FLT_OR_DBL Q_c_rem;
  FLT_OR_DBL Q_cH_rem;
  FLT_OR_DBL Q_cI_rem;
  FLT_OR_DBL Q_cM_rem;
  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
};
#endif
};

/**
 *  @brief  Add Dynamic Programming (DP) matrices (allocate memory)
 *
 *  This function adds DP matrices of a specific type to the provided
 *  #vrna_fold_compound_t, such that successive DP recursion can be applied.
 *  The function caller has to specify which type of DP matrix is requested,
 *  see #vrna_mx_type_e, and what kind of recursive algorithm will be applied
 *  later on, using the parameters type, and options, respectively. For the
 *  latter, Minimum free energy (MFE), and Partition function (PF)
 *  computations are distinguished. A third option that may be passed
 *  is #VRNA_OPTION_HYBRID, indicating that auxiliary DP arrays are
 *  required for RNA-RNA interaction prediction.
 *
 *  @note Usually, there is no need to call this function, since
 *  the constructors of #vrna_fold_compound_t are handling all the DP
 *  matrix memory allocation.
 *
 *  @see vrna_mx_mfe_add(), vrna_mx_pf_add(), vrna_fold_compound(),
 *  vrna_fold_compound_comparative(), vrna_fold_compound_free(),
 *  vrna_mx_pf_free(), vrna_mx_mfe_free(), #vrna_mx_type_e,
 *  #VRNA_OPTION_MFE, #VRNA_OPTION_PF, #VRNA_OPTION_HYBRID, #VRNA_OPTION_EVAL_ONLY
 *
 *  @param  vc      The #vrna_fold_compound_t that holds pointers to the DP matrices
 *  @param  type    The type of DP matrices requested
 *  @param  options Option flags that specify the kind of DP matrices, such
 *                  as MFE or PF arrays, and auxiliary requirements
 *  @returns        1 if DP matrices were properly allocated and attached,
 *                  0 otherwise
 */
int
vrna_mx_add(vrna_fold_compound_t  *vc,
            vrna_mx_type_e        type,
            unsigned int          options);


int
vrna_mx_mfe_add(vrna_fold_compound_t  *vc,
                vrna_mx_type_e        mx_type,
                unsigned int          options);


int
vrna_mx_pf_add(vrna_fold_compound_t *vc,
               vrna_mx_type_e       mx_type,
               unsigned int         options);


int
vrna_mx_prepare(vrna_fold_compound_t  *vc,
                unsigned int          options);


/**
 *  @brief  Free memory occupied by the Minimum Free Energy (MFE) Dynamic Programming (DP) matrices
 *
 *  @see vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_fold_compound_free(), vrna_mx_pf_free()
 *
 *  @param  vc  The #vrna_fold_compound_t storing the MFE DP matrices that are to be erased from memory
 */
void
vrna_mx_mfe_free(vrna_fold_compound_t *vc);


/**
 *  @brief  Free memory occupied by the Partition Function (PF) Dynamic Programming (DP) matrices
 *
 *  @see vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_fold_compound_free(), vrna_mx_mfe_free()
 *
 *  @param  vc  The #vrna_fold_compound_t storing the PF DP matrices that are to be erased from memory
 */
void
vrna_mx_pf_free(vrna_fold_compound_t *vc);


/**
 *  @}
 */

#endif
