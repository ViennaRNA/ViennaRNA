#ifndef VIENNA_RNA_PACKAGE_PARTFUNC_GQUAD_H
#define VIENNA_RNA_PACKAGE_PARTFUNC_GQUAD_H

#include <string.h>

#include "ViennaRNA/datastructures/sparse_mx.h"
#include <ViennaRNA/fold_compound.h>

/**
 *  @file       ViennaRNA/partfunc/gquad.h
 *  @ingroup    gquads
 *  @brief      G-quadruplex Partition Function API
 */

/**
 *  @addtogroup gquad_eval
 *  @{
 */
FLT_OR_DBL
vrna_gq_int_loop_pf(vrna_fold_compound_t  *fc,
                    unsigned int          i,
                    unsigned int          j);


/**
 *  @}
 */


/**
 *  @addtogroup gquad_dp
 *  @{
 */

vrna_smx_csr_FLT_OR_DBL_t *
vrna_gq_pos_pf(vrna_fold_compound_t *fc);


/**
 *  @}
 */

/**
 *  @addtogroup gquad_parse
 *  @{
 */
void
get_gquad_pattern_pf(short            *S,
                     int              i,
                     int              j,
                     vrna_exp_param_t *pf,
                     int              *L,
                     int              l[3]);


void
vrna_get_gquad_pattern_pf(vrna_fold_compound_t  *fc,
                          unsigned int          i, 
                          unsigned int          j,
                          unsigned int          *L,
                          unsigned int          [3]);

/**
 *  @}
 */


/**
 *  @addtogroup gquad_other
 *  @{
 */
plist *
  get_plist_gquad_from_pr(short *S,
                          int gi,
                          int gj,
                          vrna_smx_csr(FLT_OR_DBL) * q_gq,
                          FLT_OR_DBL * probs,
                          FLT_OR_DBL * scale,
                          vrna_exp_param_t * pf);


vrna_ep_t *
vrna_plist_gquad_from_pr(vrna_fold_compound_t *fc,
                         int                  gi,
                         int                  gj);


plist *
  get_plist_gquad_from_pr_max(short *S,
                              int gi,
                              int gj,
                              vrna_smx_csr(FLT_OR_DBL) * q_gq,
                              FLT_OR_DBL * probs,
                              FLT_OR_DBL * scale,
                              int *L,
                              int l[3],
                              vrna_exp_param_t * pf);


vrna_ep_t *
vrna_plist_gquad_from_pr_max(vrna_fold_compound_t *fc,
                             unsigned int                  gi,
                             unsigned int                  gj,
                             unsigned int                  *Lmax,
                             unsigned int                  lmax[3]);


/**
 *  @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup gquad_deprecated
 *  @{
 */


DEPRECATED(FLT_OR_DBL * get_gquad_pf_matrix(short *S,
                                            FLT_OR_DBL * scale,
                                            vrna_exp_param_t * pf),
           "Use vrna_gq_pos_pf() instead");


DEPRECATED(FLT_OR_DBL * get_gquad_pf_matrix_comparative(unsigned int n,
                                                        short *S_cons,
                                                        short **S,
                                                        unsigned int **a2s,
                                                        FLT_OR_DBL * scale,
                                                        unsigned int n_seq,
                                                        vrna_exp_param_t * pf),
           "Use vrna_gq_pos_pf() instead");

/**
 * @}
 */


#endif

#endif
