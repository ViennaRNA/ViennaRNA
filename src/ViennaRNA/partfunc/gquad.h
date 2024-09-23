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
