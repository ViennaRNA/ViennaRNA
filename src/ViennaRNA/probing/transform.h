#ifndef VIENNA_RNA_PACKAGE_PROBING_TRANSFORM_H
#define VIENNA_RNA_PACKAGE_PROBING_TRANSFORM_H

#include <ViennaRNA/fold_compound.h>

/**
 *  @file     ViennaRNA/probing/transform.h
 *  @ingroup  probing_data
 *  @brief    This module provides function to transform probing data
 */

/**
 *  @addtogroup probing_data
 *  @{
 */

typedef   double (*vrna_probing_transform_f) (double,
                                              void *);


#define VRNA_REACTIVITY_TRANS_DEFAULT                             0U
#define VRNA_REACTIVITY_TRANS_IDEN                                1U
#define VRNA_REACTIVITY_TRANS_NEG_IGNORE                          2U
#define VRNA_REACTIVITY_TRANS_NEG_ZERO                            3U
#define VRNA_REACTIVITY_TRANS_LOG1P                               4U


#define VRNA_REACTIVITY_MISSING                                   -999.


vrna_probing_transform_f
vrna_reactivity_trans_default(unsigned int flag);


vrna_probing_transform_f
vrna_reactivity_trans_method(unsigned int flag);


FLT_OR_DBL *
vrna_reactivity_transform(unsigned int n,
                          const double *reactivity,
                          vrna_probing_transform_f transform_cb,
                          void *options);

/**
 *  @}
 */

#endif
