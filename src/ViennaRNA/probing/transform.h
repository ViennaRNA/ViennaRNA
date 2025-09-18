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


#define VRNA_TRANSFORM_BIN_OPTION_PROJECT               (1 << 0)
#define VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_UPPERBOUND  (1 << 1)
#define VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_LOWERBOUND  (1 << 2)

#define VRNA_TRANSFORM_LM_OPTION_LOG                    (1 << 0)


vrna_probing_transform_f
vrna_data_transform_method_bin(const double         **thresholds,
                               unsigned int         thresholds_num,
                               double               oolb_value,
                               double               ooub_value,
                               unsigned char        options,
                               void                 **transform_options_p,
                               vrna_auxdata_free_f  *transform_options_free);


PUBLIC vrna_probing_transform_f
vrna_data_transform_method_lm(double               slope,
                              double               intercept,
                              double               oob_value,
                              unsigned char        options,
                              void                 **transform_options_p,
                              vrna_auxdata_free_f  *transform_options_free);


/**
 *  @}
 */

#endif
