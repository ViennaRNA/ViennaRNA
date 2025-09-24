#ifndef VIENNA_RNA_PACKAGE_PROBING_TRANSFORM_H
#define VIENNA_RNA_PACKAGE_PROBING_TRANSFORM_H

#include <ViennaRNA/fold_compound.h>

/**
 *  @file     ViennaRNA/probing/transform.h
 *  @ingroup  probing_data
 *  @brief    This module provides function to transform linear data
 */

/**
 *  @addtogroup probing_data
 *  @{
 */

typedef   void                  *vrna_data_lin_trans_opt_t;
typedef   vrna_auxdata_free_f   vrna_data_lin_trans_opt_free_f;
typedef   double                (*vrna_data_lin_trans_f) (double,
                                                          vrna_data_lin_trans_opt_t);


double *
vrna_data_lin_transform(const double              *data,
                        size_t                    data_size,
                        vrna_data_lin_trans_f     transform_cb,
                        vrna_data_lin_trans_opt_t transform_opt);


#define VRNA_TRANSFORM_BIN_OPTION_PROJECT               (1 << 0)
#define VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_UPPERBOUND  (1 << 1)
#define VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_LOWERBOUND  (1 << 2)
#define VRNA_TRANSFORM_BIN_OPTION_DEFAULT               0

#define VRNA_TRANSFORM_LM_OPTION_LOG                    (1 << 0)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_LOW        (1 << 1)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_HIGH       (1 << 2)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_LOW        (1 << 3)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_HIGH       (1 << 4)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE            (VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_LOW | VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_HIGH)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET            (VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_LOW | VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_HIGH)
#define VRNA_TRANSFORM_LM_OPTION_CLIP                   (VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_LOW | VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_HIGH | VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_LOW | VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_HIGH)
#define VRNA_TRANSFORM_LM_OPTION_DEFAULT                (VRNA_TRANSFORM_LM_OPTION_CLIP)


vrna_data_lin_trans_f
vrna_data_transform_method_bin(double                         (*thresholds)[2],
                               unsigned int                   thresholds_num,
                               double                         oolb_value,
                               double                         ooub_value,
                               unsigned char                  options,
                               vrna_data_lin_trans_opt_t      *transform_options_p,
                               vrna_data_lin_trans_opt_free_f *transform_options_free);


vrna_data_lin_trans_f
vrna_data_transform_method_lm(double                          slope,
                              double                          intercept,
                              double                          domain[4],
                              double                          oob_value,
                              unsigned char                   options,
                              vrna_data_lin_trans_opt_t       *transform_options_p,
                              vrna_data_lin_trans_opt_free_f  *transform_options_free);


vrna_data_lin_trans_f
vrna_data_transform_method_log(double                         value_shift,
                               double                         oob_value,
                               vrna_data_lin_trans_opt_t      *transform_options_p,
                               vrna_data_lin_trans_opt_free_f *transform_options_free);


/**
 *  @}
 */

#endif
