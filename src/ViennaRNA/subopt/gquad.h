#ifndef VIENNA_RNA_PACKAGE_SUBOPT_GQUAD_H
#define VIENNA_RNA_PACKAGE_SUBOPT_GQUAD_H

#include "ViennaRNA/datastructures/array.h"
#include <ViennaRNA/fold_compound.h>

/**
 *  @file       ViennaRNA/subopt/gquad.h
 *  @ingroup    gquads
 *  @brief      G-quadruplex suboptimals
 */

/**
 *  @addtogroup gquad_eval
 *  @{
 */

vrna_array(int)
vrna_gq_int_loop_subopt(vrna_fold_compound_t * fc,
                        unsigned int i,
                        unsigned int j,
                        vrna_array(int) * p_p,
                        vrna_array(int) * q_p,
                        int threshold);
/**
 *  @}
 */

#endif
