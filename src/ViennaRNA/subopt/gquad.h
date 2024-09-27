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
vrna_gq_int_loop_subopt(vrna_fold_compound_t      *fc,
                        unsigned int              i,
                        unsigned int              j,
                        vrna_array(unsigned int)  *p_p,
                        vrna_array(unsigned int)  *q_p,
                        int                       threshold);

void
get_gquad_pattern_exhaustive(short        *S,
                             unsigned int i,
                             unsigned int j,
                             vrna_param_t *P,
                             unsigned int *L,
                             unsigned int *l,
                             int          threshold);



unsigned int
get_gquad_count(short         *S,
                unsigned int  i,
                unsigned int  j);

/**
 *  @}
 */

#endif
