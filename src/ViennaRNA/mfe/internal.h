#ifndef VIENNA_RNA_PACKAGE_MFE_INTERNAL_H
#define VIENNA_RNA_PACKAGE_MFE_INTERNAL_H

#include <math.h>

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>

/**
 *  @file     ViennaRNA/mfe/internal.h
 *  @ingroup  mfe
 *  @brief    Energy evaluation of interior loops for MFE computations
 */

/**
 *  @addtogroup   eval_loops_int
 *  @{
 */


/**
 *  @name Minimum Free Energy API
 *  @{
 */

int
vrna_E_int_loop(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j);


int
vrna_E_ext_int_loop(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   *ip,
                    int                   *iq);


int
vrna_E_stack(vrna_fold_compound_t *fc,
             int                  i,
             int                  j);


/* End basic interface */
/**@}*/


/**
 * @}
 */


#endif
