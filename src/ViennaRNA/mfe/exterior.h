#ifndef VIENNA_RNA_PACKAGE_MFE_EXTERIOR_H
#define VIENNA_RNA_PACKAGE_MFE_EXTERIOR_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>
/**
 *  @file     ViennaRNA/mfe/exterior.h
 *  @ingroup  mfe
 *  @brief    Energy evaluation of exterior loops for MFE computations
 */

/**
 *  @addtogroup   eval_loops_ext
 *  @{
 */


/**
 *  @name Minimum Free Energy API
 *  @{
 */


int
vrna_E_ext_loop_5(vrna_fold_compound_t *fc);


int
vrna_E_ext_loop_3(vrna_fold_compound_t  *fc,
                  int                   i);


/* End basic interface */
/**@}*/


/**
 * @}
 */

#endif
