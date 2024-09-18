#ifndef VIENNA_RNA_PACKAGE_MFE_MULTIBRANCH_H
#define VIENNA_RNA_PACKAGE_MFE_MULTIBRANCH_H

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>

/**
 *  @file     ViennaRNA/mfe/multibranch.h
 *  @ingroup  mfe
 *  @brief    Energy evaluation of multibranch loops for MFE
 */

/**
 *  @addtogroup  eval_loops_mb
 *  @{
 */


/**
 *  @name Minimum Free Energy API
 *  @{
 */


/**
 *  @brief Evaluate energy of multi branch loop helices stacking onto closing pair (i,j)
 *
 *  Computes total free energy for coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1)
 */
int
vrna_E_mb_loop_stack(vrna_fold_compound_t *fc,
                     int                  i,
                     int                  j);


int
vrna_E_mb_loop_fast(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   *dmli1,
                    int                   *dmli2);


int
E_ml_rightmost_stem(int                   i,
                    int                   j,
                    vrna_fold_compound_t  *fc);


int
vrna_E_ml_stems_fast(vrna_fold_compound_t *fc,
                     int                  i,
                     int                  j,
                     int                  *fmi,
                     int                  *dmli);


/* End basic interface */
/**@}*/


/**
 * @}
 */

#endif
