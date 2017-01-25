#ifndef VIENNA_RNA_PACKAGE_WALK_H_
#define VIENNA_RNA_PACKAGE_WALK_H_

/**
 *  @file     walk.h
 *  @ingroup  walks
 *  @brief    Methods to compute walks through the energy landscape of an RNA sequence
 */

/**
 *  @addtogroup   gradient_walks
 *  @brief Methods to obtain different features of gradient walks (final structure, single moves etc.)
 *
 *  @{
 *  @ingroup  gradient_walks
 */

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/neighbor.h>


/**
 * @brief Computes a steepest descent walk. The result is written to the input structure as pairtable.
 * (use vrna_ptable_copy(pt) if you want to keep the initial structure)
 * The path will be written to the move list.
 */
void
vrna_walk_gradient(vrna_fold_compound_t *vc,
                   short                *ptStartAndResultStructure,
                   int                  shifts);


move *
vrna_walk_gradient_path(vrna_fold_compound_t  *vc,
                        short                 *ptStartAndResultStructure,
                        int                   shifts);


void
vrna_walk_random(vrna_fold_compound_t *vc,
                 short                *ptStartAndResultStructure,
                 unsigned int         steps,
                 int                  shifts);


move *
vrna_walk_random_path(vrna_fold_compound_t  *vc,
                      short                 *ptStartAndResultStructure,
                      unsigned int          steps,
                      int                   shifts);


/**
 *  @}
 */
#endif /* VIENNA_RNA_PACKAGE_WALK_H_ */
