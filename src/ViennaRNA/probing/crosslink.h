#ifndef VIENNA_RNA_PACKAGE_PROBING_XLINK_H
#define VIENNA_RNA_PACKAGE_PROBING_XLINK_H

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/probing/basic.h>

/**
 *  @file     ViennaRNA/probing/crosslink.h
 *  @ingroup probing_data
 *  @brief  This module provides function to incorporate RNA structure
 *          probing data derived from cross-linking experiments into the folding recursions
 *          by means of soft constraints and grammar extensions
 */


/**
 *  @ingroup probing_data
 */
int
vrna_sc_xl_prepare(vrna_fold_compound_t *fc);

int
vrna_mfe_multibranch_loop_xlink(vrna_fold_compound_t         *fc,
                                unsigned int                 i,
                                unsigned int                 j,
                                struct vrna_mx_mfe_aux_ml_s  *helpers);

void
vrna_sc_xl_update_mfe(vrna_fold_compound_t *fc,
                      unsigned int         i);

#endif
