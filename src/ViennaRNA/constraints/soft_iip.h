#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_IIP_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_IIP_H

#include <ViennaRNA/fold_compound.h>
#include "ViennaRNA/mfe/multibranch.h"
#include "ViennaRNA/partfunc/multibranch.h"
#include <ViennaRNA/probabilities/basic.h>
#include <ViennaRNA/probing/basic.h>

/**
 *  @file     ViennaRNA/constraints/soft_iip.h
 *  @ingroup constraints
 *  @brief  This module provides function to incorporate RNA structure
 *          probing data derived from cross-linking experiments into the folding recursions
 *          by means of soft constraints and grammar extensions
 */

typedef struct vrna_crosslink_data_s  *vrna_sc_iip_t;


/**
 *  @ingroup probing_data
 */
int
vrna_sc_iip(vrna_fold_compound_t *fc);


int
vrna_sc_add_iip(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          j,
                unsigned int          k,
                unsigned int          l,
                double                energy,
                unsigned int          options);

int
vrna_sc_add_iip_file(vrna_fold_compound_t  *fc,
                     const char            *filename);


int
vrna_sc_iip_prepare(vrna_fold_compound_t *fc);


int
vrna_mfe_multibranch_loop_xlink(vrna_fold_compound_t        *fc,
                                unsigned int                i,
                                unsigned int                j,
                                struct vrna_mx_mfe_aux_ml_s *helpers);


unsigned int
vrna_bt_multibranch_loop_xlink(vrna_fold_compound_t *fc,
                               unsigned int         i,
                               unsigned int         j,
                               int                  en,
                               vrna_bps_t           bp_stack,
                               vrna_bts_t           bt_stack);


void
vrna_sc_iip_update_mfe(vrna_fold_compound_t  *fc,
                      unsigned int          i);


FLT_OR_DBL
vrna_pf_multibranch_loop_xlink(vrna_fold_compound_t *fc,
                               unsigned int         i,
                               unsigned int         j,
                               vrna_mx_pf_aux_ml_t  helpers);


void
vrna_sc_iip_update_pf(vrna_fold_compound_t *fc,
                     unsigned int         j);


vrna_eqm_probs_component_t
vrna_sc_iip_get_outside_cb(vrna_fold_compound_t  *fc);

#endif
