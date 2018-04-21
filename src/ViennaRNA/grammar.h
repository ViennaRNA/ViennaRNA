#ifndef VIENNA_RNA_PACKAGE_GRAMMAR_H
#define VIENNA_RNA_PACKAGE_GRAMMAR_H

/**
 *  @file     grammar.h
 *  @ingroup  grammar
 *  @brief    Implementations for the RNA folding grammar
 */

/**
 * @addtogroup grammar
 * @{
 *    @brief  The RNA folding grammar as implemented in RNAlib
*/

#include <ViennaRNA/data_structures.h>

typedef int (vrna_callback_gr_rule)(vrna_fold_compound_t *vc, int i, int j, void *data);
typedef FLT_OR_DBL (vrna_callback_gr_rule_exp)(vrna_fold_compound_t *vc, int i, int j, void *data);

typedef void (vrna_callback_gr_free_auxdata)(void *data);

typedef struct vrna_gr_aux_s  vrna_gr_aux_t;

struct vrna_gr_aux_s {

  vrna_callback_gr_rule *cb_aux_f;
  vrna_callback_gr_rule *cb_aux_c;
  vrna_callback_gr_rule *cb_aux_m;
  vrna_callback_gr_rule *cb_aux_m1;
  vrna_callback_gr_rule *cb_aux;

  vrna_callback_gr_rule_exp *cb_aux_exp_f;
  vrna_callback_gr_rule_exp *cb_aux_exp_c;
  vrna_callback_gr_rule_exp *cb_aux_exp_m;
  vrna_callback_gr_rule_exp *cb_aux_exp_m1;
  vrna_callback_gr_rule_exp *cb_aux_exp;

  void                          *auxdata;
  vrna_callback_gr_free_auxdata *free_auxdata;
};

/**
 * @}
 */

/**
 *  @addtogroup domains
 *  @brief This module covers simple and straight-forward extensions to the RNA folding grammar.
*/

#endif
