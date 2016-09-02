#ifndef VIENNA_RNA_PACKAGE_GRAMMAR_H
#define VIENNA_RNA_PACKAGE_GRAMMAR_H

/**
 *  @file     grammar.h
 *  @ingroup  grammar
 *  @brief    Implementations for the RNA folding grammar
 */

#include <ViennaRNA/data_structures.h>

typedef void (vrna_callback_gr_rule_aux)(vrna_fold_compound_t *vc, int i, int j, void *data);

typedef void (vrna_callback_gr_free_auxdata)(void *data);

typedef struct vrna_gr_aux_s  vrna_gr_aux_t;

struct vrna_gr_aux_s {

  vrna_callback_gr_rule_aux     *cb_aux_f;
  vrna_callback_gr_rule_aux     *cb_aux_c;
  vrna_callback_gr_rule_aux     *cb_aux_m;
  vrna_callback_gr_rule_aux     *cb_aux_m1;
  vrna_callback_gr_rule_aux     *cb_aux;

  void                          *auxdata;
  vrna_callback_gr_free_auxdata *free_auxdata;
};


#endif
