#ifndef VRNA_INTERN_GRAMMAR_DAT_H
#define VRNA_INTERN_GRAMMAR_DAT_H

#include <ViennaRNA/utils/basic.h>

#include "ViennaRNA/grammar/basic.h"
#include "ViennaRNA/grammar/mfe.h"
#include "ViennaRNA/grammar/partfunc.h"

typedef struct {
  vrna_gr_inside_f        cb;
  vrna_gr_outside_f       cb_bt;
  void                    *data;
  vrna_auxdata_prepare_f  prepare;
  vrna_auxdata_free_f     release;
} grammar_rule_wrap_t;


typedef struct {
  vrna_gr_inside_exp_f    cb;
  vrna_gr_outside_exp_f   cb_out;
  void                    *data;
  vrna_auxdata_prepare_f  prepare;
  vrna_auxdata_free_f     release;
} grammar_rule_exp_wrap_t;


struct vrna_gr_aux_s {
  vrna_array(grammar_rule_wrap_t) f;
  vrna_array(grammar_rule_wrap_t) c;
  vrna_array(grammar_rule_wrap_t) m;
  vrna_array(grammar_rule_wrap_t) m1;
  vrna_array(grammar_rule_wrap_t) m2;
  vrna_array(grammar_rule_wrap_t) aux;

  vrna_array(grammar_rule_exp_wrap_t) exp_f;
  vrna_array(grammar_rule_exp_wrap_t) exp_c;
  vrna_array(grammar_rule_exp_wrap_t) exp_m;
  vrna_array(grammar_rule_exp_wrap_t) exp_m1;
  vrna_array(grammar_rule_exp_wrap_t) exp_m2;
  vrna_array(grammar_rule_exp_wrap_t) exp_aux;

  vrna_array(vrna_recursion_status_f) cbs_status;
  vrna_array(void *)                  datas;
  vrna_array(vrna_auxdata_prepare_f)  prepare_datas;
  vrna_array(vrna_auxdata_free_f)     free_datas;

  vrna_gr_serialize_bp_f  serialize_bp;
  void                    *serialize_bp_data;
  vrna_auxdata_prepare_f  serialize_bp_prepare_data;
  vrna_auxdata_free_f     serialize_bp_free_data;
};


PRIVATE void
init_aux_grammar(vrna_fold_compound_t *fc) VRNA_UNUSED;


PRIVATE void
init_aux_grammar(vrna_fold_compound_t *fc)
{
  fc->aux_grammar = (struct vrna_gr_aux_s *)vrna_alloc(sizeof(struct vrna_gr_aux_s));

  if (fc->aux_grammar) {
    vrna_array_init(fc->aux_grammar->f);
    vrna_array_init(fc->aux_grammar->c);
    vrna_array_init(fc->aux_grammar->m);
    vrna_array_init(fc->aux_grammar->m1);
    vrna_array_init(fc->aux_grammar->m2);
    vrna_array_init(fc->aux_grammar->aux);

    vrna_array_init(fc->aux_grammar->exp_f);
    vrna_array_init(fc->aux_grammar->exp_c);
    vrna_array_init(fc->aux_grammar->exp_m);
    vrna_array_init(fc->aux_grammar->exp_m1);
    vrna_array_init(fc->aux_grammar->exp_m2);
    vrna_array_init(fc->aux_grammar->exp_aux);

    vrna_array_init(fc->aux_grammar->cbs_status);
    vrna_array_init(fc->aux_grammar->datas);
    vrna_array_init(fc->aux_grammar->free_datas);
    vrna_array_init(fc->aux_grammar->prepare_datas);

    fc->aux_grammar->serialize_bp               = NULL;
    fc->aux_grammar->serialize_bp_data          = NULL;
    fc->aux_grammar->serialize_bp_prepare_data  = NULL;
    fc->aux_grammar->serialize_bp_free_data     = NULL;
  }
}

#endif
