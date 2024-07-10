#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/grammar.h"

#include "ViennaRNA/grammar.inc"

PRIVATE void
add_aux_grammar(vrna_fold_compound_t *fc);


PUBLIC int
vrna_gr_prepare(vrna_fold_compound_t  *fc,
                unsigned int          options)
{
  int ret = 1;

  if ((fc) &&
      (fc->aux_grammar) &&
      (fc->aux_grammar->prepare_data))
    ret &= fc->aux_grammar->prepare_data(fc,
                                         fc->aux_grammar->data,
                                         options,
                                         NULL);

  return ret;
}


PUBLIC unsigned int
vrna_gr_register_bt_flag(vrna_fold_compound_t *fc)
{
  unsigned int flag = 0;

  if ((fc) &&
      (fc->aux_grammar)) {
    flag = 1 + 
           VRNA_MX_FLAG_MAX + 
           vrna_array_size(fc->aux_grammar->bt_flags);

    vrna_array_append(fc->aux_grammar->bt_flags,
                      flag);
  }

  return flag;
}


PUBLIC int
vrna_gr_set_aux_f(vrna_fold_compound_t  *fc,
                  vrna_grammar_rule_f   cb,
                  vrna_grammar_bt_f     cb_bt)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_f = cb;
    fc->aux_grammar->cb_bt_f  = cb_bt;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_exp_f(vrna_fold_compound_t    *fc,
                      vrna_grammar_rule_f_exp cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_exp_f = cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_c(vrna_fold_compound_t  *fc,
                  vrna_grammar_rule_f   cb,
                  vrna_grammar_bt_f     cb_bt)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_c = cb;
    fc->aux_grammar->cb_bt_c  = cb_bt;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_exp_c(vrna_fold_compound_t    *fc,
                      vrna_grammar_rule_f_exp cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_exp_c = cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_m(vrna_fold_compound_t  *fc,
                  vrna_grammar_rule_f   cb,
                  vrna_grammar_bt_f     cb_bt)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_m = cb;
    fc->aux_grammar->cb_bt_m  = cb_bt;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_exp_m(vrna_fold_compound_t    *fc,
                      vrna_grammar_rule_f_exp cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_exp_m = cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_m1(vrna_fold_compound_t *fc,
                   vrna_grammar_rule_f  cb,
                  vrna_grammar_bt_f     cb_bt)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_m1 = cb;
    fc->aux_grammar->cb_bt_m1  = cb_bt;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_exp_m1(vrna_fold_compound_t     *fc,
                       vrna_grammar_rule_f_exp  cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_exp_m1 = cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux(vrna_fold_compound_t    *fc,
                vrna_grammar_rule_f_aux cb,
                  vrna_grammar_bt_f     cb_bt)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux     = cb;
    fc->aux_grammar->cb_bt_aux  = cb_bt;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_exp(vrna_fold_compound_t        *fc,
                    vrna_grammar_rule_f_aux_exp cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_exp = cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_data(vrna_fold_compound_t   *fc,
                 void                   *data,
                 vrna_auxdata_prepare_f prepare_cb,
                 vrna_auxdata_free_f    free_cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->data         = data;
    fc->aux_grammar->free_data    = free_cb;
    fc->aux_grammar->prepare_data = prepare_cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_cond(vrna_fold_compound_t     *fc,
                 vrna_recursion_status_f  cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_proc = cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_reset(vrna_fold_compound_t *fc)
{
  int ret = 0;

  if ((fc) && (fc->aux_grammar)) {
    if (fc->aux_grammar->free_data)
      fc->aux_grammar->free_data(fc->aux_grammar->data);

    vrna_array_free(fc->aux_grammar->bt_flags);
    free(fc->aux_grammar);
    fc->aux_grammar = NULL;
  }

  return ret;
}


PRIVATE void
add_aux_grammar(vrna_fold_compound_t *fc)
{
  fc->aux_grammar = (struct vrna_gr_aux_s *)vrna_alloc(sizeof(struct vrna_gr_aux_s));

  fc->aux_grammar->cb_proc = NULL;

  fc->aux_grammar->cb_aux     = NULL;
  fc->aux_grammar->cb_aux_f   = NULL;
  fc->aux_grammar->cb_aux_c   = NULL;
  fc->aux_grammar->cb_aux_m   = NULL;
  fc->aux_grammar->cb_aux_m1  = NULL;

  fc->aux_grammar->cb_bt_f    = NULL;
  fc->aux_grammar->cb_bt_c    = NULL;
  fc->aux_grammar->cb_bt_m    = NULL;
  fc->aux_grammar->cb_bt_m1   = NULL;
  fc->aux_grammar->cb_bt_aux  = NULL;

  fc->aux_grammar->cb_aux_exp     = NULL;
  fc->aux_grammar->cb_aux_exp_f   = NULL;
  fc->aux_grammar->cb_aux_exp_c   = NULL;
  fc->aux_grammar->cb_aux_exp_m   = NULL;
  fc->aux_grammar->cb_aux_exp_m1  = NULL;

  fc->aux_grammar->data         = NULL;
  fc->aux_grammar->free_data    = NULL;
  fc->aux_grammar->prepare_data = NULL;

  vrna_array_init(fc->aux_grammar->bt_flags);
}
