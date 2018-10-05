#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/grammar.h"


PRIVATE void
add_aux_grammar(vrna_fold_compound_t *fc);


PUBLIC int
vrna_gr_set_aux_f(vrna_fold_compound_t  *fc,
                  vrna_callback_gr_rule *cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_f = cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_exp_f(vrna_fold_compound_t      *fc,
                      vrna_callback_gr_rule_exp *cb)
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
                  vrna_callback_gr_rule *cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_c = cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_exp_c(vrna_fold_compound_t      *fc,
                      vrna_callback_gr_rule_exp *cb)
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
                  vrna_callback_gr_rule *cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_m = cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_exp_m(vrna_fold_compound_t      *fc,
                      vrna_callback_gr_rule_exp *cb)
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
vrna_gr_set_aux_m1(vrna_fold_compound_t   *fc,
                   vrna_callback_gr_rule  *cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_m1 = cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_exp_m1(vrna_fold_compound_t       *fc,
                       vrna_callback_gr_rule_exp  *cb)
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
vrna_gr_set_aux(vrna_fold_compound_t  *fc,
                vrna_callback_gr_rule_aux *cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux = cb;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_aux_exp(vrna_fold_compound_t      *fc,
                    vrna_callback_gr_rule_aux_exp *cb)
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
vrna_gr_set_data(vrna_fold_compound_t       *fc,
                 void                       *data,
                 vrna_callback_gr_free_data *free_data)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->data       = data;
    fc->aux_grammar->free_data  = free_data;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_gr_set_cond(vrna_fold_compound_t   *fc,
                 vrna_callback_gr_cond  *cb)
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

    free(fc->aux_grammar);
    fc->aux_grammar = NULL;
  }

  return ret;
}


PRIVATE void
add_aux_grammar(vrna_fold_compound_t *fc)
{
  fc->aux_grammar = (struct vrna_gr_aux_s *)vrna_alloc(sizeof(struct vrna_gr_aux_s *));

  fc->aux_grammar->cb_proc = NULL;

  fc->aux_grammar->cb_aux     = NULL;
  fc->aux_grammar->cb_aux_f   = NULL;
  fc->aux_grammar->cb_aux_c   = NULL;
  fc->aux_grammar->cb_aux_m   = NULL;
  fc->aux_grammar->cb_aux_m1  = NULL;

  fc->aux_grammar->cb_aux_exp     = NULL;
  fc->aux_grammar->cb_aux_exp_f   = NULL;
  fc->aux_grammar->cb_aux_exp_c   = NULL;
  fc->aux_grammar->cb_aux_exp_m   = NULL;
  fc->aux_grammar->cb_aux_exp_m1  = NULL;

  fc->aux_grammar->data       = NULL;
  fc->aux_grammar->free_data  = NULL;
}
