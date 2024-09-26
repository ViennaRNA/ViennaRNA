#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/grammar/partfunc.h"

#include "ViennaRNA/intern/grammar_dat.h"

PRIVATE void
add_rule_exp(vrna_array(grammar_rule_exp_wrap_t)  where,
             vrna_gr_inside_exp_f                 cb,
             vrna_gr_outside_exp_f                cb_out,
             void                                 *data,
             vrna_auxdata_prepare_f               data_prepare,
             vrna_auxdata_free_f                  data_release)
{
  grammar_rule_exp_wrap_t r = {
    .cb       = cb,
    .cb_out   = cb_out,
    .data     = data,
    .prepare  = data_prepare,
    .release  = data_release
  };

  vrna_array_append(where, r);
}


PUBLIC unsigned int
vrna_gr_add_aux_exp_f(vrna_fold_compound_t    *fc,
                      vrna_gr_inside_exp_f    cb,
                      vrna_gr_outside_exp_f   cb_out,
                      void                    *data,
                      vrna_auxdata_prepare_f  data_prepare,
                      vrna_auxdata_free_f     data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_out)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule_exp(fc->aux_grammar->exp_f,
                 cb,
                 cb_out,
                 data,
                 data_prepare,
                 data_release);

    ret = vrna_array_size(fc->aux_grammar->exp_f);
  }

  return ret;
}


PUBLIC unsigned int
vrna_gr_add_aux_exp_c(vrna_fold_compound_t    *fc,
                      vrna_gr_inside_exp_f    cb,
                      vrna_gr_outside_exp_f   cb_out,
                      void                    *data,
                      vrna_auxdata_prepare_f  data_prepare,
                      vrna_auxdata_free_f     data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_out)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule_exp(fc->aux_grammar->exp_c,
                 cb,
                 cb_out,
                 data,
                 data_prepare,
                 data_release);

    ret = vrna_array_size(fc->aux_grammar->exp_c);
  }

  return ret;
}


PUBLIC unsigned int
vrna_gr_add_aux_exp_m(vrna_fold_compound_t    *fc,
                      vrna_gr_inside_exp_f    cb,
                      vrna_gr_outside_exp_f   cb_out,
                      void                    *data,
                      vrna_auxdata_prepare_f  data_prepare,
                      vrna_auxdata_free_f     data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_out)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule_exp(fc->aux_grammar->exp_m,
                 cb,
                 cb_out,
                 data,
                 data_prepare,
                 data_release);

    ret = vrna_array_size(fc->aux_grammar->exp_m);
  }

  return ret;
}


PUBLIC unsigned int
vrna_gr_add_aux_exp_m1(vrna_fold_compound_t   *fc,
                       vrna_gr_inside_exp_f   cb,
                       vrna_gr_outside_exp_f  cb_out,
                       void                   *data,
                       vrna_auxdata_prepare_f data_prepare,
                       vrna_auxdata_free_f    data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_out)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule_exp(fc->aux_grammar->exp_m1,
                 cb,
                 cb_out,
                 data,
                 data_prepare,
                 data_release);

    ret = vrna_array_size(fc->aux_grammar->exp_m1);
  }

  return ret;
}


PUBLIC unsigned int
vrna_gr_add_aux_exp_m2(vrna_fold_compound_t   *fc,
                       vrna_gr_inside_exp_f   cb,
                       vrna_gr_outside_exp_f  cb_out,
                       void                   *data,
                       vrna_auxdata_prepare_f data_prepare,
                       vrna_auxdata_free_f    data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_out)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule_exp(fc->aux_grammar->exp_m2,
                 cb,
                 cb_out,
                 data,
                 data_prepare,
                 data_release);

    ret = vrna_array_size(fc->aux_grammar->exp_m2);
  }

  return ret;
}


PUBLIC unsigned int
vrna_gr_add_aux_exp_aux(vrna_fold_compound_t    *fc,
                        vrna_gr_inside_exp_f    cb,
                        vrna_gr_outside_exp_f   cb_out,
                        void                    *data,
                        vrna_auxdata_prepare_f  data_prepare,
                        vrna_auxdata_free_f     data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_out)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule_exp(fc->aux_grammar->exp_aux,
                 cb,
                 cb_out,
                 data,
                 data_prepare,
                 data_release);

    ret = 1 +
          VRNA_MX_FLAG_MAX +
          vrna_array_size(fc->aux_grammar->exp_aux);
  }

  return ret;
}
