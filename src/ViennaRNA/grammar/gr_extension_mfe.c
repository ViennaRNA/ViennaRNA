#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/grammar/mfe.h"

#include "ViennaRNA/intern/grammar_dat.h"

PRIVATE void
add_rule(vrna_array(grammar_rule_wrap_t)  where,
         vrna_gr_inside_f                 cb,
         vrna_gr_outside_f                cb_bt,
         void                             *data,
         vrna_auxdata_prepare_f           data_prepare,
         vrna_auxdata_free_f              data_release)
{
  grammar_rule_wrap_t r = {
    .cb       = cb,
    .cb_bt    = cb_bt,
    .data     = data,
    .prepare  = data_prepare,
    .release  = data_release
  };

  vrna_array_append(where, r);
}


PUBLIC unsigned int
vrna_gr_add_aux_f(vrna_fold_compound_t    *fc,
                  vrna_gr_inside_f        cb,
                  vrna_gr_outside_f       cb_bt,
                  void                    *data,
                  vrna_auxdata_prepare_f  data_prepare,
                  vrna_auxdata_free_f     data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_bt)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule(fc->aux_grammar->f,
             cb,
             cb_bt,
             data,
             data_prepare,
             data_release);

    ret = vrna_array_size(fc->aux_grammar->f);
  }

  return ret;
}


PUBLIC unsigned int
vrna_gr_add_aux_c(vrna_fold_compound_t    *fc,
                  vrna_gr_inside_f        cb,
                  vrna_gr_outside_f       cb_bt,
                  void                    *data,
                  vrna_auxdata_prepare_f  data_prepare,
                  vrna_auxdata_free_f     data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_bt)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule(fc->aux_grammar->c,
             cb,
             cb_bt,
             data,
             data_prepare,
             data_release);

    ret = vrna_array_size(fc->aux_grammar->c);
  }

  return ret;
}


PUBLIC unsigned int
vrna_gr_add_aux_m(vrna_fold_compound_t    *fc,
                  vrna_gr_inside_f        cb,
                  vrna_gr_outside_f       cb_bt,
                  void                    *data,
                  vrna_auxdata_prepare_f  data_prepare,
                  vrna_auxdata_free_f     data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_bt)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule(fc->aux_grammar->m,
             cb,
             cb_bt,
             data,
             data_prepare,
             data_release);

    ret = vrna_array_size(fc->aux_grammar->m);
  }

  return ret;
}


PUBLIC unsigned int
vrna_gr_add_aux_m1(vrna_fold_compound_t   *fc,
                   vrna_gr_inside_f       cb,
                   vrna_gr_outside_f      cb_bt,
                   void                   *data,
                   vrna_auxdata_prepare_f data_prepare,
                   vrna_auxdata_free_f    data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_bt)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule(fc->aux_grammar->m1,
             cb,
             cb_bt,
             data,
             data_prepare,
             data_release);

    ret = vrna_array_size(fc->aux_grammar->m1);
  }

  return ret;
}


PUBLIC unsigned int
vrna_gr_add_aux_m2(vrna_fold_compound_t   *fc,
                   vrna_gr_inside_f       cb,
                   vrna_gr_outside_f      cb_bt,
                   void                   *data,
                   vrna_auxdata_prepare_f data_prepare,
                   vrna_auxdata_free_f    data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_bt)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule(fc->aux_grammar->m2,
             cb,
             cb_bt,
             data,
             data_prepare,
             data_release);

    ret = vrna_array_size(fc->aux_grammar->m2);
  }

  return ret;
}


PUBLIC unsigned int
vrna_gr_add_aux(vrna_fold_compound_t    *fc,
                vrna_gr_inside_f        cb,
                vrna_gr_outside_f       cb_bt,
                void                    *data,
                vrna_auxdata_prepare_f  data_prepare,
                vrna_auxdata_free_f     data_release)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb || cb_bt)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    add_rule(fc->aux_grammar->aux,
             cb,
             cb_bt,
             data,
             data_prepare,
             data_release);

    ret = 1 +
          VRNA_MX_FLAG_MAX +
          vrna_array_size(fc->aux_grammar->aux);
  }

  return ret;
}
