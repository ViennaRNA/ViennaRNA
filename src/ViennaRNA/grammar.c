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
init_aux_grammar(vrna_fold_compound_t *fc);


PUBLIC unsigned int
vrna_gr_prepare(vrna_fold_compound_t  *fc,
                unsigned int          options)
{
  unsigned int ret = 1;

  if ((fc) &&
      (fc->aux_grammar)) {
    /* prepare status data */
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->datas); i++)
      if (fc->aux_grammar->prepare_datas[i])
        ret &= fc->aux_grammar->prepare_datas[i](fc,
                                                 fc->aux_grammar->datas[i],
                                                 options,
                                                 NULL);

    /* prepare f data */
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->f); i++)
      if (fc->aux_grammar->f[i].prepare)
        ret &= fc->aux_grammar->f[i].prepare(fc,
                                             fc->aux_grammar->f[i].data,
                                             options,
                                             NULL);

    /* prepare c data */
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->c); i++)
      if (fc->aux_grammar->c[i].prepare)
        ret &= fc->aux_grammar->c[i].prepare(fc,
                                             fc->aux_grammar->c[i].data,
                                             options,
                                             NULL);

    /* prepare m data */
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->m); i++)
      if (fc->aux_grammar->m[i].prepare)
        ret &= fc->aux_grammar->m[i].prepare(fc,
                                             fc->aux_grammar->m[i].data,
                                             options,
                                             NULL);

    /* prepare m1 data */
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->f); i++)
      if (fc->aux_grammar->m1[i].prepare)
        ret &= fc->aux_grammar->m1[i].prepare(fc,
                                              fc->aux_grammar->m1[i].data,
                                              options,
                                              NULL);

    /* prepare aux data */
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->aux); i++)
      if (fc->aux_grammar->aux[i].prepare)
        ret &= fc->aux_grammar->aux[i].prepare(fc,
                                               fc->aux_grammar->aux[i].data,
                                               options,
                                               NULL);
  }

  return ret;
}


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


PUBLIC unsigned int
vrna_gr_add_status(vrna_fold_compound_t     *fc,
                   vrna_recursion_status_f  cb,
                   void                     *data,
                   vrna_auxdata_prepare_f   prepare_cb,
                   vrna_auxdata_free_f      free_cb)
{
  unsigned int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    vrna_array_append(fc->aux_grammar->cbs_status, cb);
    vrna_array_append(fc->aux_grammar->datas, data);
    vrna_array_append(fc->aux_grammar->prepare_datas, prepare_cb);
    vrna_array_append(fc->aux_grammar->free_datas, free_cb);

    ret = vrna_array_size(fc->aux_grammar->cbs_status);
  }

  return ret;
}


PUBLIC void
vrna_gr_reset(vrna_fold_compound_t *fc)
{
  if ((fc) &&
      (fc->aux_grammar)) {
    /* release memory for MFE aux grammars */
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->f); i++)
      if (fc->aux_grammar->f[i].release)
        fc->aux_grammar->f[i].release(fc->aux_grammar->f[i].data);

    vrna_array_free(fc->aux_grammar->f);

    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->c); i++)
      if (fc->aux_grammar->c[i].release)
        fc->aux_grammar->c[i].release(fc->aux_grammar->c[i].data);

    vrna_array_free(fc->aux_grammar->c);

    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->m); i++)
      if (fc->aux_grammar->m[i].release)
        fc->aux_grammar->m[i].release(fc->aux_grammar->m[i].data);

    vrna_array_free(fc->aux_grammar->m);

    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->m1); i++)
      if (fc->aux_grammar->m1[i].release)
        fc->aux_grammar->m1[i].release(fc->aux_grammar->m1[i].data);

    vrna_array_free(fc->aux_grammar->m1);

    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->aux); i++)
      if (fc->aux_grammar->aux[i].release)
        fc->aux_grammar->aux[i].release(fc->aux_grammar->aux[i].data);

    vrna_array_free(fc->aux_grammar->aux);

    /* release memory for partition function aux grammars */
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->exp_f); i++)
      if (fc->aux_grammar->exp_f[i].release)
        fc->aux_grammar->exp_f[i].release(fc->aux_grammar->exp_f[i].data);

    vrna_array_free(fc->aux_grammar->exp_f);

    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->exp_c); i++)
      if (fc->aux_grammar->exp_c[i].release)
        fc->aux_grammar->exp_c[i].release(fc->aux_grammar->exp_c[i].data);

    vrna_array_free(fc->aux_grammar->exp_c);

    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->exp_m); i++)
      if (fc->aux_grammar->exp_m[i].release)
        fc->aux_grammar->exp_m[i].release(fc->aux_grammar->exp_m[i].data);

    vrna_array_free(fc->aux_grammar->exp_m);

    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->exp_m1); i++)
      if (fc->aux_grammar->exp_m1[i].release)
        fc->aux_grammar->exp_m1[i].release(fc->aux_grammar->exp_m1[i].data);

    vrna_array_free(fc->aux_grammar->exp_m1);

    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->exp_aux); i++)
      if (fc->aux_grammar->exp_aux[i].release)
        fc->aux_grammar->exp_aux[i].release(fc->aux_grammar->exp_aux[i].data);

    vrna_array_free(fc->aux_grammar->exp_aux);

    /* release memory for aux grammar status messages */
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->datas); i++)
      if (fc->aux_grammar->free_datas[i])
        fc->aux_grammar->free_datas[i](fc->aux_grammar->datas[i]);

    vrna_array_free(fc->aux_grammar->datas);
    vrna_array_free(fc->aux_grammar->prepare_datas);
    vrna_array_free(fc->aux_grammar->free_datas);
    vrna_array_free(fc->aux_grammar->cbs_status);

    /* release memory for aux grammar structure itself */
    free(fc->aux_grammar);
    fc->aux_grammar = NULL;
  }
}


PRIVATE void
init_aux_grammar(vrna_fold_compound_t *fc)
{
  fc->aux_grammar = (struct vrna_gr_aux_s *)vrna_alloc(sizeof(struct vrna_gr_aux_s));

  if (fc->aux_grammar) {
    vrna_array_init(fc->aux_grammar->aux);
    vrna_array_init(fc->aux_grammar->f);
    vrna_array_init(fc->aux_grammar->c);
    vrna_array_init(fc->aux_grammar->m);
    vrna_array_init(fc->aux_grammar->m1);

    vrna_array_init(fc->aux_grammar->exp_f);
    vrna_array_init(fc->aux_grammar->exp_c);
    vrna_array_init(fc->aux_grammar->exp_m);
    vrna_array_init(fc->aux_grammar->exp_m1);
    vrna_array_init(fc->aux_grammar->exp_aux);

    vrna_array_init(fc->aux_grammar->cbs_status);
    vrna_array_init(fc->aux_grammar->datas);
    vrna_array_init(fc->aux_grammar->free_datas);
    vrna_array_init(fc->aux_grammar->prepare_datas);
  }
}
