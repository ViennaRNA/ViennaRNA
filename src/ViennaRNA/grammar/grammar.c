#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/grammar/basic.h"

#include "ViennaRNA/intern/grammar_dat.h"

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
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->m1); i++)
      if (fc->aux_grammar->m1[i].prepare)
        ret &= fc->aux_grammar->m1[i].prepare(fc,
                                              fc->aux_grammar->m1[i].data,
                                              options,
                                              NULL);

    /* prepare m2 data */
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->m2); i++)
      if (fc->aux_grammar->m2[i].prepare)
        ret &= fc->aux_grammar->m2[i].prepare(fc,
                                              fc->aux_grammar->m2[i].data,
                                              options,
                                              NULL);

    /* prepare aux data */
    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->aux); i++)
      if (fc->aux_grammar->aux[i].prepare)
        ret &= fc->aux_grammar->aux[i].prepare(fc,
                                               fc->aux_grammar->aux[i].data,
                                               options,
                                               NULL);

    if ((fc->aux_grammar->serialize_bp) &&
        (fc->aux_grammar->serialize_bp_prepare_data))
      ret &= fc->aux_grammar->serialize_bp_prepare_data(fc,
                                                        fc->aux_grammar->serialize_bp_data,
                                                         options,
                                                         NULL);
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

  if ((fc) &&
      (cb)) {
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


PUBLIC unsigned int
vrna_gr_set_serialize_bp(vrna_fold_compound_t   *fc,
                         vrna_gr_serialize_bp_f cb,
                         void                   *data,
                         vrna_auxdata_prepare_f prepare_cb,
                         vrna_auxdata_free_f    free_cb)
{
  unsigned int ret = 0;

  if ((fc) &&
      (cb)) {
    if (!fc->aux_grammar)
      init_aux_grammar(fc);

    fc->aux_grammar->serialize_bp = cb;
    fc->aux_grammar->serialize_bp_data = data;
    fc->aux_grammar->serialize_bp_prepare_data = prepare_cb;
    fc->aux_grammar->serialize_bp_free_data = free_cb;
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

    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->m2); i++)
      if (fc->aux_grammar->m2[i].release)
        fc->aux_grammar->m2[i].release(fc->aux_grammar->m2[i].data);

    vrna_array_free(fc->aux_grammar->m2);

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

    for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->exp_m2); i++)
      if (fc->aux_grammar->exp_m2[i].release)
        fc->aux_grammar->exp_m2[i].release(fc->aux_grammar->exp_m2[i].data);

    vrna_array_free(fc->aux_grammar->exp_m2);

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

    if (fc->aux_grammar->serialize_bp_free_data)
      fc->aux_grammar->serialize_bp_free_data(fc->aux_grammar->serialize_bp_data);

    /* release memory for aux grammar structure itself */
    free(fc->aux_grammar);
    fc->aux_grammar = NULL;
  }
}
