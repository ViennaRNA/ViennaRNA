#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/constraints/soft.h"


#ifndef INLINE
# ifdef __GNUC__
#   define INLINE inline
# else
#   define INLINE
# endif
#endif

/*
 #################################
 # Type definitions              #
 #################################
 */

typedef struct {
  vrna_array(vrna_sc_direct_f)        cbs;
  vrna_array(vrna_sc_exp_direct_f)    cbs_exp;
  vrna_array(void *)                  data;
  vrna_array(void *)                  data_exp;
  vrna_array(vrna_auxdata_prepare_f)  prepare_data;
  vrna_array(vrna_auxdata_free_f)     free_data;
} sc_cb_container_t;


typedef struct {
  vrna_sc_direct_f  cb;
  void              *data;
} sc_cb_en_wrap;


typedef struct {
  vrna_fold_compound_t  *fc;
  sc_cb_container_t     data[VRNA_DECOMP_TYPES_MAX];
} sc_multi_s;


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE INLINE int
sc_collect(int            i,
           int            j,
           int            k,
           int            l,
           unsigned char  d,
           void           *data);


PRIVATE INLINE FLT_OR_DBL
sc_exp_collect(int            i,
               int            j,
               int            k,
               int            l,
               unsigned char  d,
               void           *data);


PRIVATE void
sc_multi_free(void *data);


PRIVATE int
sc_multi_prepare(vrna_fold_compound_t *fc,
                 void                 *data,
                 unsigned int         event);


PRIVATE INLINE FLT_OR_DBL
cb_exp_default(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               int                  k,
               int                  l,
               void                 *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC size_t
vrna_sc_multi_cb_add(vrna_fold_compound_t   *fc,
                     vrna_sc_direct_f       cb,
                     vrna_sc_exp_direct_f   cb_exp,
                     void                   *data,
                     vrna_auxdata_prepare_f prepare_cb,
                     vrna_auxdata_free_f    free_cb,
                     unsigned int           d)
{
  vrna_sc_t         *sc;
  sc_cb_container_t *data_multi;
  sc_multi_s        *multi_s;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_SINGLE) &&
      (cb) &&
      (d) &&
      (d < VRNA_DECOMP_TYPES_MAX)) {
    if (!fc->sc)
      vrna_sc_init(fc);

    sc      = fc->sc;
    multi_s = NULL;

    if (sc->f != &sc_collect) {
      multi_s = (sc_multi_s *)vrna_alloc(sizeof(sc_multi_s));
      memset(&(multi_s->data[0]), 0, sizeof(sc_cb_container_t) * VRNA_DECOMP_TYPES_MAX);
      multi_s->fc = fc;

      vrna_sc_add_auxdata(fc, multi_s, &sc_multi_prepare, &sc_multi_free);
      vrna_sc_add_f(fc, &sc_collect);
      vrna_sc_add_exp_f(fc, &sc_exp_collect);
    } else {
      multi_s = fc->sc->data;
    }

    if (multi_s) {
      if (!multi_s->data[d].cbs) {
        vrna_array_init(multi_s->data[d].cbs);
        vrna_array_init(multi_s->data[d].cbs_exp);
        vrna_array_init(multi_s->data[d].data);
        vrna_array_init(multi_s->data[d].data_exp);
        vrna_array_init(multi_s->data[d].prepare_data);
        vrna_array_init(multi_s->data[d].free_data);
      }

      data_multi = &(multi_s->data[d]);
      vrna_array_append(data_multi->cbs, cb);
      vrna_array_append(data_multi->data, data);
      vrna_array_append(data_multi->prepare_data, prepare_cb);
      vrna_array_append(data_multi->free_data, free_cb);

      if (cb_exp) {
        vrna_array_append(data_multi->cbs_exp, cb_exp);
        vrna_array_append(data_multi->data_exp, data);
      } else {
        sc_cb_en_wrap *wrapper = (sc_cb_en_wrap *)vrna_alloc(sizeof(sc_cb_en_wrap));
        wrapper->cb   = cb;
        wrapper->data = data;
        vrna_array_append(data_multi->cbs_exp, &cb_exp_default);
        vrna_array_append(data_multi->data_exp, wrapper);
      }

      return vrna_array_size(data_multi->cbs);
    }
  }

  return 0;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE int
sc_collect(int            i,
           int            j,
           int            k,
           int            l,
           unsigned char  d,
           void           *data)
{
  int         e     = 0;
  sc_multi_s  *msc  = (sc_multi_s *)data;

  if (msc->data[d].cbs) {
    vrna_fold_compound_t  *fc     = msc->fc;
    vrna_sc_direct_f      *cbs    = msc->data[d].cbs;
    void                  **data  = msc->data[d].data;
    size_t                stop    = vrna_array_size(cbs);

    for (size_t c = 0; c < stop; c++)
      e += cbs[c](fc, i, j, k, l, data[c]);
  }

  return e;
}


PRIVATE INLINE FLT_OR_DBL
sc_exp_collect(int            i,
               int            j,
               int            k,
               int            l,
               unsigned char  d,
               void           *data)
{
  FLT_OR_DBL  q     = 1.0;
  sc_multi_s  *msc  = (sc_multi_s *)data;

  if (msc->data[d].cbs_exp) {
    vrna_fold_compound_t  *fc       = msc->fc;
    vrna_sc_exp_direct_f  *cbs_exp  = msc->data[d].cbs_exp;
    void                  **data    = msc->data[d].data_exp;
    size_t                stop      = vrna_array_size(cbs_exp);

    for (size_t c = 0; c < stop; c++)
      q *= cbs_exp[c](fc, i, j, k, l, data[c]);
  }

  return q;
}


PRIVATE INLINE FLT_OR_DBL
cb_exp_default(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               int                  k,
               int                  l,
               void                 *data)
{
  sc_cb_en_wrap *wrapper  = (sc_cb_en_wrap *)data;
  FLT_OR_DBL    kT        = fc->exp_params->kT;
  int           en        = wrapper->cb(fc, i, j, k, l, wrapper->data) * 10.;

  return (FLT_OR_DBL)exp(-(FLT_OR_DBL)en / kT);
}


PRIVATE void
sc_multi_free(void *data)
{
  if (data) {
    sc_multi_s *msc = (sc_multi_s *)data;

    /* go through all distinguished loop types */
    for (unsigned char d = 1; d < VRNA_DECOMP_TYPES_MAX; d++) {
      if (msc->data[d].cbs) {
        /* go through all auxiliary data for current loop type callbacks and release memory if required */
        for (size_t c = 0; c < vrna_array_size(msc->data[d].data); c++)
          if (msc->data[d].free_data[c])
            msc->data[d].free_data[c](msc->data[d].data[c]);

        /* release wrapper memory for default partition function callbacks */
        for (size_t c = 0; c < vrna_array_size(msc->data[d].cbs_exp); c++)
          if (msc->data[d].cbs_exp[c] == &cb_exp_default)
            free(msc->data[d].data_exp[c]);

        /* release memory of callback-, data-, prepare_data-, free_data-arrays for current loop type */
        vrna_array_free(msc->data[d].cbs);
        vrna_array_free(msc->data[d].cbs_exp);
        vrna_array_free(msc->data[d].data);
        vrna_array_free(msc->data[d].data_exp);
        vrna_array_free(msc->data[d].prepare_data);
        vrna_array_free(msc->data[d].free_data);
      }
    }

    free(msc);
  }
}


PRIVATE int
sc_multi_prepare(vrna_fold_compound_t *fc,
                 void                 *data,
                 unsigned int         event)
{
  int ret = 0;

  if (data) {
    sc_multi_s *msc = (sc_multi_s *)data;

    /* go through all distinguished loop types */
    for (unsigned char d = 1; d < VRNA_DECOMP_TYPES_MAX; d++) {
      if (msc->data[d].cbs) {
        for (size_t c = 0; c < vrna_array_size(msc->data[d].data); c++)
          if (msc->data[d].prepare_data[c])
            ret |= msc->data[d].prepare_data[c](fc, msc->data[d].data[c], event);
      }
    }
  }

  return ret;
}
