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

typedef vrna_callback_sc_direct *sc_cb_array_t;
typedef vrna_callback_sc_exp_direct *sc_cb_exp_array_t;
typedef void *sc_cb_data_array_t;
typedef vrna_callback_free_auxdata *sc_cb_data_free_array_t;

typedef struct {
  sc_cb_array_t           *cbs;
  sc_cb_exp_array_t       *cbs_exp;
  sc_cb_data_array_t      *data;
  sc_cb_data_array_t      *data_exp;
  sc_cb_data_free_array_t *free_data;
} sc_cb_container_t;


typedef struct {
  vrna_callback_sc_direct *cb;
  void                    *data;
} sc_cb_en_wrap;

typedef int (sc_cb_multi_intern)(vrna_fold_compound_t *fc,
                                 int                  i,
                                 int                  j,
                                 int                  k,
                                 int                  l,
                                 sc_cb_container_t    *data);


typedef struct {
  vrna_fold_compound_t  *fc;
  sc_cb_container_t     data[VRNA_DECOMP_TYPES_MAX];
} sc_multi_s;


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE sc_cb_array_t *
sc_cb_array_append(sc_cb_array_t            *container,
                   vrna_callback_sc_direct  *callback);


PRIVATE sc_cb_exp_array_t *
sc_cb_exp_array_append(sc_cb_exp_array_t            *container,
                       vrna_callback_sc_exp_direct  *callback);


PRIVATE sc_cb_data_array_t *
sc_cb_data_array_append(sc_cb_data_array_t  *container,
                        void                *data);


PRIVATE sc_cb_data_free_array_t *
sc_cb_free_data_array_append(sc_cb_data_free_array_t    *container,
                             vrna_callback_free_auxdata *free_data);


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
vrna_sc_multi_cb_add(vrna_fold_compound_t         *fc,
                     vrna_callback_sc_direct      *cb,
                     vrna_callback_sc_exp_direct  *cb_exp,
                     void                         *data,
                     vrna_callback_free_auxdata   *free_data,
                     unsigned int                 d)
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

      vrna_sc_add_data(fc, multi_s, &sc_multi_free);
      vrna_sc_add_f(fc, &sc_collect);
      vrna_sc_add_exp_f(fc, &sc_exp_collect);
    } else {
      multi_s = fc->sc->data;
    }

    if (multi_s) {
      if (!multi_s->data[d].cbs) {
        multi_s->data[d].cbs        = vrna_array_init(sizeof(sc_cb_array_t), 1);
        multi_s->data[d].cbs_exp    = vrna_array_init(sizeof(sc_cb_exp_array_t), 1);
        multi_s->data[d].data       = vrna_array_init(sizeof(void *), 1);
        multi_s->data[d].data_exp   = vrna_array_init(sizeof(void *), 1);
        multi_s->data[d].free_data  = vrna_array_init(sizeof(vrna_callback_free_auxdata *), 1);
      }

      data_multi            = &(multi_s->data[d]);
      data_multi->cbs       = sc_cb_array_append(data_multi->cbs, cb);
      data_multi->data      = sc_cb_data_array_append(data_multi->data, data);
      data_multi->free_data = sc_cb_free_data_array_append(data_multi->free_data, free_data);

      if (cb_exp) {
        data_multi->cbs_exp   = sc_cb_exp_array_append(data_multi->cbs_exp, cb_exp);
        data_multi->data_exp  = sc_cb_data_array_append(data_multi->data_exp, data);
      } else {
        sc_cb_en_wrap *wrapper = (sc_cb_en_wrap *)vrna_alloc(sizeof(sc_cb_en_wrap));
        wrapper->cb           = cb;
        wrapper->data         = data;
        data_multi->cbs_exp   = sc_cb_exp_array_append(data_multi->cbs_exp, &cb_exp_default);
        data_multi->data_exp  = sc_cb_data_array_append(data_multi->data_exp, wrapper);
      }

      return VRNA_ARRAY_SIZE(data_multi->cbs);
    }
  }

  return 0;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE sc_cb_array_t *
sc_cb_array_append(sc_cb_array_t            *container,
                   vrna_callback_sc_direct  *callback)
{
  if (VRNA_ARRAY_RESIZE_REQUIRED(container, 1))
    container = vrna_array_grow(container, 1);

  container[VRNA_ARRAY_SIZE(container)++] = callback;

  return container;
}


PRIVATE sc_cb_exp_array_t *
sc_cb_exp_array_append(sc_cb_exp_array_t            *container,
                       vrna_callback_sc_exp_direct  *callback)
{
  if (VRNA_ARRAY_RESIZE_REQUIRED(container, 1))
    container = vrna_array_grow(container, 1);

  container[VRNA_ARRAY_SIZE(container)++] = callback;

  return container;
}


PRIVATE sc_cb_data_array_t *
sc_cb_data_array_append(sc_cb_data_array_t  *container,
                        void                *data)
{
  if (VRNA_ARRAY_RESIZE_REQUIRED(container, 1))
    container = vrna_array_grow(container, 1);

  container[VRNA_ARRAY_SIZE(container)++] = data;

  return container;
}


PRIVATE sc_cb_data_free_array_t *
sc_cb_free_data_array_append(sc_cb_data_free_array_t    *container,
                             vrna_callback_free_auxdata *free_data)
{
  if (VRNA_ARRAY_RESIZE_REQUIRED(container, 1))
    container = vrna_array_grow(container, 1);

  container[VRNA_ARRAY_SIZE(container)++] = free_data;

  return container;
}


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
    sc_cb_array_t         *cbs    = msc->data[d].cbs;
    void                  **data  = msc->data[d].data;
    size_t                stop    = VRNA_ARRAY_SIZE(cbs);

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
    sc_cb_exp_array_t     *cbs_exp  = msc->data[d].cbs_exp;
    void                  **data    = msc->data[d].data_exp;
    size_t                stop      = VRNA_ARRAY_SIZE(cbs_exp);

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
  sc_multi_s *msc = (sc_multi_s *)data;

  /* go through all distinguished loop types */
  for (unsigned char d = 1; d < VRNA_DECOMP_TYPES_MAX; d++) {
    if (msc->data[d].cbs) {
      /* go through all auxiliary data for current loop type callbacks and release memory if required */
      for (size_t c = 0; c < VRNA_ARRAY_SIZE(msc->data[d].data); c++)
        if (msc->data[d].free_data[c])
          msc->data[d].free_data[c](msc->data[d].data[c]);

      /* release wrapper memory for default partition function callbacks */
      for (size_t c = 0; c < VRNA_ARRAY_SIZE(msc->data[d].cbs_exp); c++)
        if (msc->data[d].cbs_exp[c] == &cb_exp_default)
          free(msc->data[d].data_exp[c]);

      /* release memory of callback-, data-, free_data-arrays for current loop type */
      vrna_array_free((void *)msc->data[d].cbs);
      vrna_array_free((void *)msc->data[d].cbs_exp);
      vrna_array_free((void *)msc->data[d].data);
      vrna_array_free((void *)msc->data[d].data_exp);
      vrna_array_free((void *)msc->data[d].free_data);
    }
  }
}
