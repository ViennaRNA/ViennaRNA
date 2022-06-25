#include <stdlib.h>
#include <string.h>

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

typedef vrna_sc_cb_multi            * sc_cb_multi_ptr_t;
typedef void                        * sc_cb_data_array_t;
typedef vrna_callback_free_auxdata  * sc_cb_data_free_array_t;

typedef struct {
  sc_cb_multi_ptr_t           *cbs;
  sc_cb_data_array_t      *datas;
  sc_cb_data_free_array_t *free_data;
} sc_cb_container_t;


typedef int (sc_cb_multi_intern)(vrna_fold_compound_t  *fc,
                                     int                   i,
                               int                   j,
                               int                   k,
                               int                   l,
                               sc_cb_container_t       *data);


typedef struct {
  vrna_fold_compound_t  *fc;
  sc_cb_multi_intern    *cb[VRNA_DECOMP_TYPES_MAX];
  sc_cb_container_t     data[VRNA_DECOMP_TYPES_MAX];
} sc_multi_s;





PRIVATE sc_cb_multi_ptr_t *
sc_cb_array_append(sc_cb_multi_ptr_t *container,
                   vrna_sc_cb_multi *callback) {

  if (VRNA_ARRAY_RESIZE_REQUIRED(container, 1))
    container = vrna_array_grow(container, 1);

  container[VRNA_ARRAY_SIZE(container)++] = callback;

  return container;
}


PRIVATE sc_cb_data_array_t *
sc_cb_data_array_append(sc_cb_data_array_t *container,
                        void               *data) {
  if (VRNA_ARRAY_RESIZE_REQUIRED(container, 1))
    container = vrna_array_grow(container, 1);

  container[VRNA_ARRAY_SIZE(container)++] = data;

  return container;
}

PRIVATE sc_cb_data_free_array_t *
sc_cb_free_data_array_append(sc_cb_data_free_array_t    *container,
                             vrna_callback_free_auxdata *free_data) {

  if (VRNA_ARRAY_RESIZE_REQUIRED(container, 1))
    container = vrna_array_grow(container, 1);

  container[VRNA_ARRAY_SIZE(container)++] = free_data;

  return container;
}


PRIVATE int
sc_multi_cb_main(int i,
                 int j,
                 int k,
                 int l,
                 unsigned char d,
                 void *data)
{
  sc_multi_s  *msc = (sc_multi_s *)data;

  if (msc->cb[d])
    return msc->cb[d](msc->fc, i, j, k, l, &(msc->data[d]));

  return 0;
}


PRIVATE INLINE int
sc_multi_collect(vrna_fold_compound_t *fc,
                 int                  i,
                 int                  j,
                 int                  k,
                 int                  l,
                 sc_cb_container_t      *container)
{
  size_t  stop = VRNA_ARRAY_SIZE(container->cbs);
  int e = 0;

  for (size_t c = 0; c < stop; c++)
    e += container->cbs[c](fc, i, j, k, l, container->datas[c]);

  return e;
}


PRIVATE void
sc_multi_free(void *data)
{
  sc_multi_s *msc = (sc_multi_s *)data;

  /* go through all distinguished loop types */
  for (unsigned char d = 1; d < VRNA_DECOMP_TYPES_MAX; d++) {
    if (msc->cb[d]) {
      /* go through all auxiliary data for current loop type callbacks and release memory if required */
      for (size_t c = 0; c < VRNA_ARRAY_SIZE(msc->data[d].datas); c++)
        if (msc->data[d].free_data[c])
          msc->data[d].free_data[c](msc->data[d].datas[c]);

      /* release memory of callback-, data-, free_data-arrays for currecnt loop type */
      vrna_array_free((void *)msc->data[d].cbs);
      vrna_array_free((void *)msc->data[d].datas);
      vrna_array_free((void *)msc->data[d].free_data);
    }
  }
}


PUBLIC size_t
vrna_sc_multi_cb_add(vrna_fold_compound_t       *fc,
                     vrna_sc_cb_multi           *cb,
                     void                       *data,
                     vrna_callback_free_auxdata *free_data,
                     unsigned int               decomp_type)
{
  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_SINGLE) &&
      (cb) &&
      (decomp_type) &&
      (decomp_type < VRNA_DECOMP_TYPES_MAX)) {

    if (!fc->sc)
      vrna_sc_init(fc);

    vrna_sc_t *sc = fc->sc;

    sc_multi_s *multi_s = NULL;

    if (sc->f != &sc_multi_cb_main) {
      multi_s = (sc_multi_s *)vrna_alloc(sizeof(sc_multi_s));
      memset(multi_s->cb, 0, sizeof(sc_cb_multi_intern) * VRNA_DECOMP_TYPES_MAX);
      memset(&(multi_s->data[0]), 0, sizeof(sc_cb_container_t) * VRNA_DECOMP_TYPES_MAX);
      multi_s->fc = fc;

      vrna_sc_add_data(fc, multi_s, &sc_multi_free);
      vrna_sc_add_f(fc, sc_multi_cb_main);
    } else {
      multi_s = fc->sc->data;
    }

    if (multi_s) {
      sc_cb_container_t *data_multi;

      if (!multi_s->cb[decomp_type]) {
        multi_s->cb[decomp_type]              = &sc_multi_collect;
        multi_s->data[decomp_type].cbs        = vrna_array_init(sizeof(sc_cb_multi_ptr_t), 1);
        multi_s->data[decomp_type].datas      = vrna_array_init(sizeof(void *), 1);
        multi_s->data[decomp_type].free_data  = vrna_array_init(sizeof(vrna_callback_free_auxdata *), 1);
      }

      data_multi            = &(multi_s->data[decomp_type]);
      data_multi->cbs       = sc_cb_array_append(data_multi->cbs, cb);
      data_multi->datas     = sc_cb_data_array_append(data_multi->datas, data);
      data_multi->free_data = sc_cb_free_data_array_append(data_multi->free_data, free_data);

      return VRNA_ARRAY_SIZE(data_multi->cbs);
    }
  }

  return 0;
}

