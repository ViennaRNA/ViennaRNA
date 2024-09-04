#include <stdlib.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/sparse_mx.h"


#define VRNA_SMX_CSR_DEFINE_INIT( TYPE )                                            \
  vrna_smx_csr_##TYPE##_t *                                                         \
  vrna_smx_csr_##TYPE##_init(unsigned int n)                                        \
  {                                                                                 \
    vrna_smx_csr_##TYPE##_t *mx = NULL;                                             \
    if (n > 0) {                                                                    \
      mx = (vrna_smx_csr_##TYPE##_t *)vrna_alloc(sizeof(vrna_smx_csr_##TYPE##_t));  \
      mx->dirty = (unsigned char)1;                                                 \
      vrna_array_init(mx->v);                                                       \
      vrna_array_init(mx->col_idx);                                                 \
      vrna_array_init(mx->row_idx);                                                 \
      vrna_array_set_capacity(mx->row_idx, n + 1);                                  \
    }                                                                               \
    return mx;                                                                      \
  }


#define VRNA_SMX_CSR_DEFINE_FREE( TYPE )                  \
  void                                                    \
  vrna_smx_csr_##TYPE##_free(vrna_smx_csr_##TYPE##_t *mx) \
  {                                                       \
    if (mx) {                                             \
      vrna_array_free(mx->v);                             \
      vrna_array_free(mx->col_idx);                       \
      vrna_array_free(mx->row_idx);                       \
      free(mx);                                           \
    }                                                     \
  }


/*  note, here we assume data to be stored arrives in an ordered way,
 *  in particular row-wise (low to high) and columns low to high, in
 *  that order!
 */
#define VRNA_SMX_CSR_DEFINE_INSERT( TYPE )      \
  void  \
  vrna_smx_csr_##TYPE##_insert(vrna_smx_csr_##TYPE##_t  *mx,  \
                               unsigned int             i,    \
                               unsigned int             j,    \
                               TYPE                     e)    \
  {                                                           \
    vrna_array_append(mx->v, e);                              \
    vrna_array_append(mx->col_idx, j);                        \
    mx->row_idx[i + 1]++;                                     \
  }


#define VRNA_SMX_CSR_DEFINE_GET( TYPE )   \
  TYPE  \
  vrna_smx_csr_##TYPE##_get(vrna_smx_csr_##TYPE##_t *mx,        \
                            unsigned int            i,          \
                            unsigned int            j,          \
                            TYPE                    default_v)  \
  {                                                             \
    unsigned int  p, d, s;                                      \
    if (mx->dirty) {                                            \
      for (s = 1; s < vrna_array_capacity(mx->row_idx); s++)    \
        mx->row_idx[s] += mx->row_idx[s - 1];                   \
      mx->dirty = (unsigned char)0;                             \
    }                                                           \
    s = mx->row_idx[i];                                         \
    d = mx->row_idx[i + 1];                                     \
    if (s < d) {                                                \
      d -= s;                                                   \
      for (p = 0; p < d; p++)                                   \
        if (mx->col_idx[s + p] == j)                            \
          return mx->v[s + p];                                  \
    }                                                           \
    return default_v;                                           \
  }

#define VRNA_SMX_CSR_DEFINE_GET_SIZE( TYPE )                    \
  size_t                                                        \
  vrna_smx_csr_##TYPE##_get_size(vrna_smx_csr_##TYPE##_t  *mx)  \
  {                                                             \
    if(mx)                                                      \
      return vrna_array_size(mx->v);                            \
                                                                \
    return 0;                                                   \
  }

#define VRNA_SMX_CSR_DEFINE_GET_ENTRY( TYPE )   \
  TYPE  \
  vrna_smx_csr_##TYPE##_get_entry(vrna_smx_csr_##TYPE##_t *mx,        \
                                  size_t                  pos,        \
                                  unsigned int            *i,         \
                                  unsigned int            *j,         \
                                  TYPE                    default_v)  \
  {                                                                   \
    unsigned int  d, s;                                               \
    if ((mx) &&                                                       \
        (pos < vrna_array_size(mx->v)) &&                             \
        (i) && (j)) {                                                 \
      if (mx->dirty) {                                                \
        for (s = 1; s < vrna_array_capacity(mx->row_idx); s++)        \
          mx->row_idx[s] += mx->row_idx[s - 1];                       \
        mx->dirty = (unsigned char)0;                                 \
      }                                                               \
      /* get j-position */                                            \
      *j = mx->col_idx[pos];                                          \
      /* go through row_idx to find i position */                     \
      for (d = 1; d < vrna_array_capacity(mx->row_idx); d++) {        \
        /* found i-position */                                        \
        if (mx->row_idx[d] > pos) {                                   \
          *i = d - 1;                                                 \
          return mx->v[pos];                                          \
        }                                                             \
      }                                                               \
    }                                                                 \
    return default_v;                                                 \
  }


#define VRNA_SMX_CSR_DEFINE_ALL( TYPE )  \
        VRNA_SMX_CSR_DEFINE_INIT(TYPE)   \
        VRNA_SMX_CSR_DEFINE_FREE(TYPE)   \
        VRNA_SMX_CSR_DEFINE_INSERT(TYPE) \
        VRNA_SMX_CSR_DEFINE_GET(TYPE) \
        VRNA_SMX_CSR_DEFINE_GET_SIZE(TYPE) \
        VRNA_SMX_CSR_DEFINE_GET_ENTRY(TYPE)

VRNA_SMX_CSR_DEFINE_ALL(int)
VRNA_SMX_CSR_DEFINE_ALL(float)
VRNA_SMX_CSR_DEFINE_ALL(double)
VRNA_SMX_CSR_DEFINE_ALL(FLT_OR_DBL)
