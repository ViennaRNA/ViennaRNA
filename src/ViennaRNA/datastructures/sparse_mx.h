#ifndef   VIENNA_RNA_PACKAGE_SPARSE_MATRIX_H
#define   VIENNA_RNA_PACKAGE_SPARSE_MATRIX_H

#include <stdlib.h>
#include <ViennaRNA/datastructures/array.h>


#define vrna_smx_csr( TYPE ) \
  vrna_smx_csr_##TYPE##_t

#define VRNA_SMX_CSR_DECLARE( TYPE )      \
  typedef struct vrna_smx_csr_##TYPE##_s  \
  {                                       \
    unsigned char    dirty;               \
    vrna_array(TYPE) v;                   \
    vrna_array(unsigned int) col_idx;     \
    vrna_array(unsigned int) row_idx;     \
  } vrna_smx_csr_##TYPE##_t;


#define VRNA_SMX_CSR_DECLARE_INIT( TYPE )     \
  vrna_smx_csr_##TYPE##_t *                   \
  vrna_smx_csr_##TYPE##_init(unsigned int n);


#define VRNA_SMX_CSR_DECLARE_FREE( TYPE )                   \
  void                                                      \
  vrna_smx_csr_##TYPE##_free(vrna_smx_csr_##TYPE##_t *mx);


/*  note, here we assume data to be stored arrives in an ordered way,
 *  in particular row-wise (low to high) and columns low to high, in
 *  that order!
 */
#define VRNA_SMX_CSR_DECLARE_INSERT( TYPE )      \
  void  \
  vrna_smx_csr_##TYPE##_insert(vrna_smx_csr_##TYPE##_t  *mx,  \
                               unsigned int             i,    \
                               unsigned int             j,    \
                               TYPE                     e);


#define VRNA_SMX_CSR_DECLARE_GET( TYPE )   \
  TYPE  \
  vrna_smx_csr_##TYPE##_get(vrna_smx_csr_##TYPE##_t *mx,        \
                            unsigned int            i,          \
                            unsigned int            j,          \
                            TYPE                    default_v);


#define VRNA_SMX_CSR_DECLARE_GET_SIZE( TYPE )                   \
  size_t                                                        \
  vrna_smx_csr_##TYPE##_get_size(vrna_smx_csr_##TYPE##_t  *mx);



#define VRNA_SMX_CSR_DECLARE_GET_ENTRY( TYPE )                   \
  TYPE  \
  vrna_smx_csr_##TYPE##_get_entry(vrna_smx_csr_##TYPE##_t *mx,        \
                                  size_t                  pos,        \
                                  unsigned int            *i,         \
                                  unsigned int            *j,         \
                                  TYPE                    default_v);


#define VRNA_SMX_CSR_DECLARE_ALL( TYPE )      \
        VRNA_SMX_CSR_DECLARE(TYPE)            \
        VRNA_SMX_CSR_DECLARE_INIT(TYPE)       \
        VRNA_SMX_CSR_DECLARE_FREE(TYPE)       \
        VRNA_SMX_CSR_DECLARE_INSERT(TYPE)     \
        VRNA_SMX_CSR_DECLARE_GET(TYPE)        \
        VRNA_SMX_CSR_DECLARE_GET_SIZE(TYPE)   \
        VRNA_SMX_CSR_DECLARE_GET_ENTRY(TYPE)


#ifndef VRNA_DISABLE_C11_FEATURES

#define GENERIC_SMX_CSR_PTR_FN_DEFINE(SMX_CSR_PTR, FN_NAME) \
  _Generic((SMX_CSR_PTR), \
      vrna_smx_csr_int_t *: (vrna_smx_csr_int ## _ ## FN_NAME), \
      vrna_smx_csr_float_t *: (vrna_smx_csr_float ## _ ## FN_NAME), \
      vrna_smx_csr_double_t *: (vrna_smx_csr_double ## _ ## FN_NAME), \
      vrna_smx_csr_FLT_OR_DBL_t *: (vrna_smx_csr_FLT_OR_DBL ## _ ## FN_NAME))

#define vrna_smx_csr_insert(X, i, j, v) \
  GENERIC_SMX_CSR_PTR_FN_DEFINE((X), insert)((X), i, j, v)

#define vrna_smx_csr_get(X, i, j, d) \
  GENERIC_SMX_CSR_PTR_FN_DEFINE((X), get)((X), i, j, d)

#define vrna_smx_csr_get_size(X) \
  GENERIC_SMX_CSR_PTR_FN_DEFINE((X), get_size)((X))

#define vrna_smx_csr_get_entry(X, p, i, j, d) \
  GENERIC_SMX_CSR_PTR_FN_DEFINE((X), get_entry)((X), p, i, j, d)

#define vrna_smx_csr_free(X) \
  GENERIC_SMX_CSR_PTR_FN_DEFINE((X), free)((X))

#endif

/* Below follows a list of sparse matrix declarations for different data types */
VRNA_SMX_CSR_DECLARE_ALL(int)
VRNA_SMX_CSR_DECLARE_ALL(float)
VRNA_SMX_CSR_DECLARE_ALL(double)
VRNA_SMX_CSR_DECLARE_ALL(FLT_OR_DBL)


#endif
