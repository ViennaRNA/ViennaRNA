#ifndef VIENNA_RNA_PACKAGE_GQUAD_MFE_H
#define VIENNA_RNA_PACKAGE_GQUAD_MFE_H

#include <string.h>

#include <ViennaRNA/datastructures/basic.h>
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/datastructures/sparse_mx.h"
#include <ViennaRNA/utils/log.h>
#include <ViennaRNA/sequences/alphabet.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>

#ifndef INLINE
#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif
#endif

#ifdef VRNA_WARN_DEPRECATED
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

/**
 *  @file       ViennaRNA/mfe/gquad.h
 *  @ingroup    gquads
 *  @brief      G-quadruplexes MFE API
 */


/**
 *  @addtogroup gquad_eval
 *  @{
 */

int
vrna_mfe_gquad_internal_loop(vrna_fold_compound_t *fc,
                             unsigned int         i,
                             unsigned int         j);


/**
 *  @}
 */


/**
 *  @addtogroup gquad_dp
 *  @{
 */


/**
 *  @brief  Get G-Quadruplexes (MFE)
 *
 *  This function yields a sparse, two-dimensional matrix @f$ G @f@ that at
 *  position @f$ G(i,j) @f$ stores the minimum free energy of a G-Quadruplex
 *  that starts at @f$ i @f$ and ends at @f$ j @f$.
 *
 *  @param  fc    The fold_compound
 *  @return       A sparse matrix with all MFE G-Quadruplexes
 */
vrna_smx_csr_int_t *
vrna_mfe_gquad_mx(vrna_fold_compound_t *fc);


int **
get_gquad_L_matrix(short        *S,
                   int          start,
                   int          maxdist,
                   int          n,
                   int          **g,
                   vrna_param_t *P);


void
vrna_gquad_mx_local_update(vrna_fold_compound_t *fc,
                           int                  start);


/**
 *  @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup gquad_deprecated
 *  @{
 */


/**
 *  @brief Get a triangular matrix prefilled with minimum free energy
 *  contributions of G-quadruplexes.
 *
 *  At each position ij in the matrix, the minimum free energy of any
 *  G-quadruplex delimited by i and j is stored. If no G-quadruplex formation
 *  is possible, the matrix element is set to INF.
 *  Access the elements in the matrix via matrix[indx[j]+i]. To get
 *  the integer array indx see get_jindx().
 *
 *  @see get_jindx(), encode_sequence()
 *
 *  @param S  The encoded sequence
 *  @param P  A pointer to the data structure containing the precomputed energy contributions
 *  @return   A pointer to the G-quadruplex contribution matrix
 */
DEPRECATED(int *
           get_gquad_matrix(short *S,
                            vrna_param_t * P),
           "Use vrna_mfe_gquad_mx() instead");


DEPRECATED(int *
           get_gquad_ali_matrix(unsigned int n,
                                short *S_cons,
                                short **S,
                                unsigned int **a2s,
                                int n_seq,
                                vrna_param_t * P),
           "Use vrna_mfe_gquad_mx() instead");


DEPRECATED(int
           E_GQuad_IntLoop_L_comparative(int          i,
                                         int          j,
                                         unsigned int *tt,
                                         short        *S_cons,
                                         short        **S5,
                                         short        **S3,
                                         unsigned int **a2s,
                                         int          **ggg,
                                         int          n_seq,
                                         vrna_param_t *P),
           "Use vrna_mfe_gquad_internal_loop() instead");


DEPRECATED(int
           E_GQuad_IntLoop_L(int          i,
                             int          j,
                             int          type,
                             short        *S,
                             int          **ggg,
                             int          maxdist,
                             vrna_param_t *P),
           "Use vrna_mfe_gquad_internal_loop() instead");

/**
 * @}
 */


#endif

#endif
