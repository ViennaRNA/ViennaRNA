#ifndef VIENNA_RNA_PACKAGE_EVAL_GQUAD_H
#define VIENNA_RNA_PACKAGE_EVAL_GQUAD_H

#include <string.h>

#include <ViennaRNA/datastructures/basic.h>
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/datastructures/sparse_mx.h"
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
 *  @file       ViennaRNA/eval/gquad.h
 *  @ingroup    gquads
 *  @brief      G-quadruplex Energy Evaluation
 */

/**
 *  @addtogroup gquad_eval
 *  @{
 */

int
vrna_E_gquad(unsigned int L,
             unsigned int l[3],
             vrna_param_t * P);


FLT_OR_DBL
vrna_exp_E_gquad(unsigned int L,
                 unsigned int l[3],
                 vrna_exp_param_t * pf);


void
vrna_E_consensus_gquad(unsigned int L,
                       unsigned int l[3],
                       unsigned int position,
                       unsigned int length,
                       unsigned int n_seq,
                       const short **S,
                       const unsigned int **a2s,
                       vrna_param_t * P,
                       int en[2]);


FLT_OR_DBL
vrna_exp_E_consensus_gquad(unsigned int L,
                           unsigned int l[3],
                           vrna_exp_param_t * pf,
                           unsigned int position,
                           unsigned int length,
                           unsigned int n_seq,
                           const short **S,
                           const unsigned int **a2s);

/**
 *  @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup gquad_deprecated
 *  @{
 */
DEPRECATED(int
           E_gquad(int          L,
                   int          l[3],
                   vrna_param_t *P),
           "Use vrna_E_gquad() instead");


DEPRECATED(FLT_OR_DBL
           exp_E_gquad(int              L,
                       int              l[3],
                       vrna_exp_param_t *pf),
           "Use vrna_exp_E_gquad() instead");


DEPRECATED(void
           E_gquad_ali_en(int           i,
                          int           L,
                          int           l[3],
                          const short   **S,
                          unsigned int  **a2s,
                          unsigned int  n_seq,
                          vrna_param_t  *P,
                          int           en[2]),
           "Use vrna_E_consensus_gquad() instead");


DEPRECATED(FLT_OR_DBL
           exp_E_gquad_ali(int              i,
                           int              L,
                           int              l[3],
                           short            **S,
                           unsigned int     **a2s,
                           int              n_seq,
                           vrna_exp_param_t *pf),
           "Use vrna_exp_E_consensus_gquad() instead");


/**
 * @}
 */


#endif

#endif
