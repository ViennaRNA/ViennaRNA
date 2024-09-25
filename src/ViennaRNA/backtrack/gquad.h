#ifndef VIENNA_RNA_PACKAGE_BT_GQUAD_H
#define VIENNA_RNA_PACKAGE_BT_GQUAD_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>

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
 *  @addtogroup gquad_backtrack
 *  @{
 */
int
vrna_bt_gquad(vrna_fold_compound_t *fc,
              unsigned int         i,
              unsigned int         j,
              unsigned int         *L,
              unsigned int         l[3]);

int
vrna_bt_gquad_mfe(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j,
                  vrna_bps_t            bp_stack);


int
vrna_bt_gquad_internal(vrna_fold_compound_t  *fc,
                       unsigned int          i,
                       unsigned int          j,
                       int                   en,
                       vrna_bps_t            bp_stack,
                       vrna_bts_t            bt_stack);


/**
 * @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup gquad_deprecated
 *  @{
 */
DEPRECATED(int
           vrna_BT_gquad_mfe(vrna_fold_compound_t *fc,
                             int                  i,
                             int                  j,
                             vrna_bp_stack_t      *bp_stack,
                             unsigned int         *stack_count),
           "Use vrna_bt_gquad_mfe() instead");


DEPRECATED(int
           vrna_BT_gquad_int(vrna_fold_compound_t *fc,
                             int                  i,
                             int                  j,
                             int                  en,
                             vrna_bp_stack_t      *bp_stack,
                             unsigned int         *stack_count),
           "Use vrna_bt_gquad_internal() instead");


/**
 *  backtrack an internal loop like enclosed g-quadruplex
 *  with closing pair (i,j) with underlying Lfold matrix
 *
 *  @param c      The total contribution the loop should resemble
 *  @param i      position i of enclosing pair
 *  @param j      position j of enclosing pair
 *  @param type   base pair type of enclosing pair (must be reverse type)
 *  @param S      integer encoded sequence
 *  @param ggg    triangular matrix containing g-quadruplex contributions
 *  @param p      here the 5' position of the gquad is stored
 *  @param q      here the 3' position of the gquad is stored
 *  @param P      the datastructure containing the precalculated contibutions
 *
 *  @return       1 on success, 0 if no gquad found
 */
DEPRECATED(int
           backtrack_GQuad_IntLoop_L(int          c,
                                     int          i,
                                     int          j,
                                     int          type,
                                     short        *S,
                                     int          **ggg,
                                     int          maxdist,
                                     int          *p,
                                     int          *q,
                                     vrna_param_t *P),
           "Use vrna_bt_gquad_internal() instead");


DEPRECATED(int
           backtrack_GQuad_IntLoop_L_comparative(int          c,
                                                 int          i,
                                                 int          j,
                                                 unsigned int *type,
                                                 short        *S_cons,
                                                 short        **S5,
                                                 short        **S3,
                                                 unsigned int **a2s,
                                                 int          **ggg,
                                                 int          *p,
                                                 int          *q,
                                                 int          n_seq,
                                                 vrna_param_t *P),
           "Use vrna_bt_gquad_internal() instead");


/** @} */

#endif

#endif
