#ifndef VIENNA_RNA_PACKAGE_GQUAD_H
#define VIENNA_RNA_PACKAGE_GQUAD_H

#include <string.h>

#include <ViennaRNA/datastructures/basic.h>
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/datastructures/sparse_mx.h"
#include <ViennaRNA/utils/log.h>
#include <ViennaRNA/alphabet.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>

/* backward compatibility */
#include <ViennaRNA/mfe/gquad.h>

/* backward compatibility */
#include <ViennaRNA/eval/gquad.h>

/* backward compatibility */
#include <ViennaRNA/partfunc/gquad.h>


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
 *  @file       ViennaRNA/loops/gquad.h
 *  @ingroup    gquads
 *  @brief      G-quadruplexes
 */

#define   VRNA_GQUAD_DB_SYMBOL      '+'
#define   VRNA_GQUAD_DB_SYMBOL_END  '~'


/**
 *  @addtogroup gquad_parse
 *  @{
 */
void
get_gquad_pattern_mfe(short         *S,
                      int           i,
                      int           j,
                      vrna_param_t  *P,
                      int           *L,
                      int           l[3]);


void
get_gquad_pattern_exhaustive(short        *S,
                             int          i,
                             int          j,
                             vrna_param_t *P,
                             int          *L,
                             int          *l,
                             int          threshold);


void
get_gquad_pattern_pf(short            *S,
                     int              i,
                     int              j,
                     vrna_exp_param_t *pf,
                     int              *L,
                     int              l[3]);


void
vrna_get_gquad_pattern_pf(vrna_fold_compound_t  *fc,
                          unsigned int          i, 
                          unsigned int          j,
                          unsigned int          *L,
                          unsigned int          [3]);


plist *
get_plist_gquad_from_db(const char  *structure,
                        float       pr);


int
get_gquad_count(short *S,
                int   i,
                int   j);


int
get_gquad_layer_count(short *S,
                      int   i,
                      int   j);


void
get_gquad_pattern_mfe_ali(short         **S,
                          unsigned int  **a2s,
                          short         *S_cons,
                          int           n_seq,
                          int           i,
                          int           j,
                          vrna_param_t  *P,
                          int           *L,
                          int           l[3]);


/**
 *  @brief  Parse a G-Quadruplex from a dot-bracket structure string
 *
 *  Given a dot-bracket structure (possibly) containing gquads encoded
 *  by '+' signs (and an optional '~' end sign, find first gquad, return
 *  end position (1-based) or 0 if none found.
 *  Upon return L and l[] contain the number of stacked layers, as well as
 *  the lengths of the linker regions.
 *
 *  @note   For circular RNAs and G-Quadruplexes spanning the n,1-junction
 *          the sum of linkers and g-runs is lower than the end position.
 *          This condition can be used to check whether or not to accept
 *          a G-Quadruplex parsed from the dot-bracket string. Also note,
 *          that such n,1-junction spanning G-Quadruplexes must end with
 *          a `~` sign, to be unambigous.
 *
 *
 *  To parse a string with many gquads, call vrna_gq_parse() repeatedly e.g.
 *
 *  @code
 *  end1 = vrna_gq_parse(struc, &L, l); ... ;
 *  end2 = vrna_gq_parse(struc+end1, &L, l); end2+=end1; ... ;
 *  end3 = vrna_gq_parse(struc+end2, &L, l); end3+=end2; ... ;
 *  @endcode
 *
 *  @param  db_string   The input structure in dot-bracket notation
 *  @param  L           A pointer to an unsigned integer to store the layer (stack) size
 *  @param  l           An array of three values to store the respective linker lenghts
 *  @return             The end position of the G-Quadruplex (1-based) or 0 if not found
 */
unsigned int
vrna_gq_parse(const char *db_string,
              unsigned int *L,
              unsigned int l[3]);

void
vrna_db_insert_gq(char          *db,
                  unsigned int  i,
                  unsigned int  L,
                  unsigned int l[3],
                  unsigned int  n);

/**
 *  @}
 */


/**
 *  @addtogroup gquad_other
 *  @{
 */
plist *
  get_plist_gquad_from_pr(short *S,
                          int gi,
                          int gj,
                          vrna_smx_csr(FLT_OR_DBL) * q_gq,
                          FLT_OR_DBL * probs,
                          FLT_OR_DBL * scale,
                          vrna_exp_param_t * pf);


vrna_ep_t *
vrna_plist_gquad_from_pr(vrna_fold_compound_t *fc,
                         int                  gi,
                         int                  gj);


plist *
  get_plist_gquad_from_pr_max(short *S,
                              int gi,
                              int gj,
                              vrna_smx_csr(FLT_OR_DBL) * q_gq,
                              FLT_OR_DBL * probs,
                              FLT_OR_DBL * scale,
                              int *L,
                              int l[3],
                              vrna_exp_param_t * pf);


vrna_ep_t *
vrna_plist_gquad_from_pr_max(vrna_fold_compound_t *fc,
                             unsigned int                  gi,
                             unsigned int                  gj,
                             unsigned int                  *Lmax,
                             unsigned int                  lmax[3]);


/**
 *  @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup gquad_deprecated
 *  @{
 */


int
E_gquad(int           L,
        int           l[3],
        vrna_param_t  *P);


FLT_OR_DBL
exp_E_gquad(int               L,
            int               l[3],
            vrna_exp_param_t  *pf);


void
E_gquad_ali_en(int          i,
               int          L,
               int          l[3],
               const short  **S,
               unsigned int **a2s,
               unsigned int n_seq,
               vrna_param_t *P,
               int          en[2]);


FLT_OR_DBL
exp_E_gquad_ali(int               i,
                int               L,
                int               l[3],
                short             **S,
                unsigned int      **a2s,
                int               n_seq,
                vrna_exp_param_t  *pf);


/**
 *  @brief  Parse a G-Quadruplex from a dot-bracket structure string
 *
 *  Given a dot-bracket structure (possibly) containing gquads encoded
 *  by '+' signs, find first gquad, return end position or 0 if none found
 *  Upon return L and l[] contain the number of stacked layers, as well as
 *  the lengths of the linker regions.
 *  To parse a string with many gquads, call parse_gquad repeatedly e.g.
 *  end1 = parse_gquad(struc, &L, l); ... ;
 *  end2 = parse_gquad(struc+end1, &L, l); end2+=end1; ... ;
 *  end3 = parse_gquad(struc+end2, &L, l); end3+=end2; ... ;
 */
DEPRECATED(int
           parse_gquad(const char *struc,
                       int        *L,
                       int        l[3]),
           "Use vrna_gq_parse() instead");


DEPRECATED(FLT_OR_DBL * get_gquad_pf_matrix(short *S,
                                            FLT_OR_DBL * scale,
                                            vrna_exp_param_t * pf),
           "Use vrna_gq_pos_pf() instead");


DEPRECATED(FLT_OR_DBL * get_gquad_pf_matrix_comparative(unsigned int n,
                                                        short *S_cons,
                                                        short **S,
                                                        unsigned int **a2s,
                                                        FLT_OR_DBL * scale,
                                                        unsigned int n_seq,
                                                        vrna_exp_param_t * pf),
           "Use vrna_gq_pos_pf() instead");

/**
 * @}
 */


#endif

#endif
