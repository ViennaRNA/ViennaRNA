#ifndef VIENNA_RNA_PACKAGE_GQUAD_H
#define VIENNA_RNA_PACKAGE_GQUAD_H

#include <ViennaRNA/datastructures/basic.h>
#include "ViennaRNA/datastructures/sparse_mx.h"
#include <ViennaRNA/alphabet.h>
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
 *  @file       gquad.h
 *  @ingroup    paired_modules
 *  @brief      G-quadruplexes
 */

/**
 *  @addtogroup gquads
 *  @{
 *
 *  @brief Various functions related to G-quadruplex computations
 */


int         E_gquad(int           L,
                    int           l[3],
                    vrna_param_t  *P);


FLT_OR_DBL  exp_E_gquad(int               L,
                        int               l[3],
                        vrna_exp_param_t  *pf);


void E_gquad_ali_en(int           i,
                    int           L,
                    int           l[3],
                    const short   **S,
                    unsigned int  **a2s,
                    unsigned int  n_seq,
                    vrna_param_t  *P,
                    int           en[2]);



vrna_smx_csr(int) *
vrna_gq_pos_mfe(vrna_fold_compound_t *fc);

vrna_smx_csr(FLT_OR_DBL) *
vrna_gq_pos_pf(vrna_fold_compound_t *fc);


int **get_gquad_L_matrix(short        *S,
                         int          start,
                         int          maxdist,
                         int          n,
                         int          **g,
                         vrna_param_t *P);


void
vrna_gquad_mx_local_update(vrna_fold_compound_t *fc,
                           int                  start);


void get_gquad_pattern_mfe(short        *S,
                           int          i,
                           int          j,
                           vrna_param_t *P,
                           int          *L,
                           int          l[3]);


void
get_gquad_pattern_exhaustive(short        *S,
                             int          i,
                             int          j,
                             vrna_param_t *P,
                             int          *L,
                             int          *l,
                             int          threshold);


void get_gquad_pattern_pf(short             *S,
                          int               i,
                          int               j,
                          vrna_exp_param_t  *pf,
                          int               *L,
                          int               l[3]);


plist *get_plist_gquad_from_pr(short            *S,
                               int              gi,
                               int              gj,
                               vrna_smx_csr(FLT_OR_DBL)  *q_gq,
                               FLT_OR_DBL       *probs,
                               FLT_OR_DBL       *scale,
                               vrna_exp_param_t *pf);

vrna_ep_t *
vrna_plist_gquad_from_pr(vrna_fold_compound_t *fc,
                         int               gi,
                         int               gj);


plist *get_plist_gquad_from_pr_max(short            *S,
                                   int              gi,
                                   int              gj,
                                   vrna_smx_csr(FLT_OR_DBL)  *q_gq,
                                   FLT_OR_DBL       *probs,
                                   FLT_OR_DBL       *scale,
                                   int              *L,
                                   int              l[3],
                                   vrna_exp_param_t *pf);


plist *get_plist_gquad_from_db(const char *structure,
                               float      pr);


vrna_ep_t *
vrna_plist_gquad_from_pr_max(vrna_fold_compound_t *fc,
                             int                  gi,
                             int                  gj,
                             int                  *Lmax,
                             int                  lmax[3]);


int         get_gquad_count(short *S,
                            int   i,
                            int   j);


int         get_gquad_layer_count(short *S,
                                  int   i,
                                  int   j);


void get_gquad_pattern_mfe_ali(short        **S,
                               unsigned int **a2s,
                               short        *S_cons,
                               int          n_seq,
                               int          i,
                               int          j,
                               vrna_param_t *P,
                               int          *L,
                               int          l[3]);


/**
 *  given a dot-bracket structure (possibly) containing gquads encoded
 *  by '+' signs, find first gquad, return end position or 0 if none found
 *  Upon return L and l[] contain the number of stacked layers, as well as
 *  the lengths of the linker regions.
 *  To parse a string with many gquads, call parse_gquad repeatedly e.g.
 *  end1 = parse_gquad(struc, &L, l); ... ;
 *  end2 = parse_gquad(struc+end1, &L, l); end2+=end1; ... ;
 *  end3 = parse_gquad(struc+end2, &L, l); end3+=end2; ... ;
 */
int parse_gquad(const char  *struc,
                int         *L,
                int         l[3]);




INLINE PRIVATE int backtrack_GQuad_IntLoop_L(int          c,
                                             int          i,
                                             int          j,
                                             int          type,
                                             short        *S,
                                             int          **ggg,
                                             int          maxdist,
                                             int          *p,
                                             int          *q,
                                             vrna_param_t *P);


PRIVATE INLINE int
vrna_BT_gquad_int(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  int                   en,
                  vrna_bp_stack_t       *bp_stack,
                  unsigned int          *stack_count);


PRIVATE INLINE int
vrna_bt_gquad_int(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  int                   en,
                  vrna_bps_t            bp_stack);


PRIVATE INLINE int
vrna_bt_gquad_mfe(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  vrna_bps_t            bp_stack)
{
  /*
   * here we do some fancy stuff to backtrace the stacksize and linker lengths
   * of the g-quadruplex that should reside within position i,j
   */
  short         *S;
  int           l[3], L, a, n_seq;
  vrna_param_t  *P;

  if (fc) {
    P = fc->params;
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        S = fc->sequence_encoding2;
        L = -1;

        get_gquad_pattern_mfe(S, i, j, P, &L, l);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        n_seq = fc->n_seq;
        L     = -1;
        get_gquad_pattern_mfe_ali(fc->S, fc->a2s, fc->S_cons, n_seq, i, j, P, &L, l);
        break;
    }

    if (L != -1) {
      /* fill the G's of the quadruplex into base_pair2 */
      for (a = 0; a < L; a++) {
        vrna_bps_push(bp_stack,
                      (vrna_bp_t){
                        .i = i + a,
                        .j = i + a
                      });
        vrna_bps_push(bp_stack,
                      (vrna_bp_t){
                        .i  = i + L + l[0] + a,
                        .j    = i + L + l[0] + a
                      });
        vrna_bps_push(bp_stack,
                      (vrna_bp_t){
                        .i  = i + L + l[0] + L + l[1] + a,
                        .j    = i + L + l[0] + L + l[1] + a
                      });
        vrna_bps_push(bp_stack,
                      (vrna_bp_t){
                        .i  = i + L + l[0] + L + l[1] + L + l[2] + a,
                        .j    = i + L + l[0] + L + l[1] + L + l[2] + a
                      });
      }

      return 1;
    } else {
      return 0;
    }
  }

  return 0;
}


PRIVATE INLINE int
vrna_BT_gquad_mfe(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  vrna_bp_stack_t       *bp_stack,
                  unsigned int          *stack_count)
{
  int r = 0;

  if ((fc) &&
      (bp_stack) &&
      (stack_count)) {
    vrna_bps_t bps = vrna_bps_init(4);
    r = vrna_bt_gquad_mfe(fc, i, j, bps);

    while (vrna_bps_size(bps) > 0) {
      vrna_bp_t bp = vrna_bps_pop(bps);
      bp_stack[++(*stack_count)].i = bp.i;
      bp_stack[*stack_count].j = bp.j;
    }

    vrna_bps_free(bps);
  }

  return r;
}


PRIVATE INLINE int
vrna_bt_gquad_int(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  int                   en,
                  vrna_bps_t            bp_stack)
{
  int           energy, e_gq, dangles, p, q, maxl, minl, c0, l1, u1, u2;
  unsigned char type;
  unsigned int  **a2s, n_seq, s;
  char          *ptype;
  short         si, sj, *S, *S1, **SS, **S5, **S3;
  vrna_smx_csr(int) *c_gq;

  vrna_param_t  *P;
  vrna_md_t     *md;

  n_seq   = fc->n_seq;
  S       = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;
  S1      = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : fc->S_cons;
  SS      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s     = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  c_gq    = fc->matrices->c_gq;
  P       = fc->params;
  md      = &(P->model_details);
  dangles = md->dangles;
  si      = S1[i + 1];
  sj      = S1[j - 1];
  energy  = 0;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      type    = vrna_get_ptype_md(S[i], S[j], md);
      if (dangles == 2)
        energy += P->mismatchI[type][si][sj];

      if (type > 2)
        energy += P->TerminalAU;
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      for (s = 0; s < n_seq; s++) {
        type = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
        if (md->dangles == 2)
          energy += P->mismatchI[type][S3[s][i]][S5[s][j]];

        if (type > 2)
          energy += P->TerminalAU;
      }
      break;

    default:
      return INF;
  }

  p = i + 1;
  if (S1[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(c0, minl);
      c0    = j - 3;
      maxl  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(c0, maxl);
      for (q = minl; q < maxl; q++) {
        if (S1[q] != 3)
          continue;
#ifndef VRNA_DISABLE_C11_FEATURES
        e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
        e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif

        if (e_gq != INF) {
          c0 = energy + e_gq;

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              c0  += P->internal_loop[j - q - 1];
              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                u1  = a2s[s][j - 1] - a2s[s][q];
                c0 += P->internal_loop[u1];
              }
              break;
          }

          if (en == c0)
            return vrna_bt_gquad_mfe(fc, p, q, bp_stack);
        }
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S1[p] != 3)
      continue;

    minl  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(c0, minl);
    c0    = j - 1;
    maxl  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxl  = MIN2(c0, maxl);
    for (q = minl; q < maxl; q++) {
      if (S1[q] != 3)
        continue;

#ifndef VRNA_DISABLE_C11_FEATURES
      e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
      e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif

      if (e_gq != INF) {
        c0  = energy + e_gq;

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            c0 += P->internal_loop[l1 + j - q - 1];
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              u1  = a2s[s][p - 1] - a2s[s][i];
              u2  = a2s[s][j - 1] - a2s[s][q];
              c0 += P->internal_loop[u1 + u2];
            }
            break;
        }

        if (en == c0)
          return vrna_bt_gquad_mfe(fc, p, q, bp_stack);
      }
    }
  }

  q = j - 1;
  if (S1[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S1[p] != 3)
        continue;

#ifndef VRNA_DISABLE_C11_FEATURES
      e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
      e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#endif

      if (e_gq != INF) {
        c0  = energy + e_gq;

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            c0  += P->internal_loop[l1];
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              u1  = a2s[s][p - 1] - a2s[s][i];
              c0 += P->internal_loop[u1];
            }
            break;
        }

        if (en == c0)
          return vrna_bt_gquad_mfe(fc, p, q, bp_stack);
      }
    }

  return 0;
}


PRIVATE INLINE int
vrna_BT_gquad_int(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  int                   en,
                  vrna_bp_stack_t       *bp_stack,
                  unsigned int          *stack_count)
{
  int r = 0;

  if ((fc) &&
      (bp_stack) &&
      (stack_count)) {
    vrna_bps_t bps = vrna_bps_init(4);
    r = vrna_bt_gquad_int(fc, i, j, en, bps);
    while (vrna_bps_size(bps) > 0) {
      vrna_bp_t bp = vrna_bps_pop(bps);
      bp_stack[++(*stack_count)].i = bp.i;
      bp_stack[*stack_count].j = bp.j;
    }
    
    vrna_bps_free(bps);
  }

  return r;
}


/**
 *  backtrack an interior loop like enclosed g-quadruplex
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
INLINE PRIVATE int
backtrack_GQuad_IntLoop_L(int           c,
                          int           i,
                          int           j,
                          int           type,
                          short         *S,
                          int           **ggg,
                          int           maxdist,
                          int           *p,
                          int           *q,
                          vrna_param_t  *P)
{
  int   energy, dangles, k, l, maxl, minl, c0, l1;
  short si, sj;

  dangles = P->model_details.dangles;
  si      = S[i + 1];
  sj      = S[j - 1];
  energy  = 0;

  if (dangles == 2)
    energy += P->mismatchI[type][si][sj];

  if (type > 2)
    energy += P->TerminalAU;

  k = i + 1;
  if (S[k] == 3) {
    if (k < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - i + k - MAXLOOP - 2;
      c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(c0, minl);
      c0    = j - 3;
      maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(c0, maxl);
      for (l = minl; l < maxl; l++) {
        if (S[l] != 3)
          continue;

        if (c == energy + ggg[k][l - k] + P->internal_loop[j - l - 1]) {
          *p  = k;
          *q  = l;
          return 1;
        }
      }
    }
  }

  for (k = i + 2;
       k < j - VRNA_GQUAD_MIN_BOX_SIZE;
       k++) {
    l1 = k - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S[k] != 3)
      continue;

    minl  = j - i + k - MAXLOOP - 2;
    c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(c0, minl);
    c0    = j - 1;
    maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxl  = MIN2(c0, maxl);
    for (l = minl; l < maxl; l++) {
      if (S[l] != 3)
        continue;

      if (c == energy + ggg[k][l - k] + P->internal_loop[l1 + j - l - 1]) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }
  }

  l = j - 1;
  if (S[l] == 3)
    for (k = i + 4;
         k < j - VRNA_GQUAD_MIN_BOX_SIZE;
         k++) {
      l1 = k - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S[k] != 3)
        continue;

      if (c == energy + ggg[k][l - k] + P->internal_loop[l1]) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }

  return 0;
}


INLINE PRIVATE int
backtrack_GQuad_IntLoop_L_comparative(int           c,
                                      int           i,
                                      int           j,
                                      unsigned int  *type,
                                      short         *S_cons,
                                      short         **S5,
                                      short         **S3,
                                      unsigned int  **a2s,
                                      int           **ggg,
                                      int           *p,
                                      int           *q,
                                      int           n_seq,
                                      vrna_param_t  *P)
{
  /*
   * The case that is handled here actually resembles something like
   * an interior loop where the enclosing base pair is of regular
   * kind and the enclosed pair is not a canonical one but a g-quadruplex
   * that should then be decomposed further...
   */
  int mm, dangle_model, k, l, maxl, minl, c0, l1, ss, tt, eee, u1, u2;

  dangle_model = P->model_details.dangles;

  mm = 0;
  for (ss = 0; ss < n_seq; ss++) {
    tt = type[ss];

    if (dangle_model == 2)
      mm += P->mismatchI[tt][S3[ss][i]][S5[ss][j]];

    if (tt > 2)
      mm += P->TerminalAU;
  }

  for (k = i + 2;
       k < j - VRNA_GQUAD_MIN_BOX_SIZE;
       k++) {
    if (S_cons[k] != 3)
      continue;

    l1 = k - i - 1;
    if (l1 > MAXLOOP)
      break;

    minl  = j - i + k - MAXLOOP - 2;
    c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(c0, minl);
    c0    = j - 1;
    maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxl  = MIN2(c0, maxl);
    for (l = minl; l < maxl; l++) {
      if (S_cons[l] != 3)
        continue;

      eee = 0;

      for (ss = 0; ss < n_seq; ss++) {
        u1  = a2s[ss][k - 1] - a2s[ss][i];
        u2  = a2s[ss][j - 1] - a2s[ss][l];
        eee += P->internal_loop[u1 + u2];
      }

      c0 = mm +
           ggg[k][l - k] +
           eee;

      if (c == c0) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }
  }
  k = i + 1;
  if (S_cons[k] == 3) {
    if (k < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - i + k - MAXLOOP - 2;
      c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(c0, minl);
      c0    = j - 3;
      maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(c0, maxl);
      for (l = minl; l < maxl; l++) {
        if (S_cons[l] != 3)
          continue;

        eee = 0;

        for (ss = 0; ss < n_seq; ss++) {
          u1  = a2s[ss][j - 1] - a2s[ss][l];
          eee += P->internal_loop[u1];
        }

        if (c == mm + ggg[k][l - k] + eee) {
          *p  = k;
          *q  = l;
          return 1;
        }
      }
    }
  }

  l = j - 1;
  if (S_cons[l] == 3) {
    for (k = i + 4; k < j - VRNA_GQUAD_MIN_BOX_SIZE; k++) {
      l1 = k - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S_cons[k] != 3)
        continue;

      eee = 0;

      for (ss = 0; ss < n_seq; ss++) {
        u1  = a2s[ss][k - 1] - a2s[ss][i];
        eee += P->internal_loop[u1];
      }

      if (c == mm + ggg[k][l - k] + eee) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }
  }

  return 0;
}


PRIVATE INLINE
int
vrna_E_gq_intLoop(vrna_fold_compound_t *fc,
                  int           i,
                  int           j)
{
  unsigned int      type, s, n_seq, **a2s;
  int               energy, ge, e_gq, dangles, p, q, l1, minq, maxq, c0, u1, u2;
  short             *S, *S1, si, sj, **SS, **S5, **S3;
  vrna_param_t      *P;
  vrna_md_t         *md;
  vrna_smx_csr(int) *c_gq;

  n_seq   = fc->n_seq;
  S       = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;
  S1      = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : fc->S_cons;
  SS      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s     = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  c_gq    = fc->matrices->c_gq;
  P       = fc->params;
  md      = &(P->model_details);
  dangles = md->dangles;
  si      = S1[i + 1];
  sj      = S1[j - 1];
  energy  = 0;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      type    = vrna_get_ptype_md(S[i], S[j], md);
      if (dangles == 2)
        energy += P->mismatchI[type][si][sj];

      if (type > 2)
        energy += P->TerminalAU;
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      for (s = 0; s < n_seq; s++) {
        type = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
        if (md->dangles == 2)
          energy += P->mismatchI[type][S3[s][i]][S5[s][j]];

        if (type > 2)
          energy += P->TerminalAU;
      }
      break;

    default:
      return INF;
  }

  ge = INF;

  p = i + 1;
  if (S1[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 3;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S1[q] != 3)
          continue;

#ifndef VRNA_DISABLE_C11_FEATURES
        e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
        e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif

        if (e_gq != INF) {
          c0 = energy + e_gq;

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              c0  += P->internal_loop[j - q - 1];
              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                u1  = a2s[s][j - 1] - a2s[s][q];
                c0 += P->internal_loop[u1];
              }

              break;
          }

          ge  = MIN2(ge, c0);
        }
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S1[p] != 3)
      continue;

    minq  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minq  = MAX2(c0, minq);
    c0    = j - 1;
    maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxq  = MIN2(c0, maxq);
    for (q = minq; q < maxq; q++) {
      if (S1[q] != 3)
        continue;

#ifndef VRNA_DISABLE_C11_FEATURES
      e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
      e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif

      if (e_gq != INF) {
        c0  = energy + e_gq;

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            c0 += P->internal_loop[l1 + j - q - 1];
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              u1  = a2s[s][p - 1] - a2s[s][i];
              u2  = a2s[s][j - 1] - a2s[s][q];
              c0 += P->internal_loop[u1 + u2];
            }
            break;
        }

        ge  = MIN2(ge, c0);
      }
    }
  }

  q = j - 1;
  if (S1[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S1[p] != 3)
        continue;

#ifndef VRNA_DISABLE_C11_FEATURES
      e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
      e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#endif

      if (e_gq != INF) {
        c0  = energy + e_gq;

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            c0  += P->internal_loop[l1];
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              u1  = a2s[s][p - 1] - a2s[s][i];
              c0 += P->internal_loop[u1];
            }
            break;
        }

        ge  = MIN2(ge, c0);
      }
    }

  return ge;
}


PRIVATE INLINE
int
E_GQuad_IntLoop_L_comparative(int           i,
                              int           j,
                              unsigned int  *tt,
                              short         *S_cons,
                              short         **S5,
                              short         **S3,
                              unsigned int  **a2s,
                              int           **ggg,
                              int           n_seq,
                              vrna_param_t  *P)
{
  unsigned int  type;
  int           eee, energy, ge, p, q, l1, u1, u2, minq, maxq, c0, s;
  vrna_md_t     *md;

  md      = &(P->model_details);
  energy  = 0;

  for (s = 0; s < n_seq; s++) {
    type = tt[s];
    if (md->dangles == 2)
      energy += P->mismatchI[type][S3[s][i]][S5[s][j]];

    if (type > 2)
      energy += P->TerminalAU;
  }

  ge = INF;

  p = i + 1;
  if (S_cons[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 3;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S_cons[q] != 3)
          continue;

        eee = 0;

        for (s = 0; s < n_seq; s++) {
          u1  = a2s[s][j - 1] - a2s[s][q];
          eee += P->internal_loop[u1];
        }

        c0 = energy +
             ggg[p][q - p] +
             eee;
        ge = MIN2(ge, c0);
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S_cons[p] != 3)
      continue;

    minq  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minq  = MAX2(c0, minq);
    c0    = j - 1;
    maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxq  = MIN2(c0, maxq);
    for (q = minq; q < maxq; q++) {
      if (S_cons[q] != 3)
        continue;

      eee = 0;

      for (s = 0; s < n_seq; s++) {
        u1  = a2s[s][p - 1] - a2s[s][i];
        u2  = a2s[s][j - 1] - a2s[s][q];
        eee += P->internal_loop[u1 + u2];
      }

      c0 = energy +
           ggg[p][q - p] +
           eee;
      ge = MIN2(ge, c0);
    }
  }

  q = j - 1;
  if (S_cons[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S_cons[p] != 3)
        continue;

      eee = 0;

      for (s = 0; s < n_seq; s++) {
        u1  = a2s[s][p - 1] - a2s[s][i];
        eee += P->internal_loop[u1];
      }

      c0 = energy +
           ggg[p][q - p] +
           eee;
      ge = MIN2(ge, c0);
    }

  return ge;
}


PRIVATE INLINE
int *
vrna_E_gq_intLoop_exhaustive(vrna_fold_compound_t *fc,
                            int          i,
                            int          j,
                           int          **p_p,
                           int          **q_p,
                           int          threshold)
{
  int   type;
  int   energy, *ge, e_gq, dangles, p, q, l1, minq, maxq, c0;
  short *S, *S1, si, sj;
  int   cnt = 0;

  vrna_param_t  *P = fc->params;
  vrna_md_t     *md = &(P->model_details);

  vrna_smx_csr(int) *c_gq = fc->matrices->c_gq;

  S       = fc->sequence_encoding2;
  S1       = fc->sequence_encoding;
  type    = vrna_get_ptype_md(S[i], S[j], md);
  dangles = md->dangles;
  si      = S[i + 1];
  sj      = S[j - 1];
  energy  = 0;

  if (dangles == 2)
    energy += P->mismatchI[type][si][sj];

  if (type > 2)
    energy += P->TerminalAU;

  /* guess how many gquads are possible in interval [i+1,j-1] */
  *p_p  = (int *)vrna_alloc(sizeof(int) * 256);
  *q_p  = (int *)vrna_alloc(sizeof(int) * 256);
  ge    = (int *)vrna_alloc(sizeof(int) * 256);

  p = i + 1;
  if (S[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 3;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S[q] != 3)
          continue;

#ifndef VRNA_DISABLE_C11_FEATURES
        e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
        e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif

        if (e_gq != INF) {
          c0 = energy + e_gq + P->internal_loop[j - q - 1];

          if (c0 <= threshold) {
            ge[cnt]       = energy + P->internal_loop[j - q - 1];
            (*p_p)[cnt]   = p;
            (*q_p)[cnt++] = q;
          }
        }
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S[p] != 3)
      continue;

    minq  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minq  = MAX2(c0, minq);
    c0    = j - 1;
    maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxq  = MIN2(c0, maxq);
    for (q = minq; q < maxq; q++) {
      if (S[q] != 3)
        continue;

#ifndef VRNA_DISABLE_C11_FEATURES
      e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
      e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif

      if (e_gq != INF) {
        c0 = energy + e_gq + P->internal_loop[l1 + j - q - 1];

        if (c0 <= threshold) {
          ge[cnt]       = energy + P->internal_loop[l1 + j - q - 1];
          (*p_p)[cnt]   = p;
          (*q_p)[cnt++] = q;
        }
      }
    }
  }

  q = j - 1;
  if (S[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S[p] != 3)
        continue;

#ifndef VRNA_DISABLE_C11_FEATURES
      e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
      e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif

      if (e_gq != INF) {
        c0 = energy + e_gq + P->internal_loop[l1];
        if (c0 <= threshold) {
          ge[cnt]       = energy + P->internal_loop[l1];
          (*p_p)[cnt]   = p;
          (*q_p)[cnt++] = q;
        }
      }
    }

  (*p_p)[cnt] = -1;

  return ge;
}


PRIVATE INLINE
int
E_GQuad_IntLoop_L(int           i,
                  int           j,
                  int           type,
                  short         *S,
                  int           **ggg,
                  int           maxdist,
                  vrna_param_t  *P)
{
  int   energy, ge, dangles, p, q, l1, minq, maxq, c0;
  short si, sj;

  dangles = P->model_details.dangles;
  si      = S[i + 1];
  sj      = S[j - 1];
  energy  = 0;

  if (dangles == 2)
    energy += P->mismatchI[type][si][sj];

  if (type > 2)
    energy += P->TerminalAU;

  ge = INF;

  p = i + 1;
  if (S[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 3;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S[q] != 3)
          continue;

        c0  = energy + ggg[p][q - p] + P->internal_loop[j - q - 1];
        ge  = MIN2(ge, c0);
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S[p] != 3)
      continue;

    minq  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minq  = MAX2(c0, minq);
    c0    = j - 1;
    maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxq  = MIN2(c0, maxq);
    for (q = minq; q < maxq; q++) {
      if (S[q] != 3)
        continue;

      c0  = energy + ggg[p][q - p] + P->internal_loop[l1 + j - q - 1];
      ge  = MIN2(ge, c0);
    }
  }

  q = j - 1;
  if (S[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S[p] != 3)
        continue;

      c0  = energy + ggg[p][q - p] + P->internal_loop[l1];
      ge  = MIN2(ge, c0);
    }

  return ge;
}


PRIVATE INLINE
FLT_OR_DBL
vrna_exp_E_gq_intLoop(vrna_fold_compound_t *fc,
                      int               i,
                      int               j)
{
  unsigned int      type, s, n_seq, **a2s;
  int               k, l, minl, maxl, u, u1, u2, r;
  FLT_OR_DBL        q, qe, q_g;
  double            *expintern;
  short             *S, *S1, si, sj, **SS, **S5, **S3;
  FLT_OR_DBL        *scale;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_smx_csr(FLT_OR_DBL) *q_gq;

  n_seq     = fc->n_seq;                                                               
  S         = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;       
  S1        = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : fc->S_cons;  
  SS        = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;                        
  S5        = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;                       
  S3        = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;                       
  a2s       = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;                      
  q_gq      = fc->exp_matrices->q_gq;                                                  
  scale     = fc->exp_matrices->scale;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  dangles   = md->dangles;
  si        = S1[i + 1];
  sj        = S1[j - 1];

  qe = 1.;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      type    = vrna_get_ptype_md(S[i], S[j], md);
      if (dangles == 2)
        qe *=  (FLT_OR_DBL)pf_params->expmismatchI[type][si][sj];

      if (type > 2)
        qe *= (FLT_OR_DBL)pf_params->expTermAU;
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      for (s = 0; s < n_seq; s++) {
        type = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
        if (md->dangles == 2)
          qe *= (FLT_OR_DBL)pf_params->expmismatchI[type][S3[s][i]][S5[s][j]];

        if (type > 2)
          qe *= (FLT_OR_DBL)pf_params->expTermAU;
      }
      break;

    default:
      return 0.;
  }

  expintern = &(pf_params->expinternal[0]);
  q         = 0;
  k         = i + 1;

  if (S1[k] == 3) {
    if (k < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - MAXLOOP - 1;
      u     = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(u, minl);
      u     = j - 3;
      maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(u, maxl);
      for (l = minl; l < maxl; l++) {
        if (S1[l] != 3)
          continue;

#ifndef VRNA_DISABLE_C11_FEATURES
        q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
        q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

        if (q_g != 0.) {
          q_g *= qe * scale[j - l + 1];

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              q_g *= (FLT_OR_DBL)expintern[j - l - 1];
              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                u1  = a2s[s][j - 1] - a2s[s][l];
                q_g *= (FLT_OR_DBL)expintern[u1];
              }

              break;
          }

          q += q_g;
        }

      }
    }
  }

  for (k = i + 2;
       k <= j - VRNA_GQUAD_MIN_BOX_SIZE;
       k++) {
    u = k - i - 1;
    if (u > MAXLOOP)
      break;

    if (S1[k] != 3)
      continue;

    minl  = j - i + k - MAXLOOP - 2;
    r     = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(r, minl);
    maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    r     = j - 1;
    maxl  = MIN2(r, maxl);
    for (l = minl; l < maxl; l++) {
      if (S1[l] != 3)
        continue;

#ifndef VRNA_DISABLE_C11_FEATURES
      q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
      q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

      if (q_g != 0.) {
        q_g *= qe * scale[u + j - l + 1];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            q_g *= (FLT_OR_DBL)expintern[u + j - l - 1];
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              u1  = a2s[s][k - 1] - a2s[s][i];
              u2  = a2s[s][j - 1] - a2s[s][l];
              q_g *= (FLT_OR_DBL)expintern[u1 + u2];
            }
            break;
        }

        q += q_g;
      }
    }
  }

  l = j - 1;
  if (S1[l] == 3)
    for (k = i + 4; k <= j - VRNA_GQUAD_MIN_BOX_SIZE; k++) {
      u = k - i - 1;
      if (u > MAXLOOP)
        break;

      if (S1[k] != 3)
        continue;

#ifndef VRNA_DISABLE_C11_FEATURES
      q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
      q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

      if (q_g != 0.) {
        q_g *= qe * scale[u + 2];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            q_g *= (FLT_OR_DBL)expintern[u];
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              u1  = a2s[s][k - 1] - a2s[s][i];
              q_g *= (FLT_OR_DBL)expintern[u1];
            }
            break;
        }

        q += q_g;
      }
    }

  return q;
}


/**
 * @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

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
DEPRECATED(int *get_gquad_matrix(short *S, vrna_param_t *P),
           "Use vrna_gq_pos_mfe() instead");

DEPRECATED(int *get_gquad_ali_matrix(unsigned int  n,
                          short         *S_cons,
                          short         **S,
                          unsigned int  **a2s,
                          int           n_seq,
                          vrna_param_t  *P),
           "Use vrna_gq_pos_mfe() instead");


DEPRECATED(FLT_OR_DBL *get_gquad_pf_matrix(short             *S,
                                FLT_OR_DBL        *scale,
                                vrna_exp_param_t  *pf),
           "Use vrna_gq_pos_pf() instead");


DEPRECATED(FLT_OR_DBL *get_gquad_pf_matrix_comparative(unsigned int  n,
                                            short             *S_cons,
                                            short             **S,
                                            unsigned int      **a2s,
                                            FLT_OR_DBL        *scale,
                                            unsigned int      n_seq,
                                            vrna_exp_param_t  *pf),
           "Use vrna_gq_pos_pf() instead");


#endif

#endif
