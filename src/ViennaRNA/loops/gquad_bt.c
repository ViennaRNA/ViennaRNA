/*
 * gquad.c
 *
 * Ronny Lorenz 2012
 *
 * ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/params/constants.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/loops/gquad.h"

#include "ViennaRNA/loops/gquad_intern.h"


#ifndef INLINE
#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif
#endif


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/**
 *  MFE callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_mfe(int   i,
          int   L,
          int   *l,
          void  *data,
          void  *P,
          void  *NA,
          void  *NA2);


PRIVATE
void
gquad_mfe_pos(int   i,
              int   L,
              int   *l,
              void  *data,
              void  *P,
              void  *Lmfe,
              void  *lmfe);


PRIVATE void
gquad_mfe_ali_pos(int   i,
                  int   L,
                  int   *l,
                  void  *data,
                  void  *helper,
                  void  *Lmfe,
                  void  *lmfe);


PRIVATE
void
gquad_pos_exhaustive(int  i,
                     int  L,
                     int  *l,
                     void *data,
                     void *P,
                     void *Lex,
                     void *lex);


/**
 * Partition function callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_pf(int  i,
         int  L,
         int  *l,
         void *data,
         void *P,
         void *NA,
         void *NA2);


PRIVATE void
gquad_pf_ali(int  i,
             int  L,
             int  *l,
             void *data,
             void *helper,
             void *NA,
             void *NA2);


/**
 * Partition function callback for process_gquad_enumeration()
 * in contrast to gquad_pf() it stores the stack size L and
 * the linker lengths l[3] of the g-quadruplex that dominates
 * the interval [i,j]
 * (FLT_OR_DBL *)data must be 0. on entry
 */
PRIVATE
void
gquad_pf_pos(int  i,
             int  L,
             int  *l,
             void *data,
             void *pf,
             void *Lmax,
             void *lmax);


PRIVATE void
gquad_pf_pos_ali(int  i,
                 int  L,
                 int  *l,
                 void *data,
                 void *helper,
                 void *NA1,
                 void *NA2);


/**
 * MFE (alifold) callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_mfe_ali(int   i,
              int   L,
              int   *l,
              void  *data,
              void  *helper,
              void  *NA,
              void  *NA2);


/**
 * MFE (alifold) callback for process_gquad_enumeration()
 * with seperation of free energy and penalty contribution
 */
PRIVATE
void
gquad_mfe_ali_en(int  i,
                 int  L,
                 int  *l,
                 void *data,
                 void *helper,
                 void *NA,
                 void *NA2);


PRIVATE
void
gquad_interact(int  i,
               int  L,
               int  *l,
               void *data,
               void *pf,
               void *index,
               void *NA2);


PRIVATE
void
gquad_interact_ali(int  i,
                   int  L,
                   int  *l,
                   void *data,
                   void *index,
                   void *helper,
                   void *NA);


PRIVATE
void
gquad_count(int   i,
            int   L,
            int   *l,
            void  *data,
            void  *NA,
            void  *NA2,
            void  *NA3);


PRIVATE
void
gquad_count_layers(int  i,
                   int  L,
                   int  *l,
                   void *data,
                   void *NA,
                   void *NA2,
                   void *NA3);


/* other useful static functions */

PRIVATE int
E_gquad_consensus(int                 L,
                  int                 l[3],
                  unsigned int        position,
                  unsigned int        length,
                  unsigned int        n_seq,
                  const unsigned int  **a2s,
                  vrna_param_t        *P);


PRIVATE
int
E_gquad_ali_penalty(int           L,
                    int           l[3],
                    unsigned int  i,
                    unsigned int  length,
                    unsigned int  n_seq,
                    const short   **S,
                    vrna_param_t  *P);


PRIVATE
FLT_OR_DBL
exp_E_gquad_ali_penalty(int               L,
                        int               l[3],
                        unsigned int      i,
                        unsigned int      length,
                        unsigned int      n_seq,
                        const short       **S,
                        vrna_exp_param_t  *P);


PRIVATE void
count_gquad_layer_mismatches(int          L,
                             int          l[3],
                             unsigned int i,
                             unsigned int n,
                             unsigned int n_seq,
                             const short  **S,
                             unsigned int mm[2]);


PRIVATE void
gquad_mfe_ali_pos(int   i,
                  int   L,
                  int   *l,
                  void  *data,
                  void  *P,
                  void  *Lmfe,
                  void  *lmfe);


/*
 #########################################
 # BEGIN OF PUBLIC FUNCTION DEFINITIONS  #
 #      (all available in RNAlib)        #
 #########################################
 */

PUBLIC int
vrna_bt_gquad_mfe(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  vrna_bps_t            bp_stack)
{
  /*
   * here we do some fancy stuff to backtrace the stacksize and linker lengths
   * of the g-quadruplex that should reside within position i,j
   */
  short         *S_enc, *S_tmp;
  unsigned int  n, n2;
  int           l[3], L, a, n_seq;
  vrna_param_t  *P;

  if (fc) {
    n     = fc->length;
    P     = fc->params;
    L     = -1;
    S_tmp = NULL;

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      S_enc = fc->S_cons;
      n_seq = fc->n_seq;
    } else {
      S_enc = fc->sequence_encoding2;
      n_seq = 1;
    }

    if (P->model_details.circ) {
      n2 = MIN2(n, VRNA_GQUAD_MAX_BOX_SIZE) - 1;

      S_tmp = (short *)vrna_alloc(sizeof(short) * (n + n2 + 1));
      memcpy(S_tmp, S_enc, sizeof(short) * (n + 1));
      memcpy(S_tmp + (n + 1), S_enc + 1, sizeof(short) * n2);
      S_tmp[0]  = n + n2;
      S_enc     = S_tmp;
      if (j < i)
        j += n;
    }

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      get_gquad_pattern_mfe_ali(fc->S, fc->a2s, fc->S_cons, n_seq, i, j, P, &L, l);
    } else {
      get_gquad_pattern_mfe(S_enc, i, j, P, &L, l);
    }

    if (L != -1) {
      /* fill the G's of the quadruplex into bp_stack */
      for (a = 0; a < L; a++) {
        int p1, p2, p3, p4;
        p1  = i + a;
        p2  = p1 + L + l[0];
        p3  = p2 + L + l[1];
        p4  = p3 + L + l[2];
        if (p1 > n) {
          p1  = ((p1 - 1) % n) + 1;
          p2  = ((p2 - 1) % n) + 1;
          p3  = ((p3 - 1) % n) + 1;
          p4  = ((p4 - 1) % n) + 1;
        } else if (p2 > n) {
          p2  = ((p2 - 1) % n) + 1;
          p3  = ((p3 - 1) % n) + 1;
          p4  = ((p4 - 1) % n) + 1;
        } else if (p3 > n) {
          p3  = ((p3 - 1) % n) + 1;
          p4  = ((p4 - 1) % n) + 1;
        } else if (p4 > n) {
          p4 = ((p4 - 1) % n) + 1;
        }

        vrna_bps_push(bp_stack,
                      (vrna_bp_t){
          .i  = p1,
          .j  = p1
        });
        vrna_bps_push(bp_stack,
                      (vrna_bp_t){
          .i  = p2,
          .j  = p2
        });
        vrna_bps_push(bp_stack,
                      (vrna_bp_t){
          .i  = p3,
          .j  = p3
        });
        vrna_bps_push(bp_stack,
                      (vrna_bp_t){
          .i  = p4,
          .j  = p4
        });
      }
      free(S_tmp);
      return 1;
    } else {
      free(S_tmp);
      return 0;
    }
  }

  return 0;
}


PUBLIC int
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
      bp_stack[++(*stack_count)].i  = bp.i;
      bp_stack[*stack_count].j      = bp.j;
    }

    vrna_bps_free(bps);
  }

  return r;
}


PUBLIC int
vrna_bt_gquad_int(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j,
                  int                   en,
                  vrna_bps_t            bp_stack,
                  vrna_bts_t            bt_stack)
{
  unsigned char type;
  short         si, sj, *S, *S1, **SS, **S5, **S3;
  int           energy, e_gq, dangles, c0;
  unsigned int  **a2s, n_seq, s, p, q, l1, u1, u2, maxl, minl;

  vrna_smx_csr(int) * c_gq;

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
      type = vrna_get_ptype_md(S[i], S[j], md);
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
    if (p + VRNA_GQUAD_MIN_BOX_SIZE < j) {
      minl = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      if (minl + 1 + MAXLOOP < j)
        minl = j - MAXLOOP - 1;

      maxl = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      if (maxl + 3 >= j)
        maxl = j - 3;

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
              c0 += P->internal_loop[j - q - 1];
              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                u1  = a2s[s][j - 1] - a2s[s][q];
                c0  += P->internal_loop[u1];
              }
              break;
          }

          if (en == c0) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = p,
                            .j = q,
                            .ml = VRNA_MX_FLAG_G}));

            return 1;
          }
        }
      }
    }
  }

  for (p = i + 2;
       p + VRNA_GQUAD_MIN_BOX_SIZE < j;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S1[p] != 3)
      continue;

    minl = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    if (minl + 1 + MAXLOOP - l1 < j)
      minl = j - MAXLOOP + l1 - 1;

    maxl = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    if (maxl >= j)
      maxl = j - 1;

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
            c0 += P->internal_loop[l1 + j - q - 1];
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              u1  = a2s[s][p - 1] - a2s[s][i];
              u2  = a2s[s][j - 1] - a2s[s][q];
              c0  += P->internal_loop[u1 + u2];
            }
            break;
        }

        if (en == c0) {
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = p,
                          .j = q,
                          .ml = VRNA_MX_FLAG_G}));

          return 1;
        }
      }
    }
  }

  q = j - 1;
  if (S1[q] == 3)
    for (p = i + 4;
         p + VRNA_GQUAD_MIN_BOX_SIZE - 1 < j;
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
        c0 = energy + e_gq;

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            c0 += P->internal_loop[l1];
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              u1  = a2s[s][p - 1] - a2s[s][i];
              c0  += P->internal_loop[u1];
            }
            break;
        }

        if (en == c0) {
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = p,
                          .j = q,
                          .ml = VRNA_MX_FLAG_G}));

          return 1;
        }
      }
    }

  return 0;
}


PUBLIC int
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
    vrna_bts_t bts = vrna_bts_init(0);

    r = vrna_bt_gquad_int(fc, i, j, en, bps, bts);

    while (vrna_bts_size(bts) > 0) {
      vrna_sect_t s= vrna_bts_pop(bts);
      r = vrna_bt_gquad_mfe(fc, s.i, s.j, bps);
    }

    while (vrna_bps_size(bps) > 0) {
      vrna_bp_t bp = vrna_bps_pop(bps);
      bp_stack[++(*stack_count)].i  = bp.i;
      bp_stack[*stack_count].j      = bp.j;
    }

    vrna_bps_free(bps);
    vrna_bts_free(bts);
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
int
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


int
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


/*
 #########################################
 # BEGIN OF PRIVATE FUNCTION DEFINITIONS #
 #          (internal use only)          #
 #########################################
 */
PRIVATE INLINE int
aln_linker_length(unsigned int        start,
                  unsigned int        end,
                  unsigned int        n,
                  const unsigned int  *a2ss)
{
  if (start <= end) {
    return a2ss[end] - a2ss[start - 1];
  } else {
    return a2ss[n] - a2ss[start - 1] + a2ss[end];
  }
}


/* compute (individual) lengths of the unpaired linker sequences */
/*  note here, that we might have a GQ spanning the n,1 junction,
 *  so we first need to transform the linker start and end
 *  positions accordingly
 */
PRIVATE INLINE void
aln_linker_positions(int          L,
                     int          l[3],
                     unsigned int position,
                     unsigned int length,
                     unsigned int starts[3],
                     unsigned int ends[3])
{
  if ((length > 0) &&
      (position + 4 * L + l[0] + l[1] + l[2] >= length)) {
    starts[0] = (position + L - 1) % (length) + 1;
    ends[0]   = (position + L + l[0] - 1 - 1) % (length) + 1;
    starts[1] = (position + 2 * L + l[0] - 1) % (length) + 1;
    ends[1]   = (position + 2 * L + l[0] + l[1] - 1 - 1) % (length) + 1;
    starts[2] = (position + 3 * L + l[0] + l[1] - 1) % (length) + 1;
    ends[2]   = (position + 3 * L + l[0] + l[1] + l[2] - 1 - 1) % (length) + 1;
  } else {
    starts[0] = position + L;
    ends[0]   = starts[0] + l[0] - 1;
    starts[1] = ends[0] + L + 1;
    ends[1]   = starts[1] + l[1] - 1;
    starts[2] = ends[1] + L + 1;
    ends[2]   = starts[2] + l[2] - 1;
  }
}


/* retrieve a set of sequence coordinates for the Gs involved
 * in a layer (1-based) of a GQ with stack size L and linker
 * lengths l starting at position i. The GQ may cross the n,1
 * junction so the total length of the sequence (alignment) has
 * to be passed through variable n
 */
PRIVATE void
gq_layer_pos(int          L,
             int          l[3],
             unsigned int layer,
             unsigned int i,
             unsigned int n,
             unsigned int layer_pos[4])
{
  if ((n > 0) &&
      (i + 4 * L + l[0] + l[1] + l[2] >= n)) {
    layer_pos[0]  = (i + layer - 1 - 1) % (n) + 1;
    layer_pos[1]  = (i + layer + L + l[0] - 1 - 1) % (n) + 1;
    layer_pos[2]  = (i + layer + 2 * L + l[0] + l[1] - 1 - 1) % (n) + 1;
    layer_pos[3]  = (i + layer + 3 * L + l[0] + l[1] + l[2] - 1 - 1) % (n) + 1;
  } else {
    layer_pos[0]  = i + layer - 1;
    layer_pos[1]  = i + layer + L + l[0] - 1;
    layer_pos[2]  = i + layer + 2 * L + l[0] + l[1] - 1;
    layer_pos[3]  = i + layer + 3 * L + l[0] + l[1] + l[2] - 1;
  }
}


PRIVATE int
E_gquad_consensus(int                 L,
                  int                 l[3],
                  unsigned int        position,
                  unsigned int        length,
                  unsigned int        n_seq,
                  const unsigned int  **a2s,
                  vrna_param_t        *P)
{
  unsigned int  l_start[3], l_end[3], s, u1, u2, u3;
  int           e;


  e = 0;

  aln_linker_positions(L, l, position, length, l_start, l_end);

  for (s = 0; s < n_seq; s++) {
    u1  = aln_linker_length(l_start[0], l_end[0], length, a2s[s]);
    u2  = aln_linker_length(l_start[1], l_end[1], length, a2s[s]);
    u3  = aln_linker_length(l_start[2], l_end[2], length, a2s[s]);

    e += P->gquad[L][u1 + u2 + u3];
  }

  return e;
}


PRIVATE FLT_OR_DBL
exp_E_gquad_consensus(int                 L,
                      int                 l[3],
                      unsigned int        position,
                      unsigned int        length,
                      unsigned int        n_seq,
                      const unsigned int  **a2s,
                      vrna_exp_param_t    *pf)
{
  unsigned int  l_start[3], l_end[3], s, u1, u2, u3;
  FLT_OR_DBL    q;


  q = 1.;

  aln_linker_positions(L, l, position, length, l_start, l_end);

  for (s = 0; s < n_seq; s++) {
    u1  = aln_linker_length(l_start[0], l_end[0], length, a2s[s]);
    u2  = aln_linker_length(l_start[1], l_end[1], length, a2s[s]);
    u3  = aln_linker_length(l_start[2], l_end[2], length, a2s[s]);

    q *= pf->expgquad[L][u1 + u2 + u3];
  }

  return q;
}


PRIVATE int
E_gquad_ali_penalty(int           L,
                    int           l[3],
                    unsigned int  i,
                    unsigned int  length,
                    unsigned int  n_seq,
                    const short   **S,
                    vrna_param_t  *P)
{
  unsigned int mm[2];

  count_gquad_layer_mismatches(L, l, i, length, n_seq, S, mm);

  if (mm[1] > P->gquadLayerMismatchMax)
    return INF;
  else
    return P->gquadLayerMismatch * mm[0];
}


PRIVATE FLT_OR_DBL
exp_E_gquad_ali_penalty(int               L,
                        int               l[3],
                        unsigned int      i,
                        unsigned int      n,
                        unsigned int      n_seq,
                        const short       **S,
                        vrna_exp_param_t  *pf)
{
  unsigned int mm[2];

  count_gquad_layer_mismatches(L, l, i, n, n_seq, S, mm);

  if (mm[1] > pf->gquadLayerMismatchMax)
    return (FLT_OR_DBL)0.;
  else
    return (FLT_OR_DBL)pow(pf->expgquadLayerMismatch, (double)mm[0]);
}


PRIVATE void
count_gquad_layer_mismatches(int          L,
                             int          l[3],
                             unsigned int i,
                             unsigned int n,
                             unsigned int n_seq,
                             const short  **S,
                             unsigned int mm[2])
{
  unsigned int  s, layer_pos[4], k, ld, mismatch;
  int           cnt;

  mm[0] = mm[1] = 0;


  /* check for compatibility in the alignment */
  for (s = 0; s < n_seq; s++) {
    mismatch = 0;

    /* check bottom layer */
    gq_layer_pos(L, l, 1, i, n, layer_pos);

    if (S[s][layer_pos[0]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[1]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[2]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[3]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    }

    /* check top layer */
    gq_layer_pos(L, l, L, i, n, layer_pos);

    if (S[s][layer_pos[0]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[1]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[2]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[3]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    }

    /* check inner layers */
    for (cnt = 2; cnt < L; cnt++) {
      gq_layer_pos(L, l, cnt, i, n, layer_pos);

      if (S[s][layer_pos[0]] != 3) {
        /* add 2x penalty for missing inner layer */
        mismatch += 2;
      } else if (S[s][layer_pos[1]] != 3) {
        /* add 2x penalty for missing inner layer */
        mismatch += 2;
      } else if (S[s][layer_pos[2]] != 3) {
        /* add 2x penalty for missing inner layer */
        mismatch += 2;
      } else if (S[s][layer_pos[3]] != 3) {
        /* add 2x penalty for missing inner layer */
        mismatch += 2;
      }
    }

    mm[0] += mismatch;

    if (mismatch >= 2 * (L - 1))
      mm[1]++;
  }
}


PRIVATE void
gquad_mfe(int   i,
          int   L,
          int   *l,
          void  *data,
          void  *P,
          void  *NA,
          void  *NA2)
{
  int cc = ((vrna_param_t *)P)->gquad[L][l[0] + l[1] + l[2]];

  if (cc < *((int *)data))
    *((int *)data) = cc;
}


PRIVATE void
gquad_mfe_pos(int   i,
              int   L,
              int   *l,
              void  *data,
              void  *P,
              void  *Lmfe,
              void  *lmfe)
{
  int cc = ((vrna_param_t *)P)->gquad[L][l[0] + l[1] + l[2]];

  if (cc < *((int *)data)) {
    *((int *)data)        = cc;
    *((int *)Lmfe)        = L;
    *((int *)lmfe)        = l[0];
    *(((int *)lmfe) + 1)  = l[1];
    *(((int *)lmfe) + 2)  = l[2];
  }
}


PRIVATE void
gquad_mfe_ali_pos(int   i,
                  int   L,
                  int   *l,
                  void  *data,
                  void  *helper,
                  void  *Lmfe,
                  void  *lmfe)
{
  int cc = INF;

  gquad_mfe_ali(i, L, l, (void *)&cc, helper, NULL, NULL);

  if (cc < *((int *)data)) {
    *((int *)data)        = cc;
    *((int *)Lmfe)        = L;
    *((int *)lmfe)        = l[0];
    *(((int *)lmfe) + 1)  = l[1];
    *(((int *)lmfe) + 2)  = l[2];
  }
}


PRIVATE
void
gquad_pos_exhaustive(int  i,
                     int  L,
                     int  *l,
                     void *data,
                     void *P,
                     void *Lex,
                     void *lex)
{
  int cnt;
  int cc = ((vrna_param_t *)P)->gquad[L][l[0] + l[1] + l[2]];

  if (cc <= *((int *)data)) {
    /*  since Lex is an array of L values and lex an
     *  array of l triples we need to find out where
     *  the current gquad position is to be stored...
     * the below implementation might be slow but we
     * still use it for now
     */
    for (cnt = 0; ((int *)Lex)[cnt] != -1; cnt++);

    *((int *)Lex + cnt)             = L;
    *((int *)Lex + cnt + 1)         = -1;
    *(((int *)lex) + (3 * cnt) + 0) = l[0];
    *(((int *)lex) + (3 * cnt) + 1) = l[1];
    *(((int *)lex) + (3 * cnt) + 2) = l[2];
  }
}


PRIVATE
void
gquad_count(int   i,
            int   L,
            int   *l,
            void  *data,
            void  *NA,
            void  *NA2,
            void  *NA3)
{
  *((int *)data) += 1;
}


PRIVATE
void
gquad_count_layers(int  i,
                   int  L,
                   int  *l,
                   void *data,
                   void *NA,
                   void *NA2,
                   void *NA3)
{
  *((int *)data) += L;
}


PRIVATE void
gquad_pf(int  i,
         int  L,
         int  *l,
         void *data,
         void *pf,
         void *NA,
         void *NA2)
{
  *((FLT_OR_DBL *)data) += ((vrna_exp_param_t *)pf)->expgquad[L][l[0] + l[1] + l[2]];
}


PRIVATE void
gquad_pf_ali(int  i,
             int  L,
             int  *l,
             void *data,
             void *helper,
             void *NA,
             void *NA2)
{
  const short             **S;
  const unsigned int      **a2s;
  unsigned int            n, s, n_seq;
  int                     u1, u2, u3;
  FLT_OR_DBL              penalty;
  vrna_exp_param_t        *pf;
  struct gquad_ali_helper *gq_help;

  gq_help = (struct gquad_ali_helper *)helper;
  S       = gq_help->S;
  a2s     = gq_help->a2s;
  n       = gq_help->length;
  n_seq   = gq_help->n_seq;
  pf      = gq_help->pf;
  penalty = exp_E_gquad_ali_penalty(L, l, (unsigned int)i, n, n_seq, S, pf);

  if (penalty != 0.)
    *((FLT_OR_DBL *)data) += penalty *
                             exp_E_gquad_consensus(L, l, (unsigned int)i, n, n_seq, a2s, pf);
}


PRIVATE void
gquad_pf_pos(int  i,
             int  L,
             int  *l,
             void *data,
             void *pf,
             void *Lmax,
             void *lmax)
{
  FLT_OR_DBL gq = 0.;

  gquad_pf(i, L, l, (void *)&gq, pf, NULL, NULL);

  if (gq > *((FLT_OR_DBL *)data)) {
    *((FLT_OR_DBL *)data) = gq;
    *((int *)Lmax)        = L;
    *((int *)lmax)        = l[0];
    *(((int *)lmax) + 1)  = l[1];
    *(((int *)lmax) + 2)  = l[2];
  }
}


PRIVATE void
gquad_pf_pos_ali(int  i,
                 int  L,
                 int  *l,
                 void *data,
                 void *helper,
                 void *NA,
                 void *NA2)
{
  FLT_OR_DBL              gq        = 0.;
  struct gquad_ali_helper *gq_help  = (struct gquad_ali_helper *)helper;

  gquad_pf_ali(i, L, l, (void *)&gq, helper, NULL, NULL);

  if (gq > *((FLT_OR_DBL *)data)) {
    *((FLT_OR_DBL *)data) = gq;
    gq_help->L            = L;
    gq_help->l[0]         = l[0];
    gq_help->l[1]         = l[1];
    gq_help->l[2]         = l[2];
  }
}


PRIVATE void
gquad_mfe_ali(int   i,
              int   L,
              int   *l,
              void  *data,
              void  *helper,
              void  *NA,
              void  *NA2)
{
  int en[2], cc;

  en[0] = en[1] = INF;

  CHECK_GQUAD(L, l, return );

  gquad_mfe_ali_en(i, L, l, (void *)(&(en[0])), helper, NULL, NULL);
  if (en[1] != INF) {
    cc = en[0] + en[1];
    if (cc < *((int *)data))
      *((int *)data) = cc;
  }
}


PRIVATE void
gquad_mfe_ali_en(int  i,
                 int  L,
                 int  *l,
                 void *data,
                 void *helper,
                 void *NA,
                 void *NA2)
{
  const short             **S;
  const unsigned int      **a2s;
  unsigned int            s, n_seq, n;
  int                     en[2], cc, dd, u1, u2, u3;
  vrna_param_t            *P;
  struct gquad_ali_helper *gq_help;

  gq_help = (struct gquad_ali_helper *)helper;
  S       = gq_help->S;
  a2s     = gq_help->a2s;
  n       = gq_help->length;
  n_seq   = gq_help->n_seq;
  P       = gq_help->P;

  en[0] = E_gquad_consensus(L, l, i, n, n_seq, (const unsigned int **)a2s, P);
  en[1] = E_gquad_ali_penalty(L, l, i, n, n_seq, S, P);

  if (en[1] != INF) {
    cc  = en[0] + en[1];
    dd  = ((int *)data)[0] + ((int *)data)[1];
    if (cc < dd) {
      ((int *)data)[0]  = en[0];
      ((int *)data)[1]  = en[1];
    }
  }
}


PRIVATE void
gquad_interact(int  i,
               int  L,
               int  *l,
               void *data,
               void *pf,
               void *index,
               void *NA2)
{
  int         x, *idx;
  FLT_OR_DBL  gq, *pp;

  idx = (int *)index;
  pp  = (FLT_OR_DBL *)data;
  gq  = exp_E_gquad(L, l, (vrna_exp_param_t *)pf);

  for (x = 0; x < L; x++) {
    pp[idx[i + x] - (i + x + 3 * L + l[0] + l[1] + l[2])]                       += gq;
    pp[idx[i + x] - (i + x + L + l[0])]                                         += gq;
    pp[idx[i + x + L + l[0]] - (i + x + 2 * L + l[0] + l[1])]                   += gq;
    pp[idx[i + x + 2 * L + l[0] + l[1]] - (i + x + 3 * L + l[0] + l[1] + l[2])] += gq;
  }
}


PRIVATE void
gquad_interact_ali(int  i,
                   int  L,
                   int  *l,
                   void *data,
                   void *index,
                   void *helper,
                   void *NA)
{
  int         x, *idx, bad;
  FLT_OR_DBL  gq, *pp;

  idx = (int *)index;
  pp  = (FLT_OR_DBL *)data;
  bad = 0;

  CHECK_GQUAD(L, l, bad = 1);

  gq = 0.;

  if (!bad) {
    gquad_pf_ali(i, L, l,
                 (void *)(&gq),
                 helper,
                 NULL,
                 NULL);
  }

  for (x = 0; x < L; x++) {
    pp[idx[i + x] - (i + x + 3 * L + l[0] + l[1] + l[2])]                       += gq;
    pp[idx[i + x] - (i + x + L + l[0])]                                         += gq;
    pp[idx[i + x + L + l[0]] - (i + x + 2 * L + l[0] + l[1])]                   += gq;
    pp[idx[i + x + 2 * L + l[0] + l[1]] - (i + x + 3 * L + l[0] + l[1] + l[2])] += gq;
  }
}

