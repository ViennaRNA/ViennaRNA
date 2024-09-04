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
gquad_mfe(unsigned int   i,
          unsigned int   L,
          unsigned int   *l,
          void  *data,
          void  *P,
          void  *NA,
          void  *NA2);


PRIVATE
void
gquad_mfe_pos(unsigned int   i,
              unsigned int   L,
              unsigned int   *l,
              void  *data,
              void  *P,
              void  *Lmfe,
              void  *lmfe);


PRIVATE void
gquad_mfe_ali_pos(unsigned int   i,
                  unsigned int   L,
                  unsigned int   *l,
                  void  *data,
                  void  *helper,
                  void  *Lmfe,
                  void  *lmfe);


PRIVATE
void
gquad_pos_exhaustive(unsigned int  i,
                     unsigned int  L,
                     unsigned int  *l,
                     void *data,
                     void *P,
                     void *Lex,
                     void *lex);


/**
 * Partition function callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_pf(unsigned int  i,
         unsigned int  L,
         unsigned int  *l,
         void *data,
         void *P,
         void *NA,
         void *NA2);


PRIVATE void
gquad_pf_ali(unsigned int  i,
             unsigned int  L,
             unsigned int  *l,
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
gquad_pf_pos(unsigned int  i,
             unsigned int  L,
             unsigned int  *l,
             void *data,
             void *pf,
             void *Lmax,
             void *lmax);


PRIVATE void
gquad_pf_pos_ali(unsigned int  i,
                 unsigned int  L,
                 unsigned int  *l,
                 void *data,
                 void *helper,
                 void *NA1,
                 void *NA2);


/**
 * MFE (alifold) callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_mfe_ali(unsigned int   i,
              unsigned int   L,
              unsigned int   *l,
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
gquad_mfe_ali_en(unsigned int  i,
                 unsigned int  L,
                 unsigned int  *l,
                 void *data,
                 void *helper,
                 void *NA,
                 void *NA2);


PRIVATE
void
gquad_interact(unsigned int  i,
               unsigned int  L,
               unsigned int  *l,
               void *data,
               void *pf,
               void *index,
               void *NA2);


PRIVATE
void
gquad_interact_ali(unsigned int  i,
                   unsigned int  L,
                   unsigned int  *l,
                   void *data,
                   void *index,
                   void *helper,
                   void *NA);


PRIVATE
void
gquad_count(unsigned int   i,
            unsigned int   L,
            unsigned int   *l,
            void  *data,
            void  *NA,
            void  *NA2,
            void  *NA3);


PRIVATE
void
gquad_count_layers(unsigned int  i,
                   unsigned int  L,
                   unsigned int  *l,
                   void *data,
                   void *NA,
                   void *NA2,
                   void *NA3);


/* other useful static functions */

PRIVATE int
E_gquad_consensus(unsigned int                 L,
                  unsigned int                 l[3],
                  unsigned int        position,
                  unsigned int        length,
                  unsigned int        n_seq,
                  const unsigned int  **a2s,
                  vrna_param_t        *P);


PRIVATE
int
E_gquad_ali_penalty(unsigned int           L,
                    unsigned int           l[3],
                    unsigned int  i,
                    unsigned int  length,
                    unsigned int  n_seq,
                    const short   **S,
                    vrna_param_t  *P);


PRIVATE
FLT_OR_DBL
exp_E_gquad_ali_penalty(unsigned int               L,
                        unsigned int               l[3],
                        unsigned int      i,
                        unsigned int      length,
                        unsigned int      n_seq,
                        const short       **S,
                        vrna_exp_param_t  *P);


PRIVATE void
count_gquad_layer_mismatches(unsigned int          L,
                             unsigned int          l[3],
                             unsigned int i,
                             unsigned int n,
                             unsigned int n_seq,
                             const short  **S,
                             unsigned int mm[2]);


PRIVATE void
gquad_mfe_ali_pos(unsigned int   i,
                  unsigned int   L,
                  unsigned int   *l,
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

/* 1. Energy evaluation */


/********************************
 * Here are the G-quadruplex energy
 * contribution functions
 *********************************/
PUBLIC int
vrna_gq_energy(unsigned int          L,
               unsigned int          l[3],
               vrna_param_t *P)
{
  int c = INF;

  if (P) {
    CHECK_GQUAD(L, l, return c);

    gquad_mfe(0, L, l,
              (void *)(&c),
              (void *)P,
              NULL,
              NULL);
  }

  return c;
}


PUBLIC FLT_OR_DBL
vrna_gq_exp_energy(unsigned int              L,
                   unsigned int              l[3],
                   vrna_exp_param_t *pf)
{
  FLT_OR_DBL q = 0.;

  if (pf) {
    CHECK_GQUAD(L, l, return q);

    gquad_pf(0, L, l,
             (void *)(&q),
             (void *)pf,
             NULL,
             NULL);
  }

  return q;
}


PUBLIC void
vrna_gq_consensus_energy(unsigned int          L,
                         unsigned int          l[3],
                         unsigned int position,
                         unsigned int length,
                         unsigned int n_seq,
                         const short  **S,
                         unsigned int **a2s,
                         vrna_param_t *P,
                         int          en[2])
{
  unsigned int  s;
  int           ee, ee2, u1, u2, u3;

  en[0] = en[1] = INF;

  if (P) {
    CHECK_GQUAD(L, l, return );

    /* compute actual quadruplex contribution for subalignment */
    ee = E_gquad_consensus(L, l, position, length, n_seq, (const unsigned int **)a2s, P);

    /* get penalty from incompatible sequences in alignment */
    ee2 = E_gquad_ali_penalty(L, l, position, length, n_seq, S, P);

    /* assign return values */
    if (ee2 != INF) {
      en[0] = ee;
      en[1] = ee2;
    }
  }
}


PUBLIC FLT_OR_DBL
vrna_gq_consensus_exp_energy(unsigned int                L,
                             unsigned int                l[3],
                             vrna_exp_param_t   *pf,
                             unsigned int       position,
                             unsigned int       length,
                             unsigned int       n_seq,
                             const short        **S,
                             const unsigned int **a2s)
{
  FLT_OR_DBL q = 0.;

  if ((pf) &&
      (S) &&
      (a2s)) {
    CHECK_GQUAD(L, l, return q);

    struct gquad_ali_helper gq_help = {
      .S      = S,
      .a2s    = a2s,
      .length = length,
      .n_seq  = n_seq,
      .pf     = pf
    };

    gquad_pf_ali(position, L, l,
                 (void *)(&q),
                 (void *)&gq_help,
                 NULL,
                 NULL);
  }

  return q;
}


PUBLIC int
vrna_gq_int_loop_mfe(vrna_fold_compound_t *fc,
                     unsigned int         i,
                     unsigned int         j)
{
  unsigned int  type, s, n_seq, **a2s, p, q, u1, u2, minq, maxq, l1;
  int           energy, ge, e_gq, dangles, c0;
  short         *S, *S1, si, sj, **SS, **S5, **S3;
  vrna_param_t  *P;
  vrna_md_t     *md;

  vrna_smx_csr(int) * c_gq;

  ge = INF;

  if ((fc) &&
      (i > 0) &&
      (i + VRNA_GQUAD_MIN_BOX_SIZE < j)) {
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
        return ge;
    }


    p = i + 1;
    if (S1[p] == 3) {
      if (p + VRNA_GQUAD_MIN_BOX_SIZE < j) {
        minq = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
        if (minq + 1 + MAXLOOP < j)
          minq = j - MAXLOOP - 1;

        maxq = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
        if (maxq + 3 > j)
          maxq = j - 3;

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
                c0 += P->internal_loop[j - q - 1];
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++) {
                  u1  = a2s[s][j - 1] - a2s[s][q];
                  c0  += P->internal_loop[u1];
                }

                break;
            }

            ge = MIN2(ge, c0);
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

      minq = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      if (minq + 1 + MAXLOOP - l1 < j)
        minq = j - MAXLOOP + l1 - 1;

      maxq = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      if (maxq >= j)
        maxq = j - 1;

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

          ge = MIN2(ge, c0);
        }
      }
    }

    q = j - 1;
    if (S1[q] == 3)
      for (p = (i + 4 + VRNA_GQUAD_MAX_BOX_SIZE - 1 < q) ? q - VRNA_GQUAD_MAX_BOX_SIZE + 1 : i + 4;
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

          ge = MIN2(ge, c0);
        }
      }
  }

  return ge;
}


vrna_array(int)
vrna_gq_int_loop_subopt(vrna_fold_compound_t * fc,
                        unsigned int i,
                        unsigned int j,
                        vrna_array(int) * ps,
                        vrna_array(int) * qs,
                        int threshold){
  unsigned int   type, p, q, l1, minq, maxq;
  int   energy, *ge, e_gq, dangles, c0;
  short *S, *S1, si, sj;

  ge  = NULL;
  *ps = NULL;
  *qs = NULL;

  if ((fc) &&
      (ps) &&
      (qs)) {
    vrna_param_t  *P  = fc->params;
    vrna_md_t     *md = &(P->model_details);

    vrna_smx_csr(int) * c_gq = fc->matrices->c_gq;

    S       = fc->sequence_encoding2;
    S1      = fc->sequence_encoding;
    type    = vrna_get_ptype_md(S[i], S[j], md);
    dangles = md->dangles;
    si      = S1[i + 1];
    sj      = S1[j - 1];
    energy  = 0;

    if (dangles == 2)
      energy += P->mismatchI[type][si][sj];

    if (type > 2)
      energy += P->TerminalAU;

    vrna_array_init(*ps);
    vrna_array_init(*qs);
    vrna_array_init(ge);

    p = i + 1;
    if (S[p] == 3) {
      if (p + VRNA_GQUAD_MIN_BOX_SIZE < j) {
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
              vrna_array_append(ge, energy + P->internal_loop[j - q - 1]);
              vrna_array_append(*ps, p);
              vrna_array_append(*qs, q);
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
            vrna_array_append(ge, energy + P->internal_loop[l1 + j - q - 1]);
            vrna_array_append(*ps, p);
            vrna_array_append(*qs, q);
          }
        }
      }
    }

    q = j - 1;
    if (S[q] == 3)
      for (p = i + 4;
           p + VRNA_GQUAD_MIN_BOX_SIZE < j;
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
            vrna_array_append(ge, energy + P->internal_loop[l1]);
            vrna_array_append(*ps, p);
            vrna_array_append(*qs, q);
          }
        }
      }
  }

  return ge;
}


PUBLIC int
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


PUBLIC int
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


PUBLIC FLT_OR_DBL
vrna_gq_int_loop_pf(vrna_fold_compound_t  *fc,
                    unsigned int          i,
                    unsigned int          j)
{
  short             *S, *S1, si, sj, **SS, **S5, **S3;
  unsigned int      type, s, n_seq, k, l, minl, maxl, u, u1, u2, **a2s;
  FLT_OR_DBL        q, qe, q_g, *scale;
  double            *expintern;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;

  vrna_smx_csr(FLT_OR_DBL) * q_gq;

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

  qe  = 1.;
  q   = 0.;

  if ((fc) &&
      (i > 0) &&
      (i + VRNA_GQUAD_MIN_BOX_SIZE < j)) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type = vrna_get_ptype_md(S[i], S[j], md);
        if (dangles == 2)
          qe *= (FLT_OR_DBL)pf_params->expmismatchI[type][si][sj];

        if (type > 2)
          qe *= (FLT_OR_DBL)pf_params->expTermAU;

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++) {
          type = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
          if (dangles == 2)
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
      if (k + VRNA_GQUAD_MIN_BOX_SIZE < j) {
        minl = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
        if (minl + 1 + MAXLOOP < j)
          minl = j - MAXLOOP - 1;

        maxl = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
        if (maxl + 3 > j)
          maxl = j - 3;

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
         k + VRNA_GQUAD_MIN_BOX_SIZE < j;
         k++) {
      u = k - i - 1;
      if (u > MAXLOOP)
        break;

      if (S1[k] != 3)
        continue;

      minl = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      if (minl + 1 + MAXLOOP - u < j)
        minl = j - MAXLOOP + u - 1;

      maxl = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      if (maxl >= j)
        maxl = j - 1;

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
      for (k = (i + 4 + VRNA_GQUAD_MAX_BOX_SIZE - 1 < l) ? l - VRNA_GQUAD_MAX_BOX_SIZE + 1 : i + 4;
           k + VRNA_GQUAD_MIN_BOX_SIZE - 1 < j;
           k++) {
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
  }

  return q;
}


/* 4. Parsing */

PUBLIC plist *
get_plist_gquad_from_db(const char  *structure,
                        float       pr)
{
  unsigned int  x, L, l[3], n, size, ge, ee, actual_size, gb;
  plist         *pl;

  actual_size = 0;
  ge          = 0;
  n           = 2;
  size        = strlen(structure);
  pl          = (plist *)vrna_alloc(n * size * sizeof(plist));

  while ((ee = vrna_gq_parse(structure + ge, &L, l)) > 0) {
    ge  += ee;

    if (4 * L + l[0] + l[1] + l[2] > ee) {
      gb = size + ge - L * 4 - l[0] - l[1] - l[2] + 1;
    } else {
      gb = ge - L * 4 - l[0] - l[1] - l[2] + 1;
    }

    /* add pseudo-base pair enclosing gquad */
    if (actual_size >= n * size - 5) {
      n   *= 2;
      pl  = (plist *)vrna_realloc(pl, n * size * sizeof(plist));
    }

    pl[actual_size].i       = gb;
    pl[actual_size].j       = ge;
    pl[actual_size].p       = pr;
    pl[actual_size++].type  = VRNA_PLIST_TYPE_GQUAD;

    for (x = 0; x < L; x++) {
      if (actual_size >= n * size - 5) {
        n   *= 2;
        pl  = (plist *)vrna_realloc(pl, n * size * sizeof(plist));
      }

      pl[actual_size].i       = (gb + x - 1) % (size) + 1;
      pl[actual_size].j       = (ge + x - L + 1 - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;

      pl[actual_size].i       = (gb + x - 1) % (size) + 1;
      pl[actual_size].j       = (gb + x + l[0] + L - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;

      pl[actual_size].i       = (gb + x + l[0] + L - 1) % (size) + 1;
      pl[actual_size].j       = (ge + x - 2 * L - l[2] + 1 - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;

      pl[actual_size].i       = (ge + x - 2 * L - l[2] + 1 - 1) % (size) + 1;
      pl[actual_size].j       = (ge + x - L + 1 - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;
    }
  }

  pl[actual_size].i   = pl[actual_size].j = 0;
  pl[actual_size++].p = 0;
  pl                  = (plist *)vrna_realloc(pl, actual_size * sizeof(plist));
  return pl;
}


PUBLIC void
get_gquad_pattern_mfe(short         *S,
                      int           i,
                      int           j,
                      vrna_param_t  *P,
                      int           *L,
                      int           l[3])
{
  unsigned int *gg = get_g_islands_sub(S, (unsigned int)i, (unsigned int)j);
  int c   = INF;

  process_gquad_enumeration(gg, i, j,
                            &gquad_mfe_pos,
                            (void *)(&c),
                            (void *)P,
                            (void *)L,
                            (void *)l);

  gg += i - 1;
  free(gg);
}


PUBLIC void
get_gquad_pattern_mfe_ali(short         **S,
                          unsigned int  **a2s,
                          short         *S_cons,
                          int           n_seq,
                          int           i,
                          int           j,
                          vrna_param_t  *P,
                          int           *L,
                          int           l[3])
{
  int           mfe;
  unsigned int *gg, LL, ll[3];


  gg  = get_g_islands_sub(S_cons, (unsigned int)i, (unsigned int)j);
  mfe = INF;

  struct gquad_ali_helper gq_help = {
    .S      = (const short **)S,
    .a2s    = (const unsigned int **)a2s,
    .n_seq  = (unsigned int)n_seq,
    .P      = P
  };

  process_gquad_enumeration(gg, i, j,
                            &gquad_mfe_ali_pos,
                            (void *)(&mfe),
                            (void *)&gq_help,
                            (void *)&LL,
                            (void *)ll);

  gg += (unsigned int)i - 1;

  *L = LL;
  l[0] = ll[0];
  l[1] = ll[1];
  l[2] = ll[2];
  free(gg);
}


PUBLIC void
get_gquad_pattern_exhaustive(short        *S,
                             int          i,
                             int          j,
                             vrna_param_t *P,
                             int          *L,
                             int          *l,
                             int          threshold)
{
  unsigned int *gg, LL, ll[3];

  gg = get_g_islands_sub(S, (unsigned int)i, (unsigned int)j);

  process_gquad_enumeration(gg, i, j,
                            &gquad_pos_exhaustive,
                            (void *)(&threshold),
                            (void *)P,
                            (void *)&LL,
                            (void *)ll);

  gg += (unsigned int)i - 1;

  *L = LL;
  l[0] = ll[0];
  l[1] = ll[1];
  l[2] = ll[2];
  free(gg);
}


PUBLIC void
get_gquad_pattern_pf(short            *S,
                     int              i,
                     int              j,
                     vrna_exp_param_t *pf,
                     int              *L,
                     int              l[3])
{
  unsigned int *gg, LL, ll[3];

  gg = get_g_islands_sub(S, (unsigned int)i, (unsigned int)j);

  FLT_OR_DBL  q   = 0.;

  process_gquad_enumeration(gg, i, j,
                            &gquad_pf_pos,
                            (void *)(&q),
                            (void *)pf,
                            (void *)&LL,
                            (void *)ll);

  gg += (unsigned int)i - 1;

  *L = LL;
  l[0] = ll[0];
  l[1] = ll[1];
  l[2] = ll[2];
  free(gg);
}


PUBLIC void
vrna_get_gquad_pattern_pf(vrna_fold_compound_t  *fc,
                          unsigned int          i,
                          unsigned int          j,
                          unsigned int          *L,
                          unsigned int          l[3])
{
  short             *S_enc, *S_tmp;
  unsigned int      n, n2, *gg;
  FLT_OR_DBL        q   = 0.;
  vrna_exp_param_t  *pf_params;
  void              *data;
  void              ( *process_f )(unsigned int,
                                   unsigned int,
                                   unsigned int *,
                                   void *,
                                   void *,
                                   void *,
                                   void *);
  struct gquad_ali_helper tmp = {
    0
  };

  n         = fc->length;
  n2        = 0;
  pf_params = fc->exp_params;
  S_tmp     = NULL;

  *L   = 0;
  l[0] = l[1] = l[2] = 0;

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      S_enc = fc->S_cons;
      struct gquad_ali_helper gq_help = {
        .S      = (const short **)fc->S,
        .a2s    = (const unsigned int **)fc->a2s,
        .length = fc->length,
        .n_seq  = fc->n_seq,
        .pf     = pf_params,
        .L      = 0,
        .l      = &(l[0])
      };
      tmp       = gq_help;
      data      = (void *)&tmp;
      process_f = &gquad_pf_pos_ali;
      break;

    default:
      S_enc     = fc->sequence_encoding2;
      data      = (void *)pf_params;
      process_f = &gquad_pf_pos;
      break;
  }

  if ((pf_params->model_details.circ) &&
      (j < i)) {
    j += n;
    /* G-Quadruplex wraps around the n,1 junction */
    n2 = MIN2(n, VRNA_GQUAD_MAX_BOX_SIZE) - 1;
    S_tmp = (short *)vrna_alloc(sizeof(short) * (n + n2 + 1));
    memcpy(S_tmp, S_enc, sizeof(short) * (n + 1));
    memcpy(S_tmp + (n + 1), S_enc + 1, sizeof(short) * n2);
    S_tmp[0]  = n + n2;
    S_enc     = S_tmp;
    n         += n2;
  }

  gg = get_g_islands_sub(S_enc, i, j);

  process_gquad_enumeration(gg, i, j,
                            process_f,
                            (void *)(&q),
                            data,
                            (void *)L,
                            (void *)&(l[0]));

  if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
    *L = tmp.L;

  gg += i - 1;
  free(gg);
  free(S_tmp);
}


PUBLIC plist *
get_plist_gquad_from_pr(short                     *S,
                        int                       gi,
                        int                       gj,
                        vrna_smx_csr(FLT_OR_DBL)  *q_gq,
                        FLT_OR_DBL                *probs,
                        FLT_OR_DBL                *scale,
                        vrna_exp_param_t          *pf)
{
  int L, l[3];

  return get_plist_gquad_from_pr_max(S, gi, gj, q_gq, probs, scale, &L, l, pf);
}


PUBLIC vrna_ep_t *
vrna_plist_gquad_from_pr(vrna_fold_compound_t *fc,
                         int                  gi,
                         int                  gj)
{
  unsigned int L, l[3];

  return vrna_plist_gquad_from_pr_max(fc, (unsigned int)gi, (unsigned int)gj, &L, l);
}


PUBLIC plist *
get_plist_gquad_from_pr_max(short                     *S,
                            int                       gi,
                            int                       gj,
                            vrna_smx_csr(FLT_OR_DBL)  *q_gq,
                            FLT_OR_DBL                *probs,
                            FLT_OR_DBL                *scale,
                            int                       *Lmax,
                            int                       lmax[3],
                            vrna_exp_param_t          *pf)
{
  unsigned int  n, size, *gg, counter, i, j, L, l[3];
  int         *my_index;
  FLT_OR_DBL  pp, *tempprobs;
  plist       *pl;

  n         = S[0];
  size      = (n * (n + 1)) / 2 + 2;
  tempprobs = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  pl        = (plist *)vrna_alloc((S[0] * S[0]) * sizeof(plist));
  gg        = get_g_islands_sub(S, (unsigned int)gi, (unsigned int)gj);
  counter   = 0;
  my_index  = vrna_idx_row_wise(n);

  process_gquad_enumeration(gg, gi, gj,
                            &gquad_interact,
                            (void *)tempprobs,
                            (void *)pf,
                            (void *)my_index,
                            NULL);

  pp = 0.;
  process_gquad_enumeration(gg, gi, gj,
                            &gquad_pf_pos,
                            (void *)(&pp),
                            (void *)pf,
                            (void *)&L,
                            (void *)l);

#ifndef VRNA_DISABLE_C11_FEATURES
  pp = probs[my_index[gi] - gj] * scale[gj - gi + 1] /
       vrna_smx_csr_get(q_gq, gi, gj, 0.);
#else
  pp = probs[my_index[gi] - gj] * scale[gj - gi + 1] /
       vrna_smx_csr_FLT_OR_DBL_get(q_gq, gi, gj, 0.);
#endif
  for (i = gi; i < gj; i++) {
    for (j = i; j <= gj; j++) {
      if (tempprobs[my_index[i] - j] > 0.) {
        pl[counter].i = i;
        pl[counter].j = j;
        pl[counter].p = pp *
                        tempprobs[my_index[i] - j];
        pl[counter++].type = VRNA_PLIST_TYPE_TRIPLE;
      }
    }
  }
  pl[counter].i   = pl[counter].j = 0;
  pl[counter++].p = 0.;
  /* shrink memory to actual size needed */
  pl = (plist *)vrna_realloc(pl, counter * sizeof(plist));

  gg += gi - 1;
  free(gg);
  free(my_index);
  free(tempprobs);

  *Lmax = L;
  lmax[0] = l[0];
  lmax[1] = l[1];
  lmax[2] = l[2];
  return pl;
}


PUBLIC plist *
vrna_plist_gquad_from_pr_max(vrna_fold_compound_t *fc,
                             unsigned int                  gi,
                             unsigned int                  gj,
                             unsigned int                  *Lmax,
                             unsigned int                  lmax[3])
{
  short             *S_enc, *S_tmp;
  unsigned int      n, n2, **a2s, real_j, *gg;
  int               size, counter, i, j, *my_index;
  FLT_OR_DBL        pp, *tempprobs, *probs, *scale;
  plist             *pl;
  vrna_exp_param_t  *pf;
  vrna_smx_csr(FLT_OR_DBL) * q_gq;
  struct gquad_ali_helper tmp = {
    0
  };

  n         = fc->length;
  n2        = 0;
  pf        = fc->exp_params;
  q_gq      = fc->exp_matrices->q_gq;
  probs     = fc->exp_matrices->probs;
  scale     = fc->exp_matrices->scale;
  counter   = 0;
  real_j    = gj;
  S_tmp     = NULL;

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      S_enc = fc->S_cons;
      struct gquad_ali_helper gq_help = {
        .S      = (const short **)fc->S,
        .a2s    = (const unsigned int **)fc->a2s,
        .length = fc->length,
        .n_seq  = fc->n_seq,
        .pf     = pf,
        .L      = *Lmax,
        .l      = &(lmax[0])
      };
      tmp       = gq_help;
      break;

    default:
      S_enc     = fc->sequence_encoding2;
      break;
  }

  if ((pf->model_details.circ) &&
      (gi > gj)) {
    n2 = MIN2(n, VRNA_GQUAD_MAX_BOX_SIZE) - 1;

    gj += n;

    S_tmp = (short *)vrna_alloc(sizeof(short) * (n + n2 + 1));
    memcpy(S_tmp, S_enc, sizeof(short) * (n + 1));
    memcpy(S_tmp + (n + 1), S_enc + 1, sizeof(short) * n2);
    S_tmp[0]  = n + n2;
    S_enc     = S_tmp;
    n         += n2;
  }

  size      = (n * (n + 1)) / 2 + 2;
  tempprobs = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  pl        = (plist *)vrna_alloc((n * n) * sizeof(plist));
  my_index  = vrna_idx_row_wise(n);

  gg        = get_g_islands_sub(S_enc, gi, gj);

  pp = 0.;

  if (fc->type == VRNA_FC_TYPE_SINGLE) {
    process_gquad_enumeration(gg, gi, gj,
                              &gquad_interact,
                              (void *)tempprobs,
                              (void *)pf,
                              (void *)my_index,
                              NULL);

    process_gquad_enumeration(gg, gi, gj,
                              &gquad_pf_pos,
                              (void *)(&pp),
                              (void *)pf,
                              (void *)Lmax,
                              (void *)lmax);
  } else {
    process_gquad_enumeration(gg, gi, gj,
                              &gquad_interact_ali,
                              (void *)tempprobs,
                              (void *)my_index,
                              (void *)&(tmp),
                              NULL);

    process_gquad_enumeration(gg, gi, gj,
                              &gquad_pf_pos_ali,
                              (void *)(&pp),
                              (void *)&tmp,
                              NULL,
                              NULL);
    *Lmax = tmp.L;
  }

  if (gj > real_j) {
    vrna_smx_csr(FLT_OR_DBL) *p_gq = fc->exp_matrices->p_gq;
    pp = scale[real_j + n - n2 - gi + 1] *
#ifndef VRNA_DISABLE_C11_FEATURES
         vrna_smx_csr_get(p_gq, gi, real_j, 0.);
#else
         vrna_smx_csr_FLT_OR_DBL_get(p_gq, gi, real_j, 0.);
#endif
  } else {
    pp = probs[my_index[gi] - gj] *
       scale[gj - gi + 1];
  }
  
  pp  /= 
#ifndef VRNA_DISABLE_C11_FEATURES
       vrna_smx_csr_get(q_gq, gi, real_j, 0.);
#else
       vrna_smx_csr_FLT_OR_DBL_get(q_gq, gi, real_j, 0.);
#endif

  n -= n2;

  for (i = gi; i < gj; i++) {
    for (j = i; j <= gj; j++) {
      if (tempprobs[my_index[i] - j] > 0.) {
        pl[counter].i = (i - 1) % (n) + 1;
        pl[counter].j = (j - 1) % (n) + 1;
        pl[counter].p = pp *
                        tempprobs[my_index[i] - j];
        pl[counter++].type = VRNA_PLIST_TYPE_TRIPLE;
      }
    }
  }
  pl[counter].i   = pl[counter].j = 0;
  pl[counter++].p = 0.;
  /* shrink memory to actual size needed */
  pl = (plist *)vrna_realloc(pl, counter * sizeof(plist));

  gg += gi - 1;
  free(gg);
  free(my_index);
  free(S_tmp);
  free(tempprobs);
  return pl;
}


PUBLIC int
get_gquad_count(short *S,
                int   i,
                int   j)
{
  unsigned int p, q, *gg;
  int           counter;

  gg = get_g_islands_sub(S, (unsigned int)i, (unsigned int)j);
  counter = 0;

  FOR_EACH_GQUAD(p, q, (unsigned int)i, (unsigned int)j)
  process_gquad_enumeration(gg, p, q,
                            &gquad_count,
                            (void *)(&counter),
                            NULL,
                            NULL,
                            NULL);

  gg += (unsigned int)i - 1;
  free(gg);
  return counter;
}


PUBLIC int
get_gquad_layer_count(short *S,
                      int   i,
                      int   j)
{
  unsigned int  p, q, *gg;
  int           counter;

  gg = get_g_islands_sub(S, (unsigned int)i, (unsigned int)j);
  counter = 0;

  FOR_EACH_GQUAD(p, q, (unsigned int)i, (unsigned int)j)
  process_gquad_enumeration(gg, p, q,
                            &gquad_count_layers,
                            (void *)(&counter),
                            NULL,
                            NULL,
                            NULL);

  gg += (unsigned int)i - 1;
  free(gg);
  return counter;
}


PUBLIC unsigned int
vrna_gq_parse(const char    *db_string,
              unsigned int  *L,
              unsigned int  l[3])
{
  unsigned int i, n, end, start, stop, pre, tetrad, stacks, LL[4];

  end = 0;

  if ((db_string == NULL) ||
      (L == NULL))
    return 0;

  /*
   *  scan along the dot-bracket string to identify the first character that
   *  indicates a part of a G-Quadruplex
   */

  n       = strlen(db_string);
  LL[0]   = LL[1] = LL[2] = LL[3] = 0;
  stop    = 0;
  *L      = 0;
  l[0]    = l[1] = l[2] = 0;
  tetrad  = 0;

  for (i = 0;
       (db_string[i]) &&
       (db_string[i] != VRNA_GQUAD_DB_SYMBOL) &&
       (db_string[i] != VRNA_GQUAD_DB_SYMBOL_END);
       i++);

  if (db_string[i]) {
    pre = i;

    /* we've encountered some piece that suspicially looks like a G-Quadruplex */
    if (db_string[i] == VRNA_GQUAD_DB_SYMBOL_END) {
      stop  = 1;
      LL[0] = 1;
    }

    while (stop == 0) {
      start = i;
      while (db_string[++i] == VRNA_GQUAD_DB_SYMBOL)
        /* stop consuming gquad symbols if we already know the number of stacks */
        if ((tetrad > 1) &&
            (i - start == LL[tetrad - 1]))
          break;

      end         = i;
      LL[tetrad]  = end - start;

      if (db_string[i] == VRNA_GQUAD_DB_SYMBOL_END)
        LL[tetrad]++;

      if ((tetrad > 1) &&
          (LL[tetrad] != LL[tetrad - 1])) {
        vrna_log_debug("unequal stack sizes (%u vs. %u) in G-quadruplex",
                       LL[tetrad - 1],
                       LL[tetrad]);
        return 0;
      }

      /* check whether we stopped due to linker or end symbol */
      if (db_string[i] == VRNA_GQUAD_DB_SYMBOL_END) {
        stop = i + 1;
        break;
      }

      if (tetrad == 3)
        break; /* no more linkers */

      /*
       * next character must be linker!
       * count number of linker nucleotides
       */
      while (db_string[i++] == '.');
      l[tetrad] = i - end - 1;

      if (l[tetrad] == 0) {
        vrna_log_debug("zero-length linker in G-Quadruplex");
        return 0;
      } else {
        i--;
      }

      if (db_string[i] == '\0') {
        return 0;
      }

      tetrad++;
    }

    if (stop > 0) {
      if (tetrad < 3) {
        /* move currently parsed lengths to correct positions */
        unsigned int s;
        for (s = 0; s <= tetrad; s++)
          LL[3 - s] = LL[tetrad - s];

        for (; s <= 3; s++)
          LL[3 - s] = 0;

        for (s = 0; s < tetrad; s++)
          l[2 - s] = l[tetrad - s - 1];

        l[2 - s++] = pre;

        for (; s < 3; s++)
          l[2 - s] = 0;
      }

      /* only continue parsing if we did not find a complete G-Quadruplex yet */
      if (LL[0] != LL[3]) {
        i = n - 1;
        /* check whether the end of the structure is a continuation of a g-run or a linker */
        if (db_string[i] == '.') {
          while (db_string[--i] == '.');
          l[2 - tetrad] += n - i - 1;
          tetrad++;
        } else if (pre > 0) {
          tetrad++;
        }

        while (tetrad < 4) {
          start = i;

          while (db_string[--i] == VRNA_GQUAD_DB_SYMBOL) {
            /* stop consuming gquad symbols if we already know the number of stacks */
            if ((tetrad == 3) &&
                (start - i + LL[3 - tetrad] >= LL[3]))
              break;
          }

          end             = i;
          LL[3 - tetrad]  += start - end;

          if ((tetrad > 0) &&
              (LL[3 - tetrad] != LL[4 - tetrad])) {
            vrna_log_debug("unequal stack sizes (%u vs. %u) in G-quadruplex",
                           LL[3 - tetrad],
                           LL[4 - tetrad]);
            return 0;
          }

          if (tetrad == 3)
            break; /* no more linkers */

          /*
           * next character must be linker!
           * count number of linker nucleotides
           */
          while (db_string[i--] == '.');
          l[2 - tetrad] = end - i - 1;

          if (l[2 - tetrad] == 0) {
            vrna_log_debug("zero-length linker in G-Quadruplex");
            return 0;
          } else {
            i++;
          }

          if (i <= stop + 1) {
            return 0;
          }

          tetrad++;
        }
      }

      end = stop;
    }
  }

  /* last sanity check */
  if ((LL[0] == LL[1]) &&
      (LL[1] == LL[2]) &&
      (LL[2] == LL[3])) {
    *L = LL[0];
  } else {
    vrna_log_debug("Unequal G-quadruplex stack tetrad lengths (%u, %u, %u, %d)",
                   LL[0],
                   LL[1],
                   LL[2],
                   LL[3]);
    return 0;
  }

  return end;
}


PUBLIC void
vrna_db_insert_gq(char          *db,
                  unsigned int  i,
                  unsigned int  L,
                  unsigned int l[3],
                  unsigned int  n)
{
  if (db) {
    if (n == 0)
      n = strlen(db);

    if (4 * L + l[0] + l[1] + l[2] <= n) {
      unsigned int j, ll;
      for (ll = 0; ll < L; ll++) {
        j     = (i + ll - 1) % (n);
        db[j] = VRNA_GQUAD_DB_SYMBOL;
      }

      i += L + l[0];
      for (ll = 0; ll < L; ll++) {
        j     = (i + ll - 1) % (n);
        db[j] = VRNA_GQUAD_DB_SYMBOL;
      }

      i += L + l[1];
      for (ll = 0; ll < L; ll++) {
        j     = (i + ll - 1) % (n);
        db[j] = VRNA_GQUAD_DB_SYMBOL;
      }

      i += L + l[2];
      for (ll = 0; ll < L - 1; ll++) {
        j     = (i + ll - 1) % (n);
        db[j] = VRNA_GQUAD_DB_SYMBOL;
      }
      j     = (i + ll - 1) % (n);
      db[j] = VRNA_GQUAD_DB_SYMBOL_END;
    }
  }
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
aln_linker_positions(unsigned int          L,
                     unsigned int          l[3],
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
gq_layer_pos(unsigned int          L,
             unsigned int          l[3],
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
E_gquad_consensus(unsigned int                 L,
                  unsigned int                 l[3],
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
exp_E_gquad_consensus(unsigned int                 L,
                      unsigned int                 l[3],
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
E_gquad_ali_penalty(unsigned int           L,
                    unsigned int           l[3],
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
exp_E_gquad_ali_penalty(unsigned int               L,
                        unsigned int               l[3],
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
count_gquad_layer_mismatches(unsigned int          L,
                             unsigned int          l[3],
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
gquad_mfe(unsigned int   i,
          unsigned int   L,
          unsigned int   *l,
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
gquad_mfe_pos(unsigned int   i,
              unsigned int   L,
              unsigned int   *l,
              void  *data,
              void  *P,
              void  *Lmfe,
              void  *lmfe)
{
  int cc = ((vrna_param_t *)P)->gquad[L][l[0] + l[1] + l[2]];

  if (cc < *((int *)data)) {
    *((int *)data)        = cc;
    *((unsigned int *)Lmfe)        = L;
    *((unsigned int *)lmfe)        = l[0];
    *(((unsigned int *)lmfe) + 1)  = l[1];
    *(((unsigned int *)lmfe) + 2)  = l[2];
  }
}


PRIVATE void
gquad_mfe_ali_pos(unsigned int   i,
                  unsigned int   L,
                  unsigned int   *l,
                  void  *data,
                  void  *helper,
                  void  *Lmfe,
                  void  *lmfe)
{
  int cc = INF;

  gquad_mfe_ali(i, L, l, (void *)&cc, helper, NULL, NULL);

  if (cc < *((int *)data)) {
    *((int *)data)        = cc;
    *((unsigned int *)Lmfe)        = L;
    *((unsigned int *)lmfe)        = l[0];
    *(((unsigned int *)lmfe) + 1)  = l[1];
    *(((unsigned int *)lmfe) + 2)  = l[2];
  }
}


PRIVATE
void
gquad_pos_exhaustive(unsigned int  i,
                     unsigned int  L,
                     unsigned int  *l,
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
    for (cnt = 0; ((unsigned int *)Lex)[cnt] != 0; cnt++);

    *((unsigned int *)Lex + cnt)             = L;
    *((unsigned int *)Lex + cnt + 1)         = 0;
    *(((unsigned int *)lex) + (3 * cnt) + 0) = l[0];
    *(((unsigned int *)lex) + (3 * cnt) + 1) = l[1];
    *(((unsigned int *)lex) + (3 * cnt) + 2) = l[2];
  }
}


PRIVATE
void
gquad_count(unsigned int   i,
            unsigned int   L,
            unsigned int   *l,
            void  *data,
            void  *NA,
            void  *NA2,
            void  *NA3)
{
  *((int *)data) += 1;
}


PRIVATE
void
gquad_count_layers(unsigned int  i,
                   unsigned int  L,
                   unsigned int  *l,
                   void *data,
                   void *NA,
                   void *NA2,
                   void *NA3)
{
  *((unsigned int *)data) += L;
}


PRIVATE void
gquad_pf(unsigned int  i,
         unsigned int  L,
         unsigned int  *l,
         void *data,
         void *pf,
         void *NA,
         void *NA2)
{
  *((FLT_OR_DBL *)data) += ((vrna_exp_param_t *)pf)->expgquad[L][l[0] + l[1] + l[2]];
}


PRIVATE void
gquad_pf_ali(unsigned int  i,
             unsigned int  L,
             unsigned int  *l,
             void *data,
             void *helper,
             void *NA,
             void *NA2)
{
  const short             **S;
  const unsigned int      **a2s;
  unsigned int            n, s, n_seq;
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
                             exp_E_gquad_consensus(L, l, i, n, n_seq, a2s, pf);
}


PRIVATE void
gquad_pf_pos(unsigned int  i,
             unsigned int  L,
             unsigned int  *l,
             void *data,
             void *pf,
             void *Lmax,
             void *lmax)
{
  FLT_OR_DBL gq = 0.;

  gquad_pf(i, L, l, (void *)&gq, pf, NULL, NULL);

  if (gq > *((FLT_OR_DBL *)data)) {
    *((FLT_OR_DBL *)data) = gq;
    *((unsigned int *)Lmax)        = L;
    *((unsigned int *)lmax)        = l[0];
    *(((unsigned int *)lmax) + 1)  = l[1];
    *(((unsigned int *)lmax) + 2)  = l[2];
  }
}


PRIVATE void
gquad_pf_pos_ali(unsigned int  i,
                 unsigned int  L,
                 unsigned int  *l,
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
gquad_mfe_ali(unsigned int   i,
              unsigned int   L,
              unsigned int   *l,
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
gquad_mfe_ali_en(unsigned int  i,
                 unsigned int  L,
                 unsigned int  *l,
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
gquad_interact(unsigned int  i,
               unsigned int  L,
               unsigned int  *l,
               void *data,
               void *pf,
               void *index,
               void *NA2)
{
  int         x, *idx, LL, ll[3];
  FLT_OR_DBL  gq, *pp;

  idx = (int *)index;
  pp  = (FLT_OR_DBL *)data;
  LL  = (int)L;
  ll[0] = (int)l[0];
  ll[1] = (int)l[1];
  ll[2] = (int)l[2];

  gq  = exp_E_gquad(LL, ll, (vrna_exp_param_t *)pf);

  for (x = 0; x < L; x++) {
    pp[idx[i + x] - (i + x + 3 * L + l[0] + l[1] + l[2])]                       += gq;
    pp[idx[i + x] - (i + x + L + l[0])]                                         += gq;
    pp[idx[i + x + L + l[0]] - (i + x + 2 * L + l[0] + l[1])]                   += gq;
    pp[idx[i + x + 2 * L + l[0] + l[1]] - (i + x + 3 * L + l[0] + l[1] + l[2])] += gq;
  }
}


PRIVATE void
gquad_interact_ali(unsigned int  i,
                   unsigned int  L,
                   unsigned int  *l,
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


/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC int
parse_gquad(const char  *struc,
            int         *L,
            int         l[3])
{
  unsigned int ret, LL, ll[3];

  ret = vrna_gq_parse(struc, &LL, ll);

  if (ret) {
    *L    = (int)LL;
    l[0]  = (int)ll[0];
    l[1]  = (int)ll[1];
    l[2]  = (int)ll[2];
  }

  return ret;
}


PUBLIC int
E_gquad(int           L,
        int           l[3],
        vrna_param_t  *P)
{
  unsigned int LL, ll[3];

  LL = (unsigned int)L;
  ll[0] = (unsigned int)l[0];
  ll[1] = (unsigned int)l[1];
  ll[2] = (unsigned int)l[2];

  return vrna_gq_energy(LL, ll, P);
}


PUBLIC FLT_OR_DBL
exp_E_gquad(int               L,
            int               l[3],
            vrna_exp_param_t  *pf)
{
  unsigned int LL, ll[3];

  LL = (unsigned int)L;
  ll[0] = (unsigned int)l[0];
  ll[1] = (unsigned int)l[1];
  ll[2] = (unsigned int)l[2];

  return vrna_gq_exp_energy(LL, ll, pf);
}


PUBLIC void
E_gquad_ali_en(int          i,
               int          L,
               int          l[3],
               const short  **S,
               unsigned int **a2s,
               unsigned int n_seq,
               vrna_param_t *P,
               int          en[2])
{
  unsigned int LL, ll[3];

  LL = (unsigned int)L;
  ll[0] = (unsigned int)l[0];
  ll[1] = (unsigned int)l[1];
  ll[2] = (unsigned int)l[2];

  vrna_gq_consensus_energy(LL, ll, (unsigned int)i, 0, n_seq, S, a2s, P, en);
}


PUBLIC FLT_OR_DBL
exp_E_gquad_ali(int               i,
                int               L,
                int               l[3],
                short             **S,
                unsigned int      **a2s,
                int               n_seq,
                vrna_exp_param_t  *pf)
{
  unsigned int LL, ll[3];

  LL = (unsigned int)L;
  ll[0] = (unsigned int)l[0];
  ll[1] = (unsigned int)l[1];
  ll[2] = (unsigned int)l[2];

  return vrna_gq_consensus_exp_energy(LL,
                                      ll,
                                      pf,
                                      (unsigned int)i,
                                      0,
                                      n_seq,
                                      (const short **)S,
                                      (const unsigned int **)a2s);
}


#endif
