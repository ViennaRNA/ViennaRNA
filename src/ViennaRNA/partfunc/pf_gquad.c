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

#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/params/constants.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/partfunc/gquad.h"

#include "ViennaRNA/intern/gquad_helpers.h"


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


PRIVATE
FLT_OR_DBL
exp_E_gquad_ali_penalty(unsigned int               L,
                        unsigned int               l[3],
                        unsigned int      i,
                        unsigned int      length,
                        unsigned int      n_seq,
                        const short       **S,
                        vrna_exp_param_t  *P);


PRIVATE FLT_OR_DBL
exp_E_gquad_consensus(unsigned int                 L,
                      unsigned int                 l[3],
                      unsigned int        position,
                      unsigned int        length,
                      unsigned int        n_seq,
                      const unsigned int  **a2s,
                      vrna_exp_param_t    *pf);

PRIVATE void
count_gquad_layer_mismatches(unsigned int          L,
                             unsigned int          l[3],
                             unsigned int i,
                             unsigned int n,
                             unsigned int n_seq,
                             const short  **S,
                             unsigned int mm[2]);


PRIVATE INLINE void
aln_linker_positions(unsigned int          L,
                     unsigned int          l[3],
                     unsigned int position,
                     unsigned int length,
                     unsigned int starts[3],
                     unsigned int ends[3]);

PRIVATE void
gq_layer_pos(unsigned int          L,
             unsigned int          l[3],
             unsigned int layer,
             unsigned int i,
             unsigned int n,
             unsigned int layer_pos[4]);


PRIVATE INLINE int
aln_linker_length(unsigned int        start,
                  unsigned int        end,
                  unsigned int        n,
                  const unsigned int  *a2ss);

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


/*
 #########################################
 # BEGIN OF PUBLIC FUNCTION DEFINITIONS  #
 #      (all available in RNAlib)        #
 #########################################
 */

PUBLIC
vrna_smx_csr(FLT_OR_DBL) *
vrna_gq_pos_pf(vrna_fold_compound_t * fc){
  vrna_smx_csr(FLT_OR_DBL) * q_gq = NULL;

  if (fc) {
    unsigned int      i, j, n, n2, *gg;
    short             *S_tmp, *S_enc;
    FLT_OR_DBL        q, *scale;
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
    q_gq      = vrna_smx_csr_FLT_OR_DBL_init(n + 1);
    scale     = fc->exp_matrices->scale;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        S_enc     = fc->sequence_encoding2;
        data      = (void *)pf_params;
        process_f = &gquad_pf;
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        S_enc = fc->S_cons;
        struct gquad_ali_helper gq_help = {
          .S      = (const short **)fc->S,
          .a2s    = (const unsigned int **)fc->a2s,
          .length = fc->length,
          .n_seq  = fc->n_seq,
          .pf     = pf_params
        };
        tmp       = gq_help;
        data      = (void *)&tmp;
        process_f = &gquad_pf_ali;
        break;

      default:
        return NULL;
    }

    if (pf_params->model_details.circ) {
      n2 = MIN2(n, VRNA_GQUAD_MAX_BOX_SIZE) - 1;

      S_tmp = (short *)vrna_alloc(sizeof(short) * (n + n2 + 1));
      memcpy(S_tmp, S_enc, sizeof(short) * (n + 1));
      memcpy(S_tmp + (n + 1), S_enc + 1, sizeof(short) * n2);
      S_tmp[0]  = n + n2;
      S_enc     = S_tmp;
      n         += n2;
    }

    gg = get_g_islands(S_enc);

    FOR_EACH_GQUAD_INC(i, j, 1, n) {
      q = 0.;

      if (i > n - n2)
        break;

      process_gquad_enumeration(gg, i, j,
                                process_f,
                                (void *)(&q),
                                data,
                                NULL,
                                NULL);
      if ((q != 0.) &&
          (j - i + 1 <= n - n2)) {
#ifndef VRNA_DISABLE_C11_FEATURES
        vrna_smx_csr_insert(q_gq, i, (j - 1) % (n - n2) + 1, q * scale[j - i + 1]);

#else
        vrna_smx_csr_FLT_OR_DBL_insert(q_gq, i, (j - 1) % (n - n2) + 1, q * scale[j - i + 1]);
#endif
      }
    }

    free(S_tmp);
    free(gg);
  }

  return q_gq;
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
  unsigned int      n, n2, real_j, *gg;
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


/*
 #########################################
 # BEGIN OF PRIVATE FUNCTION DEFINITIONS #
 #          (internal use only)          #
 #########################################
 */
PRIVATE void
gquad_pf(unsigned int  i VRNA_UNUSED,
         unsigned int  L,
         unsigned int  *l,
         void *data,
         void *pf,
         void *NA VRNA_UNUSED,
         void *NA2 VRNA_UNUSED)
{
  *((FLT_OR_DBL *)data) += ((vrna_exp_param_t *)pf)->expgquad[L][l[0] + l[1] + l[2]];
}


PRIVATE void
gquad_pf_ali(unsigned int  i,
             unsigned int  L,
             unsigned int  *l,
             void *data,
             void *helper,
             void *NA VRNA_UNUSED,
             void *NA2 VRNA_UNUSED)
{
  const short             **S;
  const unsigned int      **a2s;
  unsigned int            n, n_seq;
  FLT_OR_DBL              penalty;
  vrna_exp_param_t        *pf;
  struct gquad_ali_helper *gq_help;

  gq_help = (struct gquad_ali_helper *)helper;
  S       = gq_help->S;
  a2s     = gq_help->a2s;
  n       = gq_help->length;
  n_seq   = gq_help->n_seq;
  pf      = gq_help->pf;
  penalty = exp_E_gquad_ali_penalty(L, l, i, n, n_seq, S, pf);

  if (penalty != 0.)
    *((FLT_OR_DBL *)data) += penalty *
                             exp_E_gquad_consensus(L, l, i, n, n_seq, a2s, pf);
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


PRIVATE void
count_gquad_layer_mismatches(unsigned int          L,
                             unsigned int          l[3],
                             unsigned int i,
                             unsigned int n,
                             unsigned int n_seq,
                             const short  **S,
                             unsigned int mm[2])
{
  unsigned int  s, layer_pos[4], mismatch, cnt;

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
                 void *NA VRNA_UNUSED,
                 void *NA2 VRNA_UNUSED)
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
gquad_interact(unsigned int  i,
               unsigned int  L,
               unsigned int  *l,
               void *data,
               void *pf,
               void *index,
               void *NA2 VRNA_UNUSED)
{
  unsigned int  x;
  int           *idx;
  FLT_OR_DBL    gq, *pp;

  idx = (int *)index;
  pp  = (FLT_OR_DBL *)data;

  gq  = vrna_exp_E_gquad(L, l, (vrna_exp_param_t *)pf);

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
                   void *NA VRNA_UNUSED)
{
  unsigned int  x;
  int           *idx, bad;
  FLT_OR_DBL    gq, *pp;

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

PUBLIC FLT_OR_DBL *
get_gquad_pf_matrix(short             *S,
                    FLT_OR_DBL        *scale,
                    vrna_exp_param_t  *pf)
{
  unsigned int  n, size, *gg, i, j;
  int           *my_index;
  FLT_OR_DBL    *data;


  n         = S[0];
  size      = (n * (n + 1)) / 2 + 2;
  data      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  gg        = get_g_islands(S);
  my_index  = vrna_idx_row_wise(n);

  FOR_EACH_GQUAD(i, j, 1, n){
    process_gquad_enumeration(gg, i, j,
                              &gquad_pf,
                              (void *)(&(data[my_index[i] - j])),
                              (void *)pf,
                              NULL,
                              NULL);
    data[my_index[i] - j] *= scale[j - i + 1];
  }

  free(my_index);
  free(gg);
  return data;
}


PUBLIC FLT_OR_DBL *
get_gquad_pf_matrix_comparative(unsigned int      n,
                                short             *S_cons,
                                short             **S,
                                unsigned int      **a2s,
                                FLT_OR_DBL        *scale,
                                unsigned int      n_seq,
                                vrna_exp_param_t  *pf)
{
  unsigned int         size, *gg, i, j;
  int           *my_index;
  FLT_OR_DBL  *data;


  size      = (n * (n + 1)) / 2 + 2;
  data      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  gg        = get_g_islands(S_cons);
  my_index  = vrna_idx_row_wise(n);

  struct gquad_ali_helper gq_help = {
    .S      = (const short **)S,
    .a2s    = (const unsigned int **)a2s,
    .length = n,
    .n_seq  = n_seq,
    .pf     = pf
  };

  FOR_EACH_GQUAD(i, j, 1, n){
    process_gquad_enumeration(gg, i, j,
                              &gquad_pf_ali,
                              (void *)(&(data[my_index[i] - j])),
                              (void *)&gq_help,
                              NULL,
                              NULL);
    data[my_index[i] - j] *= scale[j - i + 1];
  }

  free(my_index);
  free(gg);
  return data;
}


#endif

