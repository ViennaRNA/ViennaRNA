#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/mfe/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/internal.h"


#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/internal_hc.inc"
#include "ViennaRNA/constraints/internal_sc.inc"

#include "ViennaRNA/mfe/internal.h"


typedef struct {
  eval_hc               hc_wrapper;
  struct hc_int_def_dat hc_dat;
  struct sc_int_dat     sc_wrapper;

  unsigned int          *tt;
} helper_data_t;


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
mfe_internal_loop(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j);


PRIVATE int
mfe_internal_loop_ext(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      unsigned int          *ip,
                      unsigned int          *iq);


PRIVATE helper_data_t *
get_intloop_helpers(vrna_fold_compound_t  *fc,
                    unsigned int          i,
                    unsigned int          j);


PRIVATE void
free_intloop_helpers(helper_data_t *h);


PRIVATE int
mfe_stacks(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           helper_data_t        *helpers);


PRIVATE int
mfe_bulges(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           helper_data_t        *helpers);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_mfe_internal(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j)
{
  if (fc)
    return mfe_internal_loop(fc, i, j);

  return INF;
}


PUBLIC int
vrna_mfe_internal_ext(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      unsigned int          *ip,
                      unsigned int          *iq)
{
  if (fc)
    return mfe_internal_loop_ext(fc, i, j, ip, iq);

  return INF;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE helper_data_t *
get_intloop_helpers(vrna_fold_compound_t  *fc,
                    unsigned int          i,
                    unsigned int          j)
{
  helper_data_t *h = (helper_data_t *)vrna_alloc(sizeof(helper_data_t));

  /* init hard constraints wrapper */
  h->hc_wrapper = prepare_hc_int_def(fc, &(h->hc_dat));

  /* init soft constraints wrapper */
  init_sc_int(fc, &(h->sc_wrapper));


  if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
    vrna_md_t *md = &(fc->params->model_details);
    short     **S = fc->S;

    h->tt = (unsigned int *)vrna_alloc(sizeof(unsigned int) * fc->n_seq);

    for (unsigned int s = 0; s < fc->n_seq; s++)
      h->tt[s] = vrna_get_ptype_md(S[s][i], S[s][j], md);
  } else {
    h->tt = NULL;
  }

  return h;
}


PRIVATE void
free_intloop_helpers(helper_data_t *h)
{
  free_sc_int(&(h->sc_wrapper));
  free(h->tt);
  free(h);
}


PRIVATE int
mfe_stacks(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           helper_data_t        *helpers)
{
  unsigned char sliding_window, hc_decompose, *hc_mx, **hc_mx_local, eval;
  char          *ptype, **ptype_local;
  short         *S, **SS, **S5, **S3;
  unsigned int  k, l, *sn, n_seq, s, n, type, type2;
  int           e, eee, *idx, ij, *c, *rtype, **c_local, kl;
  vrna_param_t  *P;
  vrna_md_t     *md;

  e = INF;

  n               = fc->length;
  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  sn              = fc->strand_number;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  idx             = fc->jindx;
  ij              = (sliding_window) ? 0 : idx[j] + i;
  hc_mx           = (sliding_window) ? NULL : fc->hc->mx;
  hc_mx_local     = (sliding_window) ? fc->hc->matrix_local : NULL;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     =
    (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? fc->ptype_local : NULL) : NULL;
  S       = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3      = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  c       = (sliding_window) ? NULL : fc->matrices->c;
  c_local = (sliding_window) ? fc->matrices->c_local : NULL;
  P       = fc->params;
  md      = &(P->model_details);
  rtype   = &(md->rtype[0]);

  k = i + 1;
  l = j - 1;

  if (k < l) {
    /* first, determine whether we are allowed to form this stack */
    hc_decompose  = (sliding_window) ? hc_mx_local[i][j - i] : hc_mx[n * i + j];
    eval          = hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP;
    hc_decompose  = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[n * k + l];
    eval          = (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) ? eval : (unsigned char)0;
    eval          = (helpers->hc_wrapper(i, j, k, l, &(helpers->hc_dat))) ? eval : (unsigned char)0;

    if ((sn[i] != sn[k]) ||
        (sn[l] != sn[j]))
      eval = (unsigned char)0;

    if (eval) {
      type = 0;

      if (fc->type == VRNA_FC_TYPE_SINGLE)
        type = sliding_window ?
               vrna_get_ptype_window(i, j, ptype_local) :
               vrna_get_ptype(ij, ptype);

      kl  = (sliding_window) ? 0 : idx[l] + k;
      eee = (sliding_window) ? c_local[k][l - k] : c[kl];

      if (eee != INF) {
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type2 = sliding_window ?
                    rtype[vrna_get_ptype_window(k, l, ptype_local)] :
                    rtype[vrna_get_ptype(kl, ptype)];

            eee += vrna_E_internal(0, 0, type, type2, S[i + 1], S[j - 1], S[i], S[j], P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type2 = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
              eee   += vrna_E_internal(0,
                                       0,
                                       helpers->tt[s],
                                       type2,
                                       S3[s][i],
                                       S5[s][j],
                                       S5[s][k],
                                       S3[s][l],
                                       P);
            }

            break;
        }
        if (helpers->sc_wrapper.pair)
          eee += helpers->sc_wrapper.pair(i, j, k, l, &(helpers->sc_wrapper));

        e = MIN2(e, eee);
      }
    }
  }

  return e;
}


PRIVATE int
mfe_bulges(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           helper_data_t        *helpers)
{
  unsigned char sliding_window, hc_decompose, *hc_mx, **hc_mx_local;
  char          *ptype, **ptype_local;
  short         *S, **SS, **S5, **S3;
  unsigned int  *sn, **a2s, n_seq, s, n, *hc_up, with_ud, noclose,
                type, type2, has_nick, k, l, last_k, first_l, u1, u2,
                noGUclosure, u1_local;
  int           e, eee, *idx, ij, kl, *c, *rtype, **c_local;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;

  e = INF;

  n               = fc->length;
  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  sn              = fc->strand_number;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  idx             = fc->jindx;
  ij              = (sliding_window) ? 0 : idx[j] + i;
  hc_mx           = (sliding_window) ? NULL : fc->hc->mx;
  hc_mx_local     = (sliding_window) ? fc->hc->matrix_local : NULL;
  hc_up           = fc->hc->up_int;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     =
    (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? fc->ptype_local : NULL) : NULL;
  S           = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s         = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  c           = (sliding_window) ? NULL : fc->matrices->c;
  c_local     = (sliding_window) ? fc->matrices->c_local : NULL;
  P           = fc->params;
  md          = &(P->model_details);
  rtype       = &(md->rtype[0]);
  domains_up  = fc->domains_up;
  with_ud     = ((domains_up) && (domains_up->energy_cb)) ? 1 : 0;

  hc_decompose = (sliding_window) ? hc_mx_local[i][j - i] : hc_mx[n * i + j];

  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    has_nick    = sn[i] != sn[j] ? 1 : 0;
    noGUclosure = md->noGUclosure;
    type        = 0;

    if (fc->type == VRNA_FC_TYPE_SINGLE)
      type = sliding_window ?
             vrna_get_ptype_window(i, j, ptype_local) :
             vrna_get_ptype(ij, ptype);

    noclose = ((noGUclosure) && (type == 3 || type == 4)) ? 1 : 0;

    if (!noclose) {
      /* only proceed if the enclosing pair is allowed */

      /* handle bulges in 5' side */
      l = j - 1;
      if (l > i + 2) {
        last_k = l - 1;

        if (last_k > i + 1 + MAXLOOP)
          last_k = i + 1 + MAXLOOP;

        if (last_k > i + 1 + hc_up[i + 1])
          last_k = i + 1 + hc_up[i + 1];

        u1 = 1;

        k   = i + 2;
        kl  = (sliding_window) ? 0 : idx[l] + k;

        hc_mx += n * l;

        for (; k <= last_k; k++, u1++, kl++) {
          hc_decompose = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[k];

          if ((hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (helpers->hc_wrapper(i, j, k, l, &(helpers->hc_dat)))) {
            eee = (sliding_window) ? c_local[k][l - k] : c[kl];

            if (eee < INF) {
              switch (fc->type) {
                case VRNA_FC_TYPE_SINGLE:
                  type2 = sliding_window ?
                          rtype[vrna_get_ptype_window(k, l, ptype_local)] :
                          rtype[vrna_get_ptype(kl, ptype)];

                  if ((noGUclosure) && (type2 == 3 || type2 == 4))
                    continue;

                  if ((has_nick) && ((sn[i] != sn[k]) || (sn[j - 1] != sn[j]))) {
                    eee = INF;
                  } else {
                    eee += vrna_E_internal(u1,
                                           0,
                                           type,
                                           type2,
                                           S[i + 1],
                                           S[j - 1],
                                           S[k - 1],
                                           S[l + 1],
                                           P);
                  }

                  break;

                case VRNA_FC_TYPE_COMPARATIVE:
                  for (s = 0; s < n_seq; s++) {
                    u1_local  = a2s[s][k - 1] - a2s[s][i];
                    type2     = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    eee       +=
                      vrna_E_internal(u1_local, 0, helpers->tt[s], type2, S3[s][i], S5[s][j],
                                      S5[s][k], S3[s][l],
                                      P);
                  }

                  break;
              }

              if (helpers->sc_wrapper.pair)
                eee += helpers->sc_wrapper.pair(i, j, k, l, &(helpers->sc_wrapper));

              e = MIN2(e, eee);

              if (with_ud) {
                eee += domains_up->energy_cb(fc,
                                             i + 1, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                             domains_up->data);
                e = MIN2(e, eee);
              }
            }
          }
        }

        hc_mx -= n * l;
      }

      /* handle bulges in 3' side */
      k = i + 1;
      if (k + 2 < j) {
        first_l = k + 1;
        if (first_l + MAXLOOP + 1 < j)
          first_l = j - 1 - MAXLOOP;

        u2    = 1;
        hc_mx += n * k;

        for (l = j - 2; l >= first_l; l--, u2++) {
          if (u2 > hc_up[l + 1])
            break;

          kl            = (sliding_window) ? 0 : idx[l] + k;
          hc_decompose  = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[l];

          if ((hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (helpers->hc_wrapper(i, j, k, l, &(helpers->hc_dat)))) {
            eee = (sliding_window) ? c_local[k][l - k] : c[kl];

            if (eee < INF) {
              switch (fc->type) {
                case VRNA_FC_TYPE_SINGLE:
                  type2 = sliding_window ?
                          rtype[vrna_get_ptype_window(k, l, ptype_local)] :
                          rtype[vrna_get_ptype(kl, ptype)];

                  if ((noGUclosure) && (type2 == 3 || type2 == 4))
                    continue;

                  if ((has_nick) && ((sn[i] != sn[i + 1]) || (sn[j] != sn[l]))) {
                    eee = INF;
                  } else {
                    eee += vrna_E_internal(0,
                                           u2,
                                           type,
                                           type2,
                                           S[i + 1],
                                           S[j - 1],
                                           S[k - 1],
                                           S[l + 1],
                                           P);
                  }

                  break;

                case VRNA_FC_TYPE_COMPARATIVE:
                  for (s = 0; s < n_seq; s++) {
                    unsigned int u2_local = a2s[s][j - 1] - a2s[s][l];
                    type2 = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    eee   +=
                      vrna_E_internal(0, u2_local, helpers->tt[s], type2, S3[s][i], S5[s][j],
                                      S5[s][k], S3[s][l],
                                      P);
                  }

                  break;
              }

              if (helpers->sc_wrapper.pair)
                eee += helpers->sc_wrapper.pair(i, j, k, l, &(helpers->sc_wrapper));

              e = MIN2(e, eee);

              if (with_ud) {
                eee += domains_up->energy_cb(fc,
                                             l + 1, j - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                             domains_up->data);
                e = MIN2(e, eee);
              }
            }
          }
        }

        hc_mx -= n * k;
      }
    }
  }

  return e;
}


PRIVATE int
mfe_internal_loop(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j)
{
  unsigned char sliding_window, hc_decompose, *hc_mx, **hc_mx_local;
  char          *ptype, **ptype_local;
  short         *S, **SS, **S5, **S3;
  unsigned int  *sn, **a2s, n_seq, s, n, *hc_up, type, type2, has_nick,
                k, l, last_k, first_l, u1, u2, noGUclosure, with_ud,
                with_gquad, noclose, u1_local, u2_local;
  int           e, eee, e3, e5, *idx, ij, kl, *c, *rtype, **c_local;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;
  helper_data_t *helpers;

  helpers = get_intloop_helpers(fc, i, j);

  e = INF;

  n               = fc->length;
  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  sn              = fc->strand_number;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  idx             = fc->jindx;
  ij              = (sliding_window) ? 0 : idx[j] + i;
  hc_mx           = (sliding_window) ? NULL : fc->hc->mx;
  hc_mx_local     = (sliding_window) ? fc->hc->matrix_local : NULL;
  hc_up           = fc->hc->up_int;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     =
    (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? fc->ptype_local : NULL) : NULL;
  S           = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s         = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  c           = (sliding_window) ? NULL : fc->matrices->c;
  c_local     = (sliding_window) ? fc->matrices->c_local : NULL;
  P           = fc->params;
  md          = &(P->model_details);
  rtype       = &(md->rtype[0]);
  domains_up  = fc->domains_up;
  with_ud     = ((domains_up) && (domains_up->energy_cb)) ? 1 : 0;
  with_gquad  = md->gquad;

  hc_decompose = (sliding_window) ? hc_mx_local[i][j - i] : hc_mx[n * i + j];

  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    has_nick    = sn[i] != sn[j] ? 1 : 0;
    noGUclosure = md->noGUclosure;
    type        = 0;

    if (fc->type == VRNA_FC_TYPE_SINGLE)
      type = sliding_window ?
             vrna_get_ptype_window(i, j, ptype_local) :
             vrna_get_ptype(ij, ptype);

    noclose = ((noGUclosure) && (type == 3 || type == 4)) ? 1 : 0;

    e = MIN2(e, mfe_stacks(fc, i, j, helpers));

    e = MIN2(e, mfe_bulges(fc, i, j, helpers));

    if (!noclose) {
      /* only proceed if the enclosing pair is allowed */

      /* last but not least, all other internal loops */
      first_l = i + 2 + 1;
      if (first_l + MAXLOOP + 1 < j)
        first_l = j - 1 - MAXLOOP;

      u2 = 1;
      for (l = j - 2; l >= first_l; l--, u2++) {
        if (u2 > hc_up[l + 1])
          break;

        last_k = l - 1;

        if (last_k + u2 > i + 1 + MAXLOOP)
          last_k = i + 1 + MAXLOOP - u2;

        if (last_k > i + 1 + hc_up[i + 1])
          last_k = i + 1 + hc_up[i + 1];

        u1  = 1;
        k   = i + 2;
        kl  = (sliding_window) ? 0 : idx[l] + k;

        hc_mx += n * l;

        for (; k <= last_k; k++, u1++, kl++) {
          hc_decompose = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[k];

          if ((hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (helpers->hc_wrapper(i, j, k, l, &(helpers->hc_dat)))) {
            eee = (sliding_window) ? c_local[k][l - k] : c[kl];

            if (eee < INF) {
              switch (fc->type) {
                case VRNA_FC_TYPE_SINGLE:
                  type2 = sliding_window ?
                          rtype[vrna_get_ptype_window(k, l, ptype_local)] :
                          rtype[vrna_get_ptype(kl, ptype)];

                  if ((noGUclosure) && (type2 == 3 || type2 == 4))
                    continue;

                  if ((has_nick) && ((sn[i] != sn[k]) || (sn[j] != sn[l]))) {
                    eee = INF;
                  } else {
                    eee +=
                      vrna_E_internal(u1, u2, type, type2, S[i + 1], S[j - 1], S[k - 1], S[l + 1],
                                      P);
                  }

                  break;

                case VRNA_FC_TYPE_COMPARATIVE:
                  for (s = 0; s < n_seq; s++) {
                    u1_local  = a2s[s][k - 1] - a2s[s][i];
                    u2_local  = a2s[s][j - 1] - a2s[s][l];
                    type2     = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    eee       += vrna_E_internal(u1_local,
                                                 u2_local,
                                                 helpers->tt[s],
                                                 type2,
                                                 S3[s][i],
                                                 S5[s][j],
                                                 S5[s][k],
                                                 S3[s][l],
                                                 P);
                  }

                  break;
              }

              if (helpers->sc_wrapper.pair)
                eee += helpers->sc_wrapper.pair(i, j, k, l, &(helpers->sc_wrapper));

              e = MIN2(e, eee);

              if (with_ud) {
                e5 = domains_up->energy_cb(fc,
                                           i + 1, k - 1,
                                           VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                           domains_up->data);
                e3 = domains_up->energy_cb(fc,
                                           l + 1, j - 1,
                                           VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                           domains_up->data);

                e = MIN2(e, eee + e5);
                e = MIN2(e, eee + e3);
                e = MIN2(e, eee + e5 + e3);
              }
            }
          }
        }

        hc_mx -= n * l;
      }

      if (with_gquad) {
        /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
        eee = vrna_mfe_gquad_internal_loop(fc, i, j);
        e   = MIN2(e, eee);
      }
    }
  }

  free_intloop_helpers(helpers);

  return e;
}


PRIVATE int
mfe_internal_loop_ext(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      unsigned int          *ip,
                      unsigned int          *iq)
{
  int                   energy, e, *indx, *c;
  unsigned char         *hc, eval_loop;
  unsigned int          n, q, p, s, n_seq, u1, u2, qmin, *tt, *hc_up;
  short                 **SS;
  vrna_md_t             *md;
  vrna_param_t          *P;
  eval_hc               evaluate;
  struct hc_int_def_dat hc_dat_local;

  n     = fc->length;
  n_seq = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  SS    = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  indx  = fc->jindx;
  c     = fc->matrices->c;
  hc    = fc->hc->mx;
  hc_up = fc->hc->up_int;
  P     = fc->params;
  md    = &(P->model_details);
  tt    = NULL;

  e = INF;

  evaluate = prepare_hc_int_def(fc, &hc_dat_local);

  /* CONSTRAINED internal LOOP start */
  if (hc[n * i + j] & (VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) {
    /* prepare necessary variables */
    if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      tt = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);

      for (s = 0; s < n_seq; s++)
        tt[s] = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
    }

    for (p = j + 1; p < n; p++) {
      u1 = p - j - 1;
      if (u1 + i - 1 > MAXLOOP)
        break;

      if (hc_up[j + 1] < u1)
        break;

      qmin = p + 1;

      if (u1 + i + n > p + MAXLOOP + 2)
        qmin = u1 + i + n - MAXLOOP - 1;

      for (q = n; q >= qmin; q--) {
        u2 = i - 1 + n - q;
        if (hc_up[q + 1] < u2)
          break;

        if (u1 + u2 > MAXLOOP)
          continue;

        int pq = indx[q] + p;

        eval_loop = hc[n * p + q] &
                    (VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC);

        if (eval_loop && evaluate(i, j, p, q, &hc_dat_local)) {
          energy = c[pq];
          if (energy < INF) {
            energy += vrna_eval_internal(fc, i, j, p, q, VRNA_EVAL_LOOP_NO_HC);

            if (energy < e) {
              e = energy;
              if ((ip != NULL) && (iq != NULL)) {
                *ip = p;
                *iq = q;
              }
            }
          }
        }
      }
    }
  }

  free(tt);

  return e;
}


/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC int
vrna_E_int_loop(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j)
{
  return vrna_mfe_internal(fc,
                           (unsigned int)i,
                           (unsigned int)j);
}


PUBLIC int
vrna_E_ext_int_loop(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   *ip,
                    int                   *iq)
{
  int e = INF;

  if (fc) {
    unsigned int p, q;
    e = vrna_mfe_internal_ext(fc,
                              (unsigned int)i,
                              (unsigned int)j,
                              &p,
                              &q);
    if (ip)
      *ip = (int)p;

    if (iq)
      *iq = (int)q;
  }

  return e;
}


#endif
