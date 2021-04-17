#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/loops/external.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/loops/internal.h"


#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "internal_hc.inc"
#include "internal_sc.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
E_internal_loop(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j);


PRIVATE int
E_ext_internal_loop(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   *ip,
                    int                   *iq);


PRIVATE int
E_stack(vrna_fold_compound_t  *fc,
        int                   i,
        int                   j);


PRIVATE INLINE int
eval_int_loop(vrna_fold_compound_t  *fc,
              int                   i,
              int                   j,
              int                   k,
              int                   l);


PRIVATE INLINE int
eval_ext_int_loop(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  int                   k,
                  int                   l);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_E_int_loop(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j)
{
  int e = INF;

  if (fc)
    e = E_internal_loop(fc, i, j);

  return e;
}


PUBLIC int
vrna_E_ext_int_loop(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   *ip,
                    int                   *iq)
{
  int e = INF;

  if (fc)
    e = E_ext_internal_loop(fc, i, j, ip, iq);

  return e;
}


PUBLIC int
vrna_E_stack(vrna_fold_compound_t *fc,
             int                  i,
             int                  j)
{
  int e = INF;

  if ((fc) && (i > 0) && (i < j) && (j - i > 3))
    e = E_stack(fc, i, j);

  return e;
}


PUBLIC int
vrna_eval_int_loop(vrna_fold_compound_t *fc,
                   int                  i,
                   int                  j,
                   int                  k,
                   int                  l)
{
  int e = INF;

  if (fc)
    e = eval_int_loop(fc, i, j, k, l);

  return e;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE int
eval_int_loop(vrna_fold_compound_t  *fc,
              int                   i,
              int                   j,
              int                   k,
              int                   l)
{
  unsigned int      *sn, *ss, n_seq, s, **a2s;
  int               e, *rtype, type, type2, with_ud;
  short             *S, *S2, **SS, **S5, **S3;
  vrna_param_t      *P;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;
  struct sc_int_dat sc_wrapper;

  n_seq       = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  P           = fc->params;
  md          = &(P->model_details);
  sn          = fc->strand_number;
  ss          = fc->strand_start;
  rtype       = &(md->rtype[0]);
  S           = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  S2          = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s         = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  domains_up  = fc->domains_up;
  with_ud     = ((domains_up) && (domains_up->energy_cb)) ? 1 : 0;
  e           = INF;

  init_sc_int(fc, &sc_wrapper);

  {
    int energy, e5, e3, u1, u2;

    energy = 0;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type  = vrna_get_ptype_md(S2[i], S2[j], md);
        type2 = vrna_get_ptype_md(S2[l], S2[k], md);

        u1  = k - i - 1;
        u2  = j - l - 1;

        if ((sn[i] == sn[k]) && (sn[l] == sn[j])) {
          /* regular interior loop */
          energy = E_IntLoop(u1, u2, type, type2, S[i + 1], S[j - 1], S[k - 1], S[l + 1], P);
        } else {
          /* interior loop like cofold structure */
          short Si, Sj;
          Si      = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;
          Sj      = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;
          energy  = E_IntLoop_Co(rtype[type], rtype[type2],
                                 i, j, k, l,
                                 ss[fc->strand_order[1]], /* serves as cut point substitute */
                                 Si, Sj,
                                 S[k - 1], S[l + 1],
                                 md->dangles,
                                 P);
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++) {
          type    = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
          type2   = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
          u1      = a2s[s][k - 1] - a2s[s][i];
          u2      = a2s[s][j - 1] - a2s[s][l];
          energy  += E_IntLoop(u1, u2, type, type2, S3[s][i], S5[s][j], S5[s][k], S3[s][l], P);
        }

        break;
    }

    /* add soft constraints */
    if (sc_wrapper.pair)
      energy += sc_wrapper.pair(i, j, k, l, &sc_wrapper);

    e = energy;

    if (with_ud) {
      e5 = e3 = 0;

      u1  = k - i - 1;
      u2  = j - l - 1;

      if (u1 > 0) {
        e5 = domains_up->energy_cb(fc,
                                   i + 1, k - 1,
                                   VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                   domains_up->data);
      }

      if (u2 > 0) {
        e3 = domains_up->energy_cb(fc,
                                   l + 1, j - 1,
                                   VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                   domains_up->data);
      }

      e = MIN2(e, energy + e5);
      e = MIN2(e, energy + e3);
      e = MIN2(e, energy + e5 + e3);
    }
  }

  free_sc_int(&sc_wrapper);

  return e;
}


PRIVATE INLINE int
eval_ext_int_loop(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  int                   k,
                  int                   l)
{
  unsigned int      n, n_seq, s, **a2s;
  int               e, type, type2, with_ud;
  short             *S, *S2, **SS, **S5, **S3;
  vrna_param_t      *P;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;
  struct sc_int_dat sc_wrapper;

  n           = fc->length;
  n_seq       = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  P           = fc->params;
  md          = &(P->model_details);
  S           = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  S2          = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s         = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  domains_up  = fc->domains_up;
  with_ud     = ((domains_up) && (domains_up->energy_cb)) ? 1 : 0;
  e           = INF;

  init_sc_int(fc, &sc_wrapper);

  {
    int energy, e5, e3, u1, u2, u3;

    energy = 0;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type  = vrna_get_ptype_md(S2[j], S2[i], md);
        type2 = vrna_get_ptype_md(S2[l], S2[k], md);

        u1  = i - 1;
        u2  = k - j - 1;
        u3  = n - l;

        /* regular interior loop */
        energy = E_IntLoop(u2, u1 + u3, type, type2, S[j + 1], S[i - 1], S[k - 1], S[l + 1], P);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++) {
          type    = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
          type2   = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
          u1      = a2s[s][i - 1];
          u2      = a2s[s][k - 1] - a2s[s][j];
          u3      = a2s[s][n] - a2s[s][l];
          energy  += E_IntLoop(u2, u1 + u3, type, type2, S3[s][j], S5[s][i], S5[s][k], S3[s][l], P);
        }

        break;
    }

    /* add soft constraints */
    if (sc_wrapper.pair_ext)
      energy += sc_wrapper.pair_ext(i, j, k, l, &sc_wrapper);

    e = energy;

    if (with_ud) {
      e5 = e3 = 0;

      u1  = i - 1;
      u2  = k - j - 1;
      u3  = n - l;

      if (u2 > 0) {
        e5 = domains_up->energy_cb(fc,
                                   j + 1, k - 1,
                                   VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                   domains_up->data);
      }

      if (u1 + u3 > 0) {
        e3 = domains_up->energy_cb(fc,
                                   l + 1, i - 1,
                                   VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                   domains_up->data);
      }

      e = MIN2(e, energy + e5);
      e = MIN2(e, energy + e3);
      e = MIN2(e, energy + e5 + e3);
    }
  }

  free_sc_int(&sc_wrapper);

  return e;
}


PRIVATE int
E_internal_loop(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j)
{
  unsigned char         sliding_window, hc_decompose, *hc_mx, **hc_mx_local;
  char                  *ptype, **ptype_local;
  short                 *S, **SS, **S5, **S3;
  unsigned int          *sn, *ss, **a2s, n_seq, s, n;
  int                   e, eee, *idx, ij, *c, *ggg, *rtype, with_ud, with_gquad, noclose,
                        *hc_up, **c_local, **ggg_local;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  struct hc_int_def_dat hc_dat_local;
  eval_hc               *evaluate;
  struct sc_int_dat     sc_wrapper;

  evaluate = prepare_hc_int_def(fc, &hc_dat_local);
  init_sc_int(fc, &sc_wrapper);

  e = INF;

  n               = fc->length;
  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  sn              = fc->strand_number;
  ss              = fc->strand_start;
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
  ggg         = (sliding_window) ? NULL : fc->matrices->ggg;
  c_local     = (sliding_window) ? fc->matrices->c_local : NULL;
  ggg_local   = (sliding_window) ? fc->matrices->ggg_local : NULL;
  P           = fc->params;
  md          = &(P->model_details);
  rtype       = &(md->rtype[0]);
  domains_up  = fc->domains_up;
  with_ud     = ((domains_up) && (domains_up->energy_cb)) ? 1 : 0;
  with_gquad  = md->gquad;

  hc_decompose = (sliding_window) ? hc_mx_local[i][j - i] : hc_mx[n * i + j];

  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    unsigned int  type, type2, has_nick, *tt;
    int           k, l, kl, last_k, first_l, u1, u2, turn, noGUclosure;

    has_nick    = sn[i] != sn[j] ? 1 : 0;
    turn        = md->min_loop_size;
    noGUclosure = md->noGUclosure;
    tt          = NULL;
    type        = 0;

    if (fc->type == VRNA_FC_TYPE_SINGLE)
      type = sliding_window ?
             vrna_get_ptype_window(i, j, ptype_local) :
             vrna_get_ptype(ij, ptype);

    noclose = ((noGUclosure) && (type == 3 || type == 4)) ? 1 : 0;

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      tt = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);
      for (s = 0; s < n_seq; s++)
        tt[s] = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
    }

    /* handle stacks separately */
    k = i + 1;
    l = j - 1;
    if (k < l) {
      kl            = (sliding_window) ? 0 : idx[l] + k;
      hc_decompose  = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[n * k + l];

      if ((hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
          (evaluate(i, j, k, l, &hc_dat_local))) {
        eee = (sliding_window) ? c_local[k][l - k] : c[kl];

        if (eee != INF) {
          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              type2 = sliding_window ?
                      rtype[vrna_get_ptype_window(k, l, ptype_local)] :
                      rtype[vrna_get_ptype(kl, ptype)];

              if ((has_nick) && ((sn[i] != sn[i + 1]) || (sn[j - 1] != sn[j]))) {
                /* interior loop like cofold structure */
                short Si, Sj;
                Si  = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;
                Sj  = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;
                eee += E_IntLoop_Co(rtype[type], rtype[type2],
                                    i, j, k, l,
                                    ss[fc->strand_order[1]], /* serves as cutpoint replacement */
                                    Si, Sj,
                                    S[i], S[j],
                                    md->dangles,
                                    P);
              } else {
                eee += E_IntLoop(0, 0, type, type2, S[i + 1], S[j - 1], S[i], S[j], P);
              }

              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                type2 = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                eee   += E_IntLoop(0, 0, tt[s], type2, S3[s][i], S5[s][j], S5[s][k], S3[s][l], P);
              }

              break;
          }

          if (sc_wrapper.pair)
            eee += sc_wrapper.pair(i, j, k, l, &sc_wrapper);

          e = MIN2(e, eee);
        }
      }
    }

    if (!noclose) {
      /* only proceed if the enclosing pair is allowed */

      /* handle bulges in 5' side */
      l = j - 1;
      if (l > i + 2) {
        last_k = l - turn - 1;

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
              (evaluate(i, j, k, l, &hc_dat_local))) {
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
                    /* interior loop like cofold structure */
                    short Si, Sj;
                    Si  = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;
                    Sj  = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;
                    eee += E_IntLoop_Co(rtype[type], rtype[type2],
                                        i, j, k, l,
                                        ss[fc->strand_order[1]],
                                        Si, Sj,
                                        S[k - 1], S[j],
                                        md->dangles,
                                        P);
                  } else {
                    eee += E_IntLoop(u1, 0, type, type2, S[i + 1], S[j - 1], S[k - 1], S[l + 1], P);
                  }

                  break;

                case VRNA_FC_TYPE_COMPARATIVE:
                  for (s = 0; s < n_seq; s++) {
                    int u1_local = a2s[s][k - 1] - a2s[s][i];
                    type2 = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    eee   +=
                      E_IntLoop(u1_local, 0, tt[s], type2, S3[s][i], S5[s][j], S5[s][k], S3[s][l],
                                P);
                  }

                  break;
              }

              if (sc_wrapper.pair)
                eee += sc_wrapper.pair(i, j, k, l, &sc_wrapper);

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
      if (k < j - 2) {
        first_l = k + turn + 1;
        if (first_l < j - 1 - MAXLOOP)
          first_l = j - 1 - MAXLOOP;

        u2    = 1;
        hc_mx += n * k;

        for (l = j - 2; l >= first_l; l--, u2++) {
          if (u2 > hc_up[l + 1])
            break;

          kl            = (sliding_window) ? 0 : idx[l] + k;
          hc_decompose  = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[l];

          if ((hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (evaluate(i, j, k, l, &hc_dat_local))) {
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
                    /* interior loop like cofold structure */
                    short Si, Sj;
                    Si  = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;
                    Sj  = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;
                    eee += E_IntLoop_Co(rtype[type], rtype[type2],
                                        i, j, k, l,
                                        ss[fc->strand_order[1]],
                                        Si, Sj,
                                        S[i], S[l + 1],
                                        md->dangles,
                                        P);
                  } else {
                    eee += E_IntLoop(0, u2, type, type2, S[i + 1], S[j - 1], S[k - 1], S[l + 1], P);
                  }

                  break;

                case VRNA_FC_TYPE_COMPARATIVE:
                  for (s = 0; s < n_seq; s++) {
                    int u2_local = a2s[s][j - 1] - a2s[s][l];
                    type2 = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    eee   +=
                      E_IntLoop(0, u2_local, tt[s], type2, S3[s][i], S5[s][j], S5[s][k], S3[s][l],
                                P);
                  }

                  break;
              }

              if (sc_wrapper.pair)
                eee += sc_wrapper.pair(i, j, k, l, &sc_wrapper);

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

      /* last but not least, all other internal loops */
      first_l = i + 2 + turn + 1;
      if (first_l < j - 1 - MAXLOOP)
        first_l = j - 1 - MAXLOOP;

      u2 = 1;
      for (l = j - 2; l >= first_l; l--, u2++) {
        if (u2 > hc_up[l + 1])
          break;

        last_k = l - turn - 1;

        if (last_k > i + 1 + MAXLOOP - u2)
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
              (evaluate(i, j, k, l, &hc_dat_local))) {
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
                    /* interior loop like cofold structure */
                    short Si, Sj;
                    Si  = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;
                    Sj  = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;
                    eee += E_IntLoop_Co(rtype[type], rtype[type2],
                                        i, j, k, l,
                                        ss[fc->strand_order[1]],
                                        Si, Sj,
                                        S[k - 1], S[l + 1],
                                        md->dangles,
                                        P);
                  } else {
                    eee +=
                      E_IntLoop(u1, u2, type, type2, S[i + 1], S[j - 1], S[k - 1], S[l + 1], P);
                  }

                  break;

                case VRNA_FC_TYPE_COMPARATIVE:
                  for (s = 0; s < n_seq; s++) {
                    int u1_local  = a2s[s][k - 1] - a2s[s][i];
                    int u2_local  = a2s[s][j - 1] - a2s[s][l];
                    type2 = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    eee   += E_IntLoop(u1_local,
                                       u2_local,
                                       tt[s],
                                       type2,
                                       S3[s][i],
                                       S5[s][j],
                                       S5[s][k],
                                       S3[s][l],
                                       P);
                  }

                  break;
              }

              if (sc_wrapper.pair)
                eee += sc_wrapper.pair(i, j, k, l, &sc_wrapper);

              e = MIN2(e, eee);

              if (with_ud) {
                int e5, e3;

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
        eee = INF;

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            if (sliding_window)
              eee = E_GQuad_IntLoop_L(i, j, type, S, ggg_local, fc->window_size, P);
            else if (sn[j] == sn[i])
              eee = E_GQuad_IntLoop(i, j, type, S, ggg, idx, P);

            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            if (sliding_window) {
              eee = E_GQuad_IntLoop_L_comparative(i,
                                                  j,
                                                  tt,
                                                  fc->S_cons,
                                                  S5,
                                                  S3,
                                                  a2s,
                                                  ggg_local,
                                                  n_seq,
                                                  P);
            } else {
              eee = E_GQuad_IntLoop_comparative(i,
                                                j,
                                                tt,
                                                fc->S_cons,
                                                S5,
                                                S3,
                                                a2s,
                                                ggg,
                                                idx,
                                                n_seq,
                                                P);
            }

            break;
        }

        e = MIN2(e, eee);
      }

      free(tt);
    }
  }

  free_sc_int(&sc_wrapper);

  return e;
}


PRIVATE int
E_ext_internal_loop(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   *ip,
                    int                   *iq)
{
  int                   q, p, e, s, u1, u2, qmin, energy,
                        n, *indx, *hc_up, *c, turn, n_seq;
  unsigned char         *hc, eval_loop;
  unsigned int          *tt;
  short                 **SS;
  vrna_md_t             *md;
  vrna_param_t          *P;
  eval_hc               *evaluate;
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
  turn  = md->min_loop_size;
  tt    = NULL;

  e = INF;

  evaluate = prepare_hc_int_def(fc, &hc_dat_local);

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc[n * i + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
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

      qmin = u1 + i - 1 + n - MAXLOOP;
      if (qmin < p + turn + 1)
        qmin = p + turn + 1;

      for (q = n; q >= qmin; q--) {
        u2 = i - 1 + n - q;
        if (hc_up[q + 1] < u2)
          break;

        if (u1 + u2 > MAXLOOP)
          continue;

        int pq = indx[q] + p;

        eval_loop = hc[n * p + q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP;

        if (eval_loop && evaluate(i, j, p, q, &hc_dat_local)) {
          energy = c[pq];
          if (energy < INF) {
            energy += eval_ext_int_loop(fc, i, j, p, q);

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


PRIVATE int
E_stack(vrna_fold_compound_t  *fc,
        int                   i,
        int                   j)
{
  unsigned char         sliding_window, hc_decompose_ij, hc_decompose_pq,
                        *hc_mx, **hc_mx_local, eval_loop;
  char                  *ptype, **ptype_local;
  short                 *S, **SS;
  unsigned int          n, *sn, *ss, type, type_2;
  int                   e, ij, pq, p, q, s, n_seq, *rtype, *indx;
  vrna_param_t          *P;
  vrna_md_t             *md;
  eval_hc               *evaluate;
  struct hc_int_def_dat hc_dat_local;
  struct sc_int_dat     sc_wrapper;

  e               = INF;
  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n               = fc->length;
  p               = i + 1;
  q               = j - 1;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  sn              = fc->strand_number;
  ss              = fc->strand_start;
  S               = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ?
                    (sliding_window ? NULL : fc->ptype) :
                    NULL;
  ptype_local = (fc->type == VRNA_FC_TYPE_SINGLE) ?
                (sliding_window ? fc->ptype_local : NULL) :
                NULL;
  P           = fc->params;
  md          = &(P->model_details);
  rtype       = &(md->rtype[0]);
  indx        = (sliding_window) ? NULL : fc->jindx;
  hc_mx       = (sliding_window) ? NULL : fc->hc->mx;
  hc_mx_local = (sliding_window) ? fc->hc->matrix_local : NULL;
  ij          = (sliding_window) ? 0 : indx[j] + i;
  pq          = (sliding_window) ? 0 : indx[q] + p;
  evaluate    = prepare_hc_int_def(fc, &hc_dat_local);

  init_sc_int(fc, &sc_wrapper);

  hc_decompose_ij = (sliding_window) ? hc_mx_local[i][j - i] : hc_mx[n * i + j];
  hc_decompose_pq = (sliding_window) ? hc_mx_local[p][q - p] : hc_mx[n * p + q];

  eval_loop = (hc_decompose_ij & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (hc_decompose_pq & VRNA_CONSTRAINT_CONTEXT_INT_LOOP);

  if (eval_loop && evaluate(i, j, p, q, &hc_dat_local)) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type = sliding_window ?
               vrna_get_ptype_window(i, j, ptype_local) :
               vrna_get_ptype(ij, ptype);
        type_2 = sliding_window ?
                 rtype[vrna_get_ptype_window(p, q, ptype_local)] :
                 rtype[vrna_get_ptype(pq, ptype)];

        if ((sn[p] == sn[i]) && (sn[j] == sn[q])) {
          /* regular stack */
          e = P->stack[type][type_2];
        } else {
          /* stack like cofold structure */
          short si, sj;
          si  = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;
          sj  = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;
          e   = E_IntLoop_Co(rtype[type], rtype[type_2],
                             i, j, p, q,
                             ss[fc->strand_order[1]],
                             si, sj,
                             S[p - 1], S[q + 1],
                             md->dangles,
                             P);
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        e = 0;
        for (s = 0; s < n_seq; s++) {
          type    = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
          type_2  = vrna_get_ptype_md(SS[s][q], SS[s][p], md);  /* q,p not p,q! */
          e       += P->stack[type][type_2];
        }

        break;
    }

    /* add soft constraints */
    if (sc_wrapper.pair)
      e += sc_wrapper.pair(i, j, p, q, &sc_wrapper);
  }

  free_sc_int(&sc_wrapper);

  return e;
}
