#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/exterior_loops.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/interior_loops.h"


#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "interior_loops.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE FLT_OR_DBL
exp_E_int_loop(vrna_fold_compound_t *vc,
               int                  i,
               int                  j);


PRIVATE FLT_OR_DBL
exp_E_int_loop_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j);


PRIVATE FLT_OR_DBL
exp_E_int_loop_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j);


PRIVATE FLT_OR_DBL
exp_E_interior_loop(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC FLT_OR_DBL
vrna_exp_E_int_loop(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j)
{
  FLT_OR_DBL q = 0.;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          q = exp_E_int_loop_window(vc, i, j);
        else
          q = exp_E_int_loop(vc, i, j);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        q = exp_E_int_loop_comparative(vc, i, j);
        break;
    }
  }

  return q;
}


PRIVATE FLT_OR_DBL
exp_E_int_loop(vrna_fold_compound_t *vc,
               int                  i,
               int                  j)
{
  unsigned char             type, type_2;
  char                      *ptype;
  unsigned char             *hc, eval_loop;
  short                     *S1, S_i1, S_j1;
  unsigned int              *sn;
  int                       k, l, u1, u2, kl, maxk, minl, *rtype, noGUclosure,
                            no_close, *my_iindx, *jindx, *hc_up, ij,
                            with_gquad, turn;
  FLT_OR_DBL                qbt1, q_temp, *qb, *G, *scale;
  vrna_sc_t                 *sc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct  default_data      hc_dat_local;

  ptype       = vc->ptype;
  S1          = vc->sequence_encoding;
  S_i1        = S1[i + 1];
  S_j1        = S1[j - 1];
  my_iindx    = vc->iindx;
  jindx       = vc->jindx;
  hc          = vc->hc->matrix;
  hc_up       = vc->hc->up_int;
  sc          = vc->sc;
  sn          = vc->strand_number;
  pf_params   = vc->exp_params;
  ij          = jindx[j] + i;
  md          = &(pf_params->model_details);
  with_gquad  = md->gquad;
  turn        = md->min_loop_size;
  qb          = vc->exp_matrices->qb;
  G           = vc->exp_matrices->G;
  scale       = vc->exp_matrices->scale;
  domains_up  = vc->domains_up;
  qbt1        = 0.;
  evaluate    = prepare_hc_default(vc, &hc_dat_local);

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    type        = (unsigned char)ptype[ij];
    rtype       = &(md->rtype[0]);
    noGUclosure = md->noGUclosure;
    no_close    = (((type == 3) || (type == 4)) && noGUclosure);
    maxk        = i + MAXLOOP + 1;
    maxk        = MIN2(maxk, j - turn - 2);
    maxk        = MIN2(maxk, i + 1 + hc_up[i + 1]);

    if (type == 0)
      type = 7;

    for (k = i + 1; k <= maxk; k++) {
      if (sn[k] != sn[i])
        break;

      u1 = k - i - 1;

      minl  = MAX2(k + turn + 1, j - 1 - MAXLOOP + u1);
      kl    = my_iindx[k] - j + 1;

      for (u2 = 0, l = j - 1; l >= minl; l--, kl++, u2++) {
        if (hc_up[l + 1] < u2)
          break;

        eval_loop =
          (hc[jindx[l] + k] &
           VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) ? (unsigned char)1 : (unsigned char)0;

        /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
        if (eval_loop && evaluate(i, j, k, l, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
          if (sn[j] != sn[l])
            break;

          type_2 = rtype[(unsigned char)ptype[jindx[l] + k]];

          if (type_2 == 0)
            type_2 = 7;

          q_temp = qb[kl] *
                   scale[u1 + u2 + 2] *
                   exp_E_IntLoop(u1, u2, type, type_2, S_i1, S_j1, S1[k - 1], S1[l + 1], pf_params);

          /* soft constraints */
          if (sc) {
            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i + 1][u1] *
                        sc->exp_energy_up[l + 1][u2];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);

            if (sc->exp_energy_stack) {
              if ((i + 1 == k) && (j - 1 == l)) {
                q_temp *= sc->exp_energy_stack[i] *
                          sc->exp_energy_stack[k] *
                          sc->exp_energy_stack[l] *
                          sc->exp_energy_stack[j];
              }
            }
          }

          qbt1 += q_temp;

          /* unstructured domains */
          if (domains_up && domains_up->exp_energy_cb) {
            FLT_OR_DBL qq5, qq3;

            qq5 = qq3 = 0.;

            if (u1 > 0) {
              qq5 = domains_up->exp_energy_cb(vc,
                                              i + 1, k - 1,
                                              VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                              domains_up->data);
            }

            if (u2 > 0) {
              qq3 = domains_up->exp_energy_cb(vc,
                                              l + 1, j - 1,
                                              VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                              domains_up->data);
            }

            qbt1  += q_temp * qq5;        /* only motifs in 5' part */
            qbt1  += q_temp * qq3;        /* only motifs in 3' part */
            qbt1  += q_temp * qq5 * qq3;  /* motifs in both parts */
          }
        }
      }
    }

    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      if ((!no_close) && (sn[j] == sn[i]))
        qbt1 += exp_E_GQuad_IntLoop(i, j, type, S1, G, scale, my_iindx, pf_params);
    }

    if (sc && sc->exp_energy_bp)
      qbt1 *= sc->exp_energy_bp[my_iindx[i] - j];
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_int_loop_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j)
{
  unsigned char             type, type_2;
  char                      **ptype;
  unsigned char             **hc, eval_loop;
  short                     *S1, S_i1, S_j1;
  unsigned int              *sn;
  int                       k, l, u1, u2, maxk, minl, *rtype, noGUclosure,
                            no_close, *my_iindx, *jindx, *hc_up,
                            with_gquad, turn;
  FLT_OR_DBL                qbt1, q_temp, **qb, **G, *scale;
  vrna_sc_t                 *sc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct  default_data      hc_dat_local;

  ptype       = vc->ptype_local;
  S1          = vc->sequence_encoding;
  S_i1        = S1[i + 1];
  S_j1        = S1[j - 1];
  jindx       = vc->jindx;
  hc          = vc->hc->matrix_local;
  hc_up       = vc->hc->up_int;
  sc          = vc->sc;
  sn          = vc->strand_number;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  with_gquad  = md->gquad;
  turn        = md->min_loop_size;
  qb          = vc->exp_matrices->qb_local;
  G           = vc->exp_matrices->G_local;
  scale       = vc->exp_matrices->scale;
  domains_up  = vc->domains_up;
  qbt1        = 0.;
  evaluate    = prepare_hc_default(vc, &hc_dat_local);

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc[i][j - i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    type        = ptype[i][j];
    rtype       = &(md->rtype[0]);
    noGUclosure = md->noGUclosure;
    no_close    = (((type == 3) || (type == 4)) && noGUclosure);
    maxk        = i + MAXLOOP + 1;
    maxk        = MIN2(maxk, j - turn - 2);
    maxk        = MIN2(maxk, i + 1 + hc_up[i + 1]);

    if (type == 0)
      type = 7;

    for (k = i + 1; k <= maxk; k++) {
      if (sn[k] != sn[i])
        break;

      u1 = k - i - 1;

      minl = MAX2(k + turn + 1, j - 1 - MAXLOOP + u1);

      for (u2 = 0, l = j - 1; l >= minl; l--, u2++) {
        if (hc_up[l + 1] < u2)
          break;

        eval_loop =
          (hc[k][l - k] &
           VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) ? (unsigned char)1 : (unsigned char)0;

        /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
        if (eval_loop && evaluate(i, j, k, l, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
          if (sn[j] != sn[l])
            break;

          type_2 = rtype[(unsigned char)ptype[k][l]];

          if (type_2 == 0)
            type_2 = 7;

          q_temp = qb[k][l] *
                   scale[u1 + u2 + 2] *
                   exp_E_IntLoop(u1, u2, type, type_2, S_i1, S_j1, S1[k - 1], S1[l + 1], pf_params);

          /* soft constraints */
          if (sc) {
            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i + 1][u1] *
                        sc->exp_energy_up[l + 1][u2];

            if (sc->exp_energy_bp_local)
              q_temp *= sc->exp_energy_bp_local[i][j - i];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);

            if (sc->exp_energy_stack) {
              if ((u1 == 0) && (u2 == 0)) {
                q_temp *= sc->exp_energy_stack[i] *
                          sc->exp_energy_stack[k] *
                          sc->exp_energy_stack[l] *
                          sc->exp_energy_stack[j];
              }
            }
          }

          qbt1 += q_temp;

          /* unstructured domains */
          if (domains_up && domains_up->exp_energy_cb) {
            FLT_OR_DBL qq5, qq3;

            qq5 = qq3 = 0.;

            if (u1 > 0) {
              qq5 = domains_up->exp_energy_cb(vc,
                                              i + 1, k - 1,
                                              VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                              domains_up->data);
            }

            if (u2 > 0) {
              qq3 = domains_up->exp_energy_cb(vc,
                                              l + 1, j - 1,
                                              VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                              domains_up->data);
            }

            qbt1  += q_temp * qq5;        /* only motifs in 5' part */
            qbt1  += q_temp * qq3;        /* only motifs in 3' part */
            qbt1  += q_temp * qq5 * qq3;  /* motifs in both parts */
          }
        }
      }
    }

#if 0
    /* no G-Quadruplexes for sliding-window partition function yet! */
    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      if ((!no_close) && (sn[j] == sn[i]))
        qbt1 += exp_E_GQuad_IntLoop(i, j, type, S1, G, scale, my_iindx, pf_params);
    }

#endif
  }

  return qbt1;
}


PUBLIC FLT_OR_DBL
vrna_exp_E_interior_loop(vrna_fold_compound_t *vc,
                         int                  i,
                         int                  j,
                         int                  k,
                         int                  l)
{
  if (vc)
    return exp_E_interior_loop(vc, i, j, k, l);

  return 0.;
}


PRIVATE FLT_OR_DBL
exp_E_interior_loop(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l)
{
  unsigned char             type, type_2;
  char                      *ptype;
  unsigned char             *hc, eval_loop;
  short                     *S1, S_i1, S_j1;
  unsigned int              *sn;
  int                       u1, u2, *rtype, *my_iindx, *jindx, *hc_up, ij;
  FLT_OR_DBL                qbt1, q_temp, *scale;
  vrna_sc_t                 *sc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  ptype       = vc->ptype;
  S1          = vc->sequence_encoding;
  S_i1        = S1[i + 1];
  S_j1        = S1[j - 1];
  my_iindx    = vc->iindx;
  jindx       = vc->jindx;
  hc          = vc->hc->matrix;
  hc_up       = vc->hc->up_int;
  sc          = vc->sc;
  pf_params   = vc->exp_params;
  ij          = jindx[j] + i;
  sn          = vc->strand_number;
  md          = &(pf_params->model_details);
  scale       = vc->exp_matrices->scale;
  domains_up  = vc->domains_up;
  qbt1        = 0.;
  u1          = k - i - 1;
  u2          = j - l - 1;

  if ((sn[k] != sn[i]) || (sn[j] != sn[l]))
    return qbt1;

  if (hc_up[l + 1] < u2)
    return qbt1;

  if (hc_up[i + 1] < u1)
    return qbt1;

  evaluate = prepare_hc_default(vc, &hc_dat_local);

  /* CONSTRAINED INTERIOR LOOP start */
  eval_loop =
    ((hc[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) &&
     (hc[jindx[l] + k] &
      VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) ? (unsigned char)1 : (unsigned char)0;

  /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
  if (eval_loop && evaluate(i, j, k, l, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
    type    = (unsigned char)ptype[ij];
    rtype   = &(md->rtype[0]);
    type    = (unsigned char)ptype[ij];
    type_2  = rtype[(unsigned char)ptype[jindx[l] + k]];

    if (type == 0)
      type = 7;

    if (type_2 == 0)
      type_2 = 7;

    q_temp = exp_E_IntLoop(u1, u2, type, type_2, S_i1, S_j1, S1[k - 1], S1[l + 1], pf_params) *
             scale[u1 + u2 + 2];

    /* soft constraints */
    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[i + 1][u1] *
                  sc->exp_energy_up[l + 1][u2];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);

      if (sc->exp_energy_stack) {
        if ((i + 1 == k) && (j - 1 == l)) {
          q_temp *= sc->exp_energy_stack[i] *
                    sc->exp_energy_stack[k] *
                    sc->exp_energy_stack[l] *
                    sc->exp_energy_stack[j];
        }
      }

      if (sc->exp_energy_bp)
        q_temp *= sc->exp_energy_bp[my_iindx[i] - j];
    }

    qbt1 += q_temp;

    /* unstructured domains */
    if (domains_up && domains_up->exp_energy_cb) {
      FLT_OR_DBL qq5, qq3;

      qq5 = qq3 = 0.;

      if (u1 > 0) {
        qq5 = domains_up->exp_energy_cb(vc,
                                        i + 1, k - 1,
                                        VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                        domains_up->data);
      }

      if (u2 > 0) {
        qq3 = domains_up->exp_energy_cb(vc,
                                        l + 1, j - 1,
                                        VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                        domains_up->data);
      }

      qbt1  += q_temp * qq5;        /* only motifs in 5' part */
      qbt1  += q_temp * qq3;        /* only motigs in 3' part */
      qbt1  += q_temp * qq5 * qq3;  /* motifs in both parts */
    }
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_int_loop_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j)
{
  unsigned char             type_2;
  unsigned char             *hc, eval_loop;
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  int                       n_seq, s, ij, jij, k, l, u1, u2, kl, maxk, minl, *types,
                            turn, with_gquad, *hc_up, *jindx, *my_iindx;
  FLT_OR_DBL                qbt1, *qb, *scale, qloop;
  vrna_sc_t                 **scs;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  types       = NULL;
  my_iindx    = vc->iindx;
  jindx       = vc->jindx;
  hc          = vc->hc->matrix;
  hc_up       = vc->hc->up_int;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  with_gquad  = md->gquad;
  turn        = md->min_loop_size;
  qb          = vc->exp_matrices->qb;
  scale       = vc->exp_matrices->scale;
  qbt1        = 0.;
  jij         = jindx[j] + i;
  ij          = my_iindx[i] - j;
  evaluate    = prepare_hc_default(vc, &hc_dat_local);

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc[jij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    S     = vc->S;
    S5    = vc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
    S3    = vc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
    a2s   = vc->a2s;
    scs   = vc->scs;
    n_seq = vc->n_seq;
    types = (int *)vrna_alloc(sizeof(int) * n_seq);

    for (s = 0; s < n_seq; s++)
      types[s] = get_pair_type(S[s][i], S[s][j], md);

    /* prepare necessary variables */
    maxk  = i + MAXLOOP + 1;
    maxk  = MIN2(maxk, j - turn - 2);
    maxk  = MIN2(maxk, i + 1 + hc_up[i + 1]);

    for (k = i + 1; k <= maxk; k++) {
      u1 = k - i - 1;

      minl  = MAX2(k + turn + 1, j - 1 - MAXLOOP + u1);
      kl    = my_iindx[k] - j + 1;

      for (l = j - 1; l >= minl; l--, kl++, u2++) {
        if (hc_up[l + 1] < j - l - 1)
          break;

        eval_loop =
          (hc[jindx[l] + k] &
           VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) ? (unsigned char)1 : (unsigned char)0;

        /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
        if (eval_loop && evaluate(i, j, k, l, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
          qloop = 1.;

          for (s = 0; s < n_seq; s++) {
            u1      = a2s[s][k - 1] - a2s[s][i];
            u2      = a2s[s][j - 1] - a2s[s][l];
            type_2  = get_pair_type(S[s][l], S[s][k], md);

            qloop *= exp_E_IntLoop(u1, u2,
                                   types[s], type_2, S3[s][i],
                                   S5[s][j], S5[s][k], S3[s][l],
                                   pf_params
                                   );
          }

          if (scs) {
            for (s = 0; s < n_seq; s++) {
              if (scs[s]) {
                u1  = a2s[s][k - 1] - a2s[s][i];
                u2  = a2s[s][j - 1] - a2s[s][l];

                if (scs[s]->exp_energy_up)
                  qloop *= scs[s]->exp_energy_up[a2s[s][i] + 1][u1] *
                           scs[s]->exp_energy_up[a2s[s][l] + 1][u2];

                if (scs[s]->exp_energy_stack) {
                  if (u1 + u2 == 0) {
                    if (S[s][i] && S[s][j] && S[s][k] && S[s][l]) {
                      /* don't allow gaps in stack */
                      qloop *= scs[s]->exp_energy_stack[i] *
                               scs[s]->exp_energy_stack[k] *
                               scs[s]->exp_energy_stack[l] *
                               scs[s]->exp_energy_stack[j];
                    }
                  }
                }
              }
            }
          }

          qbt1 += qb[my_iindx[k] - l] *
                  qloop *
                  scale[k - i + j - l];
        }
      }
    }

    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      /* not implemented yet! */
    }

    if (scs) {
      for (s = 0; s < n_seq; s++)
        if (scs[s] && scs[s]->exp_energy_bp)
          qbt1 *= scs[s]->exp_energy_bp[ij];
    }
  }

  /* cleanup */
  free(types);

  return qbt1;
}
