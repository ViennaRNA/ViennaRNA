#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/exterior_loops.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "exterior_loops.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE FLT_OR_DBL
exp_E_ext_fast(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               vrna_mx_pf_aux_el_t  *aux_mx);


PRIVATE FLT_OR_DBL
exp_E_ext_fast_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j,
                      vrna_mx_pf_aux_el_t   *aux_mx);


PRIVATE FLT_OR_DBL
exp_E_ext_fast_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j,
                           vrna_mx_pf_aux_el_t  *aux_mx);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC FLT_OR_DBL
exp_E_Stem(int              type,
           int              si1,
           int              sj1,
           int              extLoop,
           vrna_exp_param_t *P)
{
  double  energy  = 1.0;
  double  d5      = (si1 >= 0) ? P->expdangle5[type][si1] : 1.;
  double  d3      = (sj1 >= 0) ? P->expdangle3[type][sj1] : 1.;

  if (si1 >= 0 && sj1 >= 0)
    energy = (extLoop) ? P->expmismatchExt[type][si1][sj1] : P->expmismatchM[type][si1][sj1];
  else
    energy = d5 * d3;

  if (type > 2)
    energy *= P->expTermAU;

  if (!extLoop)
    energy *= P->expMLintern[type];

  return (FLT_OR_DBL)energy;
}


PUBLIC FLT_OR_DBL
exp_E_ExtLoop(int               type,
              int               si1,
              int               sj1,
              vrna_exp_param_t  *P)
{
  double energy = 1.0;

  if (si1 >= 0 && sj1 >= 0)
    energy = P->expmismatchExt[type][si1][sj1];
  else if (si1 >= 0)
    energy = P->expdangle5[type][si1];
  else if (sj1 >= 0)
    energy = P->expdangle3[type][sj1];

  if (type > 2)
    energy *= P->expTermAU;

  return (FLT_OR_DBL)energy;
}


PUBLIC vrna_mx_pf_aux_el_t *
vrna_exp_E_ext_fast_init(vrna_fold_compound_t *vc)
{
  vrna_mx_pf_aux_el_t *aux_mx = NULL;

  if (vc) {
    unsigned int              u, s;
    int                       i, j, max_j, d, n, turn, ij, *idx, *iidx, *hc_up;
    FLT_OR_DBL                *q, **q_local, *scale;
    vrna_callback_hc_evaluate *evaluate;
    struct default_data       hc_dat_local;

    n     = (int)vc->length;
    idx   = vc->jindx;
    iidx  = vc->iindx;
    turn  = vc->exp_params->model_details.min_loop_size;
    scale = vc->exp_matrices->scale;
    hc_up = vc->hc->up_ext;

    if (vc->hc->type == VRNA_HC_WINDOW) {
      hc_dat_local.mx_window  = vc->hc->matrix_local;
      hc_dat_local.hc_up      = hc_up;
      hc_dat_local.sn         = vc->strand_number;

      if (vc->hc->f) {
        evaluate            = &hc_default_user_window;
        hc_dat_local.hc_f   = vc->hc->f;
        hc_dat_local.hc_dat = vc->hc->data;
      } else {
        evaluate = &hc_default_window;
      }
    } else {
      hc_dat_local.idx    = idx;
      hc_dat_local.mx     = vc->hc->matrix;
      hc_dat_local.hc_up  = hc_up;
      hc_dat_local.sn     = vc->strand_number;

      if (vc->hc->f) {
        evaluate            = &hc_default_user;
        hc_dat_local.hc_f   = vc->hc->f;
        hc_dat_local.hc_dat = vc->hc->data;
      } else {
        evaluate = &hc_default;
      }
    }

    /* allocate memory for helper arrays */
    aux_mx            = (vrna_mx_pf_aux_el_t *)vrna_alloc(sizeof(vrna_mx_pf_aux_el_t));
    aux_mx->qq        = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qq1       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqu_size  = 0;
    aux_mx->qqu       = NULL;

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      vrna_sc_t *sc         = vc->sc;
      vrna_ud_t *domains_up = vc->domains_up;
      int       with_ud     = (domains_up && domains_up->exp_energy_cb);

      /* pre-processing ligand binding production rule(s) and auxiliary memory */
      if (with_ud) {
        int ud_max_size = 0;
        for (u = 0; u < domains_up->uniq_motif_count; u++)
          if (ud_max_size < domains_up->uniq_motif_size[u])
            ud_max_size = domains_up->uniq_motif_size[u];

        aux_mx->qqu_size  = ud_max_size;
        aux_mx->qqu       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (ud_max_size + 1));

        for (u = 0; u <= ud_max_size; u++)
          aux_mx->qqu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
      }

      if (vc->hc->type == VRNA_HC_WINDOW) {
        q_local = vc->exp_matrices->q_local;
        max_j   = MIN2(turn + 1, vc->window_size);
        max_j   = MIN2(max_j, n);
        for (j = 1; j <= max_j; j++)
          for (i = 1; i <= j; i++) {
            if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
              q_local[i][j] = scale[(j - i + 1)];
              if (sc) {
                if (sc->exp_energy_up)
                  q_local[i][j] *= sc->exp_energy_up[i][j - i + 1];

                if (sc->exp_f)
                  q_local[i][j] *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
              }

              if (with_ud) {
                q_local[i][j] += q_local[i][j] *
                                 domains_up->exp_energy_cb(vc,
                                                           i, j,
                                                           VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                                           domains_up->data);
              }
            } else {
              q_local[i][j] = 0;
            }
          }
      } else {
        q = vc->exp_matrices->q;
        for (d = 0; d <= turn; d++)
          for (i = 1; i <= n - d; i++) {
            j   = i + d;
            ij  = iidx[i] - j;

            if (j > n)
              continue;

            if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
              q[ij] = scale[d + 1];

              if (sc) {
                if (sc->exp_energy_up)
                  q[ij] *= sc->exp_energy_up[i][d + 1];

                if (sc->exp_f)
                  q[ij] *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
              }

              if (with_ud) {
                q[ij] += q[ij] *
                         domains_up->exp_energy_cb(vc,
                                                   i, j,
                                                   VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                                   domains_up->data);
              }
            } else {
              q[ij] = 0.;
            }
          }
      }
    } else if (vc->type == VRNA_FC_TYPE_COMPARATIVE) {
      vrna_sc_t     **scs = vc->scs;
      unsigned int  **a2s = vc->a2s;
      q = vc->exp_matrices->q;
      for (d = 0; d <= turn; d++)
        for (i = 1; i <= n - d; i++) {
          j   = i + d;
          ij  = iidx[i] - j;
          if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
            q[ij] = scale[d + 1];

            if (scs) {
              for (s = 0; s < vc->n_seq; s++)
                if (scs[s]) {
                  u = d + 1 /* a2s[s][j] - a2s[s][i] + 1 */;
                  if (scs[s]->exp_energy_up)
                    q[ij] *= scs[s]->exp_energy_up[a2s[s][i]][u];
                }
            }
          } else {
            q[ij] = 0.;
          }
        }
    }
  }

  return aux_mx;
}


PUBLIC void
vrna_exp_E_ext_fast_rotate(vrna_fold_compound_t *vc,
                           vrna_mx_pf_aux_el_t  *aux_mx)
{
  if (vc && aux_mx) {
    int         u;
    FLT_OR_DBL  *tmp;

    tmp         = aux_mx->qq1;
    aux_mx->qq1 = aux_mx->qq;
    aux_mx->qq  = tmp;

    /* rotate auxiliary arrays for unstructured domains */
    if (aux_mx->qqu) {
      tmp = aux_mx->qqu[aux_mx->qqu_size];
      for (u = aux_mx->qqu_size; u > 0; u--)
        aux_mx->qqu[u] = aux_mx->qqu[u - 1];
      aux_mx->qqu[0] = tmp;
    }
  }
}


PUBLIC void
vrna_exp_E_ext_fast_free(vrna_fold_compound_t *vc,
                         vrna_mx_pf_aux_el_t  *aux_mx)
{
  if (vc && aux_mx) {
    int u;

    free(aux_mx->qq);
    free(aux_mx->qq1);

    if (aux_mx->qqu) {
      for (u = 0; u <= aux_mx->qqu_size; u++)
        free(aux_mx->qqu[u]);

      free(aux_mx->qqu);
    }

    free(aux_mx);
  }
}


PUBLIC FLT_OR_DBL
vrna_exp_E_ext_fast(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    vrna_mx_pf_aux_el_t   *aux_mx)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return exp_E_ext_fast_window(vc, i, j, aux_mx);
        else
          return exp_E_ext_fast(vc, i, j, aux_mx);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return exp_E_ext_fast_comparative(vc, i, j, aux_mx);
        break;

      default:
        vrna_message_warning("vrna_exp_E_ext_fast@exterior_loops.c: Unknown fold_compound type");
        return 0.;
        break;
    }
  } else {
    return 0.;
  }
}


PRIVATE FLT_OR_DBL
exp_E_ext_fast(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               vrna_mx_pf_aux_el_t  *aux_mx)
{
  short                     *S1, *S2;
  int                       n, *iidx, k, ij, kl, with_ud, u, circular, with_gquad, type;
  FLT_OR_DBL                qbt1, *q, *qb, *qq, *qq1, **qqu, q_temp, *scale, q_temp2, *G;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n                   = (int)vc->length;
  iidx                = vc->iindx;
  ij                  = iidx[i] - j;
  qq                  = aux_mx->qq;
  qq1                 = aux_mx->qq1;
  qqu                 = aux_mx->qqu;
  q                   = vc->exp_matrices->q;
  qb                  = vc->exp_matrices->qb;
  G                   = vc->exp_matrices->G;
  scale               = vc->exp_matrices->scale;
  pf_params           = vc->exp_params;
  md                  = &(pf_params->model_details);
  hc                  = vc->hc;
  sc                  = vc->sc;
  domains_up          = vc->domains_up;
  circular            = md->circ;
  with_gquad          = md->gquad;
  with_ud             = (domains_up && domains_up->exp_energy_cb);
  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ext;
  hc_dat_local.sn     = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  qbt1 = 0.;

  /* all exterior loop parts [i, j] with exactly one stem (i, u) i < u < j */
  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
    q_temp = qq1[i] * scale[1];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[j][1];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
    }

    if (with_ud) {
      int cnt;
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if (j - u >= i) {
          if (evaluate(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            q_temp2 = qqu[u][i] *
                      domains_up->exp_energy_cb(vc,
                                                j - u + 1,
                                                j,
                                                VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                domains_up->data) *
                      scale[u];

            if (sc) {
              if (sc->exp_energy_up)
                q_temp2 *= sc->exp_energy_up[j - u + 1][u];

              if (sc->exp_f)
                q_temp2 *= sc->exp_f(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
            }

            q_temp += q_temp2;
          }
        }
      }
    }

    qbt1 += q_temp;
  }

  /* exterior loop part with stem (i, j) */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
    S1      = vc->sequence_encoding;
    S2      = vc->sequence_encoding2;
    type    = get_pair_type_md(S2[i], S2[j], md);
    q_temp  = qb[ij] *
              exp_E_ExtLoop(type,
                            ((i > 1) || circular) ? S1[i - 1] : -1,
                            ((j < n) || circular) ? S1[j + 1] : -1,
                            pf_params);

    if (sc)
      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

    qbt1 += q_temp;
  }

  if (with_gquad)
    qbt1 += G[ij];

  qq[i] = qbt1;

  if (with_ud)
    qqu[0][i] = qbt1;

  /* the entire stretch [i,j] is unpaired */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
    u       = j - i + 1;
    q_temp  = scale[u];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[i][u];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
    }

    qbt1 += q_temp;

    if (with_ud) {
      qbt1 += q_temp *
              domains_up->exp_energy_cb(vc,
                                        i, j,
                                        VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                        domains_up->data);
    }
  }

  kl = iidx[i] - j + 1;
  if (sc && sc->exp_f) {
    for (k = j; k > i; k--, kl++) {
      q_temp = q[kl] *
               qq[k] *
               sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, sc->data);

      qbt1 += q_temp;
    }
  } else {
    for (k = j; k > i; k--, kl++)
      qbt1 += q[kl] *
              qq[k];
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_ext_fast_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j,
                      vrna_mx_pf_aux_el_t   *aux_mx)
{
  short                     *S1;
  int                       n, k, with_ud, u, circular, with_gquad, type, winSize, turn;
  FLT_OR_DBL                qbt1, **q, **qb, *qq, *qq1, **qqu, q_temp, *scale, q_temp2, **G;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n           = (int)vc->length;
  winSize     = vc->window_size;
  qq          = aux_mx->qq;
  qq1         = aux_mx->qq1;
  qqu         = aux_mx->qqu;
  q           = vc->exp_matrices->q_local;
  qb          = vc->exp_matrices->qb_local;
  G           = vc->exp_matrices->G_local;
  scale       = vc->exp_matrices->scale;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  hc          = vc->hc;
  sc          = vc->sc;
  domains_up  = vc->domains_up;
  circular    = md->circ;
  with_gquad  = md->gquad;
  turn        = md->min_loop_size;
  with_ud     = (domains_up && domains_up->exp_energy_cb);

  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ext;
  hc_dat_local.sn         = vc->strand_number;

  if (vc->hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  /*
   *  init exterior loop contributions for small segments [i, j]
   *  that can only be unpaired.
   *  We do this only once for the very first segment ending at j
   */
  if (i == j - turn - 1) {
    for (k = j; k >= MAX2(1, j - turn); k--) {
      if (evaluate(k, j, k, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
        q[k][j] = scale[(j - k + 1)];
        if (sc) {
          if (sc->exp_energy_up)
            q[k][j] *= sc->exp_energy_up[k][j - k + 1];

          if (sc->exp_f)
            q[k][j] *= sc->exp_f(k, j, k, j, VRNA_DECOMP_EXT_UP, sc->data);
        }

        if (with_ud) {
          q[k][j] += q[k][j] *
                     domains_up->exp_energy_cb(vc,
                                               k, j,
                                               VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                               domains_up->data);
        }
      } else {
        q[k][j] = 0;
      }
    }
  }

  qbt1 = 0.;

  /* all exterior loop parts [i, j] with exactly one stem (i, u) i < u < j */
  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
    q_temp = qq1[i] * scale[1];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[j][1];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
    }

    if (with_ud) {
      int cnt;
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if (j - u >= i) {
          if (evaluate(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            q_temp2 = qqu[u][i] *
                      domains_up->exp_energy_cb(vc,
                                                j - u + 1,
                                                j,
                                                VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                domains_up->data) *
                      scale[u];

            if (sc) {
              if (sc->exp_energy_up)
                q_temp2 *= sc->exp_energy_up[j - u + 1][u];

              if (sc->exp_f)
                q_temp2 *= sc->exp_f(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
            }

            q_temp += q_temp2;
          }
        }
      }
    }

    qbt1 += q_temp;
  }

  /* exterior loop part with stem (i, j) */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
    S1      = vc->sequence_encoding;
    type    = get_pair_type_md(S1[i], S1[j], md);
    q_temp  = qb[i][j] *
              exp_E_ExtLoop(type,
                            ((i > 1) || circular) ? S1[i - 1] : -1,
                            ((j < n) || circular) ? S1[j + 1] : -1,
                            pf_params);

    if (sc)
      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

    qbt1 += q_temp;
  }

  if (with_gquad)
    qbt1 += G[i][j];

  qq[i] = qbt1;

  if (with_ud)
    qqu[0][i] = qbt1;

  /* the entire stretch [i,j] is unpaired */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
    u       = j - i + 1;
    q_temp  = scale[u];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[i][u];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
    }

    qbt1 += q_temp;

    if (with_ud) {
      qbt1 += q_temp *
              domains_up->exp_energy_cb(vc,
                                        i, j,
                                        VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                        domains_up->data);
    }
  }

  if (sc && sc->exp_f) {
    for (k = j; k > i; k--) {
      q_temp = q[i][k - 1] *
               qq[k] *
               sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, sc->data);

      qbt1 += q_temp;
    }
  } else {
    for (k = j; k > i; k--)
      qbt1 += q[i][k - 1] *
              qq[k];
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_ext_fast_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j,
                           vrna_mx_pf_aux_el_t  *aux_mx)
{
  int                       n, s, n_seq, *iidx, k, ij, kl, u, circular, type;
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  FLT_OR_DBL                qbt1, *q, *qb, *qq, *qq1, q_temp, *scale;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n                   = (int)vc->length;
  n_seq               = vc->n_seq;
  iidx                = vc->iindx;
  ij                  = iidx[i] - j;
  S                   = vc->S;
  S5                  = vc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
  S3                  = vc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
  a2s                 = vc->a2s;
  qq                  = aux_mx->qq;
  qq1                 = aux_mx->qq1;
  q                   = vc->exp_matrices->q;
  qb                  = vc->exp_matrices->qb;
  scale               = vc->exp_matrices->scale;
  pf_params           = vc->exp_params;
  md                  = &(pf_params->model_details);
  hc                  = vc->hc;
  scs                 = vc->scs;
  circular            = md->circ;
  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ext;
  hc_dat_local.sn     = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  qbt1 = 0.;

  /* all exterior loop parts [i, j] with exactly one stem (i, u) i < u < j */
  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
    q_temp = qq1[i] *
             scale[1];

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->exp_energy_up)
            q_temp *= scs[s]->exp_energy_up[a2s[s][j]][1];
      }
    }

    qbt1 += q_temp;
  }

  /* exterior loop part with stem (i, j) */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
    q_temp = qb[ij];

    for (s = 0; s < n_seq; s++) {
      type    = get_pair_type_md(S[s][i], S[s][j], md);
      q_temp  *= exp_E_ExtLoop(type,
                               ((i > 1) || circular) ? S5[s][i] : -1,
                               ((j < n) || circular) ? S3[s][j] : -1,
                               pf_params);
    }

    qbt1 += q_temp;
  }

  qq[i] = qbt1;

  /* the entire stretch [i,j] is unpaired */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
    u       = j - i + 1;
    q_temp  = scale[u];

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->exp_energy_up)
            q_temp *= scs[s]->exp_energy_up[a2s[s][i]][a2s[s][j] - a2s[s][i] + 1];
      }
    }

    qbt1 += q_temp;
  }

  kl = iidx[i] - j + 1;
  for (k = j; k > i; k--, kl++)
    qbt1 += q[kl] *
            qq[k];

  return qbt1;
}
