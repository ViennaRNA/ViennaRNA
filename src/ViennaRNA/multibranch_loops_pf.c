#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/exterior_loops.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/multibranch_loops.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "multibranch_loops.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   FLT_OR_DBL           *qqm1);


PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast_window(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          FLT_OR_DBL            *qqm1);


PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast_comparative(vrna_fold_compound_t *vc,
                               int                  i,
                               int                  j,
                               FLT_OR_DBL           *qqm1);


PRIVATE FLT_OR_DBL
exp_E_ml_fast(vrna_fold_compound_t  *vc,
              int                   i,
              int                   j,
              vrna_mx_pf_aux_ml_t   *aux_mx);


PRIVATE FLT_OR_DBL
exp_E_ml_fast_window(vrna_fold_compound_t *vc,
                     int                  i,
                     int                  j,
                     vrna_mx_pf_aux_ml_t  *aux_mx);


PRIVATE FLT_OR_DBL
exp_E_ml_fast_comparative(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          vrna_mx_pf_aux_ml_t   *aux_mx);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC FLT_OR_DBL
vrna_exp_E_mb_loop_fast(vrna_fold_compound_t  *vc,
                        int                   i,
                        int                   j,
                        FLT_OR_DBL            *qqm1)
{
  FLT_OR_DBL q = 0.;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          q = exp_E_mb_loop_fast_window(vc, i, j, qqm1);
        else
          q = exp_E_mb_loop_fast(vc, i, j, qqm1);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        q = exp_E_mb_loop_fast_comparative(vc, i, j, qqm1);
        break;
    }
  }

  return q;
}


PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   FLT_OR_DBL           *qqm1)
{
  char                      *ptype;
  short                     *S1;
  unsigned int              *sn;
  int                       ij, k, kl, *my_iindx, *jindx, *rtype, tt;
  FLT_OR_DBL                qbt1, temp, qqqmmm, *qm, *scale, expMLclosing;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  sc            = vc->sc;
  ptype         = vc->ptype;
  S1            = vc->sequence_encoding;
  qm            = vc->exp_matrices->qm;
  scale         = vc->exp_matrices->scale;
  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  ij            = jindx[j] + i;
  sn            = vc->strand_number;
  hc            = vc->hc;
  expMLclosing  = pf_params->expMLclosing;
  qbt1          = 0.;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.sn     = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* multiple stem loop contribution */
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML,
               &hc_dat_local) && (sn[i] == sn[i + 1]) && (sn[j - 1] == sn[j])) {
    rtype = &(md->rtype[0]);
    tt    = rtype[get_pair_type(ij, ptype)];

    qqqmmm = expMLclosing *
             exp_E_MLstem(tt, S1[j - 1], S1[i + 1], pf_params) *
             scale[2];

    temp  = 0.0;
    kl    = my_iindx[i + 1] - (i + 1);

    if (sc) {
      if (sc->exp_energy_bp)
        qqqmmm *= sc->exp_energy_bp[my_iindx[i] - j];

      if (sc->exp_f) {
        qqqmmm *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_ML, sc->data);

        if (hc->f) {
          for (k = i + 2; k <= j - 1; k++, kl--) {
            if ((sn[k - 1] == sn[k]) &&
                (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))) {
              temp += qm[kl] *
                      qqm1[k] *
                      sc->exp_f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
            }
          }
        } else {
          for (k = i + 2; k <= j - 1; k++, kl--) {
            if (sn[k - 1] == sn[k]) {
              temp += qm[kl] *
                      qqm1[k] *
                      sc->exp_f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
            }
          }
        }
      } else {
        if (hc->f) {
          for (k = i + 2; k <= j - 1; k++, kl--) {
            if ((sn[k - 1] == sn[k]) &&
                (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)))
              temp += qm[kl] *
                      qqm1[k];
          }
        } else {
          for (k = i + 2; k <= j - 1; k++, kl--) {
            if (sn[k - 1] == sn[k])
              temp += qm[kl] *
                      qqm1[k];
          }
        }
      }
    } else {
      if (hc->f) {
        for (k = i + 2; k <= j - 1; k++, kl--) {
          if ((sn[k - 1] == sn[k]) &&
              (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)))
            temp += qm[kl] *
                    qqm1[k];
        }
      } else {
        for (k = i + 2; k <= j - 1; k++, kl--) {
          if (sn[k - 1] == sn[k])
            temp += qm[kl] *
                    qqm1[k];
        }
      }
    }

    qbt1 += temp * qqqmmm;
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast_window(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          FLT_OR_DBL            *qqm1)
{
  char                      **ptype;
  short                     *S1;
  unsigned int              *sn;
  int                       ij, k, *rtype, tt;
  FLT_OR_DBL                qbt1, temp, qqqmmm, **qm, *scale, expMLclosing;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  sc            = vc->sc;
  ptype         = vc->ptype_local;
  S1            = vc->sequence_encoding;
  qm            = vc->exp_matrices->qm_local;
  scale         = vc->exp_matrices->scale;
  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  sn            = vc->strand_number;
  hc            = vc->hc;
  expMLclosing  = pf_params->expMLclosing;
  qbt1          = 0.;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.sn         = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  /* multiple stem loop contribution */
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML,
               &hc_dat_local) && (sn[i] == sn[i + 1]) && (sn[j - 1] == sn[j])) {
    rtype = &(md->rtype[0]);
    tt    = rtype[get_pair_type_window(i, j + i, ptype)];

    qqqmmm = expMLclosing *
             exp_E_MLstem(tt, S1[j - 1], S1[i + 1], pf_params) *
             scale[2];

    temp = 0.0;

    if (sc) {
      if (sc->exp_energy_bp_local)
        qqqmmm *= sc->exp_energy_bp_local[i][j - i];

      if (sc->exp_f) {
        qqqmmm *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_ML, sc->data);

        if (hc->f) {
          for (k = i + 2; k <= j - 1; k++) {
            if ((sn[k - 1] == sn[k]) &&
                (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))) {
              temp += qm[i + 1][k - 1] *
                      qqm1[k] *
                      sc->exp_f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
            }
          }
        } else {
          for (k = i + 2; k <= j - 1; k++) {
            if (sn[k - 1] == sn[k]) {
              temp += qm[i + 1][k - 1] *
                      qqm1[k] *
                      sc->exp_f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
            }
          }
        }
      } else {
        if (hc->f) {
          for (k = i + 2; k <= j - 1; k++) {
            if ((sn[k - 1] == sn[k]) &&
                (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)))
              temp += qm[i + 1][k - 1] *
                      qqm1[k];
          }
        } else {
          for (k = i + 2; k <= j - 1; k++) {
            if (sn[k - 1] == sn[k])
              temp += qm[i + 1][k - 1] *
                      qqm1[k];
          }
        }
      }
    } else {
      if (hc->f) {
        for (k = i + 2; k <= j - 1; k++) {
          if ((sn[k - 1] == sn[k]) &&
              (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)))
            temp += qm[i + 1][k - 1] *
                    qqm1[k];
        }
      } else {
        for (k = i + 2; k <= j - 1; k++) {
          if (sn[k - 1] == sn[k])
            temp += qm[i + 1][k - 1] *
                    qqm1[k];
        }
      }
    }

    qbt1 += temp * qqqmmm;
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast_comparative(vrna_fold_compound_t *vc,
                               int                  i,
                               int                  j,
                               FLT_OR_DBL           *qqm1)
{
  short                     **S, **S5, **S3;
  int                       k, kl, *my_iindx, tt, n_seq, s;
  FLT_OR_DBL                qbt1, temp, qqqmmm, *qm, *scale, expMLclosing;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  my_iindx      = vc->iindx;
  qm            = vc->exp_matrices->qm;
  scale         = vc->exp_matrices->scale;
  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  hc            = vc->hc;
  expMLclosing  = pf_params->expMLclosing;
  qbt1          = 0.;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.sn     = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* multiple stem loop contribution */
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    S     = vc->S;
    S5    = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
    S3    = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */
    scs   = vc->scs;
    n_seq = vc->n_seq;

    qqqmmm = 1.;

    for (s = 0; s < n_seq; s++) {
      tt      = get_pair_type_md(S[s][j], S[s][i], md);
      qqqmmm  *= exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params) *
                 expMLclosing;
    }

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s]) {
          if (scs[s]->exp_energy_bp)
            qqqmmm *= scs[s]->exp_energy_bp[my_iindx[i] - j];

          if (scs[s]->f)
            qqqmmm *= scs[s]->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, scs[s]->data);
        }
      }
    }

    /* multi-loop loop contribution */
    temp  = 0.;
    kl    = my_iindx[i + 1] - (i + 1);

    if (hc->f) {
      if (scs) {
        for (k = i + 2; k <= j - 1; k++, kl--) {
          if (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)) {
            qbt1 = qm[kl] * qqm1[k];
            for (s = 0; s < n_seq; s++)
              if (scs[s] && scs[s]->f)
                qbt1 *= scs[s]->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data);

            temp += qbt1;
          }
        }
      } else {
        for (k = i + 2; k <= j - 1; k++, kl--)
          if (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))
            temp += qm[kl] * qqm1[k];
      }
    } else {
      if (scs) {
        for (k = i + 2; k <= j - 1; k++, kl--) {
          qbt1 = qm[kl] * qqm1[k];
          for (s = 0; s < n_seq; s++)
            if (scs[s] && scs[s]->f)
              qbt1 *= scs[s]->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data);

          temp += qbt1;
        }
      } else {
        for (k = i + 2; k <= j - 1; k++, kl--)
          temp += qm[kl] * qqm1[k];
      }
    }

    temp *= scale[2];

    qbt1 = temp * qqqmmm;
  }

  return qbt1;
}


PUBLIC vrna_mx_pf_aux_ml_t *
vrna_exp_E_ml_fast_init(vrna_fold_compound_t *vc)
{
  vrna_mx_pf_aux_ml_t *aux_mx = NULL;

  if (vc) {
    int         i, j, d, n, u, turn, ij, *iidx;
    FLT_OR_DBL  *qm;

    n     = (int)vc->length;
    iidx  = vc->iindx;
    turn  = vc->exp_params->model_details.min_loop_size;
    qm    = vc->exp_matrices->qm;

    /* allocate memory for helper arrays */
    aux_mx            = (vrna_mx_pf_aux_ml_t *)vrna_alloc(sizeof(vrna_mx_pf_aux_ml_t));
    aux_mx->qqm       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqm1      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqmu_size = 0;
    aux_mx->qqmu      = NULL;

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      vrna_ud_t *domains_up = vc->domains_up;
      int       with_ud     = (domains_up && domains_up->exp_energy_cb);
      int       ud_max_size = 0;

      /* pre-processing ligand binding production rule(s) and auxiliary memory */
      if (with_ud) {
        for (u = 0; u < domains_up->uniq_motif_count; u++)
          if (ud_max_size < domains_up->uniq_motif_size[u])
            ud_max_size = domains_up->uniq_motif_size[u];

        aux_mx->qqmu_size = ud_max_size;
        aux_mx->qqmu      = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (ud_max_size + 1));
        for (u = 0; u <= ud_max_size; u++)
          aux_mx->qqmu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
      }
    }

    if (vc->hc->type == VRNA_HC_WINDOW) {
    } else {
      for (d = 0; d <= turn; d++)
        for (i = 1; i <= n - d; i++) {
          j   = i + d;
          ij  = iidx[i] - j;

          if (j > n)
            continue;

          qm[ij] = 0.;
        }
    }
  }

  return aux_mx;
}


PUBLIC void
vrna_exp_E_ml_fast_rotate(vrna_fold_compound_t  *vc,
                          vrna_mx_pf_aux_ml_t   *aux_mx)
{
  if (vc && aux_mx) {
    int         u;
    FLT_OR_DBL  *tmp;

    tmp           = aux_mx->qqm1;
    aux_mx->qqm1  = aux_mx->qqm;
    aux_mx->qqm   = tmp;

    /* rotate auxiliary arrays for unstructured domains */
    if (aux_mx->qqmu) {
      tmp = aux_mx->qqmu[aux_mx->qqmu_size];
      for (u = aux_mx->qqmu_size; u > 0; u--)
        aux_mx->qqmu[u] = aux_mx->qqmu[u - 1];
      aux_mx->qqmu[0] = tmp;
    }
  }
}


PUBLIC void
vrna_exp_E_ml_fast_free(vrna_fold_compound_t  *vc,
                        vrna_mx_pf_aux_ml_t   *aux_mx)
{
  if (vc && aux_mx) {
    int u;

    free(aux_mx->qqm);
    free(aux_mx->qqm1);

    if (aux_mx->qqmu) {
      for (u = 0; u <= aux_mx->qqmu_size; u++)
        free(aux_mx->qqmu[u]);

      free(aux_mx->qqmu);
    }

    free(aux_mx);
  }
}


PUBLIC FLT_OR_DBL
vrna_exp_E_ml_fast(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   vrna_mx_pf_aux_ml_t  *aux_mx)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return exp_E_ml_fast_window(vc, i, j, aux_mx);
        else
          return exp_E_ml_fast(vc, i, j, aux_mx);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return exp_E_ml_fast_comparative(vc, i, j, aux_mx);
        break;

      default:
        vrna_message_warning("vrna_exp_E_ml_fast@multibranch_loops.c: Unknown fold_compound type");
        return 0.;
        break;
    }
  } else {
    return 0.;
  }
}


PRIVATE FLT_OR_DBL
exp_E_ml_fast(vrna_fold_compound_t  *vc,
              int                   i,
              int                   j,
              vrna_mx_pf_aux_ml_t   *aux_mx)
{
  short                     *S1, *S2;
  int                       n, *iidx, k, ij, kl, maxk, ii, with_ud, u, circular, with_gquad,
                            *hc_up_ml, type;
  FLT_OR_DBL                qbt1, temp, *qm, *qb, *qqm, *qqm1, **qqmu, q_temp, q_temp2, *G,
                            *expMLbase,
                            expMLstem;
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
  qqm                 = aux_mx->qqm;
  qqm1                = aux_mx->qqm1;
  qqmu                = aux_mx->qqmu;
  qm                  = vc->exp_matrices->qm;
  qb                  = vc->exp_matrices->qb;
  G                   = vc->exp_matrices->G;
  expMLbase           = vc->exp_matrices->expMLbase;
  pf_params           = vc->exp_params;
  md                  = &(pf_params->model_details);
  hc                  = vc->hc;
  sc                  = vc->sc;
  domains_up          = vc->domains_up;
  circular            = md->circ;
  with_gquad          = md->gquad;
  with_ud             = (domains_up && domains_up->exp_energy_cb);
  hc_up_ml            = hc->up_ml;
  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.sn     = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  qbt1    = 0;
  q_temp  = 0.;

  qqm[i] = 0.;

  if (with_ud)
    qqmu[0][i] = 0.;

  if (with_gquad)
    expMLstem = exp_E_MLstem(0, -1, -1, pf_params);

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    q_temp = qqm1[i] *
             expMLbase[1];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[j][1];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
    }

    if (with_ud) {
      int cnt;
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if (j - u >= i) {
          if (evaluate(i, j, i, j - u, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
            q_temp2 = qqmu[u][i] *
                      domains_up->exp_energy_cb(vc,
                                                j - u + 1,
                                                j,
                                                VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                domains_up->data) *
                      expMLbase[u];

            if (sc) {
              if (sc->exp_energy_up)
                q_temp2 *= sc->exp_energy_up[j - u + 1][u];

              if (sc->exp_f)
                q_temp2 *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
            }

            q_temp += q_temp2;
          }
        }
      }
      qqmu[0][i] += q_temp;
    }

    qqm[i] += q_temp;
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
    S1    = vc->sequence_encoding;
    S2    = vc->sequence_encoding2;
    type  = get_pair_type_md(S2[i], S2[j], md);

    qbt1 = qb[ij] * exp_E_MLstem(type,
                                 ((i > 1) || circular) ? S1[i - 1] : -1,
                                 ((j < n) || circular) ? S1[j + 1] : -1,
                                 pf_params);
    if (sc)
      if (sc->exp_f)
        qbt1 *= sc->exp_f(i, j, i, j, VRNA_DECOMP_ML_STEM, sc->data);

    qqm[i] += qbt1;

    if (with_ud)
      qqmu[0][i] += qbt1;
  }

  if (with_gquad) {
    /*include gquads into qqm*/
    qqm[i] += G[ij] * expMLstem;

    if (with_ud)
      qqmu[0][i] += G[ij] * expMLstem;
  }

  /*
   *  construction of qm matrix containing multiple loop
   *  partition function contributions from segment i,j
   */
  temp  = 0.0;
  kl    = iidx[i] - j + 1; /* ii-k=[i,k-1] */
  if (hc->f) {
    if (sc && sc->exp_f) {
      for (k = j; k > i; k--, kl++) {
        if (hc->f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          q_temp  = qm[kl] * qqm[k];
          q_temp  *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
          temp    += q_temp;
        }
      }
    } else {
      for (k = j; k > i; k--, kl++)
        if (hc->f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))
          temp += qm[kl] * qqm[k];
    }
  } else {
    if (sc && sc->exp_f) {
      for (k = j; k > i; k--, kl++) {
        q_temp  = qm[kl] * qqm[k];
        q_temp  *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
        temp    += q_temp;
      }
    } else {
      for (k = j; k > i; k--, kl++)
        temp += qm[kl] * qqm[k];
    }
  }

  maxk  = MIN2(i + hc_up_ml[i], j);
  ii    = maxk - i; /* length of unpaired stretch */
  if (with_ud) {
    if (hc->f) {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          if (evaluate(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][ii];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

            temp  += q_temp;
            temp  += q_temp *
                     domains_up->exp_energy_cb(vc,
                                               i, k - 1,
                                               VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                               domains_up->data);
          }
        }
      } else {
        for (k = maxk; k > i; k--, ii--) {
          if (evaluate(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            temp  += q_temp;
            temp  += q_temp *
                     domains_up->exp_energy_cb(vc,
                                               i, k - 1,
                                               VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                               domains_up->data);
          }
        }
      }
    } else {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][ii];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

          temp  += q_temp;
          temp  += q_temp *
                   domains_up->exp_energy_cb(vc,
                                             i, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                             domains_up->data);
        }
      } else {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          temp  += q_temp;
          temp  += q_temp *
                   domains_up->exp_energy_cb(vc,
                                             i, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                             domains_up->data);
        }
      }
    }
  } else {
    if (hc->f) {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][ii];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

            temp += q_temp;
          }
        }
      } else {
        for (k = maxk; k > i; k--, ii--)
          if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data))
            temp += expMLbase[ii] *
                    qqm[k];
      }
    } else {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][ii];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

          temp += q_temp;
        }
      } else {
        for (k = maxk; k > i; k--, ii--)
          temp += expMLbase[ii] *
                  qqm[k];
      }
    }
  }

  return temp + qqm[i];
}


PRIVATE FLT_OR_DBL
exp_E_ml_fast_window(vrna_fold_compound_t *vc,
                     int                  i,
                     int                  j,
                     vrna_mx_pf_aux_ml_t  *aux_mx)
{
  short                     *S1;
  int                       n, k, maxk, ii, with_ud, u, circular, with_gquad,
                            *hc_up_ml, type;
  FLT_OR_DBL                qbt1, temp, **qm, **qb, *qqm, *qqm1, **qqmu, q_temp, q_temp2, **G,
                            *expMLbase,
                            expMLstem;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n           = (int)vc->length;
  qqm         = aux_mx->qqm;
  qqm1        = aux_mx->qqm1;
  qqmu        = aux_mx->qqmu;
  qm          = vc->exp_matrices->qm_local;
  qb          = vc->exp_matrices->qb_local;
  G           = vc->exp_matrices->G_local;
  expMLbase   = vc->exp_matrices->expMLbase;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  hc          = vc->hc;
  sc          = vc->sc;
  domains_up  = vc->domains_up;
  circular    = md->circ;
  with_gquad  = md->gquad;
  with_ud     = (domains_up && domains_up->exp_energy_cb);
  hc_up_ml    = hc->up_ml;


  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.sn         = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  qbt1    = 0;
  q_temp  = 0.;

  qqm[i] = 0.;

  if (with_ud)
    qqmu[0][i] = 0.;

  if (with_gquad)
    expMLstem = exp_E_MLstem(0, -1, -1, pf_params);

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    q_temp = qqm1[i] *
             expMLbase[1];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[j][1];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
    }

    if (with_ud) {
      int cnt;
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if (j - u >= i) {
          if (evaluate(i, j, i, j - u, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
            q_temp2 = qqmu[u][i] *
                      domains_up->exp_energy_cb(vc,
                                                j - u + 1,
                                                j,
                                                VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                domains_up->data) *
                      expMLbase[u];

            if (sc) {
              if (sc->exp_energy_up)
                q_temp2 *= sc->exp_energy_up[j - u + 1][u];

              if (sc->exp_f)
                q_temp2 *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
            }

            q_temp += q_temp2;
          }
        }
      }
      qqmu[0][i] += q_temp;
    }

    qqm[i] += q_temp;
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
    S1    = vc->sequence_encoding;
    type  = get_pair_type_md(S1[i], S1[j], md);

    qbt1 = qb[i][j] * exp_E_MLstem(type,
                                   ((i > 1) || circular) ? S1[i - 1] : -1,
                                   ((j < n) || circular) ? S1[j + 1] : -1,
                                   pf_params);
    if (sc)
      if (sc->exp_f)
        qbt1 *= sc->exp_f(i, j, i, j, VRNA_DECOMP_ML_STEM, sc->data);

    qqm[i] += qbt1;

    if (with_ud)
      qqmu[0][i] += qbt1;
  }

  if (with_gquad) {
    /*include gquads into qqm*/
    qqm[i] += G[i][j] * expMLstem;

    if (with_ud)
      qqmu[0][i] += G[i][j] * expMLstem;
  }

  /*
   *  construction of qm matrix containing multiple loop
   *  partition function contributions from segment i,j
   */
  temp = 0.0;
  if (hc->f) {
    if (sc && sc->exp_f) {
      for (k = j; k > i; k--) {
        if (hc->f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          q_temp  = qm[i][k - 1] * qqm[k];
          q_temp  *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
          temp    += q_temp;
        }
      }
    } else {
      for (k = j; k > i; k--)
        if (hc->f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))
          temp += qm[i][k - 1] * qqm[k];
    }
  } else {
    if (sc && sc->exp_f) {
      for (k = j; k > i; k--) {
        q_temp  = qm[i][k - 1] * qqm[k];
        q_temp  *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
        temp    += q_temp;
      }
    } else {
      for (k = j; k > i; k--)
        temp += qm[i][k - 1] * qqm[k];
    }
  }

  maxk  = MIN2(i + hc_up_ml[i], j);
  ii    = maxk - i; /* length of unpaired stretch */
  if (with_ud) {
    if (hc->f) {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          if (evaluate(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][ii];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

            temp  += q_temp;
            temp  += q_temp *
                     domains_up->exp_energy_cb(vc,
                                               i, k - 1,
                                               VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                               domains_up->data);
          }
        }
      } else {
        for (k = maxk; k > i; k--, ii--) {
          if (evaluate(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            temp  += q_temp;
            temp  += q_temp *
                     domains_up->exp_energy_cb(vc,
                                               i, k - 1,
                                               VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                               domains_up->data);
          }
        }
      }
    } else {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][ii];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

          temp  += q_temp;
          temp  += q_temp *
                   domains_up->exp_energy_cb(vc,
                                             i, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                             domains_up->data);
        }
      } else {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          temp  += q_temp;
          temp  += q_temp *
                   domains_up->exp_energy_cb(vc,
                                             i, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                             domains_up->data);
        }
      }
    }
  } else {
    if (hc->f) {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][ii];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

            temp += q_temp;
          }
        }
      } else {
        for (k = maxk; k > i; k--, ii--)
          if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data))
            temp += expMLbase[ii] *
                    qqm[k];
      }
    } else {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][ii];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

          temp += q_temp;
        }
      } else {
        for (k = maxk; k > i; k--, ii--)
          temp += expMLbase[ii] *
                  qqm[k];
      }
    }
  }

  return temp + qqm[i];
}


PRIVATE FLT_OR_DBL
exp_E_ml_fast_comparative(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          vrna_mx_pf_aux_ml_t   *aux_mx)
{
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  int                       n, s, n_seq, *iidx, k, ij, kl, maxk, ii, circular, *hc_up_ml, type;
  FLT_OR_DBL                temp, *qm, *qb, *qqm, *qqm1, q_temp, *expMLbase;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n         = (int)vc->length;
  n_seq     = vc->n_seq;
  iidx      = vc->iindx;
  ij        = iidx[i] - j;
  S         = vc->S;
  S5        = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
  S3        = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
  a2s       = vc->a2s;
  qqm       = aux_mx->qqm;
  qqm1      = aux_mx->qqm1;
  qm        = vc->exp_matrices->qm;
  qb        = vc->exp_matrices->qb;
  expMLbase = vc->exp_matrices->expMLbase;
  pf_params = vc->exp_params;
  md        = &(pf_params->model_details);
  hc        = vc->hc;
  scs       = vc->scs;
  circular  = md->circ;
  hc_up_ml  = hc->up_ml;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.sn     = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  q_temp = 0.;

  qqm[i] = 0.;


  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    q_temp = qqm1[i] *
             expMLbase[1];

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->exp_energy_up)
            q_temp *= scs[s]->exp_energy_up[a2s[s][j]][1];
      }
    }

    qqm[i] += q_temp;
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
    q_temp = qb[ij];

    for (s = 0; s < n_seq; s++) {
      type    = get_pair_type_md(S[s][i], S[s][j], md);
      q_temp  *= exp_E_MLstem(type,
                              ((i > 1) || circular) ? S5[s][i] : -1,
                              ((j < n) || circular) ? S3[s][j] : -1,
                              pf_params);
    }

    qqm[i] += q_temp;
  }

  /*
   *  construction of qm matrix containing multiple loop
   *  partition function contributions from segment i,j
   */
  temp  = 0.0;
  kl    = iidx[i] - j + 1; /* ii-k=[i,k-1] */
  if (hc->f) {
    for (k = j; k > i; k--, kl++)
      if (hc->f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))
        temp += qm[kl] * qqm[k];
  } else {
    for (k = j; k > i; k--, kl++)
      temp += qm[kl] * qqm[k];
  }

  maxk  = MIN2(i + hc_up_ml[i], j);
  ii    = maxk - i; /* length of unpaired stretch */

  if (hc->f) {
    if (scs) {
      for (k = maxk; k > i; k--, ii--) {
        if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          for (s = 0; s < n_seq; s++) {
            if (scs[s])
              if (scs[s]->exp_energy_up)
                q_temp *= scs[s]->exp_energy_up[a2s[s][i]][a2s[s][k] - a2s[s][i]];
          }
          temp += q_temp;
        }
      }
    } else {
      for (k = maxk; k > i; k--, ii--)
        if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data))
          temp += expMLbase[ii] *
                  qqm[k];
    }
  } else {
    if (scs) {
      for (k = maxk; k > i; k--, ii--) {
        q_temp = expMLbase[ii] *
                 qqm[k];

        for (s = 0; s < n_seq; s++) {
          if (scs[s])
            if (scs[s]->exp_energy_up)
              q_temp *= scs[s]->exp_energy_up[a2s[s][i]][a2s[s][k] - a2s[s][i]];
        }
        temp += q_temp;
      }
    } else {
      for (k = maxk; k > i; k--, ii--)
        temp += expMLbase[ii] *
                qqm[k];
    }
  }

  return temp + qqm[i];
}
