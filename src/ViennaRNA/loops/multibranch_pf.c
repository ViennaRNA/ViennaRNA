#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/loops/external.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/loops/multibranch.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "multibranch_hc.inc"
#include "multibranch_sc_pf.inc"

struct vrna_mx_pf_aux_ml_s {
  FLT_OR_DBL  *qqm;
  FLT_OR_DBL  *qqm1;

  int         qqmu_size;
  FLT_OR_DBL  **qqmu;
};


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast(vrna_fold_compound_t       *fc,
                   int                        i,
                   int                        j,
                   struct vrna_mx_pf_aux_ml_s *aux_mx);


PRIVATE FLT_OR_DBL
exp_E_ml_fast(vrna_fold_compound_t        *fc,
              int                         i,
              int                         j,
              struct vrna_mx_pf_aux_ml_s  *aux_mx);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC FLT_OR_DBL
vrna_exp_E_mb_loop_fast(vrna_fold_compound_t        *fc,
                        int                         i,
                        int                         j,
                        struct vrna_mx_pf_aux_ml_s  *aux_mx)
{
  FLT_OR_DBL q = 0.;

  if ((fc) && (aux_mx))
    q = exp_E_mb_loop_fast(fc, i, j, aux_mx);

  return q;
}


PUBLIC FLT_OR_DBL
vrna_exp_E_ml_fast(vrna_fold_compound_t       *fc,
                   int                        i,
                   int                        j,
                   struct vrna_mx_pf_aux_ml_s *aux_mx)
{
  FLT_OR_DBL q = 0.;

  if ((fc) && (aux_mx))
    q = exp_E_ml_fast(fc, i, j, aux_mx);

  return q;
}


PUBLIC struct vrna_mx_pf_aux_ml_s *
vrna_exp_E_ml_fast_init(vrna_fold_compound_t *fc)
{
  struct vrna_mx_pf_aux_ml_s *aux_mx = NULL;

  if (fc) {
    int         i, j, d, n, u, turn, ij, *iidx;
    FLT_OR_DBL  *qm;

    n     = (int)fc->length;
    iidx  = fc->iindx;
    turn  = fc->exp_params->model_details.min_loop_size;
    qm    = fc->exp_matrices->qm;

    /* allocate memory for helper arrays */
    aux_mx =
      (struct vrna_mx_pf_aux_ml_s *)vrna_alloc(sizeof(struct vrna_mx_pf_aux_ml_s));
    aux_mx->qqm       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqm1      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqmu_size = 0;
    aux_mx->qqmu      = NULL;

    if (fc->type == VRNA_FC_TYPE_SINGLE) {
      vrna_ud_t *domains_up = fc->domains_up;
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

    if (fc->hc->type == VRNA_HC_WINDOW) {
    } else {
      for (d = 0; d <= turn; d++)
        for (i = 1; i <= n - d; i++) {
          j   = i + d;
          ij  = iidx[i] - j;

          if (j > n)
            continue;

          qm[ij] = 0.;
        }

      if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_m)) {
        for (d = 0; d <= turn; d++)
          for (i = 1; i <= n - d; i++) {
            j   = i + d;
            ij  = iidx[i] - j;

            if (j > n)
              continue;

            qm[ij] += fc->aux_grammar->cb_aux_exp_m(fc, i, j, fc->aux_grammar->data);
          }
      }
    }
  }

  return aux_mx;
}


PUBLIC void
vrna_exp_E_ml_fast_rotate(struct vrna_mx_pf_aux_ml_s *aux_mx)
{
  if (aux_mx) {
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
vrna_exp_E_ml_fast_free(struct vrna_mx_pf_aux_ml_s *aux_mx)
{
  if (aux_mx) {
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


PUBLIC const FLT_OR_DBL *
vrna_exp_E_ml_fast_qqm(struct vrna_mx_pf_aux_ml_s *aux_mx)
{
  if (aux_mx)
    return (const FLT_OR_DBL *)aux_mx->qqm;

  return NULL;
}


PUBLIC const FLT_OR_DBL *
vrna_exp_E_ml_fast_qqm1(struct vrna_mx_pf_aux_ml_s *aux_mx)
{
  if (aux_mx)
    return (const FLT_OR_DBL *)aux_mx->qqm1;

  return NULL;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast(vrna_fold_compound_t       *fc,
                   int                        i,
                   int                        j,
                   struct vrna_mx_pf_aux_ml_s *aux_mx)
{
  unsigned char             sliding_window;
  char                      *ptype, **ptype_local;
  short                     *S1, **SS, **S5, **S3;
  unsigned int              *sn, n_seq, s, *se;
  int                       ij, k, kl, *my_iindx, *jindx, *rtype, tt;
  FLT_OR_DBL                qbt1, temp, qqqmmm, *qm, **qm_local, *scale, expMLclosing, *qqm1;
  vrna_hc_t                 *hc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;
  struct sc_wrapper_exp_ml  sc_wrapper;

  qqm1            = aux_mx->qqm1;
  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  se              = fc->strand_end;
  my_iindx        = (sliding_window) ? NULL : fc->iindx;
  jindx           = (sliding_window) ? NULL : fc->jindx;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     = (sliding_window) ? fc->ptype_local : NULL;
  S1              = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  qm              = (sliding_window) ? NULL : fc->exp_matrices->qm;
  qm_local        = (sliding_window) ? fc->exp_matrices->qm_local : NULL;
  scale           = fc->exp_matrices->scale;
  pf_params       = fc->exp_params;
  md              = &(pf_params->model_details);
  ij              = (sliding_window) ? 0 : jindx[j] + i;
  sn              = fc->strand_number;
  hc              = fc->hc;
  expMLclosing    = pf_params->expMLclosing;
  qbt1            = 0.;
  rtype           = &(md->rtype[0]);
  evaluate        = prepare_hc_default(fc, &hc_dat_local);

  init_sc_wrapper_ml(fc, &sc_wrapper);

  /* multiple stem loop contribution */
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    qqqmmm = pow(expMLclosing, (double)n_seq) *
             scale[2];

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        tt = (sliding_window) ?
             rtype[vrna_get_ptype_window(i, j + i, ptype_local)] :
             rtype[vrna_get_ptype(ij, ptype)];
        qqqmmm *= exp_E_MLstem(tt, S1[j - 1], S1[i + 1], pf_params);

        break;


      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++) {
          tt      = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
          qqqmmm  *= exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params);
        }
        break;
    }

    if (sc_wrapper.pair)
      qqqmmm *= sc_wrapper.pair(i, j, &sc_wrapper);

    FLT_OR_DBL *qqm1_tmp = qqm1;

    if (hc->f) {
      qqm1_tmp  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
      qqm1_tmp  -= i;

      for (k = i + 2; k <= j - 1; k++) {
        qqm1_tmp[k] = qqm1[k];
        if (!evaluate(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, &hc_dat_local))
          qqm1_tmp[k] = 0.;
      }
    }

    if (sc_wrapper.decomp_ml) {
      if (qqm1_tmp == qqm1) {
        qqm1_tmp  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
        qqm1_tmp  -= i;

        for (k = i + 2; k <= j - 1; k++)
          qqm1_tmp[k] = qqm1[k];
      }

      for (k = i + 2; k <= j - 1; k++)
        qqm1_tmp[k] *= sc_wrapper.decomp_ml(i + 1, j - 1, k - 1, k, &sc_wrapper);
    }

    temp = 0.0;

    /* set initial decomposition split point */
    k = i + 2;

    if (sliding_window) {
      for (; k <= j - 1; k++, kl--)
        temp += qm_local[i + 1][k - 1] *
                qqm1_tmp[k];
    } else {
      kl = my_iindx[i + 1] - (i + 1);
      /*
       *  loop over entire range but skip decompositions with in-between strand nick,
       *  this should be faster than evaluating hard constraints callback for each
       *  decomposition
       */
      while (1) {
        /* limit for-loop to last nucleotide of 5' part strand */
        int stop = MIN2(j - 1, se[sn[k - 1]]);

        for (; k <= stop; k++, kl--)
          temp += qm[kl] *
                  qqm1_tmp[k];

        k++;
        kl--;

        if (stop == j - 1)
          break;
      }
    }

    if (qqm1_tmp != qqm1) {
      qqm1_tmp += i;
      free(qqm1_tmp);
    }

    qbt1 += temp *
            qqqmmm;
  }

  free_sc_wrapper_ml(&sc_wrapper);

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_ml_fast(vrna_fold_compound_t        *fc,
              int                         i,
              int                         j,
              struct vrna_mx_pf_aux_ml_s  *aux_mx)
{
  unsigned char             sliding_window;
  short                     *S1, *S2, **SS, **S5, **S3;
  unsigned int              *sn, *ss, *se, n_seq, s;
  int                       n, *iidx, k, ij, kl, maxk, ii, with_ud, u, circular, with_gquad,
                            *hc_up_ml, type;
  FLT_OR_DBL                qbt1, temp, *qm, *qb, *qqm, *qqm1, **qqmu, q_temp, q_temp2, *G,
                            *expMLbase, **qb_local, **qm_local, **G_local;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_t                 *hc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;
  struct sc_wrapper_exp_ml  sc_wrapper;

  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n               = (int)fc->length;
  sn              = fc->strand_number;
  ss              = fc->strand_start;
  se              = fc->strand_end;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  iidx            = (sliding_window) ? NULL : fc->iindx;
  ij              = (sliding_window) ? 0 : iidx[i] - j;
  qqm             = aux_mx->qqm;
  qqm1            = aux_mx->qqm1;
  qqmu            = aux_mx->qqmu;
  qm              = (sliding_window) ? NULL : fc->exp_matrices->qm;
  qb              = (sliding_window) ? NULL : fc->exp_matrices->qb;
  G               = (sliding_window) ? NULL : fc->exp_matrices->G;
  qm_local        = (sliding_window) ? fc->exp_matrices->qm_local : NULL;
  qb_local        = (sliding_window) ? fc->exp_matrices->qb_local : NULL;
  G_local         = (sliding_window) ? fc->exp_matrices->G_local : NULL;
  expMLbase       = fc->exp_matrices->expMLbase;
  pf_params       = fc->exp_params;
  md              = &(pf_params->model_details);
  hc              = fc->hc;
  domains_up      = fc->domains_up;
  circular        = md->circ;
  with_gquad      = md->gquad;
  with_ud         = (domains_up && domains_up->exp_energy_cb);
  hc_up_ml        = hc->up_ml;
  evaluate        = prepare_hc_default(fc, &hc_dat_local);

  init_sc_wrapper_ml(fc, &sc_wrapper);

  qbt1    = 0;
  q_temp  = 0.;

  qqm[i] = 0.;

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    q_temp = qqm1[i] *
             expMLbase[1];

    if (sc_wrapper.red_ml)
      q_temp *= sc_wrapper.red_ml(i, j, i, j - 1, &sc_wrapper);

    qqm[i] += q_temp;
  }

  if (with_ud) {
    q_temp = 0.;

    int cnt;
    for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
      u = domains_up->uniq_motif_size[cnt];
      if (j - u >= i) {
        if (evaluate(i, j, i, j - u, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
          q_temp2 = qqmu[u][i] *
                    domains_up->exp_energy_cb(fc,
                                              j - u + 1,
                                              j,
                                              VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP |
                                              VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                              domains_up->data) *
                    expMLbase[u];

          if (sc_wrapper.red_ml)
            q_temp2 *= sc_wrapper.red_ml(i, j, i, j - u, &sc_wrapper);

          q_temp += q_temp2;
        }
      }
    }

    qqm[i] += q_temp;
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
    qbt1 = (sliding_window) ? qb_local[i][j] : qb[ij];

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        S1    = fc->sequence_encoding;
        S2    = fc->sequence_encoding2;
        type  = vrna_get_ptype_md(S2[i], S2[j], md);

        qbt1 *= exp_E_MLstem(type,
                             ((i > 1) || circular) ? S1[i - 1] : -1,
                             ((j < n) || circular) ? S1[j + 1] : -1,
                             pf_params);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        q_temp = 1.;
        for (s = 0; s < n_seq; s++) {
          type    = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
          q_temp  *= exp_E_MLstem(type,
                                  ((i > 1) || circular) ? S5[s][i] : -1,
                                  ((j < n) || circular) ? S3[s][j] : -1,
                                  pf_params);
        }
        qbt1 *= q_temp;
        break;
    }

    if (sc_wrapper.red_stem)
      qbt1 *= sc_wrapper.red_stem(i, j, i, j, &sc_wrapper);

    qqm[i] += qbt1;
  }

  if (with_gquad) {
    q_temp  = (sliding_window) ? G_local[i][j] : G[ij];
    qqm[i]  += q_temp *
               pow(exp_E_MLstem(0, -1, -1, pf_params), (double)n_seq);
  }

  if (with_ud)
    qqmu[0][i] = qqm[i];

  /*
   *  construction of qm matrix containing multiple loop
   *  partition function contributions from segment i,j
   */
  FLT_OR_DBL *qqm_tmp = qqm;

  /* apply hard constraints if necessary */
  if (hc->f) {
    qqm_tmp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
    qqm_tmp -= i;

    for (k = j; k > i; k--) {
      qqm_tmp[k] = qqm[k];
      if (!evaluate(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, &hc_dat_local))
        qqm_tmp[k] = 0.;
    }
  }

  /* apply soft constraints if necessary */
  if (sc_wrapper.decomp_ml) {
    if (qqm_tmp == qqm) {
      qqm_tmp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
      qqm_tmp -= i;

      for (k = j; k > i; k--)
        qqm_tmp[k] = qqm[k];
    }

    for (k = j; k > i; k--)
      qqm_tmp[k] *= sc_wrapper.decomp_ml(i, j, k - 1, k, &sc_wrapper);
  }

  /* finally, decompose segment */
  temp  = 0.0;
  k     = j;

  if (sliding_window) {
    for (; k > i; k--)
      temp += qm_local[i][k - 1] *
              qqm_tmp[k];
  } else {
    kl = iidx[i] - j + 1; /* ii-k=[i,k-1] */

    while (1) {
      /* limit for-loop to first nucleotide of 3' part strand */
      int stop = MAX2(i, ss[sn[k]]);
      for (; k > stop; k--, kl++)
        temp += qm[kl] *
                qqm_tmp[k];

      k--;
      kl++;

      if (stop == i)
        break;
    }
  }

  /* determine last nucleotide position for unpaired stretch */
  maxk = j;

  /* obey hard constraints (those that we can simply look-up) */
  if (maxk > i + hc_up_ml[i])
    maxk = i + hc_up_ml[i];

  /* obey connected components constraint (for multi-RNA folding) */
  if (maxk > se[sn[i]])
    maxk = se[sn[i]];

  if (qqm_tmp != qqm) {
    /* we've been using this helper array, so it's likely we use it again... */
    for (k = maxk; k > i; k--)
      qqm_tmp[k] = qqm[k];
  }

  /* apply hard constraints if necessary */
  if (hc->f) {
    if (qqm_tmp == qqm) {
      qqm_tmp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
      qqm_tmp -= i;

      for (k = maxk; k > i; k--)
        qqm_tmp[k] = qqm[k];
    }

    for (k = maxk; k > i; k--)
      if (!evaluate(i, j, k, j, VRNA_DECOMP_ML_ML, &hc_dat_local))
        qqm_tmp[k] = 0.;
  }

  /* apply soft constraints if necessary */
  if (sc_wrapper.red_ml) {
    if (qqm_tmp == qqm) {
      qqm_tmp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
      qqm_tmp -= i;

      for (k = maxk; k > i; k--)
        qqm_tmp[k] = qqm[k];
    }

    for (k = maxk; k > i; k--)
      qqm_tmp[k] *= sc_wrapper.red_ml(i, j, k, j, &sc_wrapper);
  }

  ii = maxk - i; /* length of unpaired stretch */

  /* finally, decompose segment */
  for (k = maxk; k > i; k--, ii--)
    temp += expMLbase[ii] *
            qqm_tmp[k];

  if (with_ud) {
    ii = maxk - i; /* length of unpaired stretch */

    for (k = maxk; k > i; k--, ii--)
      temp += expMLbase[ii] *
              qqm_tmp[k] *
              domains_up->exp_energy_cb(fc,
                                        i, k - 1,
                                        VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                        domains_up->data);
  }

  if (qqm_tmp != qqm) {
    qqm_tmp += i;
    free(qqm_tmp);
  }

  /* apply auxiliary grammar rule for multibranch loop case */
  if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_m))
    temp += fc->aux_grammar->cb_aux_exp_m(fc, i, j, fc->aux_grammar->data);

  free_sc_wrapper_ml(&sc_wrapper);

  return temp + qqm[i];
}
