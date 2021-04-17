#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/loops/external.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/loops/hairpin.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "hairpin_hc.inc"
#include "hairpin_sc_pf.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE FLT_OR_DBL
exp_eval_hp_loop(vrna_fold_compound_t *fc,
                 int                  i,
                 int                  j);


PRIVATE FLT_OR_DBL
exp_eval_ext_hp_loop(vrna_fold_compound_t *fc,
                     int                  i,
                     int                  j);


PRIVATE FLT_OR_DBL
exp_eval_hp_loop_fake(vrna_fold_compound_t  *fc,
                      int                   i,
                      int                   j);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

/**
 *  @brief High-Level function for hairpin loop energy evaluation (partition function variant)
 *
 *  @see E_hp_loop() for it's free energy counterpart
 */
PUBLIC FLT_OR_DBL
vrna_exp_E_hp_loop(vrna_fold_compound_t *fc,
                   int                  i,
                   int                  j)
{
  vrna_callback_hc_evaluate *evaluate;
  struct hc_hp_def_dat      hc_dat_local;

  if (fc->hc->type == VRNA_HC_WINDOW)
    evaluate = prepare_hc_hp_def_window(fc, &hc_dat_local);
  else
    evaluate = prepare_hc_hp_def(fc, &hc_dat_local);

  if ((i > 0) && (j > 0)) {
    if (evaluate(i, j, i, j, VRNA_DECOMP_PAIR_HP, &hc_dat_local)) {
      if (j > i)  /* linear case */
        return exp_eval_hp_loop(fc, i, j);
      else        /* circular case */
        return exp_eval_ext_hp_loop(fc, j, i);
    }
  }

  return 0.;
}


PRIVATE FLT_OR_DBL
exp_eval_hp_loop_fake(vrna_fold_compound_t  *fc,
                      int                   i,
                      int                   j)
{
  short             *S, *S2, s5, s3;
  unsigned int      *sn, *ss, *se;
  int               u, type, *iidx, *jidx;
  FLT_OR_DBL        qq, temp, *q, *scale;
  vrna_exp_param_t  *pf_params;
  vrna_sc_t         *sc;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;

  iidx        = fc->iindx;
  jidx        = fc->jindx;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  q           = fc->exp_matrices->q;
  scale       = fc->exp_matrices->scale;
  sn          = fc->strand_number;
  ss          = fc->strand_start;
  se          = fc->strand_end;
  domains_up  = fc->domains_up;

  qq = 0;

  switch (fc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      S2    = fc->sequence_encoding2;
      sc    = fc->sc;
      u     = j - i - 1;
      type  = vrna_get_ptype_md(S2[j], S2[i], md);

      temp = scale[2];

      if (u > 0) {
        /* add contribution for [i + 1, end(strand(i))] */
        if (sn[i] == sn[i + 1])
          temp *= q[iidx[i + 1] - se[sn[i]]];

        /* add contribution for [start(strand(j)), j - 1] */
        if (sn[j - 1] == sn[j])
          temp *= q[iidx[ss[sn[j]]] - (j - 1)];
      }

      s5  = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;
      s3  = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;

      temp *= vrna_exp_E_ext_stem(type, s5, s3, pf_params);

      qq += temp;

      /* add soft constraints */
      if (sc) {
        if (sc->exp_energy_up)
          qq *= sc->exp_energy_up[i + 1][u];

        if (sc->exp_energy_bp)
          qq *= sc->exp_energy_bp[jidx[j] + i];

        if (sc->exp_f)
          qq *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
      }

      if (domains_up && domains_up->exp_energy_cb) {
        /* we always consider both, bound and unbound state */
        qq += qq * domains_up->exp_energy_cb(fc,
                                             i + 1, j - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                             domains_up->data);
      }

      break;

    /* nothing */
    default:
      break;
  }

  return qq;
}


PRIVATE FLT_OR_DBL
exp_eval_hp_loop(vrna_fold_compound_t *fc,
                 int                  i,
                 int                  j)
{
  char                  **Ss;
  unsigned int          **a2s;
  short                 *S, *S2, **SS, **S5, **S3;
  unsigned int          *sn;
  int                   u, type, n_seq, s;
  FLT_OR_DBL            q, qbt1, *scale;
  vrna_exp_param_t      *P;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  struct sc_hp_exp_dat  sc_wrapper;

  P           = fc->exp_params;
  md          = &(P->model_details);
  sn          = fc->strand_number;
  scale       = fc->exp_matrices->scale;
  domains_up  = fc->domains_up;

  init_sc_hp_exp(fc, &sc_wrapper);
  q = 0.;

  if (sn[j] != sn[i])
    return exp_eval_hp_loop_fake(fc, i, j);

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      S2    = fc->sequence_encoding2;
      u     = j - i - 1;
      type  = vrna_get_ptype_md(S2[i], S2[j], md);

      if (sn[j] == sn[i]) {
        /* regular hairpin loop */
        q = exp_E_Hairpin(u, type, S[i + 1], S[j - 1], fc->sequence + i - 1, P);
      } else {
        /*
         * hairpin-like exterior loop (for cofolding)
         * this is currently handle somewhere else
         */
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      SS    = fc->S;
      S5    = fc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
      S3    = fc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
      Ss    = fc->Ss;
      a2s   = fc->a2s;
      n_seq = fc->n_seq;
      qbt1  = 1.;

      for (s = 0; s < n_seq; s++) {
        u = a2s[s][j - 1] - a2s[s][i];
        if (a2s[s][i] < 1)
          continue;

        type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
        qbt1  *= exp_E_Hairpin(u, type, S3[s][i], S5[s][j], Ss[s] + a2s[s][i] - 1, P);
      }

      q = qbt1;
      break;

    default:
      break;
  }

  /* add soft constraints */
  if (sc_wrapper.pair)
    q *= sc_wrapper.pair(i, j, &sc_wrapper);

  if (domains_up && domains_up->exp_energy_cb) {
    /* we always consider both, bound and unbound state */
    q += q * domains_up->exp_energy_cb(fc,
                                       i + 1, j - 1,
                                       VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                       domains_up->data);
  }

  q *= scale[j - i + 1];

  free_sc_hp_exp(&sc_wrapper);

  return q;
}


PRIVATE FLT_OR_DBL
exp_eval_ext_hp_loop(vrna_fold_compound_t *fc,
                     int                  i,
                     int                  j)
{
  char                  **Ss, *sequence, loopseq[10] = {
    0
  };
  unsigned int          **a2s;
  short                 *S, *S2, **SS, **S5, **S3;
  int                   u1, u2, n, type, n_seq, s, noGUclosure;
  FLT_OR_DBL            q, qbt1, *scale;
  vrna_exp_param_t      *P;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  struct sc_hp_exp_dat  sc_wrapper;

  n           = fc->length;
  P           = fc->exp_params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  scale       = fc->exp_matrices->scale;
  domains_up  = fc->domains_up;

  init_sc_hp_exp(fc, &sc_wrapper);

  q   = 0.;
  u1  = n - j;
  u2  = i - 1;

  if ((u1 + u2) < 3)
    return q;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sequence  = fc->sequence;
      S         = fc->sequence_encoding;
      S2        = fc->sequence_encoding2;
      type      = vrna_get_ptype_md(S2[j], S2[i], md);

      if (((type == 3) || (type == 4)) && noGUclosure)
        return q;

      /* get the loop sequence */
      if ((u1 + u2) < 7) {
        memcpy(loopseq, sequence + j - 1, sizeof(char) * (u1 + 1));
        memcpy(loopseq + u1 + 1, sequence, sizeof(char) * (u2 + 1));
        loopseq[u1 + u2 + 2] = '\0';
      }

      q = exp_E_Hairpin(u1 + u2, type, S[j + 1], S[i - 1], loopseq, P);

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      SS    = fc->S;
      S5    = fc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
      S3    = fc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
      Ss    = fc->Ss;
      a2s   = fc->a2s;
      n_seq = fc->n_seq;
      qbt1  = 1.;

      for (s = 0; s < n_seq; s++) {
        const int u1_local  = a2s[s][n] - a2s[s][j];
        const int u2_local  = a2s[s][i - 1];
        memset(loopseq, '\0', sizeof(loopseq));

        if ((u1_local + u2_local) < 7) {
          memcpy(loopseq, Ss[s] + a2s[s][j] - 1, sizeof(char) * (u1_local + 1));
          memcpy(loopseq + u1_local + 1, Ss[s], sizeof(char) * (u2_local + 1));
          loopseq[u1_local + u2_local + 2] = '\0';
        }

        type  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
        qbt1  *= exp_E_Hairpin(u1_local + u2_local, type, S3[s][j], S5[s][i], loopseq, P);
      }

      q = qbt1;

      break;

    default:
      break;
  }

  /* add soft constraints */
  if (sc_wrapper.pair_ext)
    q *= sc_wrapper.pair_ext(i, j, &sc_wrapper);

  if (domains_up && domains_up->exp_energy_cb) {
    /* we always consider both, bound and unbound state */
    q += q * domains_up->exp_energy_cb(fc,
                                       j + 1, i - 1,
                                       VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                       domains_up->data);
  }

  q *= scale[u1 + u2];

  free_sc_hp_exp(&sc_wrapper);

  return q;
}
