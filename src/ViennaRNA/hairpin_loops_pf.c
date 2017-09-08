#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/exterior_loops.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/hairpin_loops.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "hairpin_loops.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE FLT_OR_DBL
exp_eval_hp_loop(vrna_fold_compound_t *vc,
                 int                  i,
                 int                  j);


PRIVATE FLT_OR_DBL
exp_eval_ext_hp_loop(vrna_fold_compound_t *vc,
                     int                  i,
                     int                  j);


PRIVATE FLT_OR_DBL
exp_eval_hp_loop_fake(vrna_fold_compound_t  *vc,
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
vrna_exp_E_hp_loop(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j)
{
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  if (vc->hc->type == VRNA_HC_WINDOW)
    evaluate = prepare_hc_default_window(vc, &hc_dat_local);
  else
    evaluate = prepare_hc_default(vc, &hc_dat_local);

  if ((i > 0) && (j > 0)) {
    if (evaluate(i, j, i, j, VRNA_DECOMP_PAIR_HP, &hc_dat_local)) {
      if (j > i)  /* linear case */
        return exp_eval_hp_loop(vc, i, j);
      else        /* circular case */
        return exp_eval_ext_hp_loop(vc, j, i);
    }
  }

  return 0.;
}


PRIVATE FLT_OR_DBL
exp_eval_hp_loop_fake(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j)
{
  short             *S, *S2, s5, s3;
  unsigned int      strands, *sn, *so, *ss, *se;
  int               u, type, *iidx;
  FLT_OR_DBL        qq, temp, *q, *scale;
  vrna_exp_param_t  *pf_params;
  vrna_sc_t         *sc;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;

  iidx        = vc->iindx;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  q           = vc->exp_matrices->q;
  scale       = vc->exp_matrices->scale;
  strands     = vc->strands;
  sn          = vc->strand_number;
  so          = vc->strand_order;
  ss          = vc->strand_start;
  se          = vc->strand_end;
  domains_up  = vc->domains_up;

  qq = 0;

  switch (vc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      S2    = vc->sequence_encoding2;
      sc    = vc->sc;
      u     = j - i - 1;
      type  = get_pair_type(S2[j], S2[i], md);

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

      temp *= exp_E_ExtLoop(type, s5, s3, pf_params);

      qq += temp;

      /* add soft constraints */
      if (sc) {
        if (sc->exp_energy_up)
          qq *= sc->exp_energy_up[i + 1][u];

        if (sc->exp_energy_bp)
          qq *= sc->exp_energy_bp[iidx[i] - j];

        if (sc->exp_f)
          qq *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
      }

      if (domains_up && domains_up->exp_energy_cb) {
        /* we always consider both, bound and unbound state */
        qq += qq * domains_up->exp_energy_cb(vc,
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
exp_eval_hp_loop(vrna_fold_compound_t *vc,
                 int                  i,
                 int                  j)
{
  char              **Ss;
  unsigned int      **a2s;
  short             *S, **SS, **S5, **S3;
  unsigned int      *sn;
  int               u, type, n_seq, s, *iidx;
  FLT_OR_DBL        q, qbt1, *scale;
  vrna_exp_param_t  *P;
  vrna_sc_t         *sc, **scs;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;

  iidx        = vc->iindx;
  P           = vc->exp_params;
  md          = &(P->model_details);
  sn          = vc->strand_number;
  scale       = vc->exp_matrices->scale;
  domains_up  = vc->domains_up;

  q = 0.;

  if (sn[j] != sn[i])
    return exp_eval_hp_loop_fake(vc, i, j);

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      sc    = vc->sc;
      u     = j - i - 1;
      type  = get_pair_type(S[i], S[j], md);

      if (type == 0)
        type = 7;

      if (sn[j] == sn[i]) {
        /* regular hairpin loop */
        q = exp_E_Hairpin(u, type, S[i + 1], S[j - 1], vc->sequence + i - 1, P);
      } else {
        /* hairpin-like exterior loop (for cofolding) */
        /* this is currently handle somewhere else */
      }

      /* add soft constraints */
      if (sc) {
        if (sc->exp_energy_up)
          q *= sc->exp_energy_up[i + 1][u];

        switch (sc->type) {
          case VRNA_SC_DEFAULT:
            if (sc->exp_energy_bp)
              q *= sc->exp_energy_bp[iidx[i] - j];

            break;

          case VRNA_SC_WINDOW:
            if (sc->exp_energy_bp_local)
              q *= sc->exp_energy_bp_local[i][j - i];

            break;
        }

        if (sc->exp_f)
          q *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
      }

      q *= scale[u + 2];

      if (domains_up && domains_up->exp_energy_cb) {
        /* we always consider both, bound and unbound state */
        q += q * domains_up->exp_energy_cb(vc,
                                           i + 1, j - 1,
                                           VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                           domains_up->data);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      SS    = vc->S;
      S5    = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
      S3    = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
      Ss    = vc->Ss;
      a2s   = vc->a2s;
      scs   = vc->scs;
      n_seq = vc->n_seq;
      qbt1  = 1.;

      for (s = 0; s < n_seq; s++) {
        u = a2s[s][j - 1] - a2s[s][i];
        if (a2s[s][i] < 1)
          continue;

        char loopseq[10];
        loopseq[0] = '\0';
        if (u < 9)
          strncpy(loopseq, Ss[s] + a2s[s][i] - 1, 10);

        type = get_pair_type(SS[s][i], SS[s][j], md);

        qbt1 *= exp_E_Hairpin(u, type, S3[s][i], S5[s][j], loopseq, P);
      }

      /* add soft constraints */
      if (scs) {
        for (s = 0; s < n_seq; s++) {
          if (scs[s]) {
            u = a2s[s][j - 1] - a2s[s][i];

            if (scs[s]->exp_energy_bp)
              qbt1 *= scs[s]->exp_energy_bp[iidx[i] - j];

            if (scs[s]->exp_energy_up)
              qbt1 *= scs[s]->exp_energy_up[a2s[s][i + 1]][u];

            if (scs[s]->exp_f) {
              qbt1 *= scs[s]->exp_f(a2s[s][i],
                                    a2s[s][j],
                                    a2s[s][i],
                                    a2s[s][j],
                                    VRNA_DECOMP_PAIR_HP,
                                    scs[s]->data);
            }
          }
        }
      }

      q = qbt1 * scale[j - i + 1];
      break;

    default:
      break;
  }

  return q;
}


PRIVATE FLT_OR_DBL
exp_eval_ext_hp_loop(vrna_fold_compound_t *vc,
                     int                  i,
                     int                  j)
{
  char              **Ss, *sequence;
  unsigned int      **a2s;
  short             *S, **SS, **S5, **S3;
  int               u, u1, ij, n, type, n_seq, s, *rtype, *idx, noGUclosure;
  FLT_OR_DBL        q, qbt1, *scale;
  vrna_exp_param_t  *P;
  vrna_sc_t         *sc, **scs;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;

  n           = vc->length;
  idx         = vc->jindx;
  P           = vc->exp_params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  scale       = vc->exp_matrices->scale;
  domains_up  = vc->domains_up;
  rtype       = &(md->rtype[0]);

  q   = 0.;
  u   = n - j + i - 1;
  ij  = idx[j] + i;

  if (u < 3)
    return q;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sequence  = vc->sequence;
      S         = vc->sequence_encoding;
      sc        = vc->sc;
      type      = rtype[vc->ptype[ij]];

      if (type == 0)
        type = 7;

      if (((type == 3) || (type == 4)) && noGUclosure)
        return q;

      /* get the loop sequence */
      char loopseq[10];
      loopseq[0] = '\0';
      if (u < 7) {
        strcpy(loopseq, sequence + j - 1);
        strncat(loopseq, sequence, i);
      }

      q = exp_E_Hairpin(u, type, S[j + 1], S[i - 1], loopseq, P);

      /* add soft constraints */
      if (sc) {
        if (sc->exp_energy_up)
          q *= sc->exp_energy_up[1][i - 1] *
               sc->exp_energy_up[j + 1][n - j];

        if (sc->exp_f)
          q *= sc->exp_f(j, i, j, i, VRNA_DECOMP_PAIR_HP, sc->data);
      }

      q *= scale[u];

      if (domains_up && domains_up->exp_energy_cb) {
        /* we always consider both, bound and unbound state */
        q += q * domains_up->exp_energy_cb(vc,
                                           j + 1, i - 1,
                                           VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                           domains_up->data);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      SS    = vc->S;
      S5    = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
      S3    = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
      Ss    = vc->Ss;
      a2s   = vc->a2s;
      scs   = vc->scs;
      n_seq = vc->n_seq;
      qbt1  = 1.;

      for (s = 0; s < n_seq; s++) {
        u1 = a2s[s][i - 1] + a2s[s][n] - a2s[s][j];
        char loopseq[10];
        loopseq[0] = '\0';
        if (u1 < 7) {
          strcpy(loopseq, Ss[s] + a2s[s][j] - 1);
          strncat(loopseq, Ss[s], a2s[s][i]);
        }

        type  = get_pair_type(SS[s][j], SS[s][i], md);
        qbt1  *= exp_E_Hairpin(u1, type, S3[s][j], S5[s][i], loopseq, P);
      }
      /* add soft constraints */
      if (scs) {
        for (s = 0; s < n_seq; s++) {
          if (scs[s]) {
            if (scs[s]->exp_energy_up)
              qbt1 *= ((i > 1) ? scs[s]->exp_energy_up[a2s[s][1]][a2s[s][i] - a2s[s][1]] : 1.) *
                      ((j < n) ? scs[s]->exp_energy_up[a2s[s][j] + 1][a2s[s][n] - a2s[s][j]] : 1.);

            if (scs[s]->exp_f) {
              qbt1 *= scs[s]->exp_f(a2s[s][j],
                                    a2s[s][i],
                                    a2s[s][j],
                                    a2s[s][i],
                                    VRNA_DECOMP_PAIR_HP,
                                    scs[s]->data);
            }
          }
        }
      }

      q = qbt1 * scale[u];

      break;

    default:
      break;
  }

  return q;
}
