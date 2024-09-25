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
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/eval/hairpin.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/hairpin_hc.inc"
#include "ViennaRNA/constraints/hairpin_sc_pf.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE FLT_OR_DBL
exp_eval_hp_loop(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 unsigned int         options);


PRIVATE FLT_OR_DBL
exp_eval_ext_hp_loop(vrna_fold_compound_t *fc,
                     unsigned int         i,
                     unsigned int         j,
                     unsigned int         options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC FLT_OR_DBL
vrna_exp_E_hairpin(unsigned int     u,
                   unsigned int     type,
                   int              si1,
                   int              sj1,
                   const char       *sequence,
                   vrna_exp_param_t *P)
{
  double q, kT, salt_correction;

  q = 0.;

  if (P) {
    kT              = P->kT; /* kT in cal/mol  */
    salt_correction = 1.;

    if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
      if (u <= MAXLOOP)
        salt_correction = P->expSaltLoop[u + 1];
      else
        salt_correction =
          exp(-vrna_salt_loop_int(u + 1, P->model_details.salt, P->temperature + K0,
                                  P->model_details.backbone_length) * 10. / kT);
    }

    if (u <= 30)
      q = P->exphairpin[u];
    else
      q = P->exphairpin[30] * exp(-(P->lxc * log(u / 30.)) * 10. / kT);

    q *= salt_correction;

    if (u < 3)
      return (FLT_OR_DBL)q;         /* should only be the case when folding alignments */

    if ((sequence) &&
        (P->model_details.special_hp)) {
      if (u == 4) {
        char tl[7] = {
          0
        }, *ts;
        memcpy(tl, sequence, sizeof(char) * 6);
        tl[6] = '\0';
        if ((ts = strstr(P->Tetraloops, tl))) {
          if (type != 7)
            return (FLT_OR_DBL)(P->exptetra[(ts - P->Tetraloops) / 7] * salt_correction);
          else
            q *= P->exptetra[(ts - P->Tetraloops) / 7];
        }
      } else if (u == 6) {
        char tl[9] = {
          0
        }, *ts;
        memcpy(tl, sequence, sizeof(char) * 8);
        tl[8] = '\0';
        if ((ts = strstr(P->Hexaloops, tl)))
          return (FLT_OR_DBL)(P->exphex[(ts - P->Hexaloops) / 9] * salt_correction);
      } else if (u == 3) {
        char tl[6] = {
          0
        }, *ts;
        memcpy(tl, sequence, sizeof(char) * 5);
        tl[5] = '\0';
        if ((ts = strstr(P->Triloops, tl)))
          return (FLT_OR_DBL)(P->exptri[(ts - P->Triloops) / 6] * salt_correction);

        if (type > 2)
          return (FLT_OR_DBL)(q * P->expTermAU);
        else
          return (FLT_OR_DBL)q;
      }
    }

    if ((si1 >= 0) &&
        (sj1 >= 0))
      q *= P->expmismatchH[type][si1][sj1];
  }

  return (FLT_OR_DBL)q;
}


PUBLIC FLT_OR_DBL
vrna_exp_eval_hairpin(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      unsigned int          options)
{
  unsigned char         eval;
  vrna_hc_eval_f        evaluate;
  struct hc_hp_def_dat  hc_dat_local;

  if ((fc) &&
      (i > 0) &&
      (j > 0)) {
    /* prepare hard constraints check */
    if ((options & VRNA_EVAL_LOOP_NO_HC) ||
        (fc->hc == NULL)) {
      eval = (unsigned char)1;
    } else {
      if (fc->hc->type == VRNA_HC_WINDOW)
        evaluate = prepare_hc_hp_def_window(fc, &hc_dat_local);
      else
        evaluate = prepare_hc_hp_def(fc, &hc_dat_local);

      eval = evaluate(i, j, i, j, VRNA_DECOMP_PAIR_HP, &hc_dat_local);
    }

    /* is this base pair allowed to close a hairpin (like) loop ? */
    if (eval)
      return (i > j) ?
             exp_eval_ext_hp_loop(fc, j, i, options) :
             exp_eval_hp_loop(fc, i, j, options);
  }

  return 0.;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE FLT_OR_DBL
exp_eval_hp_loop(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 unsigned int         options)
{
  char                  **Ss;
  short                 *S, *S2, **SS, **S5, **S3;
  unsigned int          **a2s, *sn, u, type, s, n_seq, noGUclosure;
  FLT_OR_DBL            q, qbt1, *scale;
  vrna_exp_param_t      *P;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  struct sc_hp_exp_dat  sc_wrapper;

  P           = fc->exp_params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  sn          = fc->strand_number;
  scale       = fc->exp_matrices->scale;
  domains_up  = fc->domains_up;

  q = 0.;

  if (sn[j] != sn[i])
    return 0;

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      SS    = fc->S;
      S5    = fc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
      S3    = fc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
      Ss    = fc->Ss;
      a2s   = fc->a2s;
      n_seq = fc->n_seq;
      qbt1  = 1.;

      for (s = 0; s < n_seq; s++) {
        if (a2s[s][i] < 1)
          continue;

        u     = a2s[s][j - 1] - a2s[s][i];
        type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
        qbt1  *= vrna_exp_E_hairpin(u, type, S3[s][i], S5[s][j], Ss[s] + a2s[s][i] - 1, P);
      }

      q = qbt1;
      break;

    default:
      S     = fc->sequence_encoding;
      S2    = fc->sequence_encoding2;
      u     = j - i - 1;
      type  = vrna_get_ptype_md(S2[i], S2[j], md);

      if ((noGUclosure) &&
          ((type == 3) || (type == 4)) &&
          (!(options & VRNA_EVAL_LOOP_NO_HC)))
        break;

      q = vrna_exp_E_hairpin(u, type, S[i + 1], S[j - 1], fc->sequence + i - 1, P);

      break;
  }

  /* add soft constraints */
  if (!(options & VRNA_EVAL_LOOP_NO_SC)) {
    init_sc_hp_exp(fc, &sc_wrapper);

    if (sc_wrapper.pair)
      q *= sc_wrapper.pair(i, j, &sc_wrapper);

    free_sc_hp_exp(&sc_wrapper);
  }

  if (domains_up && domains_up->exp_energy_cb) {
    /* we always consider both, bound and unbound state */
    q += q * domains_up->exp_energy_cb(fc,
                                       i + 1, j - 1,
                                       VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                       domains_up->data);
  }

  q *= scale[j - i + 1];

  return q;
}


PRIVATE FLT_OR_DBL
exp_eval_ext_hp_loop(vrna_fold_compound_t *fc,
                     unsigned int         i,
                     unsigned int         j,
                     unsigned int         options)
{
  char                  **Ss, *sequence, loopseq[10] = {
    0
  };
  unsigned int          **a2s, u1, u2, n, type, n_seq, s, noGUclosure;
  short                 *S, *S2, **SS, **S5, **S3;
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

  q   = 0.;
  u1  = n - j;
  u2  = i - 1;

  if ((u1 + u2) < 3)
    return q;

  switch (fc->type) {
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
        qbt1  *= vrna_exp_E_hairpin(u1_local + u2_local, type, S3[s][j], S5[s][i], loopseq, P);
      }

      q = qbt1;

      break;

    default:
      sequence  = fc->sequence;
      S         = fc->sequence_encoding;
      S2        = fc->sequence_encoding2;
      type      = vrna_get_ptype_md(S2[j], S2[i], md);

      if ((noGUclosure) &&
          ((type == 3) || (type == 4)) &&
          (!(options & VRNA_EVAL_LOOP_NO_HC)))
        break;

      /* get the loop sequence */
      if ((u1 + u2) < 7) {
        memcpy(loopseq, sequence + j - 1, sizeof(char) * (u1 + 1));
        memcpy(loopseq + u1 + 1, sequence, sizeof(char) * (u2 + 1));
        loopseq[u1 + u2 + 2] = '\0';
      }

      q = vrna_exp_E_hairpin(u1 + u2, type, S[j + 1], S[i - 1], loopseq, P);

      break;
  }

  /* add soft constraints */
  if (!(options & VRNA_EVAL_LOOP_NO_SC)) {
    init_sc_hp_exp(fc, &sc_wrapper);

    if (sc_wrapper.pair_ext)
      q *= sc_wrapper.pair_ext(i, j, &sc_wrapper);

    free_sc_hp_exp(&sc_wrapper);
  }

  if (domains_up && domains_up->exp_energy_cb) {
    /* we always consider both, bound and unbound state */
    q += q * domains_up->exp_energy_cb(fc,
                                       j + 1, i - 1,
                                       VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                       domains_up->data);
  }

  q *= scale[u1 + u2];

  return q;
}


/*
 #####################################
 # DEPRECATED functions below        #
 #####################################
 */
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC FLT_OR_DBL
vrna_exp_E_hp_loop(vrna_fold_compound_t *fc,
                   int                  i,
                   int                  j)
{
  return vrna_exp_eval_hairpin(fc,
                               (unsigned int)i,
                               (unsigned int)j,
                               VRNA_EVAL_LOOP_DEFAULT);
}


#endif
