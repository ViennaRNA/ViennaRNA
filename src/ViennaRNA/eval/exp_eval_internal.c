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
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/internal.h"


#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/internal_hc.inc"
#include "ViennaRNA/constraints/internal_sc_pf.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE FLT_OR_DBL
exp_eval_internal(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j,
                  unsigned int          k,
                  unsigned int          l,
                  unsigned int          options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC FLT_OR_DBL
vrna_exp_E_internal(unsigned int      u1,
                    unsigned int      u2,
                    unsigned int      type,
                    unsigned int      type2,
                    int               si1,
                    int               sj1,
                    int               sp1,
                    int               sq1,
                    vrna_exp_param_t  *P)
{
  unsigned int  ul, us, backbones, no_close = 0;
  double        z, salt_stack_correction;
  double        salt_loop_correction = 1.;

  z = 0.;

  if (P) {
    no_close              = 0;
    salt_stack_correction = P->expSaltStack;
    salt_loop_correction  = 1.;

    if ((P->model_details.noGUclosure) &&
        ((type2 == 3) || (type2 == 4) || (type == 3) || (type == 4)))
      no_close = 1;

    if (u1 > u2) {
      ul  = u1;
      us  = u2;
    } else {
      ul  = u2;
      us  = u1;
    }

    if (ul == 0) /* stack */
      return (FLT_OR_DBL)P->expstack[type][type2] *
             salt_stack_correction;

    if (no_close)
      return 0.;

    if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
      /* salt correction for loop */
      backbones = ul + us + 2;
      if (backbones <= MAXLOOP + 1) {
        salt_loop_correction = P->expSaltLoop[backbones];
      } else {
        int E = vrna_salt_loop_int(backbones,
                                   P->model_details.salt,
                                   P->temperature + K0,
                                   P->model_details.backbone_length);

        salt_loop_correction = exp(-(double)E * 10. / P->kT);
      }
    }

    z = 1.;

    switch (us) {
      case 0:
        /* bulge */
        z = P->expbulge[ul];

        if (ul == 1) {
          z *= P->expstack[type][type2];
        } else {
          /* add stacking energy for 1-bulges */
          if (type > 2)
            z *= P->expTermAU;

          if (type2 > 2)
            z *= P->expTermAU;
        }

        break;

      case 1:
        if (ul == 1) {
          /* 1x1 loop */
          z = P->expint11[type][type2][si1][sj1];
        } else if (ul == 2) {
          /* 2x1 loop */
          if (u1 == 1)
            z = P->expint21[type][type2][si1][sq1][sj1];
          else
            z = P->expint21[type2][type][sq1][si1][sp1];
        } else {
          /* 1xn loop */
          z = P->expinternal[ul + us];
          z *= P->expninio[2][ul - us];
          z *= P->expmismatch1nI[type][si1][sj1] *
               P->expmismatch1nI[type2][sq1][sp1];
        }

        break;

      case 2:
        if (ul == 2) {
          /* 2x2 loop */
          z = P->expint22[type][type2][si1][sp1][sq1][sj1];
          break;
        } else if (ul == 3) {
          /* 2x3 loop */
          z = P->expinternal[5] *
              P->expninio[2][1];
          z *= P->expmismatch23I[type][si1][sj1] *
               P->expmismatch23I[type2][sq1][sp1];
          break;
        }

      /* fall through */

      default:
        /* generic internal loop */
        z = P->expinternal[ul + us];
        z *= P->expninio[2][ul - us];
        z *= P->expmismatchI[type][si1][sj1] *
             P->expmismatchI[type2][sq1][sp1];

        break;
    }

    z *= salt_loop_correction;
  }

  return (FLT_OR_DBL)z;
}


PUBLIC FLT_OR_DBL
vrna_exp_eval_internal(vrna_fold_compound_t *fc,
                       unsigned int         i,
                       unsigned int         j,
                       unsigned int         k,
                       unsigned int         l,
                       unsigned int         options)
{
  unsigned char         eval;
  eval_hc               evaluate;
  struct hc_int_def_dat hc_dat_local;

  if ((fc) &&
      (i > 0) &&
      (j > 0) &&
      (k > 0) &&
      (l > 0)) {
    /* prepare hard constraints check */
    if ((options & VRNA_EVAL_LOOP_NO_HC) ||
        (fc->hc == NULL)) {
      eval = (unsigned char)1;
    } else {
      evaluate  = prepare_hc_int_def(fc, &hc_dat_local);
      eval      = evaluate(i, j, k, l, &hc_dat_local);
    }

    if (eval)
      return exp_eval_internal(fc, i, j, k, l, options);
  }

  return 0.;
}


PRIVATE FLT_OR_DBL
exp_eval_internal(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j,
                  unsigned int          k,
                  unsigned int          l,
                  unsigned int          options)
{
  unsigned char         sliding_window, type, type2;
  char                  *ptype, **ptype_local;
  short                 *S1, **SS, **S5, **S3;
  unsigned int          u1, u2, *sn, n_seq, s, **a2s;
  int                   *rtype, *jindx;
  FLT_OR_DBL            qbt1, q_temp, *scale;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  struct sc_int_exp_dat sc_wrapper;

  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     =
    (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? fc->ptype_local : NULL) : NULL;
  S1          = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s         = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  jindx       = fc->jindx;
  pf_params   = fc->exp_params;
  sn          = fc->strand_number;
  md          = &(pf_params->model_details);
  scale       = fc->exp_matrices->scale;
  domains_up  = fc->domains_up;
  rtype       = &(md->rtype[0]);
  qbt1        = 0.;
  u1          = k - i - 1;
  u2          = j - l - 1;

  if ((sn[k] != sn[i]) ||
      (sn[j] != sn[l]))
    return qbt1;

  /* discard this configuration if (p,q) is not allowed to be enclosed pair of an internal loop */
  q_temp = 0;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      type = (sliding_window) ?
             vrna_get_ptype_window(i, j, ptype_local) :
             vrna_get_ptype(jindx[j] + i, ptype);
      type2 = (sliding_window) ?
              rtype[vrna_get_ptype_window(k, l, ptype_local)] :
              rtype[vrna_get_ptype(jindx[l] + k, ptype)];

      q_temp = vrna_exp_E_internal(u1,
                                   u2,
                                   type,
                                   type2,
                                   S1[i + 1],
                                   S1[j - 1],
                                   S1[k - 1],
                                   S1[l + 1],
                                   pf_params);

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      q_temp = 1.;

      for (s = 0; s < n_seq; s++) {
        unsigned int  u1_local  = a2s[s][k - 1] - a2s[s][i];
        unsigned int  u2_local  = a2s[s][j - 1] - a2s[s][l];
        type    = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
        type2   = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
        q_temp  *= vrna_exp_E_internal(u1_local,
                                       u2_local,
                                       type,
                                       type2,
                                       S3[s][i],
                                       S5[s][j],
                                       S5[s][k],
                                       S3[s][l],
                                       pf_params);
      }

      break;
  }

  if (q_temp > 0.) {
    if (!(options & VRNA_EVAL_LOOP_NO_SC)) {
      init_sc_int_exp(fc, &sc_wrapper);

      /* soft constraints */
      if (sc_wrapper.pair)
        q_temp *= sc_wrapper.pair(i, j, k, l, &sc_wrapper);

      free_sc_int_exp(&sc_wrapper);
    }

    /* unstructured domains */
    if (domains_up && domains_up->exp_energy_cb) {
      FLT_OR_DBL qq5, qq3;

      qq5 = qq3 = 0.;

      if (u1 > 0) {
        qq5 = domains_up->exp_energy_cb(fc,
                                        i + 1, k - 1,
                                        VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                        domains_up->data);
      }

      if (u2 > 0) {
        qq3 = domains_up->exp_energy_cb(fc,
                                        l + 1, j - 1,
                                        VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                        domains_up->data);
      }

      qbt1 += q_temp *
              qq5 *
              scale[u1 + u2 + 2];      /* only motifs in 5' part */
      qbt1 += q_temp *
              qq3 *
              scale[u1 + u2 + 2];      /* only motifs in 3' part */
      qbt1 += q_temp *
              qq5 *
              qq3 *
              scale[u1 + u2 + 2]; /* motifs in both parts */
    }

    qbt1 += q_temp *
            scale[u1 + u2 + 2];
  }

  return qbt1;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC FLT_OR_DBL
vrna_exp_E_interior_loop(vrna_fold_compound_t *fc,
                         int                  i,
                         int                  j,
                         int                  k,
                         int                  l)
{
  return vrna_exp_eval_internal(fc,
                                (unsigned int)i,
                                (unsigned int)j,
                                (unsigned int)k,
                                (unsigned int)l,
                                VRNA_EVAL_LOOP_DEFAULT);
}


#endif
