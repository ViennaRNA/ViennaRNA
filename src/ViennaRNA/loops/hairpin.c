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
#include "hairpin_sc.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
eval_hp_loop_fake(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

/**
 *  @brief  Evaluate the free energy of a hairpin loop
 *          and consider possible hard constraints
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 */
PUBLIC int
vrna_E_hp_loop(vrna_fold_compound_t *fc,
               int                  i,
               int                  j)
{
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  if (fc->hc->type == VRNA_HC_WINDOW)
    evaluate = prepare_hc_default_window(fc, &hc_dat_local);
  else
    evaluate = prepare_hc_default(fc, &hc_dat_local);

  if ((i > 0) && (j > 0)) {
    /* is this base pair allowed to close a hairpin (like) loop ? */
    if (evaluate(i, j, i, j, VRNA_DECOMP_PAIR_HP, &hc_dat_local)) {
      if (j > i)  /* linear case */
        return vrna_eval_hp_loop(fc, i, j);
      else        /* circular case */
        return vrna_eval_ext_hp_loop(fc, j, i);
    }
  }

  return INF;
}


/**
 *  @brief  Evaluate the free energy of an exterior hairpin loop
 *          and consider possible hard constraints
 */
PUBLIC int
vrna_E_ext_hp_loop(vrna_fold_compound_t *fc,
                   int                  i,
                   int                  j)
{
  return vrna_E_hp_loop(fc, j, i);
}


/**
 *  @brief Evaluate free energy of an exterior hairpin loop
 *
 *  @ingroup eval
 *
 */
PUBLIC int
vrna_eval_ext_hp_loop(vrna_fold_compound_t  *fc,
                      int                   i,
                      int                   j)
{
  char                  **Ss, loopseq[10] = {
    0
  };
  unsigned int          **a2s;
  short                 *S, *S2, **SS, **S5, **S3;
  int                   u1, u2, e, s, type, n_seq, length, noGUclosure;
  vrna_param_t          *P;
  vrna_md_t             *md;
  struct sc_wrapper_hp  sc_wrapper;

  length      = fc->length;
  P           = fc->params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  e           = INF;

  init_sc_wrapper(fc, &sc_wrapper);

  u1  = length - j;
  u2  = i - 1;

  if ((u1 + u2) < 3)
    return e;

  switch (fc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      S2    = fc->sequence_encoding2;
      type  = vrna_get_ptype_md(S2[j], S2[i], md);

      if (noGUclosure && ((type == 3) || (type == 4)))
        break;

      /* maximum special hp loop size: 6 */
      if ((u1 + u2) < 7) {
        memcpy(loopseq, fc->sequence + j - 1, sizeof(char) * (u1 + 1));
        memcpy(loopseq + u1 + 1, fc->sequence, sizeof(char) * (u2 + 1));
        loopseq[u1 + u2 + 2] = '\0';
      }

      e = E_Hairpin(u1 + u2, type, S[j + 1], S[i - 1], loopseq, P);

      break;

    /* sequence alignments */
    case  VRNA_FC_TYPE_COMPARATIVE:
      SS    = fc->S;
      S5    = fc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
      S3    = fc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
      Ss    = fc->Ss;
      a2s   = fc->a2s;
      n_seq = fc->n_seq;
      e     = 0;

      for (s = 0; s < n_seq; s++) {
        u1  = a2s[s][length] - a2s[s][j];
        u2  = a2s[s][i - 1];
        memset(loopseq, '\0', sizeof(loopseq));

        if ((u1 + u2) < 7) {
          memcpy(loopseq, Ss[s] + a2s[s][j] - 1, sizeof(char) * (u1 + 1));
          memcpy(loopseq + u1 + 1, Ss[s], sizeof(char) * (u2 + 1));
          loopseq[u1 + u2 + 2] = '\0';
        }

        if ((u1 + u2) < 3) {
          e += 600;
        } else {
          type  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
          e     += E_Hairpin(u1 + u2, type, S3[s][j], S5[s][i], loopseq, P);
        }
      }

      break;

    /* nothing */
    default:
      break;
  }

  if (e != INF)
    if (sc_wrapper.pair_ext)
      e += sc_wrapper.pair_ext(i, j, &sc_wrapper);

  free_sc_wrapper(&sc_wrapper);

  return e;
}


/**
 *  @brief Evaluate free energy of a hairpin loop
 *
 *  @ingroup eval
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @param  fc  The #vrna_fold_compound_t for the particular energy evaluation
 *  @param  i   5'-position of the base pair
 *  @param  j   3'-position of the base pair
 *  @returns    Free energy of the hairpin loop closed by @f$ (i,j) @f$ in deka-kal/mol
 */
PUBLIC int
vrna_eval_hp_loop(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j)
{
  char                  **Ss;
  unsigned int          **a2s;
  short                 *S, *S2, **SS, **S5, **S3;
  unsigned int          *sn;
  int                   u, e, s, type, n_seq, en, noGUclosure;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  struct sc_wrapper_hp  sc_wrapper;

  P           = fc->params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  sn          = fc->strand_number;
  domains_up  = fc->domains_up;
  e           = INF;

  if (sn[j] != sn[i])
    return eval_hp_loop_fake(fc, i, j);

  init_sc_wrapper(fc, &sc_wrapper);

  /* regular hairpin loop */
  switch (fc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      S2    = fc->sequence_encoding2;
      u     = j - i - 1;
      type  = vrna_get_ptype_md(S2[i], S2[j], md);

      if (noGUclosure && ((type == 3) || (type == 4)))
        break;

      e = E_Hairpin(u, type, S[i + 1], S[j - 1], fc->sequence + i - 1, P);

      break;

    /* sequence alignments */
    case  VRNA_FC_TYPE_COMPARATIVE:
      SS    = fc->S;
      S5    = fc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
      S3    = fc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
      Ss    = fc->Ss;
      a2s   = fc->a2s;
      n_seq = fc->n_seq;

      for (e = s = 0; s < n_seq; s++) {
        u = a2s[s][j - 1] - a2s[s][i];
        if (u < 3) {
          e += 600;                          /* ??? really 600 ??? */
        } else {
          type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
          e     += E_Hairpin(u, type, S3[s][i], S5[s][j], Ss[s] + (a2s[s][i - 1]), P);
        }
      }

      break;

    /* nothing */
    default:
      break;
  }

  if (e != INF) {
    if (sc_wrapper.pair)
      e += sc_wrapper.pair(i, j, &sc_wrapper);

    /* consider possible ligand binding */
    if (domains_up && domains_up->energy_cb) {
      en = domains_up->energy_cb(fc,
                                 i + 1, j - 1,
                                 VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                 domains_up->data);
      if (en != INF)
        en += e;

      e = MIN2(e, en);
    }
  }

  free_sc_wrapper(&sc_wrapper);

  return e;
}


PRIVATE int
eval_hp_loop_fake(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j)
{
  short         *S, *S2;
  unsigned int  *sn;
  int           u, e, ij, type, *idx, en, noGUclosure;
  vrna_param_t  *P;
  vrna_sc_t     *sc;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;

  idx         = fc->jindx;
  P           = fc->params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  sn          = fc->strand_number;
  domains_up  = fc->domains_up;
  e           = INF;

  switch (fc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      S2    = fc->sequence_encoding2;
      sc    = fc->sc;
      u     = j - i - 1;
      ij    = idx[j] + i;
      type  = vrna_get_ptype_md(S2[j], S2[i], md);

      if (noGUclosure && ((type == 3) || (type == 4)))
        break;

      /* hairpin-like exterior loop (for cofolding) */
      short si, sj;

      si  = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;
      sj  = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;

      if (md->dangles)
        e = vrna_E_ext_stem(type, sj, si, P);
      else
        e = vrna_E_ext_stem(type, -1, -1, P);

      /* add soft constraints */
      if (sc) {
        if (sc->energy_up)
          e += sc->energy_up[i + 1][u];

        if (sc->energy_bp)
          e += sc->energy_bp[ij];

        if (sc->f)
          e += sc->f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
      }

      /* consider possible ligand binding */
      if (domains_up && domains_up->energy_cb) {
        en = domains_up->energy_cb(fc,
                                   i + 1, j - 1,
                                   VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                   domains_up->data);
        if (en != INF)
          en += e;

        e = MIN2(e, en);
      }

      break;

    /* nothing */
    default:
      break;
  }

  return e;
}
