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

struct default_data {
  int                       n;
  int                       *idx;
  char                      *mx;
  int                       cp;
  int                       *hc_up;
  void                      *hc_dat;
  vrna_callback_hc_evaluate *hc_f;
};


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


PRIVATE char
hc_default(int  i,
           int  j,
           int  k,
           int  l,
           char d,
           void *data);


PRIVATE char
hc_default_user(int   i,
                int   j,
                int   k,
                int   l,
                char  d,
                void  *data);


PRIVATE int
eval_hp_loop_fake(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j);


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
 *  @brief  Evaluate the free energy of a hairpin loop
 *          and consider possible hard constraints
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 */
PUBLIC int
vrna_E_hp_loop(vrna_fold_compound_t *vc,
               int                  i,
               int                  j)
{
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = vc->hc->matrix;
  hc_dat_local.hc_up  = vc->hc->up_hp;
  hc_dat_local.n      = vc->length;
  hc_dat_local.cp     = vc->cutpoint;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if ((i > 0) && (j > 0)) {
    /* is this base pair allowed to close a hairpin (like) loop ? */
    if (evaluate(i, j, i, j, VRNA_DECOMP_PAIR_HP, &hc_dat_local)) {
      if (j > i)  /* linear case */
        return vrna_eval_hp_loop(vc, i, j);
      else        /* circular case */
        return vrna_eval_ext_hp_loop(vc, j, i);
    }
  }

  return INF;
}


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

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = vc->hc->matrix;
  hc_dat_local.hc_up  = vc->hc->up_hp;
  hc_dat_local.n      = vc->length;
  hc_dat_local.cp     = vc->cutpoint;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

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


/**
 *  @brief  Evaluate the free energy of an exterior hairpin loop
 *          and consider possible hard constraints
 */
PUBLIC int
vrna_E_ext_hp_loop(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j)
{
  return vrna_E_hp_loop(vc, j, i);
}


/**
 *  @brief Evaluate free energy of an exterior hairpin loop
 *
 *  @ingroup eval
 *
 */
PUBLIC int
vrna_eval_ext_hp_loop(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j)
{
  char            **Ss, loopseq[10];
  unsigned short  **a2s;
  short           *S, **SS, **S5, **S3;
  int             u, e, s, type, *types, n_seq, length;
  vrna_param_t    *P;
  vrna_sc_t       *sc, **scs;
  vrna_md_t       *md;

  length  = vc->length;
  P       = vc->params;
  md      = &(P->model_details);
  e       = INF;

  switch (vc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      sc    = vc->sc;
      u     = vc->length - j + i - 1;
      type  = md->pair[S[j]][S[i]];

      if (type == 0)
        type = 7;

      if (u < 7) {
        strcpy(loopseq, vc->sequence + j - 1);
        strncat(loopseq, vc->sequence, i);
      }

      e = E_Hairpin(u, type, S[j + 1], S[i - 1], loopseq, P);

      if (sc) {
        if (sc->energy_up)
          e += sc->energy_up[j + 1][vc->length - j]
               + sc->energy_up[1][i - 1];

        if (sc->f)
          e += sc->f(j, i, j, i, VRNA_DECOMP_PAIR_HP, sc->data);
      }
      break;

    /* sequence alignments */
    case  VRNA_FC_TYPE_COMPARATIVE:
      SS    = vc->S;
      S5    = vc->S5;                                 /*S5[s][i] holds next base 5' of i in sequence s*/
      S3    = vc->S3;                                 /*Sl[s][i] holds next base 3' of i in sequence s*/
      Ss    = vc->Ss;
      a2s   = vc->a2s;
      scs   = vc->scs;
      n_seq = vc->n_seq;
      e     = 0;
      types = (int *)vrna_alloc(sizeof(int) * n_seq);

      for (s = 0; s < n_seq; s++) {
        types[s] = md->pair[SS[s][j]][SS[s][i]];
        if (types[s] == 0) types[s] = 7;
      }

      for (s = 0; s < n_seq; s++) {
        char loopseq[10];
        u = a2s[s][length] - a2s[s][j] + a2s[s][i - 1];

        if (u < 9) {
          strcpy(loopseq, Ss[s] + a2s[s][j] - 1);
          strncat(loopseq, Ss[s], a2s[s][i]);
        }
        if (u < 3) e += 600;
        else e += E_Hairpin(u, types[s], S3[s][j], S5[s][i], loopseq, P);
      }
      if (scs) {
        for (s = 0; s < n_seq; s++) {
          if (scs[s]) {
            if (scs[s]->energy_up)
              e += ((i > 1) ? scs[s]->energy_up[1][a2s[s][i - 1]] : 0)
                   + ((j < length) ? scs[s]->energy_up[a2s[s][j + 1]][a2s[s][length] - a2s[s][j]] : 0);
            if (scs[s]->f)
              e += scs[s]->f(a2s[s][j], a2s[s][i], a2s[s][j], a2s[s][i], VRNA_DECOMP_PAIR_HP, scs[s]->data);
          }
        }
      }

      free(types);
      break;

    /* nothing */
    default:
      break;
  }

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
 *  @param  vc  The #vrna_fold_compound_t for the particular energy evaluation
 *  @param  i   5'-position of the base pair
 *  @param  j   3'-position of the base pair
 *  @returns    Free energy of the hairpin loop closed by @f$ (i,j) @f$ in deka-kal/mol
 */
PUBLIC int
vrna_eval_hp_loop(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j)
{
  char            **Ss;
  unsigned short  **a2s;
  short           *S, **SS, **S5, **S3;
  unsigned int    *sn;
  int             u, e, s, ij, type, *types, *idx, n_seq, en;
  vrna_param_t    *P;
  vrna_sc_t       *sc, **scs;
  vrna_md_t       *md;
  vrna_ud_t       *domains_up;

  idx         = vc->jindx;
  P           = vc->params;
  md          = &(P->model_details);
  sn          = vc->strand_number;
  domains_up  = vc->domains_up;
  e           = INF;

  if (sn[j] != sn[i])
    return eval_hp_loop_fake(vc, i, j);

  /* regular hairpin loop */
  switch (vc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      sc    = vc->sc;
      u     = j - i - 1;
      ij    = idx[j] + i;
      type  = md->pair[S[i]][S[j]];

      if (type == 0)
        type = 7;

      e = E_Hairpin(u, type, S[i + 1], S[j - 1], vc->sequence + i - 1, P);

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
        en = domains_up->energy_cb(vc,
                                   i + 1, j - 1,
                                   VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                   domains_up->data);
        if (en != INF)
          en += e;
        e = MIN2(e, en);
      }

      break;

    /* sequence alignments */
    case  VRNA_FC_TYPE_COMPARATIVE:
      SS    = vc->S;
      S5    = vc->S5;                                 /*S5[s][i] holds next base 5' of i in sequence s*/
      S3    = vc->S3;                                 /*Sl[s][i] holds next base 3' of i in sequence s*/
      Ss    = vc->Ss;
      a2s   = vc->a2s;
      scs   = vc->scs;
      n_seq = vc->n_seq;
      ij    = idx[j] + i;
      types = (int *)vrna_alloc(sizeof(int) * n_seq);

      for (s = 0; s < n_seq; s++) {
        types[s] = md->pair[SS[s][i]][SS[s][j]];
        if (types[s] == 0) types[s] = 7;
      }

      for (e = s = 0; s < n_seq; s++) {
        u = a2s[s][j - 1] - a2s[s][i];
        e += (u < 3) ? 600 : E_Hairpin(u, types[s], S3[s][i], S5[s][j], Ss[s] + (a2s[s][i - 1]), P);                          /* ??? really 600 ??? */
      }

      if (scs) {
        for (s = 0; s < n_seq; s++) {
          if (scs[s]) {
            u = a2s[s][j - 1] - a2s[s][i];

            if (scs[s]->energy_up)
              e += scs[s]->energy_up[a2s[s][i + 1]][u];

            if (scs[s]->energy_bp)
              e += scs[s]->energy_bp[ij];

            if (scs[s]->f)
              e += scs[s]->f(a2s[s][i], a2s[s][j], a2s[s][i], a2s[s][j], VRNA_DECOMP_PAIR_HP, scs[s]->data);
          }
        }
      }

      free(types);
      break;

    /* nothing */
    default:
      break;
  }

  return e;
}


PRIVATE int
eval_hp_loop_fake(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j)
{
  short         *S;
  unsigned int  *sn;
  int           u, e, ij, type, *idx, en;
  vrna_param_t  *P;
  vrna_sc_t     *sc;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;

  idx         = vc->jindx;
  P           = vc->params;
  md          = &(P->model_details);
  sn          = vc->strand_number;
  domains_up  = vc->domains_up;
  e           = INF;

  switch (vc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      sc    = vc->sc;
      u     = j - i - 1;
      ij    = idx[j] + i;
      type  = md->pair[S[i]][S[j]];

      if (type == 0)
        type = 7;

      /* hairpin-like exterior loop (for cofolding) */
      short si, sj;

      si  = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;
      sj  = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;

      if (md->dangles)
        e = E_ExtLoop(md->rtype[type], sj, si, P);
      else
        e = E_ExtLoop(md->rtype[type], -1, -1, P);

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
        en = domains_up->energy_cb(vc,
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


PRIVATE FLT_OR_DBL
exp_eval_hp_loop_fake(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j)
{
  short             *S, s5, s3;
  unsigned int      *sn;
  int               u, cp, type, *iidx;
  FLT_OR_DBL        qq, temp, *q, *scale;
  vrna_exp_param_t  *pf_params;
  vrna_sc_t         *sc;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;

  cp          = vc->cutpoint;
  iidx        = vc->iindx;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  q           = vc->exp_matrices->q;
  scale       = vc->exp_matrices->scale;
  sn          = vc->strand_number;
  domains_up  = vc->domains_up;

  qq = 0;

  switch (vc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      sc    = vc->sc;
      u     = j - i - 1;
      type  = md->pair[S[j]][S[i]];

      if (type == 0)
        type = 7;

      temp = q[iidx[i + 1] - (cp - 1)] * q[iidx[cp] - (j - 1)];
      if ((j == cp) && (i == cp - 1))
        temp = scale[2];
      else if (i == cp - 1)
        temp = q[iidx[cp] - (j - 1)] * scale[1];
      else if (j == cp)
        temp = q[iidx[i + 1] - (cp - 1)] * scale[1];
      if (j > cp)
        temp *= scale[1];
      if (i < cp - 1)
        temp *= scale[1];

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


/*
 *************************************
 * Partition function variants below *
 *************************************
 */
PRIVATE FLT_OR_DBL
exp_eval_hp_loop(vrna_fold_compound_t *vc,
                 int                  i,
                 int                  j)
{
  char              **Ss;
  unsigned short    **a2s;
  short             *S, **SS, **S5, **S3;
  unsigned int      *sn;
  int               u, ij, type, n_seq, s, *types, *idx, *iidx;
  FLT_OR_DBL        q, qbt1, *scale;
  vrna_exp_param_t  *P;
  vrna_sc_t         *sc, **scs;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;

  idx         = vc->jindx;
  iidx        = vc->iindx;
  P           = vc->exp_params;
  md          = &(P->model_details);
  sn          = vc->strand_number;
  scale       = vc->exp_matrices->scale;
  types       = NULL;
  domains_up  = vc->domains_up;

  q   = 0.;
  ij  = idx[j] + i;

  if (sn[j] != sn[i])
    return exp_eval_hp_loop_fake(vc, i, j);

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      sc    = vc->sc;
      u     = j - i - 1;
      type  = vc->ptype[ij];

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

        if (sc->exp_energy_bp)
          q *= sc->exp_energy_bp[iidx[i] - j];

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
      S5    = vc->S5;                                 /*S5[s][i] holds next base 5' of i in sequence s*/
      S3    = vc->S3;                                 /*Sl[s][i] holds next base 3' of i in sequence s*/
      Ss    = vc->Ss;
      a2s   = vc->a2s;
      scs   = vc->scs;
      n_seq = vc->n_seq;
      qbt1  = 1.;
      types = (int *)vrna_alloc(sizeof(int) * n_seq);

      for (s = 0; s < n_seq; s++) {
        types[s] = md->pair[SS[s][i]][SS[s][j]];
        if (types[s] == 0) types[s] = 7;
      }

      for (s = 0; s < n_seq; s++) {
        u = a2s[s][j - 1] - a2s[s][i];
        if (a2s[s][i] < 1) continue;
        char loopseq[10];
        if (u < 9)
          strncpy(loopseq, Ss[s] + a2s[s][i] - 1, 10);
        qbt1 *= exp_E_Hairpin(u, types[s], S3[s][i], S5[s][j], loopseq, P);
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

            if (scs[s]->exp_f)
              qbt1 *= scs[s]->exp_f(a2s[s][i], a2s[s][j], a2s[s][i], a2s[s][j], VRNA_DECOMP_PAIR_HP, scs[s]->data);
          }
        }
      }

      q = qbt1 * scale[j - i + 1];
      break;

    default:
      break;
  }

  free(types);
  return q;
}


PRIVATE FLT_OR_DBL
exp_eval_ext_hp_loop(vrna_fold_compound_t *vc,
                     int                  i,
                     int                  j)
{
  char              **Ss, *sequence;
  unsigned short    **a2s;
  short             *S, **SS, **S5, **S3;
  int               u, u1, ij, n, type, n_seq, s, *rtype, *types, *idx, noGUclosure;
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
  types       = NULL;
  domains_up  = vc->domains_up;
  rtype       = &(md->rtype[0]);

  q   = 0.;
  u   = n - j + i - 1;
  ij  = idx[j] + i;

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
      if (u < 7) {
        strcpy(loopseq, sequence + j - 1);
        strncat(loopseq, sequence, i);
      }

      q = exp_E_Hairpin(u, type, S[j + 1], S[i - 1], loopseq, P);

      /* add soft constraints */
      if (sc) {
        if (sc->exp_energy_up)
          q *= ((i > 1) ? sc->exp_energy_up[1][i - 1] : 1.)
               * ((j < n) ? sc->exp_energy_up[j + 1][n - j] : 1.);

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
      S5    = vc->S5;                                 /*S5[s][i] holds next base 5' of i in sequence s*/
      S3    = vc->S3;                                 /*Sl[s][i] holds next base 3' of i in sequence s*/
      Ss    = vc->Ss;
      a2s   = vc->a2s;
      scs   = vc->scs;
      n_seq = vc->n_seq;
      qbt1  = 1.;
      types = (int *)vrna_alloc(sizeof(int) * n_seq);

      for (s = 0; s < n_seq; s++) {
        types[s] = md->pair[SS[s][j]][SS[s][i]];
        if (types[s] == 0) types[s] = 7;
      }

      for (s = 0; s < n_seq; s++) {
        u1 = a2s[s][i] - 1 + a2s[s][n] - a2s[s][j];
        char loopseq[10];
        if (u1 < 7) {
          strcpy(loopseq, Ss[s] + a2s[s][j] - 1);
          strncat(loopseq, Ss[s], a2s[s][i]);
        }
        qbt1 *= exp_E_Hairpin(u1, types[s], S3[s][j], S5[s][i], loopseq, P);
      }

      /* add soft constraints */
      if (scs) {
        for (s = 0; s < n_seq; s++) {
          if (scs[s]) {
            if (scs[s]->exp_energy_up)
              qbt1 *= ((i > 1) ? scs[s]->exp_energy_up[a2s[s][1]][a2s[s][i] - a2s[s][1]] : 1.)
                      * ((j < n) ? scs[s]->exp_energy_up[a2s[s][j] + 1][a2s[s][n] - a2s[s][j]] : 1.);

            if (scs[s]->exp_f)
              qbt1 *= scs[s]->exp_f(a2s[s][j], a2s[s][i], a2s[s][j], a2s[s][i], VRNA_DECOMP_PAIR_HP, scs[s]->data);
          }
        }
      }

      q = qbt1 * scale[u];

      free(types);
      break;

    default:
      break;
  }

  return q;
}


/**
 *  @brief Backtrack a hairpin loop closed by @f$ (i,j) @f$
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 */
PUBLIC int
vrna_BT_hp_loop(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j,
                int                   en,
                vrna_bp_stack_t       *bp_stack,
                int                   *stack_count)
{
  int       e, u;
  vrna_sc_t *sc;

  sc = NULL;

  u = j - i - 1;

  if (vc->hc->up_hp[i + 1] < u)
    return 0;

  e = vrna_E_hp_loop(vc, i, j);

  if (e == en) {
    switch (vc->type) {
      case  VRNA_FC_TYPE_SINGLE:
        sc = vc->sc;
        break;

      case  VRNA_FC_TYPE_COMPARATIVE:
        if (vc->scs)
          sc = vc->scs[0];
        break;

      default:
        break;
    }

    if (sc) {
      if (sc->bt) {
        vrna_basepair_t *ptr, *aux_bps;
        aux_bps = sc->bt(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
        for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
          bp_stack[++(*stack_count)].i  = ptr->i;
          bp_stack[(*stack_count)].j    = ptr->j;
        }
        free(aux_bps);
      }
    }

    return 1;
  }

  return 0;
}


PRIVATE char
hc_default(int  i,
           int  j,
           int  k,
           int  l,
           char d,
           void *data)
{
  int                 ij, u, p, q;
  char                eval;
  struct default_data *dat = (struct default_data *)data;

  eval = (char)0;

  if (j > i) {
    /* linear case */
    p = i;
    q = j;
    u = q - p - 1;
  } else {
    /* circular case */
    p = j;
    q = i;
    u = dat->n - q + p - 1;
  }

  ij = dat->idx[q] + p;
  if (dat->mx[ij] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
    eval = (char)1;
    if (dat->hc_up[i + 1] < u)
      eval = (char)0;
  }

  return eval;
}


PRIVATE char
hc_default_user(int   i,
                int   j,
                int   k,
                int   l,
                char  d,
                void  *data)
{
  char                eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = hc_default(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (char)0;

  return eval;
}
