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
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/structures.h"

#include "ViennaRNA/intern/grammar_dat.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/exterior_hc.inc"
#include "ViennaRNA/constraints/exterior_sc.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE int
eval_exterior_stem(vrna_fold_compound_t   *fc,
                   unsigned int           i,
                   unsigned int           j,
                   unsigned int           options,
                   vrna_hc_eval_f         evaluate,
                   struct hc_ext_def_dat  *hc_dat_local);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_E_exterior_stem(unsigned int type,
                     int          n5d,
                     int          n3d,
                     vrna_param_t *p)
{
  int energy = 0;

  if (n5d >= 0 && n3d >= 0)
    energy += p->mismatchExt[type][n5d][n3d];
  else if (n5d >= 0)
    energy += p->dangle5[type][n5d];
  else if (n3d >= 0)
    energy += p->dangle3[type][n3d];

  if (type > 2)
    energy += p->TerminalAU;

  return energy;
}


PUBLIC int
vrna_E_exterior_loop(unsigned int n,
                     vrna_md_t    *md)
{
  vrna_md_t md_tmp;

  if (md == NULL) {
    vrna_md_set_default(&md_tmp);
    md = &md_tmp;
  }

  if ((md->circ) &&
      (md->circ_penalty)) {
    double kT = md->betaScale * (md->temperature + K0) * GASCONST / 1000.;            /* kT in kcal/mol */
    return (int)(100. * kT * (md->circ_alpha0 + (3. / 2.) * log((double)n)) + 0.5f);  /* return in dekacal/mol */
  } else {
    return 0;
  }
}


PUBLIC int
vrna_eval_exterior_stem(vrna_fold_compound_t  *fc,
                        unsigned int          i,
                        unsigned int          j,
                        unsigned int          options)
{
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat hc_dat_local;

  if ((fc) &&
      (i > 0) &&
      (j > 0) &&
      (i < j)) {
    evaluate = NULL;

    if (!(options & VRNA_EVAL_LOOP_NO_HC)) {
      if (fc->hc->type == VRNA_HC_WINDOW)
        evaluate = prepare_hc_ext_def_window(fc, &hc_dat_local);
      else
        evaluate = prepare_hc_ext_def(fc, &hc_dat_local);
    }

    return eval_exterior_stem(fc, i, j, options, evaluate, &hc_dat_local);
  }

  return INF;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE int
eval_exterior_stem(vrna_fold_compound_t   *fc,
                   unsigned int           i,
                   unsigned int           j,
                   unsigned int           options,
                   vrna_hc_eval_f         evaluate,
                   struct hc_ext_def_dat  *hc_dat_local)
{
  unsigned int  eval;
  char          *ptype;
  short         *S;
  unsigned int  type;
  int           ij, en, e, *idx;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_sc_t     *sc;

  S     = fc->sequence_encoding;
  idx   = fc->jindx;
  ptype = fc->ptype;
  P     = fc->params;
  md    = &(P->model_details);
  sc    = fc->sc;

  e = INF;

  if (options & VRNA_EVAL_LOOP_NO_HC) {
    eval = (unsigned char)1;
  } else {
    eval = evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, hc_dat_local);
  }

  if (eval) {
    ij    = idx[j] + i;
    type  = vrna_get_ptype(ij, ptype);

    switch (md->dangles) {
      case 2:
        e = vrna_E_exterior_stem(type, S[i - 1], S[j + 1], P);
        break;

      case 0:
      /* fall through */

      default:
        e = vrna_E_exterior_stem(type, -1, -1, P);
        break;
    }

    if (e != INF) {
      if (!(options & VRNA_EVAL_LOOP_NO_SC)) {
        if (sc)
          if (sc->f)
            e += sc->f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
      }
    }
  }

  if (md->dangles % 2) {
    if (options & VRNA_EVAL_LOOP_NO_HC) {
      eval = (unsigned char)1;
    } else {
      eval = evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local);
    }

    if (eval) {
      ij    = idx[j - 1] + i;
      type  = vrna_get_ptype(ij, ptype);
      en    = vrna_E_exterior_stem(type, -1, S[j], P);

      if (en != INF) {
        if (!(options & VRNA_EVAL_LOOP_NO_SC)) {
          if (sc)
            if (sc->f)
              en += sc->f(i, j, i, j - 1, VRNA_DECOMP_EXT_STEM, sc->data);
        }

        e = MIN2(e, en);
      }
    }

    if (options & VRNA_EVAL_LOOP_NO_HC) {
      eval = (unsigned char)1;
    } else {
      eval = evaluate(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, hc_dat_local);
    }

    if (eval) {
      ij    = idx[j] + i + 1;
      type  = vrna_get_ptype(ij, ptype);
      en    = vrna_E_exterior_stem(type, S[i], -1, P);

      if (en != INF) {
        if (!(options & VRNA_EVAL_LOOP_NO_SC)) {
          if (sc)
            if (sc->f)
              en += sc->f(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, sc->data);
        }

        e = MIN2(e, en);
      }
    }

    if (options & VRNA_EVAL_LOOP_NO_HC) {
      eval = (unsigned char)1;
    } else {
      eval = evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local);
    }

    if (eval) {
      ij    = idx[j - 1] + i + 1;
      type  = vrna_get_ptype(ij, ptype);
      en    = vrna_E_exterior_stem(type, S[i], S[j], P);

      if (en != INF) {
        if (!(options & VRNA_EVAL_LOOP_NO_SC)) {
          if (sc)
            if (sc->f)
              en += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_EXT_STEM, sc->data);
        }

        e = MIN2(e, en);
      }
    }
  }

  return e;
}


/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC int
vrna_E_ext_stem(unsigned int  type,
                int           n5d,
                int           n3d,
                vrna_param_t  *p)
{
  return vrna_E_exterior_stem(type, n5d, n3d, p);
}


PUBLIC int
E_Stem(int          type,
       int          si1,
       int          sj1,
       int          extLoop,
       vrna_param_t *P)
{
  int energy  = 0;
  int d5      = (si1 >= 0) ? P->dangle5[type][si1] : 0;
  int d3      = (sj1 >= 0) ? P->dangle3[type][sj1] : 0;

  if (type > 2)
    energy += P->TerminalAU;

  if (si1 >= 0 && sj1 >= 0)
    energy += (extLoop) ? P->mismatchExt[type][si1][sj1] : P->mismatchM[type][si1][sj1];
  else
    energy += d5 + d3;

  if (!extLoop)
    energy += P->MLintern[type];

  return energy;
}


PUBLIC int
E_ExtLoop(int           type,
          int           si1,
          int           sj1,
          vrna_param_t  *P)
{
  return vrna_E_exterior_stem((unsigned int)type, si1, sj1, P);
}


PUBLIC int
vrna_eval_ext_stem(vrna_fold_compound_t *fc,
                   int                  i,
                   int                  j)
{
  return vrna_eval_exterior_stem(fc,
                                 (unsigned int)i,
                                 (unsigned int)j,
                                 VRNA_EVAL_LOOP_DEFAULT);
}


#endif
