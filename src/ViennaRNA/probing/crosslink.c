/* Crosslinking data handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/constants.h" /* defines MINPSCORE */
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/probing/basic.h"


#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/mfe/multibranch.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/multibranch_sc.inc"

#include "ViennaRNA/probing/crosslink.h"


/* re-definition of vrna_mx_mfe_aux_ml_s to use its internals, see also mfe_multibrnch.c  */
struct vrna_mx_mfe_aux_ml_s {
  unsigned int  n;
  int           *Fmi;
  int           *DMLi;
  int           *DMLi1;
  int           *DMLi2;
};


typedef struct {
  unsigned int  a;      /* start of 5' interval */
  unsigned int  b;      /* end of 5' interval */
  unsigned int  c;      /* start of 3' interval */
  unsigned int  d;      /* end of 3' interval */
  double        score;  /* interval score */
} interact_region_s;


struct vrna_crosslink_data_s {
  vrna_array(interact_region_s)   regions;  /* list of interacting regions */
  vrna_array(interact_region_s **) overlap5;
  vrna_array(interact_region_s **) overlap3;

  int *A;
  int *B;
};



PRIVATE INLINE int
ml_pair_d2(vrna_fold_compound_t         *fc,
           unsigned int                 i,
           unsigned int                 j,
           struct vrna_mx_mfe_aux_ml_s  *helpers,
           struct sc_mb_dat             *sc_wrapper);


PRIVATE INLINE int
mfe_M2_xlink(vrna_fold_compound_t         *fc,
             unsigned int                 i,
             unsigned int                 j,
             struct vrna_mx_mfe_aux_ml_s  *helpers,
             struct sc_mb_dat             *sc_wrapper);


PRIVATE void
fill_helper_arrays(vrna_fold_compound_t *fc,
                   unsigned int         i,
                   unsigned int         j,
                   struct sc_mb_dat     *sc_wrapper);


PUBLIC int
vrna_sc_xl_prepare(vrna_fold_compound_t  *fc)
{
  if ((fc) &&
      (fc->sc) &&
      (fc->sc->xlink)) {
    /* here we assume that sc->data actually points to a vrna_crosslink_data_s */
     
  }
  
  return 0;
}



PUBLIC int
vrna_mfe_multibranch_loop_xlink(vrna_fold_compound_t         *fc,
                                unsigned int                 i,
                                unsigned int                 j,
                                struct vrna_mx_mfe_aux_ml_s  *helpers)
{
  if ((fc) &&
      (helpers)) {

    unsigned char                 connects_intervals;
    int                           dangle_model, decomp, e;
    vrna_param_t                  *P;
    vrna_md_t                     *md;
    struct sc_mb_dat              sc_wrapper;
    struct vrna_crosslink_data_s  *dat;

    dat                 = (struct vrna_crosslink_data_s *)fc->sc->data;
    connects_intervals  = 0;
    P                   = fc->params;
    md                  = &(P->model_details);
    dangle_model        = md->dangles;
    e = decomp          = INF;

    /* first, check whether (i,j) connects two intervals of a segment */
    if (dat->overlap5[i]) {
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
        interact_region_s *r = dat->overlap5[i][cnt];
        if ((r->c <= j) &&
            (r->d >= j)) {
          connects_intervals = 1;
          break;
        }
      }

      if (connects_intervals) {
        init_sc_mb(fc, &sc_wrapper);

        switch (dangle_model) {
          /* currently, only dangles = 2 implemented */
          default:
            decomp = ml_pair_d2(fc, i, j, helpers, &sc_wrapper);
            break;
        }

        free_sc_mb(&sc_wrapper);

        return decomp;
      }
    }

    /* hack using old API to prevent calling vrna_mfe_multibranch_loop_xlink in a loop */
    return vrna_E_mb_loop_fast(fc,
                               i,
                               j,
                               helpers->DMLi1,
                               helpers->DMLi2);

  }

  return INF;
}


PRIVATE INLINE int
ml_pair_d2(vrna_fold_compound_t         *fc,
           unsigned int                 i,
           unsigned int                 j,
           struct vrna_mx_mfe_aux_ml_s  *helpers,
           struct sc_mb_dat             *sc_wrapper)
{
  short         *S, *S2, **SS, **S5, **S3, si1, sj1;
  unsigned int  tt, strands, *sn, s, n_seq;
  int           e;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  e   = INF;
  hc  = fc->hc;

  if (hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc)) {

    /* get contribution for enclosed part */
    e = mfe_M2_xlink(fc, i, j, helpers, sc_wrapper);

    if (e != INF) {
      P   = fc->params;
      md  = &(P->model_details);

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          strands = fc->strands;
          sn      = fc->strand_number;
          S       = fc->sequence_encoding;
          S2      = fc->sequence_encoding2;
          tt      = vrna_get_ptype_md(S2[j], S2[i], md);

          if (md->noGUclosure && ((tt == 3) || (tt == 4)))
            return INF; /* not allowed */

          si1 = ((strands == 1) || (sn[i] == sn[i + 1])) ? S[i + 1] : -1;
          sj1 = ((strands == 1) || (sn[j - 1] == sn[j])) ? S[j - 1] : -1;

          e += vrna_E_multibranch_stem(tt, sj1, si1, P) +
               P->MLclosing;
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = fc->n_seq;
          SS    = fc->S;
          S5    = fc->S5;
          S3    = fc->S3;

          for (s = 0; s < n_seq; s++) {
            tt  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
            e   += vrna_E_multibranch_stem(tt, S5[s][j], S3[s][i], P);
          }

          e += n_seq * P->MLclosing;
          break;
      }

      if (sc_wrapper->pair)
        e += sc_wrapper->pair(i, j, sc_wrapper);
    }
  }

  return e;
}


/*  compute the MFE contribution for the enclosed part of a multibranch loop
    [i + 1, j - 1] with closing pair (i,j), given that (i,j) connects an
    interval pair
*/
PRIVATE INLINE int
mfe_M2_xlink(vrna_fold_compound_t         *fc,
             unsigned int                 i,
             unsigned int                 j,
             struct vrna_mx_mfe_aux_ml_s  *helpers,
             struct sc_mb_dat             *sc_wrapper)
{
  unsigned int                  u;
  int                           e, decomp, *fML, *fM2, *A, *B, *idx;
  vrna_param_t                  *P;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;

  dat = (struct vrna_crosslink_data_s *)fc->sc->data;
  e   = INF;
  P   = fc->params;
  hc  = fc->hc;
  idx = fc->jindx;
  fML = fc->matrices->fML;
  fM2 = fc->matrices->fM2_real;
  A   = dat->A;
  B   = dat->B;

  for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
    interact_region_s *r = dat->overlap5[i][cnt];
    if ((r->c <= j) &&
        (r->d >= j)) {

      if (hc->eval_mb(i + 1, j - 1, r->b, r->b + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
        if ((fML[idx[r->b] + i + 1] != INF) &&
            (fML[idx[j - 1] + r->b + 1] != INF)) {
          decomp =  fML[idx[r->b] + i + 1] +
                    fML[idx[j - 1] + r->b + 1] +
                    r->score;

          if (sc_wrapper->decomp_ml)
            decomp += sc_wrapper->decomp_ml(i + 1, j - 1, r->b, r->b + 1, sc_wrapper);

          e = MIN2(e, decomp);
        }
      }

      if (hc->eval_mb(i + 1, j - 1, i + 1, r->b, VRNA_DECOMP_ML_ML, hc)) {
        if (fM2[idx[r->b] + i + 1] != INF) {
          decomp =  fM2[idx[r->b] + i + 1] +
                    (j - 1 - r->b) * P->MLbase +
                    r->score;

          if (sc_wrapper->red_ml)
            decomp += sc_wrapper->red_ml(i + 1, j - 1, i + 1, r->b, sc_wrapper);

          e = MIN2(e, decomp);
        }
      }

      if (hc->eval_mb(i + 1, j - 1, r->b + 1, j - 1, VRNA_DECOMP_ML_ML, hc)) {
        if (fM2[idx[j - 1] + r->b + 1] != INF) {
          decomp =  fM2[idx[j - 1] + r->b + 1] +
                    (r->b - i) * P->MLbase +
                    r->score;

          if (sc_wrapper->red_ml)
            decomp += sc_wrapper->red_ml(i + 1, j - 1, r->b + 1, j - 1, sc_wrapper);

          e = MIN2(e, decomp);
        }
      }

      /* now, all cases where we allow for base pairs that cross the first interval boundary */
      decomp = A[u];
      e = MIN2(e, decomp);

      for (u = i + 2; u <= r->b; u++) {
        if (hc->eval_mb(i + 1, j - 1, u, j - 1, VRNA_DECOMP_ML_ML, hc)) {
          decomp = A[u];
          if (decomp != INF) {
            decomp += (u - 1 - i) * P->MLbase;

            if (sc_wrapper->red_ml)
              decomp += sc_wrapper->red_ml(i + 1, j - 1, u, j - 1, sc_wrapper);

            e = MIN2(e, decomp);
          }
        }

        if (hc->eval_mb(i + 1, j - 1, u - 1, u, VRNA_DECOMP_ML_ML_ML, hc)) {
          decomp = A[u];

          if (decomp != INF) {
            decomp += fML[idx[u - 1] + i + 1];

            if (sc_wrapper->decomp_ml)
              decomp += sc_wrapper->decomp_ml(i + 1, j - 1, u - 1, u, sc_wrapper);

            e = MIN2(e, decomp);
          }

          decomp = B[u];

          if (decomp != INF) {
            decomp += fML[idx[u - 1] + i + 1];

            if (sc_wrapper->decomp_ml)
              decomp += sc_wrapper->decomp_ml(i + 1, j - 1, u - 1, u, sc_wrapper);

            e = MIN2(e, decomp);
          }
        }
      }
    }
  }

  return e;
}



PUBLIC void
vrna_sc_xl_update_mfe(vrna_fold_compound_t *fc,
                      unsigned int         i)
{
  short                         *S, *S2, si1, sj1;
  unsigned int                  n, stem, decomp, *sn, strands, tt;
  int                           *A, *B;
  vrna_param_t                  *P;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;

  dat = (struct vrna_crosslink_data_s *)fc->sc->data;
  n   = fc->length;
  P   = fc->params;
  hc  = fc->hc;
  A   = dat->A;
  B   = dat->B;

  for (unsigned int v = i; v <= n; v++) {
    /* reset A and B */
    A[v] = B[v] = INF;

    /* got through all pairing partners u */
    for (unsigned int u = i + 1; u < v; u++) {
      if (hc->eval_mb(u, v, u, v, VRNA_DECOMP_ML_STEM, hc)) {
        strands = fc->strands;
        sn      = fc->strand_number;
        S       = fc->sequence_encoding;
        S2      = fc->sequence_encoding2;
        tt      = vrna_get_ptype_md(S2[u], S2[v], md);

        si1 = ((strands == 1) || (sn[u - 1] == sn[u])) ? S[u - 1] : -1;
        sj1 = ((strands == 1) || (sn[v] == sn[v + 1])) ? S[v + 1] : -1;

        stem = vrna_E_multibranch_stem(tt, si1, sj1, P);

        decomp = 0;

        /* go through all intervals I' that cover v */
        for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[v]); cnt++) {
          interact_region_s *r = dat->overlap3[v][cnt];

          if ((u < r->c) &&
              (u > r->b))
            decomp += r->score;
        }

        decomp += stem;

        /* fill A */

        
        /* fill B */
        
      }
    }
  }
}


PRIVATE int
xlink_hp_en(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data)
{
  int                           energy = 0;
  struct vrna_crosslink_data_s  *dat = (struct vrna_crosslink_data_s *)data;

  /* go through each interval that overlaps with position i to check whether j connects both segments */
  if (dat->overlap5[i])
    for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
      interact_region_s *r = dat->overlap5[i][cnt];
      if ((r->b < j) &&
          (r->c <= j) &&
          (r->d >= j))
      energy += r->score;
    }

  return energy;

}


PRIVATE int
xlink_il_en(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data)
{
  int                           energy = 0;
  struct vrna_crosslink_data_s  *dat = (struct vrna_crosslink_data_s *)data;

  /*  go through each interval that overlaps with position i to check whether
   *  k does not overlap with the 5' interval and j connects both segments */
  if (dat->overlap5[i])
    for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
      interact_region_s *r = dat->overlap5[i][cnt];
      if ((r->b < k) &&
          (r->c <= j) &&
          (r->d >= j))
      energy += r->score;
    }

  return energy;

}
