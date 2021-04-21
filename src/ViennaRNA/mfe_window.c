/*
 *                minimum free energy
 *                RNA secondary structure prediction
 *                with sliding window approach
 *
 *                c Ivo Hofacker, Peter Stadler, Ronny Lorenz
 *
 *                ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/constants.h" /* defines MINPSCORE */
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/ribo.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/utils/units.h"
#include "ViennaRNA/mfe_window.h"

#ifdef VRNA_WITH_SVM
#include "ViennaRNA/zscore_dat.inc"
#endif

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define MAXSECTORS                  500   /* dimension for a backtrack array */
#define INT_CLOSE_TO_UNDERFLOW(i)   ((i) <= (INT_MIN / 16))
#define UNDERFLOW_CORRECTION        (((INT_MIN / 32) / 100) * 100)

#define NONE -10000 /* score for forbidden pairs */


typedef struct {
  FILE  *output;
  int   dangle_model;
  int   csv;
} hit_data;


struct aux_arrays {
  int *cc;    /* auxilary arrays for canonical structures     */
  int *cc1;   /* auxilary arrays for canonical structures     */
  int *Fmi;   /* holds row i of fML (avoids jumps in memory)  */
  int *DMLi;  /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  int *DMLi1; /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  int *DMLi2; /*                MIN(fML[i+2,k]+fML[k+1,j])    */
};

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void
make_ptypes(vrna_fold_compound_t  *vc,
            int                   i);


PRIVATE char *
backtrack(vrna_fold_compound_t  *vc,
          int                   start,
          int                   maxdist);


PRIVATE int
fill_arrays(vrna_fold_compound_t            *vc,
            int                             *underflow,
            vrna_mfe_window_callback        *cb,
#ifdef VRNA_WITH_SVM
            vrna_mfe_window_zscore_callback *cb_z,
#endif
            void                            *data);


PRIVATE void
default_callback(int        start,
                 int        end,
                 const char *structure,
                 float      en,
                 void       *data);


PRIVATE double
cov_score(vrna_fold_compound_t  *fc,
          int                   i,
          int                   j,
          float                 **dm);


PRIVATE void
make_pscores(vrna_fold_compound_t *fc,
             int                  start,
             float                **dm);


PRIVATE void
default_callback_comparative(int        start,
                             int        end,
                             const char *structure,
                             float      en,
                             void       *data);


PRIVATE INLINE void
allocate_dp_matrices(vrna_fold_compound_t *fc);


PRIVATE INLINE void
free_dp_matrices(vrna_fold_compound_t *fc);


PRIVATE INLINE void
rotate_dp_matrices(vrna_fold_compound_t *fc,
                   int                  i);


PRIVATE INLINE void
init_constraints(vrna_fold_compound_t *fc,
                 float                **dm);


PRIVATE INLINE void
rotate_constraints(vrna_fold_compound_t *fc,
                   float                **dm,
                   int                  i);


PRIVATE INLINE struct aux_arrays *
get_aux_arrays(unsigned int maxdist);


PRIVATE INLINE void
rotate_aux_arrays(struct aux_arrays *aux,
                  unsigned int      maxdist);


PRIVATE INLINE void
free_aux_arrays(struct aux_arrays *aux);


PRIVATE INLINE int
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               struct aux_arrays    *aux);


#ifdef VRNA_WITH_SVM

PRIVATE void
default_callback_z(int        start,
                   int        end,
                   const char *structure,
                   float      en,
                   float      zscore,
                   void       *data);


PRIVATE INLINE int
want_backtrack(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               double               *z);


#endif


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC float
vrna_mfe_window(vrna_fold_compound_t  *vc,
                FILE                  *file)
{
  hit_data data;

  data.output       = (file) ? file : stdout;
  data.dangle_model = vc->params->model_details.dangles;
  data.csv          = 0; /* csv output is for backward-compatibility only */

  if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
    return vrna_mfe_window_cb(vc, &default_callback_comparative, (void *)&data);
  else
    return vrna_mfe_window_cb(vc, &default_callback, (void *)&data);
}


PUBLIC float
vrna_mfe_window_cb(vrna_fold_compound_t     *vc,
                   vrna_mfe_window_callback *cb,
                   void                     *data)
{
  int   energy, underflow, n_seq;
  float mfe_local, e_factor;

  /* keep track of how many times we were close to an integer underflow */
  underflow = 0;

  if (!vrna_fold_compound_prepare(vc, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW)) {
    vrna_message_warning("vrna_mfe_window@Lfold.c: Failed to prepare vrna_fold_compound");
    return (float)(INF / 100.);
  }

  n_seq     = (vc->type == VRNA_FC_TYPE_COMPARATIVE) ? vc->n_seq : 1;
  e_factor  = 100. * n_seq;

#ifdef VRNA_WITH_SVM
  energy = fill_arrays(vc, &underflow, cb, NULL, data);
#else
  energy = fill_arrays(vc, &underflow, cb, data);
#endif
  mfe_local = (underflow > 0) ? ((float)underflow * (float)(UNDERFLOW_CORRECTION)) / e_factor : 0.;
  mfe_local += (float)energy / e_factor;

  return mfe_local;
}


#ifdef VRNA_WITH_SVM

PUBLIC float
vrna_mfe_window_zscore(vrna_fold_compound_t *vc,
                       double               min_z,
                       FILE                 *file)
{
  hit_data data;

  data.output       = (file) ? file : stdout;
  data.dangle_model = vc->params->model_details.dangles;

  return vrna_mfe_window_zscore_cb(vc, min_z, &default_callback_z, (void *)&data);
}


PUBLIC float
vrna_mfe_window_zscore_cb(vrna_fold_compound_t            *vc,
                          double                          min_z,
                          vrna_mfe_window_zscore_callback *cb_z,
                          void                            *data)
{
  int   energy, underflow;
  float mfe_local;

  if (vc->type == VRNA_FC_TYPE_COMPARATIVE) {
    vrna_message_warning(
      "vrna_mfe_window_zscore@mfe_window.c: Comparative prediction not implemented");
    return (float)(INF / 100.);
  }

  if (!vrna_fold_compound_prepare(vc, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW)) {
    vrna_message_warning("vrna_mfe_window@Lfold.c: Failed to prepare vrna_fold_compound");
    return (float)(INF / 100.);
  }

  vrna_zsc_filter_update(vc, min_z, VRNA_ZSCORE_OPTIONS_NONE);

  /* keep track of how many times we were close to an integer underflow */
  underflow = 0;

  energy = fill_arrays(vc, &underflow, NULL, cb_z, data);

  mfe_local = (underflow > 0) ? ((float)underflow * (float)(UNDERFLOW_CORRECTION)) / 100. : 0.;
  mfe_local += (float)energy / 100.;

  return mfe_local;
}


#endif


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE void
allocate_dp_matrices(vrna_fold_compound_t *fc)
{
  int       i, j, length, maxdist, **c, **fML;
  vrna_hc_t *hc;
  vrna_sc_t *sc;

  length  = fc->length;
  maxdist = MIN2(fc->window_size, length);
  hc      = fc->hc;
  c       = fc->matrices->c_local;
  fML     = fc->matrices->fML_local;

  /* reserve additional memory for j-dimension */
  for (i = length; (i > length - maxdist - 5) && (i >= 0); i--) {
    c[i]                = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
    fML[i]              = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
    hc->matrix_local[i] = (unsigned char *)vrna_alloc(sizeof(unsigned char) * (maxdist + 5));
    if (fc->type == VRNA_FC_TYPE_SINGLE)
      fc->ptype_local[i] = vrna_alloc(sizeof(char) * (maxdist + 5));
    else if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      fc->pscore_local[i] = vrna_alloc(sizeof(int) * (maxdist + 5));
  }

  /*
   *  allocate one more entry for comparative predictions to allow for
   *  access to i - 1 when processing [i ... maxdist]. This is required
   *  for (default) hard constraints with noLP option
   */
  if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
    if (length > maxdist + 5)
      fc->pscore_local[length - maxdist - 5] = vrna_alloc(sizeof(int) * (maxdist + 5));

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;
      if (sc) {
        if (sc->energy_bp_local)
          for (i = length; (i > length - maxdist - 5) && (i >= 0); i--)
            sc->energy_bp_local[i] = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));

        if (sc->energy_up)
          for (i = length; (i > length - maxdist - 5) && (i >= 0); i--)
            sc->energy_up[i] = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));

        for (i = length; (i > length - maxdist - 5) && (i >= 0); i--)
          vrna_sc_update(fc, i, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      break;
  }

  if (fc->type == VRNA_FC_TYPE_SINGLE)
    for (j = length; j > length - maxdist - 4; j--)
      for (i = (length - maxdist - 4 > 0) ? length - maxdist - 4 : 1; i < j; i++)
        c[i][j - i] = fML[i][j - i] = INF;
  else if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
    for (j = length; j > length - maxdist - 3; j--)
      for (i = (length - maxdist - 2 > 0) ? length - maxdist - 2 : 1; i < j; i++)
        c[i][j - i] = fML[i][j - i] = INF;
}


PRIVATE INLINE void
free_dp_matrices(vrna_fold_compound_t *fc)
{
  int       i, length, maxdist, **c, **fML, **ggg, with_gquad;
  vrna_hc_t *hc;
  vrna_sc_t *sc;

  length      = fc->length;
  maxdist     = MIN2(fc->window_size, length);
  hc          = fc->hc;
  c           = fc->matrices->c_local;
  fML         = fc->matrices->fML_local;
  ggg         = fc->matrices->ggg_local;
  with_gquad  = fc->params->model_details.gquad;


  /* free additional memory for j-dimension */
  for (i = 0; (i < maxdist + 5) && (i <= length); i++) {
    if (fc->type == VRNA_FC_TYPE_SINGLE) {
      free(fc->ptype_local[i]);
      fc->ptype_local[i] = NULL;
    } else if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      free(fc->pscore_local[i]);
      fc->pscore_local[i] = NULL;
    }

    free(c[i]);
    c[i] = NULL;
    free(fML[i]);
    fML[i] = NULL;
    free(hc->matrix_local[i]);
    hc->matrix_local[i] = NULL;
  }

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;
      if (sc) {
        if (sc->energy_up) {
          for (i = 0; (i < maxdist + 5) && (i <= length); i++) {
            free(sc->energy_up[i]);
            sc->energy_up[i] = NULL;
          }
        }

        if (sc->energy_bp_local) {
          for (i = 0; (i < maxdist + 5) && (i <= length); i++) {
            free(sc->energy_bp_local[i]);
            sc->energy_bp_local[i] = NULL;
          }
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      break;
  }

  if (with_gquad) {
    for (i = 0; (i <= maxdist + 5) && (i <= length); i++)
      free(ggg[i]);
    free(ggg);
    fc->matrices->ggg_local = NULL;
  }
}


PRIVATE INLINE void
rotate_dp_matrices(vrna_fold_compound_t *fc,
                   int                  i)
{
  int       ii, maxdist, length, **c, **fML;
  vrna_hc_t *hc;
  vrna_sc_t *sc;

  length  = fc->length;
  maxdist = fc->window_size;
  c       = fc->matrices->c_local;
  fML     = fc->matrices->fML_local;
  hc      = fc->hc;

  if (i + maxdist + 4 <= length) {
    c[i - 1]                          = c[i + maxdist + 4];
    c[i + maxdist + 4]                = NULL;
    fML[i - 1]                        = fML[i + maxdist + 4];
    fML[i + maxdist + 4]              = NULL;
    hc->matrix_local[i - 1]           = hc->matrix_local[i + maxdist + 4];
    hc->matrix_local[i + maxdist + 4] = NULL;

    /* rotate base pair soft constraints */
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        sc = fc->sc;
        if (sc) {
          if (sc->energy_bp_local) {
            sc->energy_bp_local[i - 1]            = sc->energy_bp_local[i + maxdist + 4];
            sc->energy_bp_local[i + maxdist + 4]  = NULL;
          }

          if (sc->energy_up) {
            sc->energy_up[i - 1]            = sc->energy_up[i + maxdist + 4];
            sc->energy_up[i + maxdist + 4]  = NULL;
          }
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        break;
    }

    if ((fc->params->model_details.gquad) && (i > 1))
      vrna_gquad_mx_local_update(fc, i - 1);

    for (ii = 0; ii < maxdist + 5; ii++) {
      c[i - 1][ii]    = INF;
      fML[i - 1][ii]  = INF;
    }
  }
}


PRIVATE INLINE void
init_constraints(vrna_fold_compound_t *fc,
                 float                **dm)
{
  int i, length, maxdist;

  length  = fc->length;
  maxdist = fc->window_size;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      for (i = length; (i >= length - maxdist - 4) && (i > 0); i--) {
        make_ptypes(fc, i);
        vrna_hc_update(fc, i, VRNA_CONSTRAINT_WINDOW_UPDATE_3);
        vrna_sc_update(fc, i, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);
      }
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      for (i = length; (i >= length - maxdist - 4) && (i > 0); i--) {
        make_pscores(fc, i, dm);
        vrna_hc_update(fc, i, VRNA_CONSTRAINT_WINDOW_UPDATE_3);
      }

      /* for noLP option */
      if (length > maxdist + 5)
        make_pscores(fc, length - maxdist - 5, dm);

      break;
  }
}


PRIVATE INLINE void
rotate_constraints(vrna_fold_compound_t *fc,
                   float                **dm,
                   int                  i)
{
  int length, maxdist;

  length  = fc->length;
  maxdist = fc->window_size;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      if (i + maxdist + 4 <= length) {
        fc->ptype_local[i - 1]            = fc->ptype_local[i + maxdist + 4];
        fc->ptype_local[i + maxdist + 4]  = NULL;
        if (i > 1) {
          make_ptypes(fc, i - 1);
          vrna_hc_update(fc, i - 1, VRNA_CONSTRAINT_WINDOW_UPDATE_3);
          vrna_sc_update(fc, i - 1, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if (i + maxdist + 4 <= length) {
        if (i > 1) {
          fc->pscore_local[i - 2]           = fc->pscore_local[i + maxdist + 4];
          fc->pscore_local[i + maxdist + 4] = NULL;
          if (i > 2)
            make_pscores(fc, i - 2, dm);

          vrna_hc_update(fc, i - 1, VRNA_CONSTRAINT_WINDOW_UPDATE_3);
        } else if (i == 1) {
          free(fc->pscore_local[i - 1]);
          fc->pscore_local[i - 1]           = fc->pscore_local[i + maxdist + 4];
          fc->pscore_local[i + maxdist + 4] = NULL;
        }
      }

      break;
  }
}


PRIVATE int
fill_arrays(vrna_fold_compound_t            *vc,
            int                             *underflow,
            vrna_mfe_window_callback        *cb,
#ifdef VRNA_WITH_SVM
            vrna_mfe_window_zscore_callback *cb_z,
#endif
            void                            *data)
{
  /* fill "c", "fML" and "f3" arrays and return  optimal energy */

  char              *prev;
  int               i, j, length, maxdist, **c, **fML, *f3,
                    with_gquad, dangle_model, turn, n_seq,
                    prev_i, prev_j, prev_end, prev_en;
  float             **dm;
  double            e_fact;

#ifdef VRNA_WITH_SVM
  double            prevz;
  unsigned char     with_zscore;
  unsigned char     report_subsumed;
  vrna_zsc_dat_t    zsc_data;
#endif

  vrna_md_t         *md;
  struct aux_arrays *helper_arrays;

  int               olddm[7][7] = { { 0, 0, 0, 0, 0, 0, 0 },/* hamming distance between pairs */
                                    { 0, 0, 2, 2, 1, 2, 2 } /* CG */,
                                    { 0, 2, 0, 1, 2, 2, 2 } /* GC */,
                                    { 0, 2, 1, 0, 2, 1, 2 } /* GU */,
                                    { 0, 1, 2, 2, 0, 2, 1 } /* UG */,
                                    { 0, 2, 2, 1, 2, 0, 2 } /* AU */,
                                    { 0, 2, 2, 2, 1, 2, 0 } /* UA */ };

  n_seq         = (vc->type == VRNA_FC_TYPE_COMPARATIVE) ? vc->n_seq : 1;
  length        = vc->length;
  maxdist       = vc->window_size;
  md            = &(vc->params->model_details);
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  turn          = md->min_loop_size;
  do_backtrack  = 0;
  prev_i        = 0;
  prev_j        = 0;
  prev_end      = 0;
  prev          = NULL;
  prev_en       = 0;
  e_fact        = 100 * n_seq;
  dm            = NULL;
#ifdef VRNA_WITH_SVM
  prevz           = 0.;
  zsc_data        = vc->zscore_data;
  with_zscore     = (zsc_data) ? zsc_data->filter_on : 0;
  report_subsumed = (zsc_data) ? zsc_data->report_subsumed : 0;
#endif

  if (vc->type == VRNA_FC_TYPE_COMPARATIVE) {
#ifdef VRNA_WITH_SVM
    /* no z-scoring for comparative structure prediction */
    if (zsc_data) {
      vrna_zsc_filter_free(vc);
      zsc_data        = NULL;
      with_zscore     = 0;
      report_subsumed = 0;
    }
#endif

    if (md->ribo) {
      if (RibosumFile != NULL)
        dm = readribosum(RibosumFile);
      else
        dm = get_ribosum((const char **)vc->sequences, n_seq, length);
    } else {
      /*use usual matrix*/
      dm = (float **)vrna_alloc(7 * sizeof(float *));
      for (i = 0; i < 7; i++) {
        dm[i] = (float *)vrna_alloc(7 * sizeof(float));
        for (j = 0; j < 7; j++)
          dm[i][j] = (float)olddm[i][j];
      }
    }
  }

  c   = vc->matrices->c_local;
  fML = vc->matrices->fML_local;
  f3  = vc->matrices->f3_local;

  /* allocate memory for all helper arrays */
  helper_arrays = get_aux_arrays(maxdist);

  /* reserve additional memory for j-dimension */
  allocate_dp_matrices(vc);

  init_constraints(vc, dm);

  if (with_gquad)
    vrna_gquad_mx_local_update(vc, length - maxdist - 4);

  for (i = length - turn - 1; i >= 1; i--) {
    /* i,j in [1..length] */
    for (j = i + 1; j <= MIN2(i + turn, length); j++)
      c[i][j - i] = fML[i][j - i] = INF;

    for (j = i + turn + 1; j <= length && j <= i + maxdist; j++) {
      /* decompose subsegment [i, j] with pair (i, j) */
      c[i][j - i] = decompose_pair(vc, i, j, helper_arrays);

      /* decompose subsegment [i, j] that is multibranch loop part with at least one branch */
      fML[i][j - i] = vrna_E_ml_stems_fast(vc, i, j, helper_arrays->Fmi, helper_arrays->DMLi);

      if ((vc->aux_grammar) && (vc->aux_grammar->cb_aux))
        vc->aux_grammar->cb_aux(vc, i, j, vc->aux_grammar->data);
    } /* for (j...) */

    /* calculate energies of 5' and 3' fragments */
    f3[i] = vrna_E_ext_loop_3(vc, i);

    {
      char *ss = NULL;

      if (f3[i] < f3[i + 1]) {
        /*
         * instead of backtracing in the next iteration, we backtrack now
         * already. This is necessary to accomodate for change in free
         * energy due to unpaired nucleotides in the exterior loop, which
         * may happen in the case of using soft constraints
         */
        int ii, jj;
        ii  = i;
        jj  = vrna_BT_ext_loop_f3_pp(vc, &ii, maxdist);
        if (jj > 0) {
#ifdef VRNA_WITH_SVM
          double thisz = 0;
          if (want_backtrack(vc, ii, jj, &thisz)) {
#endif
          ss = backtrack(vc, ii, jj);
          if (prev) {
            if ((jj < prev_j) ||
#ifdef VRNA_WITH_SVM
                ((report_subsumed) && (prevz < thisz)) || /* yield last structure if it's z-score is higher than the current one */
#endif
                (strncmp(ss + prev_i - ii, prev, prev_j - prev_i + 1))) {
              /* ss does not contain prev */
#ifdef VRNA_WITH_SVM
              if (with_zscore)
                cb_z(prev_i, prev_end, prev, prev_en / e_fact, prevz, data);
              else
#endif
              cb(prev_i, prev_end, prev, prev_en / e_fact, data);
            }

            free(prev);
          }

          prev      = ss;
          prev_i    = ii;
          prev_j    = jj;
          prev_end  = MIN2(jj + ((dangle_model) ? 1 : 0), length);
          prev_en   = f3[ii] - f3[jj + 1];
#ifdef VRNA_WITH_SVM
          prevz = thisz;
        }

#endif
        } else if (jj == -1) {
          /* some error occured during backtracking */
          vrna_message_error("backtrack failed in short backtrack 1");
        }
      }

      if (i == 1) {
        if (prev) {
#ifdef VRNA_WITH_SVM
          if (with_zscore)
            cb_z(prev_i, prev_end, prev, prev_en / e_fact, prevz, data);
          else
#endif
          cb(prev_i, prev_end, prev, prev_en / e_fact, data);

          free(prev);
          prev = NULL;
#ifdef VRNA_WITH_SVM
        } else if ((f3[i] < 0) && (!with_zscore)) {
          /* why !with_zscore? */
#else
        } else if (f3[i] < 0) {
#endif
          int ii, jj;
          ii  = i;
          jj  = vrna_BT_ext_loop_f3_pp(vc, &ii, maxdist);
          if (jj > 0) {
#ifdef VRNA_WITH_SVM
            double thisz = 0;
            if (want_backtrack(vc, ii, jj, &thisz)) {
#endif
            ss = backtrack(vc, ii, jj);
#ifdef VRNA_WITH_SVM
            if (with_zscore)
              cb_z(ii, MIN2(jj + ((dangle_model) ? 1 : 0),
                            length), ss, (f3[1] - f3[jj + 1]) / e_fact, thisz, data);
            else
#endif
            cb(ii, MIN2(jj + ((dangle_model) ? 1 : 0),
                        length), ss, (f3[1] - f3[jj + 1]) / e_fact, data);

            free(ss);
#ifdef VRNA_WITH_SVM
          }

#endif
          } else if (jj == -1) {
            /* some error occured during backtracking */
            vrna_message_error("backtrack failed in short backtrack 2");
          }
        }
      }
    }

    /* check for values close to integer underflow */
    if (INT_CLOSE_TO_UNDERFLOW(f3[i])) {
      /* correct f3 free energies and increase underflow counter */
      int cnt;
      for (cnt = i; cnt <= MIN2(i + maxdist + 2, length); cnt++)
        f3[cnt] -= UNDERFLOW_CORRECTION;
      (*underflow)++;
    }

    rotate_aux_arrays(helper_arrays, maxdist);
    rotate_dp_matrices(vc, i);
    rotate_constraints(vc, dm, i);
  }

  /* clean up memory */
  free_aux_arrays(helper_arrays);
  free_dp_matrices(vc);

  if (dm) {
    for (i = 0; i < 7; i++)
      free(dm[i]);
    free(dm);
  }

  return f3[1];
}


#ifdef VRNA_WITH_SVM
PRIVATE INLINE int
want_backtrack(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               double               *z)
{
  int bt, *f3;

  bt  = 1; /* we want to backtrack by default */
  *z  = (double)INF;

  if ((fc->zscore_data) &&
      (fc->zscore_data->filter_on)) {
    if ((fc->zscore_data->pre_filter) &&
        (fc->zscore_data->current_i == i)) {
      *z = fc->zscore_data->current_z[j];
    } else {
      bt  = 0;
      f3  = fc->matrices->f3_local;
      *z  = vrna_zsc_compute(fc, i, j, f3[i] - f3[j + 1]);
    }

    if ((*z) <= fc->zscore_data->min_z)
      bt = 1;
  }

  return bt;
}


#endif


PRIVATE char *
backtrack(vrna_fold_compound_t  *vc,
          int                   start,
          int                   end)
{
  /*------------------------------------------------------------------
   *  trace back through the "c", "f3" and "fML" arrays to get the
   *  base pairing list. No search for equivalent structures is done.
   *  This is fast, since only few structure elements are recalculated.
   *  ------------------------------------------------------------------*/
  sect sector[MAXSECTORS];            /* backtracking sectors */
  char *structure, **ptype;
  int i, j, k, length, no_close, type, s, b, bt_type, turn,
      dangle_model, noLP, noGUclosure, **c, dangle3, ml, cij,
      **pscore, canonical, p, q, comp1, comp2, max3;
  vrna_param_t *P;
  vrna_md_t *md;
  vrna_bp_stack_t *bp_stack;

  length        = vc->length;
  ptype         = vc->ptype_local;
  pscore        = vc->pscore_local;
  P             = vc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  noLP          = md->noLP;
  noGUclosure   = md->noGUclosure;
  bt_type       = md->backtrack_type;
  turn          = md->min_loop_size;
  c             = vc->matrices->c_local;

  s         = 0;                                                                                /* depth of backtracking stack */
  b         = 0;                                                                                /* number of base pairs */
  bp_stack  = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + length / 2)));  /* add a guess of how many G's may be involved in a G quadruplex */

  sector[++s].i = start;
  sector[s].j   = MIN2(length, end);
  sector[s].ml  = (bt_type == 'M') ? 1 : ((bt_type == 'C') ? 2 : 0);

  structure = (char *)vrna_alloc((MIN2(length - start, end) + 3) * sizeof(char));

  memset(structure, '.', MIN2(length - start, end) + 1);

  dangle3 = 0;

  while (s > 0) {
    canonical = 1;     /* (i,j) closes a canonical structure */

    /* pop one element from stack */
    i   = sector[s].i;
    j   = sector[s].j;
    ml  = sector[s--].ml;  /* ml is a flag indicating if backtracking is to
                           * occur in the fML- (1) or in the f-array (0) */

    if (j < i + turn + 1)
      continue;                     /* no more pairs in this interval */

    switch (ml) {
      /* backtrack in f3 */
      case 0:
        if (vrna_BT_ext_loop_f3(vc, &i, j, &p, &q, bp_stack, &b)) {
          if (i > 0) {
            sector[++s].i = i;
            sector[s].j   = j;
            sector[s].ml  = 0;
          }

          if (p > 0) {
            if (((i == q + 2) || (dangle_model)) && (q < length))
              dangle3 = 1;

            i = p;
            j = q;
            goto repeat1;
          } else if (md->gquad) {
            /*
             * check whether last element on the bp_stack is involved in G-Quadruplex formation
             * and increase output dot-bracket string length by 1 if necessary
             */
            if ((bp_stack[b].i == bp_stack[b].j) && (bp_stack[b].i < length))
              dangle3 = 1;
          }

          continue;
        } else {
          vrna_message_error("backtracking failed in f3, segment [%d,%d]\n", i, j);
        }

        break;

      /* trace back in fML array */
      case 1:
        if (vrna_BT_mb_loop_split(vc, &i, &j, &p, &q, &comp1, &comp2, bp_stack, &b)) {
          if (i > 0) {
            sector[++s].i = i;
            sector[s].j   = j;
            sector[s].ml  = comp1;
          }

          if (p > 0) {
            sector[++s].i = p;
            sector[s].j   = q;
            sector[s].ml  = comp2;
          }

          continue;
        } else {
          vrna_message_error("backtracking failed in fML, segment [%d,%d]\n", i, j);
        }

        break;

      /* backtrack in c */
      case 2:
        bp_stack[++b].i = i;
        bp_stack[b].j   = j;
        goto repeat1;
        break;

      default:
        vrna_message_error("Backtracking failed due to unrecognized DP matrix!");
        break;
    }

repeat1:

    /*----- begin of "repeat:" -----*/
    if (canonical)
      cij = c[i][j - i];

    if (noLP) {
      if (vrna_BT_stack(vc, &i, &j, &cij, bp_stack, &b)) {
        canonical = 0;
        goto repeat1;
      }
    }

    canonical = 1;

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type = vrna_get_ptype_window(i, j, ptype);

        no_close = (((type == 3) || (type == 4)) && noGUclosure);

        if (no_close) {
          if (cij == FORBIDDEN)
            continue;
        } else {
          if (vrna_BT_hp_loop(vc, i, j, cij, bp_stack, &b))
            continue;
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        cij += pscore[i][j - i];
        if (vrna_BT_hp_loop(vc, i, j, cij, bp_stack, &b))
          continue;

        break;
    }

    if (vrna_BT_int_loop(vc, &i, &j, cij, bp_stack, &b)) {
      if (i < 0)
        continue;
      else
        goto repeat1;
    }

    /* (i.j) must close a multi-loop */
    if (vrna_BT_mb_loop(vc, &i, &j, &k, cij, &comp1, &comp2)) {
      sector[++s].i = i;
      sector[s].j   = k;
      sector[s].ml  = comp1;
      sector[++s].i = k + 1;
      sector[s].j   = j;
      sector[s].ml  = comp2;
    } else {
      vrna_message_error("backtracking failed in repeat, segment [%d,%d]\n", i, j);
    }

    /* end of repeat: --------------------------------------------------*/
  } /* end of infinite while loop */


  bp_stack[0].i = b;

  /* and now create a dot-brakcet string from the base pair stack... */
  max3 = 1;
  for (i = 1; i <= b; i++) {
    if (bp_stack[i].i == bp_stack[i].j) {
      /* Gquad bonds are marked as bp[i].i == bp[i].j */
      structure[bp_stack[i].i - start] = '+';
    } else {
      /* the following ones are regular base pairs */
      structure[bp_stack[i].i - start]  = '(';
      structure[bp_stack[i].j - start]  = ')';
    }

    if (max3 < bp_stack[i].j - start)
      max3 = bp_stack[i].j - start;
  }

  free(bp_stack);

  structure = (char *)vrna_realloc(structure,
                                   sizeof(char) * (max3 + dangle3 + 2));
  structure[max3 + dangle3 + 1] = '\0';

  return structure;
}


PRIVATE void
make_ptypes(vrna_fold_compound_t  *vc,
            int                   i)
{
  int j, k, type, n, maxdist, turn, noLP;
  short *S;
  char **ptype;
  vrna_md_t *md;

  n       = (int)vc->length;
  S       = vc->sequence_encoding2;
  ptype   = vc->ptype_local;
  maxdist = vc->window_size;
  md      = &(vc->params->model_details);
  turn    = md->min_loop_size;
  noLP    = md->noLP;

  for (k = turn + 1; k < maxdist; k++) {
    j = i + k;
    if (j > n)
      break;

    type = md->pair[S[i]][S[j]];

    if (noLP && type) {
      if (!ptype[i + 1][j - 1 - i - 1])
        if (j == n || i == 1 || (!md->pair[S[i - 1]][S[j + 1]]))
          type = 0;
    }

    ptype[i][j - i] = type;
  }
}


PRIVATE double
cov_score(vrna_fold_compound_t  *fc,
          int                   i,
          int                   j,
          float                 **dm)
{
  char **AS;
  short **S;
  int n_seq, k, l, s, type;
  double score;
  vrna_md_t *md;
  int pfreq[8] = {
    0, 0, 0, 0, 0, 0, 0, 0
  };

  n_seq = fc->n_seq;
  AS    = fc->sequences;
  S     = fc->S;
  md    = &(fc->params->model_details);

  for (s = 0; s < n_seq; s++) {
    if (S[s][i] == 0 && S[s][j] == 0) {
      type = 7;                             /* gap-gap  */
    } else {
      if ((AS[s][i] == '~') || (AS[s][j] == '~'))
        type = 7;
      else
        type = md->pair[S[s][i]][S[s][j]];
    }

    pfreq[type]++;
  }

  if (pfreq[0] * 2 + pfreq[7] > n_seq) {
    return NONE;
  } else {
    for (k = 1, score = 0.; k <= 6; k++) /* ignore pairtype 7 (gap-gap) */
      for (l = k; l <= 6; l++)
        /*
         * scores for replacements between pairtypes
         * consistent or compensatory mutations score 1 or 2
         */
        score += pfreq[k] * pfreq[l] * dm[k][l];
  }

  /* counter examples score -1, gap-gap scores -0.25   */
  return md->cv_fact * ((UNIT * score) / n_seq - md->nc_fact * UNIT * (pfreq[0] + pfreq[7] * 0.25));
}


PRIVATE void
make_pscores(vrna_fold_compound_t *fc,
             int                  i,
             float                **dm)
{
  /*
   * calculate co-variance bonus for each pair depending on
   * compensatory/consistent mutations and incompatible seqs
   * should be 0 for conserved pairs, >0 for good pairs
   */
  int n, j, **pscore, maxd, turn, noLP;
  vrna_md_t *md;

  n       = (int)fc->length;
  maxd    = fc->window_size;
  pscore  = fc->pscore_local;
  md      = &(fc->params->model_details);
  turn    = md->min_loop_size;
  noLP    = md->noLP;

  /*fill pscore[start], too close*/
  for (j = i + 1; (j < i + turn + 1) && (j <= n); j++)
    pscore[i][j - i] = NONE;
  for (j = i + turn + 1; ((j <= n) && (j <= i + maxd)); j++)
    pscore[i][j - i] = cov_score(fc, i, j, dm);

  if (noLP) {
    /* remove unwanted lonely pairs */
    int otype = 0, ntype = 0;
    for (j = i + turn; ((j < n) && (j < i + maxd)); j++) {
      if ((i > 1) && (j < n))
        otype = cov_score(fc, i - 1, j + 1, dm);

      if (i < n)
        ntype = pscore[i + 1][j - 1 - (i + 1)];
      else
        ntype = NONE;

      if ((otype < -4 * UNIT) && (ntype < -4 * UNIT)) /* worse than 2 counterex */
        pscore[i][j - i] = NONE;                      /* i.j can only form isolated pairs */
    }
  }

  if ((j - i + 1) > maxd)
    pscore[i][j - i] = NONE;
}


PRIVATE void
default_callback(int        start,
                 int        end,
                 const char *structure,
                 float      en,
                 void       *data)
{
  FILE *output      = ((hit_data *)data)->output;
  int dangle_model  = ((hit_data *)data)->dangle_model;

  if ((dangle_model == 2) && (start > 1))
    fprintf(output, ".%s (%6.2f) %4d\n", structure, en, start - 1);
  else
    fprintf(output, "%s (%6.2f) %4d\n ", structure, en, start);
}


#ifdef VRNA_WITH_SVM
PRIVATE void
default_callback_z(int        start,
                   int        end,
                   const char *structure,
                   float      en,
                   float      zscore,
                   void       *data)
{
  FILE *output      = ((hit_data *)data)->output;
  int dangle_model  = ((hit_data *)data)->dangle_model;

  if ((dangle_model == 2) && (start > 1))
    fprintf(output, ".%s (%6.2f) %4d z= %.3f\n", structure, en, start - 1, zscore);
  else
    fprintf(output, "%s (%6.2f) %4d z= %.3f\n ", structure, en, start, zscore);
}


#endif


PRIVATE void
default_callback_comparative(int        start,
                             int        end,
                             const char *structure,
                             float      en,
                             void       *data)
{
  FILE *output      = ((hit_data *)data)->output;
  int dangle_model  = ((hit_data *)data)->dangle_model;
  int csv           = ((hit_data *)data)->csv;

  if (csv == 1) {
    if ((dangle_model == 2) && (start > 1))
      fprintf(output, ".%s ,%6.2f, %4d, %4d\n", structure, en, start - 1, end);
    else
      fprintf(output, "%s ,%6.2f, %4d, %4d\n", structure, en, start, end);
  } else {
    if ((dangle_model == 2) && (start > 1))
      fprintf(output, ".%s (%6.2f) %4d - %4d\n", structure, en, start - 1, end);
    else
      fprintf(output, "%s (%6.2f) %4d - %4d\n", structure, en, start, end);
  }
}


struct block {
  vrna_fold_compound_t *fc;
  short *pt;
  unsigned long int start;
  unsigned long int end;
  unsigned long int shift;
  int energy;
  int energy_no3d;                   /* energy without 3'dangle */
  struct block *next_entry;
};


PRIVATE int
update_block(unsigned int i,
             unsigned int max_n,
             struct block *b)
{
  short *S1, *S2, d5, d3;
  unsigned int type;
  int n, i_local, j_local, ediff, dangles;
  vrna_fold_compound_t *fc;
  vrna_param_t *params;
  vrna_md_t *md;

  fc      = b->fc;
  n       = fc->length;
  S1      = fc->sequence_encoding;
  S2      = fc->sequence_encoding2;
  params  = fc->params;
  md      = &(params->model_details);
  dangles = md->dangles;

  /* re-evaluate energy for removal of base pair (i, k) */
  i_local = (int)(i - (b->start - b->shift) + 1);

  if (b->pt[i_local]) {
    j_local = (int)(b->pt[i_local]);
    /* compute energy differences due to removal of base pair (i_local, j_local) */
    ediff = vrna_eval_move_pt(fc, b->pt, -i_local, -j_local);

    /* update energy */
    b->energy += ediff;
    /* remove base pair from pair table */
    b->pt[i_local] = b->pt[j_local] = 0;

    /* since we've removed the outermost base pair, the block size
     *  decreases on the 3' end as well. At this point, we also
     *  remove any further trailing bases, if any
     */
    int end = j_local;

    do {
      b->end--;
      if (b->start == b->end)
        break;
    } while (b->pt[--end] == 0);

    if (b->end <= b->start)
      return 0;

    /* check whether we need to split the block into multiple ones */
    size_t stems        = 0;
    size_t mem_stems    = 10;
    size_t *start_stem  = (size_t *)vrna_alloc(sizeof(size_t) * mem_stems);
    size_t *end_stem    = (size_t *)vrna_alloc(sizeof(size_t) * mem_stems);

    for (size_t pos = i_local + 1; pos <= end; pos++)
      if (b->pt[pos] > pos) {
        start_stem[stems] = pos;
        end_stem[stems]   = b->pt[pos];
        stems++;
        if (stems == mem_stems) {
          mem_stems   *= 1.4;
          start_stem  = vrna_realloc(start_stem, sizeof(size_t) * mem_stems);
          end_stem    = vrna_realloc(end_stem, sizeof(size_t) * mem_stems);
        }

        pos = b->pt[pos];
      }

    if (stems > 1) {
      for (size_t k = stems - 1; k > 0; k--) {
        /* create a new block */
        struct block *new_block = (struct block *)vrna_alloc(sizeof(struct block));

        /* construct global coordinates for new block */
        new_block->start  = i + start_stem[k] - 1 - b->shift;
        new_block->end    = i + end_stem[k] - 1 - b->shift;
        new_block->shift  = (dangles == 2) ? 1 : 0;

        /* create structure pair table for new block */
        size_t l = end_stem[k] - start_stem[k] + 1;
        if (dangles == 2) {
          l++;    /* for 5' dangles */
          if (new_block->end < max_n)
            l++;  /* for 3' dangles */
        }

        new_block->pt     = (short *)vrna_alloc(sizeof(short) * (l + 1));
        new_block->pt[0]  = l;
        /* go through original structure and extract base pair positions relative to new block */
        for (size_t p = start_stem[k]; p <= end_stem[k]; p++) {
          if (b->pt[p] > p) {
            short i, j;

            i = p - start_stem[k] + 1 + new_block->shift;
            j = b->pt[p] - start_stem[k] + 1 + new_block->shift;

            new_block->pt[i]  = j;
            new_block->pt[j]  = i;

            /* remove base pair from original block */
            b->pt[b->pt[p]] = 0;
            b->pt[p]        = 0;
          }
        }

        /* create fold_compound for new block */
        char *seq_block = (char *)vrna_alloc(sizeof(char) * (l + 1));
        memcpy(seq_block, b->fc->sequence + start_stem[k] - 1 - new_block->shift,
               sizeof(char) * (l));
        new_block->fc = vrna_fold_compound(seq_block,
                                           &(b->fc->params->model_details),
                                           VRNA_OPTION_DEFAULT | VRNA_OPTION_EVAL_ONLY);
        free(seq_block);

        /* compute energy of new block */
        new_block->energy = vrna_eval_structure_pt(new_block->fc, new_block->pt);

        /* search for position where we can link the new block to */
        struct block *ptr, *ptr_prev;
        ptr_prev = ptr = b;
        do {
          ptr_prev  = ptr;
          ptr       = ptr->next_entry;
        } while ((ptr) && (ptr->start < new_block->start));

        /* insert new block */
        ptr_prev->next_entry  = new_block;
        new_block->next_entry = ptr;
      }

      /* Finally, we need to update the original block */

      /* update global end position */
      b->end = i + end_stem[0] - 1 - b->shift;

      /* update energy */
      b->energy = vrna_eval_structure_pt(b->fc, b->pt);
    }

    /* clean up memory */
    free(start_stem);
    free(end_stem);

    /* update for odd dangle models below */
    if (dangles % 1) {
    }
  } else if (dangles % 1) {
    /* position i is unpaired, so we only need to update energies if
     * position i + 1 forms a base pair due to 5' dangles
     */
    if (b->pt[i + 1]) {
      i_local = (int)i + b->start - 1 + 1;
      j_local = b->pt[i_local + 1];
      switch (dangles) {
        case 2:
          d5    = S1[i_local - 1];
          d3    = (j_local + 1 <= i_local + n - 1) ? S1[j_local + 1 - 1] : -1;
          type  = vrna_get_ptype_md(S2[i_local],
                                    S2[j_local],
                                    md);
          ediff = vrna_E_ext_stem(type, -1, d3, params) -
                  vrna_E_ext_stem(type, d5, d3, params);
          break;

        case 0:
          ediff = 0;
          break;

        default:
          ediff = 0;
          break;
      }
      /* update the energy */
      b->energy += ediff;
    }
  }

  /* increment start position and shift */
  b->start++;
  b->shift++;

  return 1;
}


PRIVATE struct block *
remove_block(struct block **before,
             struct block *block)
{
  struct block *ptr, **ptr2;

  if (*before == block)
    ptr2 = before;
  else
    ptr2 = &((*before)->next_entry);

  /* remember next entry */
  ptr = block->next_entry;

  /* cleanup memory for block we want to remove */
  vrna_fold_compound_free(block->fc);
  free(block->pt);
  /* remove block */
  free(block);

  /* link next entry */
  *ptr2 = ptr;

  return ptr;
}


PRIVATE void
truncate_blocks(unsigned long int i,
                unsigned long int n,
                struct block      **block_list)
{
  struct block *ptr_prev, *ptr;

  ptr_prev  = NULL;
  ptr       = *block_list;

  while (ptr) {
    /* 1. remove block if it was superseded by index i */
    if (ptr->end <= i) {
      ptr = remove_block((ptr_prev) ? &ptr_prev : block_list, ptr);
      continue;
    }

    /* 2. Update energy */
    if (ptr->start == i) {
      /* remove nucleotide at position i and update energies accordingly */
      /* this also splits a substructure component into consecutive
       * tokens, if necessary
       */
      if (!update_block(i, n, ptr)) {
        /* this block doesn't contain any more pairs so we remove it entirely */
        ptr = remove_block((ptr_prev) ? &ptr_prev : block_list, ptr);
        continue;
      }
    } else if (ptr->start > i) {
      break;
    }

    /* go to next block */
    ptr_prev  = ptr;
    ptr       = ptr->next_entry;
  }
}


PRIVATE struct block *
extract_Lfold_entry(FILE            *f,
                    long            line_start,
                    const char      *sequence,
                    const vrna_md_t *md)
{
  char *l, *seq_local, *structure;
  unsigned long int i, j;
  float en;
  struct block *storage;

  storage = NULL;

  if (fseek(f, line_start, SEEK_SET) != -1) {
    l         = vrna_read_line(f);
    en        = INF / 100.;
    structure = (char *)vrna_alloc(sizeof(char) * (strlen(l) + 1));

    /*  each line should consist of 3 parts,
     *  1. A dot-bracket structure
     *  2. An energy value in kcal/mol
     *  3. A starting position
     */
    if (sscanf(l, "%[.()] %*c %f %*c %lu", structure, &en, &i) == 3) {
      storage = (struct block *)vrna_alloc(sizeof(struct block));

      j = i + strlen(structure) - 1;

      seq_local = (char *)vrna_alloc(sizeof(char) * (j - i + 2));
      memcpy(seq_local, sequence + i - 1, (sizeof(char) * (j - i + 1)));

      storage->fc = vrna_fold_compound(seq_local,
                                       md,
                                       VRNA_OPTION_DEFAULT | VRNA_OPTION_EVAL_ONLY);
      storage->pt           = vrna_ptable(structure);
      storage->start        = i;
      storage->end          = j;
      storage->shift        = 0;
      storage->energy       = vrna_convert_kcal_to_dcal(en);
      storage->energy_no3d  = 0; /* energy without 3'dangle */
      storage->next_entry   = NULL;

      free(seq_local);

      /* with dangles==0 or dangles==2, we don't need the trailing and leading unpaired base(s) */
      if ((md->dangles % 1) == 0) {
        if (storage->pt[1] == 0) {
          storage->start++;
          storage->shift++;
        }

        if (storage->pt[storage->fc->length] == 0)
          storage->end--;
      }
    }

    /* we only store the pair table, and drop the dot-bracket string */
    free(structure);
    free(l);
  }

  return storage;
}


PRIVATE unsigned long int
insert_block(char         *structure,
             struct block *selected,
             int          *energy)
{
  short *pt;
  unsigned int long i, i_local, start, end, shift;
  int e;

  pt    = selected->pt;
  start = selected->start;
  end   = selected->end;
  shift = selected->shift;
  e     = selected->energy;

  /* The last block is always part of the full length MFE */
  for (i = start; i <= end; i++) {
    i_local = i - start + shift + 1;
    if (pt[i_local] > i_local) {
      structure[i - 1]                                = '(';
      structure[start - shift + pt[i_local] - 1 - 1]  = ')';
    }
  }

  *energy = *energy - e;

  return end;
}


PRIVATE void
print_block_list(struct block *block_list)
{
  struct block *ptr;
  size_t cnt = 0;

  for (ptr = block_list; ptr; ptr = ptr->next_entry, cnt++)
    printf("block %lu: en=%d, start: %lu, end: %lu, shift: %lu\n",
           cnt,
           ptr->energy,
           ptr->start,
           ptr->end,
           ptr->shift);
  printf("%lu blocks remaining\n", cnt);
  fflush(stdout);
}


PRIVATE void
append_blocks(struct block          **last_block,
              FILE                  *f,
              long                  *lines,
              size_t                *lines_left,
              vrna_fold_compound_t  *fc,
              unsigned long int     max_start)
{
  vrna_md_t *md;

  md = &(fc->params->model_details);

  while (((*lines_left) > 0) &&
         ((*last_block)->start < max_start)) {
    (*lines_left)--;

    struct block *ptr = extract_Lfold_entry(f, lines[(*lines_left)], fc->sequence, md);

    if (ptr != NULL) {
      (*last_block)->next_entry = ptr;
      *last_block               = ptr;
    }
  }
  ;
}


PUBLIC int
vrna_backtrack_window(vrna_fold_compound_t  *fc,
                      const char            *Lfold_filename,
                      long                  file_pos,
                      char                  **structure,
                      double                mfe)
{
  unsigned int n, look_ahead;
  int e, maxdist, underflows, *f3;
  int ret;
  double mfe_corr;
  FILE *f;

  ret         = 0;
  *structure  = NULL;

  if ((fc) &&
      (Lfold_filename) &&
      (structure)) {
    vrna_md_t *md;

    n           = fc->length;
    md          = &(fc->params->model_details);
    maxdist     = md->window_size;
    look_ahead  = 3 * maxdist;
    underflows  = 0;
    mfe_corr    = mfe;

    f3 = fc->matrices->f3_local;

    if (md->dangles % 2) {
      vrna_message_warning(
        "Global mfE structure backtracking not available for odd dangle models yet!");
      return ret;
    }

    /* check whether we need to adjust energies due to integer underflows in the forward recursion */
    while (vrna_convert_kcal_to_dcal(mfe_corr) < f3[1]) {
      mfe_corr -= (double)(UNDERFLOW_CORRECTION) / 100.;
      underflows++;
    }

    e = vrna_convert_kcal_to_dcal(mfe_corr);

#if 0
    if (underflows)
      printf("detected %d underflows %d vs. %d vs %6.2f\n", underflows, e, f3[1], mfe);

#endif

    e = f3[1];

    /* default to start at beginning of file */
    if (file_pos < 0)
      file_pos = 0;

    /* open Lfold output file and start parsing line positions */
    if ((f = fopen(Lfold_filename, "r"))) {
      if (fseek(f, file_pos, SEEK_SET) != -1) {
        size_t num_lines, mem_lines;
        long *lines;    /* array of file positions indicating start of lines */

        num_lines = 0;
        mem_lines = 1024;

        lines               = (long *)vrna_alloc(sizeof(long) * mem_lines);
        lines[num_lines++]  = ftell(f);

        /* 1. Fill array of file positions that coincide with start of lines */
        do {
          /* increase memory if necessary */
          if (num_lines == mem_lines) {
            mem_lines *= 1.4;
            lines     = (long *)vrna_realloc(lines, sizeof(long) * mem_lines);
          }

          /* seek to next newline char */
          while ((!feof(f)) && (fgetc(f) != '\n'));

          /* stop at end of file */
          if (feof(f))
            break;

          lines[num_lines++] = ftell(f);
        } while(1);

        /* do the actual parsing of the data */
        if (num_lines > 0) {
          num_lines--;

          size_t block_num    = 0;
          unsigned long int i = num_lines;
          size_t lines_left   = num_lines;

          /* we will use *ptr as current block pointer and save the list start in 8block_list */
          struct block *block_list, *ptr, *block_list_last;

          /* handle first block separately */
          do {
            lines_left--;

            block_list_last = block_list = ptr = extract_Lfold_entry(f,
                                                                     lines[lines_left],
                                                                     fc->sequence,
                                                                     md);
          } while ((block_list == NULL) && (lines_left > 0));

          /* start the actual backtracing procedure if we've read at least one block */
          if (block_list) {
            block_num++;
            *structure = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));
            memset(*structure, (int)'.', fc->length);

            /*  go through remaining file in reverse order
             *  to extract the relevant data
             */
            append_blocks(&block_list_last,
                          f,
                          lines,
                          &lines_left,
                          fc,
                          ptr->start + look_ahead);

            ptr = block_list;

            i = insert_block(*structure, ptr, &e);

            for (unsigned long int ii = ptr->start; ii <= i; ii++) {
              /* truncate remaining blocks */
              truncate_blocks(ii, n, &block_list);
              append_blocks(&block_list_last,
                            f,
                            lines,
                            &lines_left,
                            fc,
                            ii + look_ahead);
            }

            i++;

            /* proceed with remaining blocks */
            for (; i <= fc->length; i++) {
              unsigned long int prev_i = i;

              /* i pairs with some k, find block representing the substructure enclodes by (i,k) */
              if (e != f3[i + 1]) {
                if ((underflows) &&
                    (e == (f3[i + 1] - UNDERFLOW_CORRECTION))) {
                  underflows--;
                  e += UNDERFLOW_CORRECTION;
                } else {
                  /* go through list of blocks that start at i to check which one we can insert */
                  char found = 0;
                  for (ptr = block_list; ptr; ptr = ptr->next_entry) {
                    if (ptr->start > i)
                      break;

                    if ((ptr->start == i) &&
                        (ptr->end > i)) {
                      if (e == (ptr->energy + f3[ptr->end + 1])) {
                        /* found the block, let's insert it */
                        found = 1;

                        i = insert_block(*structure, ptr, &e);

                        break;
                      } else if ((underflows) &&
                                 (e == (ptr->energy + f3[ptr->end + 1] - UNDERFLOW_CORRECTION))) {
                        underflows--;
                        found = 1;
                        e     += UNDERFLOW_CORRECTION;
                        i     = insert_block(*structure, ptr, &e);

                        break;
                      }
                    }
                  }
                  if (!found)
                    printf("didn't find block for position %lu\n", i);
                }
              }

              /*
               *  if i is unpaired, so we do nothing except for
               *  truncating the remaining blocks
               */
              for (unsigned long int ii = prev_i; ii <= i; ii++) {
                truncate_blocks(ii, n, &block_list);
                append_blocks(&block_list_last,
                              f,
                              lines,
                              &lines_left,
                              fc,
                              ii + look_ahead);
              }
            }

            ret = 1;
          }
        }
      }

      fclose(f);
    }
  }

  return ret;
}


PRIVATE INLINE int
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               struct aux_arrays    *aux)
{
  unsigned char hc_decompose;
  int e, new_c, energy, stackEnergy, dangle_model, noLP,
      *DMLi1, *DMLi2, *cc, *cc1;

  dangle_model  = fc->params->model_details.dangles;
  noLP          = fc->params->model_details.noLP;
  hc_decompose  = fc->hc->matrix_local[i][j - i];
  DMLi1         = aux->DMLi1;
  DMLi2         = aux->DMLi2;
  cc            = aux->cc;
  cc1           = aux->cc1;
  e             = INF;

  /* do we evaluate this pair? */
  if (hc_decompose) {
    new_c = INF;

    /* check for hairpin loop */
    energy  = vrna_E_hp_loop(fc, i, j);
    new_c   = MIN2(new_c, energy);

    /* check for multibranch loops */
    energy  = vrna_E_mb_loop_fast(fc, i, j, DMLi1, DMLi2);
    new_c   = MIN2(new_c, energy);

    if (dangle_model == 3) {
      /* coaxial stacking */
      energy  = vrna_E_mb_loop_stack(fc, i, j);
      new_c   = MIN2(new_c, energy);
    }

    /* check for interior loops */
    energy  = vrna_E_int_loop(fc, i, j);
    new_c   = MIN2(new_c, energy);

    /* remember stack energy for --noLP option */
    if (noLP) {
      stackEnergy = vrna_E_stack(fc, i, j);
      new_c       = MIN2(new_c, cc1[j - 1 - (i + 1)] + stackEnergy);
      cc[j - i]   = new_c;
      if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
          (cc[j - i] != INF))
        cc[j - i] -= fc->pscore_local[i][j - i];

      e = cc1[j - 1 - (i + 1)] + stackEnergy;
    } else {
      e = new_c;
    }

    /* finally, check for auxiliary grammar rule(s) */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_c)) {
      energy  = fc->aux_grammar->cb_aux_c(fc, i, j, fc->aux_grammar->data);
      e   = MIN2(e, energy);
    }

    if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
        (e != INF))
      e -= fc->pscore_local[i][j - i];
  } /* end >> if (pair) << */

  return e;
}


PRIVATE INLINE struct aux_arrays *
get_aux_arrays(unsigned int maxdist)
{
  unsigned int j;
  struct aux_arrays *aux = (struct aux_arrays *)vrna_alloc(sizeof(struct aux_arrays));

  aux->cc     = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));   /* auxilary arrays for canonical structures     */
  aux->cc1    = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));   /* auxilary arrays for canonical structures     */
  aux->Fmi    = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));   /* holds row i of fML (avoids jumps in memory)  */
  aux->DMLi   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));   /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  aux->DMLi1  = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));   /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  aux->DMLi2  = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));   /*                MIN(fML[i+2,k]+fML[k+1,j])    */

  /* prefill helper arrays */
  for (j = 0; j < maxdist + 5; j++)
    aux->Fmi[j] = aux->DMLi[j] = aux->DMLi1[j] = aux->DMLi2[j] = INF;

  return aux;
}


PRIVATE INLINE void
rotate_aux_arrays(struct aux_arrays *aux,
                  unsigned int      maxdist)
{
  unsigned int j;
  int *FF;

  FF          = aux->DMLi2;
  aux->DMLi2  = aux->DMLi1;
  aux->DMLi1  = aux->DMLi;
  aux->DMLi   = FF;
  FF          = aux->cc1;
  aux->cc1    = aux->cc;
  aux->cc     = FF;
  for (j = 1; j < maxdist + 5; j++)
    aux->cc[j] = aux->Fmi[j] = aux->DMLi[j] = INF;
}


PRIVATE INLINE void
free_aux_arrays(struct aux_arrays *aux)
{
  free(aux->cc);
  free(aux->cc1);
  free(aux->Fmi);
  free(aux->DMLi);
  free(aux->DMLi1);
  free(aux->DMLi2);
  free(aux);
}
