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
#include <svm.h>
#include "ViennaRNA/utils/svm.h"
#endif

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define MAXSECTORS                  500   /* dimension for a backtrack array */
#define INT_CLOSE_TO_UNDERFLOW(i)   ((i) <= (INT_MIN / 16))
#define UNDERFLOW_CORRECTION        (INT_MIN / 32)

#define NONE -10000 /* score for forbidden pairs */


typedef struct {
  FILE  *output;
  int   dangle_model;
  int   csv;
} hit_data;

#ifdef VRNA_WITH_SVM
typedef struct {
  struct svm_model  *avg_model;
  struct svm_model  *sd_model;
  double            min_z;
  int               with_zsc;
} zscoring_dat;
#endif

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
            zscoring_dat                    *z_dat,
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


PRIVATE int
fill_arrays_comparative(vrna_fold_compound_t      *fc,
                        int                       *underflow,
                        vrna_mfe_window_callback  *cb,
                        void                      *data);


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


#ifdef VRNA_WITH_SVM

PRIVATE int
want_backtrack(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               zscoring_dat         *d,
               double               *z);


PRIVATE void
default_callback_z(int        start,
                   int        end,
                   const char *structure,
                   float      en,
                   float      zscore,
                   void       *data);


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
  int           energy, underflow, n_seq;
  float         mfe_local;

#ifdef VRNA_WITH_SVM
  zscoring_dat  z_dat;
#endif
  /* keep track of how many times we were close to an integer underflow */
  underflow = 0;

  if (!vrna_fold_compound_prepare(vc, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW)) {
    vrna_message_warning("vrna_mfe_window@Lfold.c: Failed to prepare vrna_fold_compound");
    return (float)(INF / 100.);
  }

  if (vc->type == VRNA_FC_TYPE_COMPARATIVE) {
    n_seq   = vc->n_seq;
    energy  = fill_arrays_comparative(vc, &underflow, cb, data) / 100.;

    mfe_local =
      (underflow > 0) ? ((float)underflow * (float)(UNDERFLOW_CORRECTION)) / (100. * n_seq) : 0.;
    mfe_local += (float)energy / (100. * n_seq);
  } else {
#ifdef VRNA_WITH_SVM
    z_dat.with_zsc  = 0;
    energy          = fill_arrays(vc, &underflow, cb, &z_dat, NULL, data);
#else
    energy = fill_arrays(vc, &underflow, cb, data);
#endif
    mfe_local = (underflow > 0) ? ((float)underflow * (float)(UNDERFLOW_CORRECTION)) / 100. : 0.;
    mfe_local += (float)energy / 100.;
  }

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
  int           energy, underflow;
  float         mfe_local;
  zscoring_dat  zsc_data;

  if (vc->type == VRNA_FC_TYPE_COMPARATIVE) {
    vrna_message_warning(
      "vrna_mfe_window_zscore@mfe_window.c: Comparative prediction not implemented");
    return (float)(INF / 100.);
  }

  if (!vrna_fold_compound_prepare(vc, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW)) {
    vrna_message_warning("vrna_mfe_window@Lfold.c: Failed to prepare vrna_fold_compound");
    return (float)(INF / 100.);
  }

  zsc_data.with_zsc   = 1;
  zsc_data.avg_model  = svm_load_model_string(avg_model_string);
  zsc_data.sd_model   = svm_load_model_string(sd_model_string);
  zsc_data.min_z      = min_z;

  /* keep track of how many times we were close to an integer underflow */
  underflow = 0;

  energy = fill_arrays(vc, &underflow, NULL, &zsc_data, cb_z, data);
  svm_free_model_content(zsc_data.avg_model);
  svm_free_model_content(zsc_data.sd_model);

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
            zscoring_dat                    *zsc_data,
            vrna_mfe_window_zscore_callback *cb_z,
#endif
            void                            *data)
{
  /* fill "c", "fML" and "f3" arrays and return  optimal energy */

  char          **ptype, *prev;
  unsigned char hc_decompose;
  int           i, j, length, energy, maxdist, **c, **fML, *f3, no_close,
                type, with_gquad, dangle_model, noLP, noGUclosure, turn,
                *cc, *cc1, *Fmi, *DMLi, *DMLi1, *DMLi2, prev_i, prev_j,
                prev_end, prev_en, new_c, stackEnergy;

#ifdef VRNA_WITH_SVM
  double        prevz;
#endif
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  length        = vc->length;
  ptype         = vc->ptype_local;
  maxdist       = vc->window_size;
  P             = vc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  noLP          = md->noLP;
  noGUclosure   = md->noGUclosure;
  turn          = md->min_loop_size;
  hc            = vc->hc;
  do_backtrack  = 0;
  prev_i        = 0;
  prev_j        = 0;
  prev_end      = 0;
  prev          = NULL;
  prev_en       = 0;
#ifdef VRNA_WITH_SVM
  prevz = 0.;
#endif

  c   = vc->matrices->c_local;
  fML = vc->matrices->fML_local;
  f3  = vc->matrices->f3_local;

  cc    = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /* linear array for calculating canonical structures */
  cc1   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /*   "     "        */
  Fmi   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /* holds row i of fML (avoids jumps in memory) */
  DMLi  = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
  DMLi1 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /*             MIN(fML[i+1,k]+fML[k+1,j])  */
  DMLi2 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /*             MIN(fML[i+2,k]+fML[k+1,j])  */

  /* reserve additional memory for j-dimension */
  allocate_dp_matrices(vc);

  init_constraints(vc, NULL);

  for (j = 0; j < maxdist + 5; j++)
    Fmi[j] = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;

  if (with_gquad)
    vrna_gquad_mx_local_update(vc, length - maxdist - 4);

  for (i = length - turn - 1; i >= 1; i--) {
    /* i,j in [1..length] */
    for (j = i + turn + 1; j <= length && j <= i + maxdist; j++) {
      hc_decompose  = hc->matrix_local[i][j - i];
      type          = vrna_get_ptype_window(i, j, ptype);

      no_close = (((type == 3) || (type == 4)) && noGUclosure);

      if (hc_decompose) {
        /* we have a pair */
        new_c       = INF;
        stackEnergy = INF;

        if (!no_close) {
          /* check for hairpin loop */
          energy  = vrna_E_hp_loop(vc, i, j);
          new_c   = MIN2(new_c, energy);

          /* check for multibranch loops */
          energy  = vrna_E_mb_loop_fast(vc, i, j, DMLi1, DMLi2);
          new_c   = MIN2(new_c, energy);
        }

        if (dangle_model == 3) {
          /* coaxial stacking */
          energy  = vrna_E_mb_loop_stack(vc, i, j);
          new_c   = MIN2(new_c, energy);
        }

        /* check for interior loops */
        energy  = vrna_E_int_loop(vc, i, j);
        new_c   = MIN2(new_c, energy);

        /* remember stack energy for --noLP option */
        if (noLP) {
          stackEnergy = vrna_E_stack(vc, i, j);
          new_c       = MIN2(new_c, cc1[j - 1 - (i + 1)] + stackEnergy);
          cc[j - i]   = new_c;
          c[i][j - i] = cc1[j - 1 - (i + 1)] + stackEnergy;
        } else {
          c[i][j - i] = new_c;
        }
      } /* end >> if (pair) << */
      else {
        c[i][j - i] = INF;
      }

      /*
       * done with c[i,j], now compute fML[i,j]
       * free ends ? -----------------------------------------
       */
      fML[i][j - i] = vrna_E_ml_stems_fast(vc, i, j, Fmi, DMLi);
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
          if (want_backtrack(vc, ii, jj, zsc_data, &thisz)) {
#endif
          ss = backtrack(vc, ii, jj);
          if (prev) {
            if ((jj < prev_j) || (strncmp(ss + prev_i - ii, prev, prev_j - prev_i + 1))) {
              /* ss does not contain prev */
#ifdef VRNA_WITH_SVM
              if (zsc_data->with_zsc)
                cb_z(prev_i, prev_end, prev, prev_en / 100., prevz, data);
              else
#endif
              cb(prev_i, prev_end, prev, prev_en / 100., data);
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
          if (zsc_data->with_zsc)
            cb_z(prev_i, prev_end, prev, prev_en / 100., prevz, data);
          else
#endif
          cb(prev_i, prev_end, prev, prev_en / 100., data);

          free(prev);
          prev = NULL;
#ifdef VRNA_WITH_SVM
        } else if ((f3[i] < 0) && (!zsc_data->with_zsc)) {
#else
        } else if (f3[i] < 0) {
#endif
          /* why !zsc? */
          int ii, jj;
          ii  = i;
          jj  = vrna_BT_ext_loop_f3_pp(vc, &ii, maxdist);
          if (jj > 0) {
#ifdef VRNA_WITH_SVM
            double thisz = 0;
            if (want_backtrack(vc, ii, jj, zsc_data, &thisz)) {
#endif
            ss = backtrack(vc, ii, jj);
#ifdef VRNA_WITH_SVM
            if (zsc_data->with_zsc)
              cb_z(ii, MIN2(jj + ((dangle_model) ? 1 : 0),
                            length), ss, (f3[1] - f3[jj + 1]) / 100., thisz, data);
            else
#endif
            cb(ii, MIN2(jj + ((dangle_model) ? 1 : 0),
                        length), ss, (f3[1] - f3[jj + 1]) / 100., data);

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

    {
      int *FF; /* rotate the auxilliary arrays */

      /* check for values close to integer underflow */
      if (INT_CLOSE_TO_UNDERFLOW(f3[i])) {
        /* correct f3 free energies and increase underflow counter */
        int cnt;
        for (cnt = i; cnt <= MIN2(i + maxdist + 2, length); cnt++)
          f3[cnt] -= UNDERFLOW_CORRECTION;
        (*underflow)++;
      }

      FF    = DMLi2;
      DMLi2 = DMLi1;
      DMLi1 = DMLi;
      DMLi  = FF;
      FF    = cc1;
      cc1   = cc;
      cc    = FF;
      for (j = 0; j < maxdist + 5; j++)
        cc[j] = Fmi[j] = DMLi[j] = INF;

      rotate_dp_matrices(vc, i);
      rotate_constraints(vc, NULL, i);
    }
  }

  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);

  free_dp_matrices(vc);

  return f3[1];
}


#ifdef VRNA_WITH_SVM
PRIVATE int
want_backtrack(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               zscoring_dat         *d,
               double               *z)
{
  int bt = 1; /* we want to backtrack by default */

  if (d->with_zsc) {
    short *S;
    int info_avg, start, end, dangle_model, length, *f3;
    double average_free_energy;
    double sd_free_energy;

    length        = vc->length;
    S             = vc->sequence_encoding2;
    f3            = vc->matrices->f3_local;
    dangle_model  = vc->params->model_details.dangles;

    bt    = 0;   /* let the z-score decide whether we backtrack or not */
    start = (dangle_model) ? MAX2(1, i - 1) : i;
    end   = (dangle_model) ? MIN2(length, j + 1) : j;

    int *AUGC = get_seq_composition(S, start, end, length);

    /*\svm*/
    average_free_energy = avg_regression(AUGC[0],
                                         AUGC[1],
                                         AUGC[2],
                                         AUGC[3],
                                         AUGC[4],
                                         d->avg_model,
                                         &info_avg);

    if (info_avg == 0) {
      double difference;
      double min_sd = minimal_sd(AUGC[0], AUGC[1], AUGC[2], AUGC[3], AUGC[4]);
      difference = (f3[i] - f3[j + 1]) / 100. - average_free_energy;

      if (difference - (d->min_z * min_sd) <= 0.0001) {
        sd_free_energy  = sd_regression(AUGC[0], AUGC[1], AUGC[2], AUGC[3], AUGC[4], d->sd_model);
        *z              = difference / sd_free_energy;
        if ((*z) <= d->min_z)
          bt = 1;
      }
    }

    free(AUGC);
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


PRIVATE int
fill_arrays_comparative(vrna_fold_compound_t      *fc,
                        int                       *underflow,
                        vrna_mfe_window_callback  *cb,
                        void                      *data)
{
  /* fill "c", "fML" and "f3" arrays and return  optimal energy */
  short **S;
  char **strings, *prev, *ss;
  int **pscore, i, j, length, energy, turn, n_seq, **c,
      **fML, *f3, *cc, *cc1, *Fmi, *DMLi, *DMLi1, *DMLi2,
      maxdist, prev_i, prev_j, prev_en, prev_end, with_gquad,
      new_c, psc, stackEnergy, dangle_model;
  float **dm;
  vrna_mx_mfe_t *matrices;
  vrna_param_t *P;
  vrna_md_t *md;
  vrna_hc_t *hc;

  int olddm[7][7] = { { 0, 0, 0, 0, 0, 0, 0 },             /* hamming distance between pairs PRIVATE needed??*/
                      { 0, 0, 2, 2, 1, 2, 2 } /* CG */,
                      { 0, 2, 0, 1, 2, 2, 2 } /* GC */,
                      { 0, 2, 1, 0, 2, 1, 2 } /* GU */,
                      { 0, 1, 2, 2, 0, 2, 1 } /* UG */,
                      { 0, 2, 2, 1, 2, 0, 2 } /* AU */,
                      { 0, 2, 2, 2, 1, 2, 0 } /* UA */ };

  do_backtrack  = 0;
  prev_i        = 0;
  prev_j        = 0;
  prev_end      = 0;
  prev_en       = INF;

  dm    = NULL;
  prev  = NULL;
  ss    = NULL;

  strings       = fc->sequences;
  S             = fc->S;
  n_seq         = fc->n_seq;
  length        = fc->length;
  maxdist       = fc->window_size;
  matrices      = fc->matrices;
  hc            = fc->hc;
  P             = fc->params;
  md            = &(P->model_details);
  turn          = md->min_loop_size;
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;

  if (md->ribo) {
    if (RibosumFile != NULL)
      dm = readribosum(RibosumFile);
    else
      dm = get_ribosum((const char **)strings, n_seq, length);
  } else {
    /*use usual matrix*/
    dm = (float **)vrna_alloc(7 * sizeof(float *));
    for (i = 0; i < 7; i++) {
      dm[i] = (float *)vrna_alloc(7 * sizeof(float));
      for (j = 0; j < 7; j++)
        dm[i][j] = (float)olddm[i][j];
    }
  }

  pscore = fc->pscore_local;      /* precomputed array of pair types */

  c   = fc->matrices->c_local;    /* energy array, given that i-j pair */
  fML = fc->matrices->fML_local;  /* multi-loop auxiliary energy array */
  f3  = fc->matrices->f3_local;   /* energy of 5' end */

  cc    = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  cc1   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  Fmi   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi  = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi1 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi2 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));

  /* reserve additional memory for j-dimension */
  allocate_dp_matrices(fc);

  init_constraints(fc, dm);

  for (j = 0; j < maxdist + 5; j++)
    Fmi[j] = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;

  if (with_gquad)
    vrna_gquad_mx_local_update(fc, length - maxdist - 4);

  for (i = length - turn - 1; i >= 1; i--) {
    /* i,j in [1..length] */
    for (j = i + 1; j <= length && j <= i + turn; j++)
      c[i][j - i] = fML[i][j - i] = INF;
    for (j = i + turn + 1; j <= length && j <= i + maxdist; j++) {
      psc = pscore[i][j - i];

      if (hc->matrix_local[i][j - i]) {
        /* we evaluate this pair */
        new_c       = INF;
        stackEnergy = INF;

        /* hairpin ----------------------------------------------*/
        energy  = vrna_E_hp_loop(fc, i, j);
        new_c   = MIN2(new_c, energy);

        /* check for multibranch loops */
        energy  = vrna_E_mb_loop_fast(fc, i, j, DMLi1, NULL);
        new_c   = MIN2(new_c, energy);

        /* check for interior loops */
        energy  = vrna_E_int_loop(fc, i, j);
        new_c   = MIN2(new_c, energy);

        /* remember stack energy for --noLP option */
        if (md->noLP) {
          stackEnergy = vrna_E_stack(fc, i, j);

          if (cc1[j - 1 - (i + 1)] != INF) {
            new_c       = MIN2(new_c, cc1[j - 1 - (i + 1)] + stackEnergy);
            c[i][j - i] = cc1[j - 1 - (i + 1)] + stackEnergy;
          } else {
            c[i][j - i] = INF;
          }

          if (new_c != INF)
            new_c -= psc;

          cc[j - i]   = new_c; /* add covariance bonnus/penalty */
        } else {
          c[i][j - i] = new_c; /* add covariance bonnus/penalty */
        }

        if (c[i][j - i] != INF)
          c[i][j - i] -= psc;
      } /* end >> if (pair) << */
      else {
        c[i][j - i] = INF;
      }

      /* done with c[i,j], now compute fML[i,j] */
      fML[i][j - i] = vrna_E_ml_stems_fast(fc, i, j, Fmi, DMLi);
    } /* for (j...) */

    /* calculate energies of 5' and 3' fragments */
    f3[i] = vrna_E_ext_loop_3(fc, i);

    if (f3[i] < f3[i + 1]) {
      /*
       * instead of backtracing in the next iteration, we backtrack now
       * already. This is necessary to accomodate for change in free
       * energy due to unpaired nucleotides in the exterior loop, which
       * may happen in the case of using soft constraints
       */
      int ii, jj;
      ii  = i;
      jj  = vrna_BT_ext_loop_f3_pp(fc, &ii, maxdist);
      if (jj > 0) {
        ss = backtrack(fc, ii, jj);
        if (prev) {
          if ((jj < prev_j) || (strncmp(ss + prev_i - ii, prev, prev_j - prev_i + 1)))
            /* ss does not contain prev */
            cb(prev_i, prev_end, prev, prev_en / (100. * n_seq), data);

          free(prev);
        }

        prev      = ss;
        prev_i    = ii;
        prev_j    = jj;
        prev_end  = MIN2(jj + ((dangle_model) ? 1 : 0), length);
        prev_en   = f3[ii] - f3[jj + 1];
      } else if (jj == -1) {
        /* some error occured during backtracking */
        vrna_message_error("backtrack failed in short backtrack 1");
      }
    }

    if (i == 1) {
      if (prev) {
        cb(prev_i, prev_end, prev, prev_en / (100. * n_seq), data);

        free(prev);
        prev = NULL;
      } else if (f3[i] < 0) {
        int ii, jj;
        ii  = i;
        jj  = vrna_BT_ext_loop_f3_pp(fc, &ii, maxdist);
        if (jj > 0) {
          ss = backtrack(fc, ii, jj);
          cb(ii, MIN2(jj + ((dangle_model) ? 1 : 0),
                      length), ss, (f3[1] - f3[jj + 1]) / (100. * n_seq), data);

          free(ss);
        } else if (jj == -1) {
          /* some error occured during backtracking */
          vrna_message_error("backtrack failed in short backtrack 2");
        }
      }
    }

    {
      int *FF; /* rotate the auxilliary arrays */

      /* check for values close to integer underflow */
      if (INT_CLOSE_TO_UNDERFLOW(f3[i])) {
        /* correct f3 free energies and increase underflow counter */
        int cnt;
        for (cnt = i; cnt <= MIN2(i + maxdist + 2, length); cnt++)
          f3[cnt] -= UNDERFLOW_CORRECTION;
        (*underflow)++;
      }

      FF    = DMLi2;
      DMLi2 = DMLi1;
      DMLi1 = DMLi;
      DMLi  = FF;
      FF    = cc1;
      cc1   = cc;
      cc    = FF;
      for (j = 0; j < maxdist + 5; j++)
        cc[j] = Fmi[j] = DMLi[j] = INF;

      rotate_dp_matrices(fc, i);
      rotate_constraints(fc, dm, i);
    }
  }

  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);

  for (i = 0; i < 7; i++)
    free(dm[i]);
  free(dm);

  free_dp_matrices(fc);

  return matrices->f3[1];
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
  vrna_fold_compound_t  *fc;
  short                 *pt;
  unsigned long int     start;
  unsigned long int     end;
  unsigned long int     shift;
  int                   energy;
  int                   energy_no3d; /* energy without 3'dangle */
  struct block          *next_entry;
};


struct local_struct {
  short               *pt;
  unsigned int        start;
  unsigned int        shift;
  unsigned int        length;
  int                 energy;
  unsigned char       valid;
  struct local_struct *next_entry;
};



PRIVATE void
update_block(unsigned int  i,
             unsigned int  max_n,
             struct block  *b)
{
  short                 *S1, *S2, d5, d3;
  unsigned int          type;
  int                   n, i_local, j_local, ediff, dangles;
  vrna_fold_compound_t  *fc;
  vrna_param_t          *params;
  vrna_md_t             *md;

  fc      = b->fc;
  n       = fc->length;
  S1      = fc->sequence_encoding;
  S2      = fc->sequence_encoding2;
  params  = fc->params;
  md      = &(params->model_details);
  dangles = md->dangles;

  /* re-evaluate energy for removal of base pair (i, k) */
  i_local = (int)(i - (b->start - b->shift) + 1);

#if DEBUG
  printf("block pre truncation: [%d:%d] (shift %d), e = %d\n", b->start, b->end, b->shift, b->energy);
#endif

  if (b->pt[i_local]) {
    j_local = (int)(b->pt[i_local]);
    /* compute energy differences due to removal of base pair (i_local, j_local) */
    ediff   = vrna_eval_move_pt(fc, b->pt, -i_local, -j_local);

    /* update energy */
    b->energy += ediff;
    /* remove base pair from pair table */
    b->pt[i_local] = b->pt[j_local] = 0;

    /* since we've removed the outermost base pair, the block size
        decreases on the 3' end as well. At this point, we also
        remove any further trailing bases, if any
    */
    int end = j_local;

    do {
      b->end--;
      if (b->start == b->end)
        break;
    } while (b->pt[--end] == 0);

    /* check whether we need to split the block into multiple ones */
    size_t  stems       = 0;
    size_t  mem_stems   = 10;
    size_t  *start_stem = (size_t *)vrna_alloc(sizeof(size_t) * mem_stems);
    size_t  *end_stem   = (size_t *)vrna_alloc(sizeof(size_t) * mem_stems);

    for (size_t pos = i_local + 1; pos <= end; pos++)
      if (b->pt[pos] > pos) {
        start_stem[stems] = pos;
        end_stem[stems]   = b->pt[pos];
        stems++;
        if (stems == mem_stems) {
          mem_stems *= 1.4;
          start_stem  = vrna_realloc(start_stem, sizeof(size_t) * mem_stems);
          end_stem    = vrna_realloc(end_stem, sizeof(size_t) * mem_stems);
        }
        pos = b->pt[pos];
      }

    if (stems > 1) {
#if DEBUG
      printf("splitting interval [%d:%d] (shift: %d) into %u blocks\n", b->start + 1, b->end, b->shift, stems);
      char *ss = vrna_db_from_ptable(b->pt);
      printf("%s\n%s (%6.2f)\n", b->fc->sequence, ss, (float)b->energy / 100.);
      free(ss);
#endif
      struct block *tmp = b->next_entry;
      for (size_t k = stems - 1; k > 0; k--) {

        /* create a new block */
        struct block *new_block = (struct block *)vrna_alloc(sizeof(struct block));

        /* construct global coordinates for new block */
        new_block->start  = i + start_stem[k] - 1 - b->shift;
        new_block->end    = i + end_stem[k] - 1 - b->shift;
        new_block->shift  = (dangles == 2) ? 1 : 0;

#if DEBUG
        printf("separating branch [%d, %d] into block [%d:%d]\n", start_stem[k], end_stem[k], new_block->start, new_block->end);
#endif

        /* create structure pair table for new block */
        size_t  l         = end_stem[k] - start_stem[k] + 1;
        if (dangles == 2) {
          l++; /* for 5' dangles */
          if (new_block->end < max_n)
            l++; /* for 3' dangles */
        }

        new_block->pt     = (short *)vrna_alloc(sizeof(short) * (l + 1));
        new_block->pt[0]  = l;
        /* go through original structure and extract base pair positions relative to new block */
        for (size_t p = start_stem[k]; p <= end_stem[k]; p++) {
          if (b->pt[p] > p) {
            short i,j;

            i = p - start_stem[k] + 1 + new_block->shift;
            j = b->pt[p] - start_stem[k] + 1 + new_block->shift;

            new_block->pt[i] = j;
            new_block->pt[j] = i;

            /* remove base pair from original block */
            b->pt[b->pt[p]] = 0;
            b->pt[p]        = 0;
          }
        }

        /* create fold_compound for new block */
        char *seq_block = (char *)vrna_alloc(sizeof(char) * (l + 1));
        memcpy(seq_block, b->fc->sequence + start_stem[k] - 1 - new_block->shift, sizeof(char) * (l));
        new_block->fc = vrna_fold_compound(seq_block, &(b->fc->params->model_details), VRNA_OPTION_DEFAULT | VRNA_OPTION_EVAL_ONLY);
        int e = vrna_eval_structure_pt(new_block->fc, new_block->pt);
#if DEBUG
        char *ss = vrna_db_from_ptable(new_block->pt);
        printf("%s\n%s (%6.2f)\n", seq_block, ss, (float)e / 100.);
        free(ss);
#endif
        free(seq_block);

        /* compute energy of new block */
        new_block->energy = e;
        
        /* link new block into list */
        new_block->next_entry = tmp;

        tmp = new_block;
      }

      /* link all newly created blocks to our main list */
      b->next_entry = tmp;

      /* Finally, we need to update the original block */

      /* update global end position */
      b->end = i + end_stem[0] - 1 - b->shift;

      /* update energy */
      b->energy = vrna_eval_structure_pt(b->fc, b->pt);
#if DEBUG
      ss = vrna_db_from_ptable(b->pt);
      printf("%s\n%s (%6.2f)\n", b->fc->sequence, ss, (float)b->energy / 100.);
      free(ss);
#endif
    }

    /* clean up memory */
    free(start_stem);
    free(end_stem);

    /* update for odd dangle models below */
    if (dangles % 1) {
    
    }
  } else if (dangles % 1) {
    /* position i is unpaired, so we only need to update energies if
       position i + 1 forms a base pair due to 5' dangles
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

#if DEBUG
  printf("block post truncation: [%d:%d] (shift %d), e = %d\n", b->start, b->end, b->shift, b->energy);
#endif
}


PRIVATE void
truncate_blocks(unsigned long i,
                unsigned long n,
                struct block  **block_list)
{
  struct block *ptr_prev, *ptr;

  ptr_prev  = NULL;
  ptr       = *block_list;

  while (ptr) {
    /* 1. remove block if it was superseded by index i */
    if (ptr->end <= i) {
      if (ptr_prev) {
        ptr_prev->next_entry = ptr->next_entry;
        vrna_fold_compound_free(ptr->fc);
        free(ptr->pt);
        free(ptr);
        ptr = ptr_prev->next_entry;
      } else {
        /* we are still processing the head of the list */
        *block_list = ptr->next_entry;
        vrna_fold_compound_free(ptr->fc);
        free(ptr->pt);
        free(ptr);
        ptr         = *block_list;
      }

      continue;
    }

    /* 2. Update energy */
    if (ptr->start == i) {
      /* remove nucleotide at position i and update energies accordingly */
      /* this also splits a substructure component into consecutive
         tokens, if necessary
      */
      update_block(i, n, ptr);
    } else if (ptr->start > i) {
      //break;
    }

    /* go to next block */
    ptr_prev  = ptr;
    ptr       = ptr->next_entry;
  }
}


PRIVATE struct block  *
extract_Lfold_entry(FILE            *f,
                    long            line_start,
                    const char      *sequence,
                    const vrna_md_t *md)
{
  char              *l, *seq_local, *structure;
  unsigned long int start, i, j;
  float             en;
  struct block      *storage;

  storage = NULL;

  if (fseek(f, line_start, SEEK_SET) != -1) {
    l         = vrna_read_line(f);
    en        = INF / 100.;
    structure = (char *)vrna_alloc(sizeof(char) * (strlen(l) + 1));

    /*  each line should consist of 3 parts,
        1. A dot-bracket structure
        2. An energy value in kcal/mol
        3. A starting position
    */
    if (sscanf(l, "%[.()] %*c %f %*c %lu", structure, &en, &i) == 3) {
      storage = (struct block *)vrna_alloc(sizeof(struct block));

      j = i + strlen(structure) - 1;

      seq_local = (char *)vrna_alloc(sizeof(char) * (j - i + 2));
      memcpy(seq_local, sequence + i - 1, (sizeof(char) * (j - i + 1)));

#if DEBUG
      printf("%s\n%s (%6.2f), start: %lu\n", seq_local, structure, en, i);
#endif
      storage->fc           = vrna_fold_compound(seq_local, md, VRNA_OPTION_DEFAULT | VRNA_OPTION_EVAL_ONLY);
      storage->pt           = vrna_ptable(structure);
      storage->start        = i;
      storage->end          = j;
      storage->shift        = 0;
      storage->energy       = vrna_convert_kcal_to_dcal(en);
      storage->energy_no3d  = 0; /* energy without 3'dangle */
      storage->next_entry   = NULL;

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


PUBLIC int
vrna_backtrack_window(vrna_fold_compound_t  *fc,
                      const char            *Lfold_filename,
                      long                  file_pos,
                      char                  **structure,
                      double                mfe)
{
  unsigned int  n;
  int           ret, min_en;
  FILE          *f;

  ret         = 0;
  *structure  = NULL;

  if ((fc) &&
      (Lfold_filename) &&
      (structure)) {
    int       maxdist;
    vrna_md_t *md;

    n       = fc->length;
    md      = &(fc->params->model_details);
    maxdist = md->window_size;
    min_en  = vrna_convert_kcal_to_dcal(mfe);

    /* default to start at beginning of file */
    if (file_pos < 0)
      file_pos = 0;

    /* open Lfold output file and start parsing line positions */
    if ((f = fopen(Lfold_filename, "r"))) {
      if (fseek(f, file_pos, SEEK_SET) != -1) {
        size_t  num_lines, mem_lines;
        long    *lines; /* array of file positions indicating start of lines */

        num_lines = 0;
        mem_lines = 1024;

        lines               = (long *)vrna_alloc(sizeof(long) * mem_lines);
        lines[num_lines++]  = ftell(f);

        /* 1. Fill array of file positions that coincide with start of lines */
        do {
          /* increase memory if necessary */
          if (num_lines == mem_lines) {
            mem_lines *= 1.4;
            lines = (long *)vrna_realloc(lines, sizeof(long) * mem_lines);
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

          size_t  block_num = 0;
          size_t  i         = num_lines;

          /* we will use *ptr as current block pointer and save the list start in 8block_list */
          struct block  *block_list, *ptr;

          /* handle first block separately */
          do {
            i--;

            block_list = ptr = extract_Lfold_entry(f, lines[i], fc->sequence, md);
          } while ((block_list == NULL) && (i > 0));

          if (block_list) {
            block_num++;

            /*  go through remaining file in reverse order
                to extract the relevant data
            */
            do {
              i--;

              ptr->next_entry = extract_Lfold_entry(f, lines[i], fc->sequence, md);

              if (ptr->next_entry != NULL) {
                block_num++;
                ptr = ptr->next_entry;
              }
            } while (i > 0);

#if DEBUG
            printf("Read %u blocks\n", block_num);
#endif

            /* start actual backtracking procedure */
            *structure = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));
            memset(*structure, (int)'.', fc->length);

            int *f3 = fc->matrices->f3_local;
            int e   = min_en;

            /* check whether we need to adjust energies due to integer underflows in the forward recursion */
            if (min_en != f3[1]) {
              printf("min_en != f3[1]! %d vs. %d\n", min_en, f3[1]);
            }

            ptr = block_list;

            /* The last block is always part of the full length MFE */
            for (i = ptr->start; i <= ptr->end; i++) {
              size_t i_local = i - ptr->start + ptr->shift + 1;
              if (ptr->pt[i_local] > i_local) {
                (*structure)[i - 1] = '(';
                (*structure)[ptr->start - ptr->shift + ptr->pt[i_local] - 1 - 1] = ')';
              }
            }

#if DEBUG
            printf("inserting block en=%d, [%lu:%lu] (shift %d) => energy: %d + %d = %d\n",
                    ptr->energy,
                    ptr->start,
                    ptr->end,
                    ptr->shift,
                    ptr->energy,
                    f3[ptr->end + 1],
                    e);
#endif

            e -= ptr->energy;
            size_t ii = ptr->end;

            
            for (i = ptr->start; i <= ii; i++)
              /* truncate remaining blocks */
              truncate_blocks((unsigned long)i, n, &block_list);

#if DEBUG
            size_t  cnt = 0;
            for (ptr = block_list; ptr; ptr = ptr->next_entry, cnt++)
              printf("block %u: en=%d, start: %lu, end: %lu, shift: %d\n", cnt, ptr->energy, ptr->start, ptr->end, ptr->shift);
            printf("%u blocks remaining\n", cnt);
#endif

            /* proceed with remaining blocks */
            for (; i <= fc->length; i++) {
#if DEBUG
              printf("e=%d, f3[%d]=%d\n", e, i, f3[i]);
#endif
              /* i pairs with some k, find block representing the substructure enclodes by (i,k) */
              if (e != f3[i + 1]) {
#if DEBUG
                printf("available:\n");
                for (ptr = block_list; ptr; ptr = ptr->next_entry) {
                  if (ptr->start == i)
                    printf("block [%d:%d] (shift %d), e=%d\n",
                    ptr->start, ptr->end, ptr->shift, ptr->energy);
                }
#endif
                /* go through list of blocks that start at i to check which one we can insert */
                for (ptr = block_list; ptr; ptr = ptr->next_entry) {
                  //if (ptr->start > i)
                  //  break;

                  if ((ptr->start == i) &&
                      (e == (ptr->energy + f3[ptr->end + 1]))) {
                    /* found the block, let's insert it */

#if DEBUG
                    printf("inserting block en=%d, [%lu:%lu] (shift %d) => energy: %d + f3[%d](%d) = %d\n",
                            ptr->energy, ptr->start, ptr->end, ptr->shift,
                            ptr->energy,
                            f3[ptr->end + 1],
                            ptr->end + 1,
                            ptr->energy + f3[ptr->end + 1]);
#endif

                    for (ii = ptr->start; ii <= ptr->end; ii++) {
                      size_t i_local = ii - ptr->start + ptr->shift + 1;
                      if (ptr->pt[i_local] > i_local) {
                        (*structure)[ii - 1] = '(';
                        (*structure)[ptr->start - ptr->shift + ptr->pt[i_local] - 1 - 1] = ')';
                      }
                    }

                    e -= ptr->energy;

                    ii = ptr->end;

                    for (size_t iii = ptr->start; iii <= ii; iii++)
                      /* truncate remaining blocks */
                      truncate_blocks((unsigned long)iii, n, &block_list);
                    
                    i = ii;
                    break;
                  }
                }
              } else {
                /* i is unpaired */
              }

              truncate_blocks((unsigned int)i, n, &block_list);
            }

            ret = 1;
          }

        }
#if 0
          /* start backtracing */
          char  *mfe_structure = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));
          memset(mfe_structure, (int)'.', fc->length);

          int *f3 = fc->matrices->f3_local;

          /* The last structure is always part of the full length MFE */
          size_t s = 0;
          for (size_t l = 1; l <= ss[s].length; l++)
            if (ss[s].pt[l] > l) {
              mfe_structure[ss[s].start + l - 1] = '(';
              mfe_structure[ss[s].start + ss[s].pt[l] - 1] = ')';
            }

          /* truncate other structures overlapping with the last one */
          size_t min_s = 1;
          for (unsigned int l = ss[s].start + 1; l < ss[s].start + ss[s].length ; l++) {
            printf("l=%u\n", l);
            for (size_t sss = min_s; sss < num_ss; sss++) {
              /* stop correction if structure starts later */
              if (ss[sss].start > l)
                break;

              /* only correct structure if it is not entirely subsumed */
              if (ss[sss].valid) {
                unsigned int i, j, i_local, j_local;
                i       = ss[sss].start;
                j       = i + ss[sss].length - 1;
                i_local = l - i + 1;
                j_local = ss[sss].pt[i_local];
                char *seq_local = (char *)vrna_alloc(sizeof(char) * (j - i + 2));
                memcpy(seq_local, fc->sequence + i - 1, (sizeof(char) * (j - i + 1)));

                printf("truncating and updating structure %u, interval [%d:%d]\n%s\n%s (%6.2f) shift:%d\n", sss, i, j, seq_local, vrna_db_from_ptable(ss[sss].pt), (float)ss[sss].energy / 100., ss[sss].shift);
                if (j <= l) {
                  ss[sss].valid = 0;
                } else if (j_local != 0) {
                  /* we will remove the base pair (i_local, j_local) */
                  /* compute energy difference of the removal step */
                  /* This requires further work to acknowledge dangle model!!!! */
                  int ediff = vrna_eval_move_pt_simple(seq_local, ss[sss].pt, -i_local, -j_local);
                  /* remove the pair */
                  ss[sss].pt[i_local] = ss[sss].pt[j_local] = 0;
                  /* update energy for the remaining substructure */
                  ss[sss].energy += ediff;
                } else if ((j > l) &&
                           (ss[sss].pt[i_local + 1] != 0)) {
                  /* i is unpaired, we will remove it's 5' dangle contribution to base pair (i + 1, p) */
                  unsigned int  type;
                  int           ediff, d5, d3;
                  
                  j_local = ss[sss].pt[i_local + 1];
                  
                  switch (dangles) {
                    case 2:
                      d5    = fc->sequence_encoding2[i + i_local - 1];
                      d3    = (j_local + 1 <= i_local + ss[sss].length - 1) ? fc->sequence_encoding2[i + j_local + 1 - 1] : -1;
                      type  = vrna_get_ptype_md(fc->sequence_encoding2[i + i_local - 1 + 1],
                                                fc->sequence_encoding2[i + j_local - 1],
                                                &(fc->params->model_details));
                      printf("pair (%d, %d) energy: %d vs. %d (%d, %d, type: %d)\n", i_local + 1, j_local, vrna_E_ext_stem(type, -1, d3, fc->params), vrna_E_ext_stem(type, d5, d3, fc->params), d5, d3, type);
                      ediff = vrna_E_ext_stem(type, -1, d3, fc->params) -
                              vrna_E_ext_stem(type, d5, d3, fc->params);
                      break;

                    case 0:
                      ediff = 0;
                      break;

                    default:
                      ediff = 0;
                      break;
                  }
                  /* update energy for the remaining substructure */
                  ss[sss].energy += ediff;
                } else if (j > l) {
                  /* nothing to do here, since we simply remove an unpaired nucleotide that contributes 0 */
                }
                ss[sss].shift++;
                printf("to\n%s (%6.2f) shift:%d\n", vrna_db_from_ptable(ss[sss].pt), (float)ss[sss].energy / 100., ss[sss].shift);
              }
            }
          }
#endif
      }
      fclose(f);
    }
  }

  return ret;
}
