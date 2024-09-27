/*
 * local pair probabilities for RNA secondary structures
 *
 * Stephan Bernhart, Ivo L Hofacker
 * Vienna RNA package
 */
/*
 * todo: compute energy z-score for each window
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/plotting/probabilities.h"
#include "ViennaRNA/partfunc/global.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/partfunc/exterior.h"
#include "ViennaRNA/partfunc/internal.h"
#include "ViennaRNA/partfunc/multibranch.h"
#include "ViennaRNA/LPfold.h"
#include "ViennaRNA/Lfold.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/partfunc/local.h"

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

typedef struct {
  int           bpp_print;  /* 1 if pairing probabilities should be written to file-handle, 0 if they are returned as vrna_ep_t */
  int           up_print;   /* 1 if unpaired probabilities should be written to file-handle, 0 if they are returned as array */

  FILE          *fp_pU;
  double        **pU;
  FLT_OR_DBL    bpp_cutoff;
  FILE          *fp_bpp;
  vrna_ep_t     *bpp;
  unsigned int  bpp_max_size;
  unsigned int  bpp_size;
  vrna_ep_t     *stack_prob;
  unsigned int  stack_prob_size;
  unsigned int  stack_prob_max_size;
} default_cb_data;

typedef struct {
  FLT_OR_DBL  *prml;
  FLT_OR_DBL  *prm_l;
  FLT_OR_DBL  *prm_l1;
  double      **pU;
  double      **pUO;
  double      **pUI;
  double      **pUM;
  double      **pUH;
} helper_arrays;

/* soft constraint contributions function (internal-loops) */
typedef FLT_OR_DBL (*sc_int)(vrna_fold_compound_t *,
                            int,
                            int,
                            int,
                            int);

/* QI5 contribution function for unpaired probability computations */
typedef void (*add_QI5)(FLT_OR_DBL **,
                       int,
                       int,
                       FLT_OR_DBL,
                       FLT_OR_DBL);

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#ifdef _OPENMP
#include <omp.h>
#endif

/* some backward compatibility stuff */
PRIVATE vrna_fold_compound_t  *backward_compat_compound = NULL;
PRIVATE int                   backward_compat           = 0;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

#endif
/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE void
alloc_helper_arrays(vrna_fold_compound_t  *vc,
                    int                   ulength,
                    helper_arrays         *aux_arrays,
                    unsigned int          options);


PRIVATE void
free_helper_arrays(vrna_fold_compound_t *vc,
                   int                  ulength,
                   helper_arrays        *aux_arrays,
                   unsigned int         options);


PRIVATE void
compute_probs(vrna_fold_compound_t        *vc,
              int                         j,
              helper_arrays               *aux_arrays,
              int                         ulength,
              vrna_probs_window_f  cb,
              void                        *data,
              unsigned int                options,
              int                         *ov);


PRIVATE void
make_ptypes(vrna_fold_compound_t  *vc,
            int                   j);


PRIVATE void
probability_correction(vrna_fold_compound_t *vc,
                       int                  i);


#if 0
PRIVATE vrna_ep_t *
get_deppp(vrna_fold_compound_t  *vc,
          vrna_ep_t             *pl,
          int                   start);


#endif

PRIVATE void
compute_pU(vrna_fold_compound_t       *vc,
           int                        k,
           int                        ulength,
           helper_arrays              *aux_arrays,
           vrna_probs_window_f cb,
           void                       *data,
           unsigned int               options);


PRIVATE FLT_OR_DBL *
compute_stack_probabilities(vrna_fold_compound_t  *vc,
                            int                   start);


PRIVATE void
return_pU(int                         size,
          int                         i,
          int                         max_size,
          helper_arrays               *aux_arrays,
          vrna_probs_window_f  cb,
          void                        *data,
          unsigned int                options);


PRIVATE void
print_bpp_callback(FLT_OR_DBL *pr,
                   int        size,
                   int        k,
                   void       *data);


PRIVATE void
store_bpp_callback(FLT_OR_DBL *pr,
                   int        size,
                   int        k,
                   void       *data);


#if 0
PRIVATE void
store_stack_prob_callback(FLT_OR_DBL  *pr,
                          int         size,
                          int         k,
                          void        *data);


#endif

PRIVATE void
print_pU_callback(double        *pU,
                  int           size,
                  int           k,
                  int           ulength,
                  unsigned int  type,
                  void          *data);


PRIVATE void
store_pU_callback(double        *pU,
                  int           size,
                  int           k,
                  int           ulength,
                  unsigned int  type,
                  void          *data);


PRIVATE void
backward_compat_callback(FLT_OR_DBL   *pr,
                         int          pr_size,
                         int          i,
                         int          max,
                         unsigned int type,
                         void         *data);


PRIVATE FLT_OR_DBL
sc_contribution(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j,
                int                   k,
                int                   l);


PRIVATE FLT_OR_DBL
sc_dummy(vrna_fold_compound_t *vc,
         int                  i,
         int                  j,
         int                  k,
         int                  l);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_ep_t *
vrna_pfl_fold(const char  *sequence,
              int         window_size,
              int         max_bp_span,
              float       cutoff)
{
  default_cb_data data;

  data.fp_pU                = NULL;
  data.pU                   = NULL;
  data.bpp_cutoff           = (FLT_OR_DBL)cutoff;
  data.fp_bpp               = NULL;
  data.bpp                  = NULL;
  data.bpp_max_size         = 0;
  data.bpp_size             = 0;
  data.stack_prob           = NULL;
  data.stack_prob_max_size  = 0;
  data.stack_prob_size      = 0;
  data.bpp_print            = 0;
  data.up_print             = 0;

  vrna_pfl_fold_cb(sequence, window_size, max_bp_span, &backward_compat_callback, (void *)&data);

  /* resize pair probability list to actual size */
  data.bpp =
    (vrna_ep_t *)vrna_realloc(data.bpp, sizeof(vrna_ep_t) * (data.bpp_size + 1));
  data.bpp[data.bpp_size].i     = 0;
  data.bpp[data.bpp_size].j     = 0;
  data.bpp[data.bpp_size].type  = VRNA_PLIST_TYPE_BASEPAIR;
  data.bpp[data.bpp_size].p     = 0;

  return data.bpp;
}


PUBLIC double **
vrna_pfl_fold_up(const char *sequence,
                 int        ulength,
                 int        window_size,
                 int        max_bp_span)
{
  unsigned int    i;
  double          **pU;
  default_cb_data data;

  pU = NULL;

  if (sequence) {
    i   = strlen(sequence);
    pU  = (double **)vrna_alloc(sizeof(double *) * (i + 2));

    data.fp_pU                = NULL;
    data.pU                   = pU;
    data.bpp_cutoff           = 0.;
    data.fp_bpp               = NULL;
    data.bpp                  = NULL;
    data.bpp_max_size         = 0;
    data.bpp_size             = 0;
    data.stack_prob           = NULL;
    data.stack_prob_max_size  = 0;
    data.stack_prob_size      = 0;
    data.bpp_print            = 0;
    data.up_print             = 0;

    vrna_pfl_fold_up_cb(sequence,
                        ulength,
                        window_size,
                        max_bp_span,
                        &backward_compat_callback,
                        (void *)&data);
  }

  return pU;
}


PRIVATE void
alloc_helper_arrays(vrna_fold_compound_t  *vc,
                    int                   ulength,
                    helper_arrays         *aux_arrays,
                    unsigned int          options)
{
  int i, n;

  n = vc->length;

  aux_arrays->pU  = NULL;
  aux_arrays->pUO = NULL;
  aux_arrays->pUH = NULL;
  aux_arrays->pUI = NULL;
  aux_arrays->pUM = NULL;

  aux_arrays->prm_l   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  aux_arrays->prm_l1  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  aux_arrays->prml    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  if ((options & VRNA_PROBS_WINDOW_UP) && (ulength > 0)) {
    aux_arrays->pU = (double **)vrna_alloc((n + 1) * sizeof(double *));
    for (i = 1; i <= n; i++)
      aux_arrays->pU[i] = (double *)vrna_alloc((MAX2(MAXLOOP, ulength) + 2) * sizeof(double));

    if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
      aux_arrays->pUO = (double **)vrna_alloc((n + 1) * sizeof(double *));
      aux_arrays->pUI = (double **)vrna_alloc((n + 1) * sizeof(double *));
      aux_arrays->pUM = (double **)vrna_alloc((n + 1) * sizeof(double *));
      aux_arrays->pUH = (double **)vrna_alloc((n + 1) * sizeof(double *));
      for (i = 1; i <= n; i++) {
        aux_arrays->pUH[i]  = (double *)vrna_alloc((MAX2(MAXLOOP, ulength) + 2) * sizeof(double));
        aux_arrays->pUI[i]  = (double *)vrna_alloc((MAX2(MAXLOOP, ulength) + 2) * sizeof(double));
        aux_arrays->pUO[i]  = (double *)vrna_alloc((MAX2(MAXLOOP, ulength) + 2) * sizeof(double));
        aux_arrays->pUM[i]  = (double *)vrna_alloc((MAX2(MAXLOOP, ulength) + 2) * sizeof(double));
      }
    }
  }
}


PRIVATE void
free_helper_arrays(vrna_fold_compound_t *vc,
                   int                  ulength,
                   helper_arrays        *aux_arrays,
                   unsigned int         options)
{
  int i, n;

  n = vc->length;

  free(aux_arrays->prm_l);
  free(aux_arrays->prm_l1);
  free(aux_arrays->prml);

  if ((options & VRNA_PROBS_WINDOW_UP) && (ulength > 0)) {
    for (i = 1; i <= n; i++)
      free(aux_arrays->pU[i]);
    free(aux_arrays->pU);

    if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
      for (i = 1; i <= n; i++) {
        free(aux_arrays->pUH[i]);
        free(aux_arrays->pUI[i]);
        free(aux_arrays->pUO[i]);
        free(aux_arrays->pUM[i]);
      }
      free(aux_arrays->pUH);
      free(aux_arrays->pUI);
      free(aux_arrays->pUO);
      free(aux_arrays->pUM);
    }
  }
}


PRIVATE void
return_pU(int                         size,
          int                         i,
          int                         max_size,
          helper_arrays               *aux_arrays,
          vrna_probs_window_f  cb,
          void                        *data,
          unsigned int                options)
{
  if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
    cb(aux_arrays->pUO[i], size, i, max_size, VRNA_PROBS_WINDOW_UP | VRNA_EXT_LOOP, data);
    cb(aux_arrays->pUH[i], size, i, max_size, VRNA_PROBS_WINDOW_UP | VRNA_HP_LOOP, data);
    cb(aux_arrays->pUI[i], size, i, max_size, VRNA_PROBS_WINDOW_UP | VRNA_INT_LOOP, data);
    cb(aux_arrays->pUM[i], size, i, max_size, VRNA_PROBS_WINDOW_UP | VRNA_MB_LOOP, data);
  } else {
    cb(aux_arrays->pU[i], size, i, max_size, VRNA_PROBS_WINDOW_UP | VRNA_ANY_LOOP, data);
  }
}


PRIVATE INLINE void
allocate_dp_matrices(vrna_fold_compound_t *vc,
                     int                  i,
                     unsigned int         options)
{
  char          **ptype;
  int           winSize;
  FLT_OR_DBL    **pR, **q, **qb, **qm, **qm2, **QI5, **qmb, **q2l;
  vrna_mx_pf_t  *mx;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  mx      = vc->exp_matrices;
  pR      = mx->pR;
  q       = mx->q_local;
  qb      = mx->qb_local;
  qm      = mx->qm_local;
  qm2     = mx->qm2_local;
  QI5     = mx->QI5;
  qmb     = mx->qmb;
  q2l     = mx->q2l;
  ptype   = vc->ptype_local;
  winSize = vc->window_size;
  hc      = vc->hc;

  /* allocate new part of arrays */
  pR[i] = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
  pR[i] -= i;
  q[i]  = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
  q[i]  -= i;
  qb[i] = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
  qb[i] -= i;
  qm[i] = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
  qm[i] -= i;
  if (options & VRNA_PROBS_WINDOW_UP) {
    qm2[i]  = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
    qm2[i]  -= i;
    QI5[i]  = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
    qmb[i]  = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
    q2l[i]  = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
  }

  hc->matrix_local[i] = (unsigned char *)vrna_alloc((winSize + 1) * sizeof(unsigned char));
  ptype[i]            = (char *)vrna_alloc((winSize + 1) * sizeof(char));
  ptype[i]            -= i;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = vc->sc;
      if (sc) {
        if (sc->exp_energy_bp_local)
          sc->exp_energy_bp_local[i] = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));

        if (sc->exp_energy_up)
          sc->exp_energy_up[i] = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));

        vrna_sc_update(vc, i, VRNA_OPTION_PF | VRNA_OPTION_WINDOW_F5);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      break;
  }
}


PRIVATE INLINE void
free_dp_matrices(vrna_fold_compound_t *vc,
                 unsigned int         options)
{
  unsigned int  i, n, winSize;
  char          **ptype;
  FLT_OR_DBL    **pR, **q, **qb, **qm, **qm2, **QI5, **qmb, **q2l;
  vrna_mx_pf_t  *mx;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  n       = vc->length;
  winSize = vc->window_size;
  mx      = vc->exp_matrices;
  pR      = mx->pR;
  q       = mx->q_local;
  qb      = mx->qb_local;
  qm      = mx->qm_local;
  ptype   = vc->ptype_local;
  hc      = vc->hc;
  sc      = vc->sc;

  i = 1;
  if (n > winSize + MAXLOOP)
    i = n - winSize - MAXLOOP;

  for (; i <= n; i++) {
    free(pR[i] + i);
    free(q[i] + i);
    free(qb[i] + i);
    free(qm[i] + i);
    pR[i] = NULL;
    q[i]  = NULL;
    qb[i] = NULL;
    qm[i] = NULL;
    if (options & VRNA_PROBS_WINDOW_UP) {
      qm2 = mx->qm2_local;
      QI5 = mx->QI5;
      qmb = mx->qmb;
      q2l = mx->q2l;
      free(qm2[i] + i);
      free(QI5[i]);
      free(qmb[i]);
      free(q2l[i]);
      qm2[i]  = NULL;
      QI5[i]  = NULL;
      qmb[i]  = NULL;
      q2l[i]  = NULL;
    }

    free(hc->matrix_local[i]);
    hc->matrix_local[i] = NULL;
    free(ptype[i] + i);
    ptype[i] = NULL;

    if (sc) {
      if (sc->exp_energy_up)
        free(sc->exp_energy_up[i]);

      if (sc->exp_energy_bp_local)
        free(sc->exp_energy_bp_local[i]);
    }
  }
}


PRIVATE INLINE void
init_dp_matrices(vrna_fold_compound_t *vc,
                 unsigned int         options)
{
  unsigned int j, max_j, winSize;

  winSize = vc->window_size;
  max_j   = MIN2(vc->length, 2 * winSize + MAXLOOP + 2);

  for (j = 1; j <= max_j; j++)
    allocate_dp_matrices(vc, j, options);
}


PRIVATE INLINE void
rotate_dp_matrices(vrna_fold_compound_t *vc,
                   int                  j,
                   unsigned int         options)
{
  size_t        i;
  char          **ptype;
  int           winSize, length;
  FLT_OR_DBL    **pR, **q, **qb, **qm, **qm2, **QI5, **qmb, **q2l;
  vrna_mx_pf_t  *mx;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  length  = vc->length;
  winSize = vc->window_size;
  mx      = vc->exp_matrices;
  pR      = mx->pR;
  q       = mx->q_local;
  qb      = mx->qb_local;
  qm      = mx->qm_local;
  ptype   = vc->ptype_local;
  hc      = vc->hc;
  sc      = vc->sc;

  if (j > 2 * winSize + MAXLOOP + 1) {
    i = (size_t)j - (2 * winSize + MAXLOOP + 1);
    /* free arrays may be faster than pointer rotation and reset to 0 values */
    free(pR[i] + i);
    free(q[i] + i);
    free(qb[i] + i);
    free(qm[i] + i);
    pR[i] = NULL;
    q[i]  = NULL;
    qb[i] = NULL;
    qm[i] = NULL;
    if (options & VRNA_PROBS_WINDOW_UP) {
      qm2 = mx->qm2_local;
      QI5 = mx->QI5;
      qmb = mx->qmb;
      q2l = mx->q2l;
      free(qm2[i] + i);
      free(QI5[i]);
      free(qmb[i]);
      free(q2l[i]);
      qm2[i]  = NULL;
      QI5[i]  = NULL;
      qmb[i]  = NULL;
      q2l[i]  = NULL;
    }

    free(hc->matrix_local[i]);
    hc->matrix_local[i] = NULL;
    free(ptype[i] + i);
    ptype[i] = NULL;

    if (sc) {
      if (sc->exp_energy_up) {
        free(sc->exp_energy_up[i]);
        sc->exp_energy_up[i] = NULL;
      }

      if (sc->exp_energy_bp_local) {
        free(sc->exp_energy_bp_local[i]);
        sc->exp_energy_bp_local[i] = NULL;
      }
    }

    if (j + 1 <= length)
      /* get arrays for next round */
      allocate_dp_matrices(vc, j + 1, options);
  }
}


PRIVATE INLINE void
init_constraints(vrna_fold_compound_t *fc,
                 unsigned int         options VRNA_UNUSED)
{
  unsigned int j, max_j, winSize;

  winSize = fc->window_size;
  max_j   = MIN2(fc->length, 2 * winSize + MAXLOOP + 2);

  for (j = 1; j <= max_j; j++) {
    make_ptypes(fc, j);
    vrna_hc_update(fc, j, VRNA_OPTION_WINDOW_F5);
    vrna_sc_update(fc, j, VRNA_OPTION_PF | VRNA_OPTION_WINDOW_F5);
  }
}


PRIVATE INLINE void
rotate_constraints(vrna_fold_compound_t *fc,
                   int                  j,
                   unsigned int         options VRNA_UNUSED)
{
  if (j + 1 <= (int)fc->length) {
    make_ptypes(fc, j + 1);
    vrna_hc_update(fc, j + 1, VRNA_OPTION_WINDOW_F5);
    vrna_sc_update(fc, j + 1, VRNA_OPTION_PF | VRNA_OPTION_WINDOW_F5);
  }
}


PUBLIC int
vrna_probs_window(vrna_fold_compound_t        *vc,
                  int                         ulength,
                  unsigned int                options,
                  vrna_probs_window_f  cb,
                  void                        *data)
{
  unsigned char       hc_decompose;
  int                 n, i, j, k, maxl, ov, winSize, pairSize, turn;
  FLT_OR_DBL          temp, Qmax, qbt1, **q, **qb, **qm, **qm2, **pR;
  double              max_real, *Fwindow;
  vrna_exp_param_t    *pf_params;
  vrna_md_t           *md;
  vrna_mx_pf_t        *matrices;
  vrna_hc_t           *hc;
  helper_arrays       aux_arrays;
  vrna_mx_pf_aux_el_t aux_mx_el;
  vrna_mx_pf_aux_ml_t aux_mx_ml;

  ov    = 0;
  Qmax  = 0;

  if ((!vc) || (!cb))
    return 0; /* failure */

  if (!vrna_fold_compound_prepare(vc, VRNA_OPTION_PF | VRNA_OPTION_WINDOW)) {
    vrna_log_warning("vrna_probs_window: "
                         "Failed to prepare vrna_fold_compound");
    return 0; /* failure */
  }

  /* here space for initializing everything */

  n         = vc->length;
  pf_params = vc->exp_params;
  md        = &(pf_params->model_details);
  matrices  = vc->exp_matrices;
  winSize   = vc->window_size;
  pairSize  = md->max_bp_span;
  turn      = md->min_loop_size;

  q   = matrices->q_local;
  qb  = matrices->qb_local;
  qm  = matrices->qm_local;
  qm2 = matrices->qm2_local;
  pR  = matrices->pR;

  hc = vc->hc;

  alloc_helper_arrays(vc, ulength, &aux_arrays, options);

  Fwindow =
    (options & VRNA_PROBS_WINDOW_PF) ? (double *)vrna_alloc(sizeof(double) * (winSize + 1)) : NULL;

  /* very short molecule ? */
  if (n < turn + 2) {
    if ((options & VRNA_PROBS_WINDOW_UP) && (ulength > 0)) {
      for (i = 1; i <= n; i++) {
        maxl = MIN2(MAX2(MAXLOOP, ulength), n);
        if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
          for (j = 0; j <= maxl; j++) {
            aux_arrays.pUO[i][j]  = 1.;
            aux_arrays.pUH[i][j]  = 0.;
            aux_arrays.pUI[i][j]  = 0.;
            aux_arrays.pUM[i][j]  = 0.;
          }
        } else {
          for (j = 0; j <= maxl; j++)
            aux_arrays.pU[i][j] = 1.;
        }

        return_pU(maxl, i, ulength, &aux_arrays, cb, data, options);
      }
    }

    free_helper_arrays(vc, ulength, &aux_arrays, options);

    return 1; /* success */
  }

  init_dp_matrices(vc, options);
  init_constraints(vc, options);

  /* init auxiliary arrays for fast exterior/multibranch loops */
  aux_mx_el = vrna_exp_E_ext_fast_init(vc);
  aux_mx_ml = vrna_exp_E_ml_fast_init(vc);

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  /* start recursions */
  for (j = 2; j <= n + winSize; j++) {
    if (j <= n) {
      vrna_exp_E_ext_fast_update(vc, j, aux_mx_el);
      for (i = j - 1; i >= MAX2(1, (j - winSize + 1)); i--) {
        hc_decompose  = hc->matrix_local[i][j - i];
        qbt1          = 0.;

        /*
         * construction of partition function of segment i,j
         * firstly that given i bound to j : qb(i,j)
         */
        if (hc_decompose) {
          /* process hairpin loop(s) */
          qbt1 += vrna_exp_eval_hairpin(vc, i, j, VRNA_EVAL_LOOP_DEFAULT);
          /* process internal loop(s) */
          qbt1 += vrna_exp_E_int_loop(vc, i, j);
          /* process multibranch loop(s) */
          qbt1 += vrna_exp_E_mb_loop_fast(vc, i, j, aux_mx_ml);
        }

        qb[i][j] = qbt1;

        /* Multibranch loop */
        qm[i][j] = vrna_exp_E_ml_fast(vc, i, j, aux_mx_ml);
        if ((options & VRNA_PROBS_WINDOW_UP) && (ulength > 0)) {
          /* new qm2 computation done here */
          const FLT_OR_DBL *qqm = vrna_exp_E_ml_fast_qqm(aux_mx_ml);
          temp = 0.0;
          for (k = i + 1; k <= j; k++)
            temp += qm[i][k - 1] *
                    qqm[k];
          qm2[i][j] = temp;
        }

        /* Exterior loop */
        q[i][j] = temp = vrna_exp_E_ext_fast(vc, i, j, aux_mx_el);

        if (temp > Qmax) {
          Qmax = temp;
          if (Qmax > max_real / 10.) {
            vrna_log_warning("vrna_probs_window: "
                                 "Q close to overflow: %d %d %g\n",
                                 i,
                                 j,
                                 temp);
          }
        }

        if (temp >= max_real) {
          vrna_log_warning("vrna_probs_window: "
                               "overflow while computing partition function for segment q[%d,%d]\n"
                               "use larger pf_scale",
                               i,
                               j);

          vrna_exp_E_ml_fast_free(aux_mx_ml);
          vrna_exp_E_ext_fast_free(aux_mx_el);
          free_helper_arrays(vc, ulength, &aux_arrays, options);

          return 0; /* failure */
        }
      } /* end for i */

      /*
       * here we return the partition function for subsegments [i...j] in terms
       * of ensemble free energies G_ij = -RT * ln(Q_ij) in kcal/mol
       */
      if (options & VRNA_PROBS_WINDOW_PF) {
        int start = MAX2(1, j - winSize + 1);
        Fwindow -= start;
        for (i = start; i <= j; i++)
          Fwindow[i] = (double)(-log(q[i][j]) - (j - i + 1) * log(pf_params->pf_scale)) *
                       pf_params->kT / 1000.0;

        cb(Fwindow, j, start, winSize, VRNA_PROBS_WINDOW_PF, data);
        Fwindow += start;
      }

      /*
       * just as a general service, I save here the free energy of the windows
       * no output is generated, however,...
       */
      if ((j >= winSize) && (options & VRNA_PROBS_WINDOW_UP)) {
        FLT_OR_DBL eee = 0.;
        eee = (FLT_OR_DBL)(-log(q[j - winSize + 1][j]) - winSize * log(pf_params->pf_scale)) *
              pf_params->kT / 1000.0;

        /* we could return this to the user via callback cb() if we were nice */

        aux_arrays.pU[j][0] = eee;
      }

      /* rotate auxiliary arrays */
      vrna_exp_E_ext_fast_rotate(aux_mx_el);
      vrna_exp_E_ml_fast_rotate(aux_mx_ml);
    }

    if (j > winSize) {
      compute_probs(vc, j, &aux_arrays, ulength, cb, data, options, &ov);

      if ((options & VRNA_PROBS_WINDOW_UP) && (j > winSize + MAXLOOP + 1))
        compute_pU(vc, j - winSize - MAXLOOP - 1, ulength, &aux_arrays, cb, data, options);

      if (j > 2 * winSize + MAXLOOP + 1) {
        int start = j - (2 * winSize + MAXLOOP + 1);
        probability_correction(vc, start);
        if (options & VRNA_PROBS_WINDOW_BPP) {
          cb(pR[start],
             MIN2(start + winSize, n),
             start,
             winSize,
             VRNA_PROBS_WINDOW_BPP,
             data);
        }

        if (options & VRNA_PROBS_WINDOW_STACKP) {
          int start = j - (2 * winSize - MAXLOOP);
          if (start > 1) {
            FLT_OR_DBL *stack_probs = compute_stack_probabilities(vc, start);
            stack_probs -= start + 1;
            cb(stack_probs,
               MIN2(n - start, pairSize),
               start,
               winSize,
               VRNA_PROBS_WINDOW_STACKP,
               data);
            stack_probs += start + 1;
            free(stack_probs);
          }
        }

        rotate_dp_matrices(vc, j, options);
        rotate_constraints(vc, j, options);
      }
    } /* end if (do_backtrack) */
  }   /* end for j */

  /* finish output */
  if (options & VRNA_PROBS_WINDOW_UP)
    for (j = MAX2(1, n - MAXLOOP); j <= n; j++)
      compute_pU(vc, j, ulength, &aux_arrays, cb, data, options);

  for (j = MAX2(n - winSize - MAXLOOP, 1); j <= n; j++) {
    probability_correction(vc, j);
    if (options & VRNA_PROBS_WINDOW_BPP) {
      cb(pR[j],
         MIN2(j + winSize, n),
         j,
         winSize,
         VRNA_PROBS_WINDOW_BPP,
         data);
    }

    if ((options & VRNA_PROBS_WINDOW_STACKP) && j < n) {
      int start = j;
      if (start > 1) {
        FLT_OR_DBL *stack_probs = compute_stack_probabilities(vc, start);
        stack_probs -= start + 1;
        cb(stack_probs,
           MIN2(n - start, pairSize),
           start,
           winSize,
           VRNA_PROBS_WINDOW_STACKP,
           data);
        stack_probs += start + 1;
        free(stack_probs);
      }
    }
  }

  if (ov > 0) {
    vrna_log_warning("vrna_probs_window: "
                         "%d overflows occurred while backtracking;\n"
                         "you might try a smaller pf_scale than %g\n",
                         ov,
                         pf_params->pf_scale);
  }

  free_dp_matrices(vc, options);
  free_helper_arrays(vc, ulength, &aux_arrays, options);

  /* free memory occupied by auxiliary arrays for fast exterior/multibranch loops */
  vrna_exp_E_ml_fast_free(aux_mx_ml);
  vrna_exp_E_ext_fast_free(aux_mx_el);

  free(Fwindow);

  return 1; /* success */
}


PRIVATE FLT_OR_DBL
sc_contribution(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j,
                int                   k,
                int                   l)
{
  FLT_OR_DBL  q;
  vrna_sc_t   *sc;

  q   = 1.;
  sc  = vc->sc;

  if (sc->exp_energy_up)
    q *= sc->exp_energy_up[i + 1][k - i - 1] *
         sc->exp_energy_up[l + 1][j - l - 1];

  if (sc->exp_energy_bp_local)
    q *= sc->exp_energy_bp_local[i][j - i];

  if ((sc->exp_energy_stack) && (i + 1 == k) && (l + 1 == j)) {
    q *= sc->exp_energy_stack[i] *
         sc->exp_energy_stack[k] *
         sc->exp_energy_stack[l] *
         sc->exp_energy_stack[j];
  }

  if (sc->exp_f)
    q *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);

  return q;
}


PRIVATE FLT_OR_DBL
sc_dummy(vrna_fold_compound_t *vc VRNA_UNUSED,
         int                  i VRNA_UNUSED,
         int                  j VRNA_UNUSED,
         int                  k VRNA_UNUSED,
         int                  l VRNA_UNUSED)
{
  return 1.;
}


PRIVATE void
add_QI5_contribution(FLT_OR_DBL **QI5,
                     int        i,
                     int        j,
                     FLT_OR_DBL q,
                     FLT_OR_DBL qkl)
{
  QI5[i][j] += q * qkl;
}


PRIVATE void
add_QI5_dummy(FLT_OR_DBL  **QI5 VRNA_UNUSED,
              int         i VRNA_UNUSED,
              int         j VRNA_UNUSED,
              FLT_OR_DBL  q VRNA_UNUSED,
              FLT_OR_DBL  qkl VRNA_UNUSED)
{
  return;
}


PRIVATE void
compute_probs(vrna_fold_compound_t  *vc,
              int                   j,
              helper_arrays         *aux_arrays,
              int                   ulength VRNA_UNUSED,
              vrna_probs_window_f   cb VRNA_UNUSED,
              void                  *data VRNA_UNUSED,
              unsigned int          options,
              int                   *ov)
{
  char              **ptype;
  short             *S1;
  int               start_i, i, k, l, n, m, winSize, type, type_2, tt, *rtype;
  FLT_OR_DBL        *prml, *prm_l, *prm_l1, **pR, **QI5, **qmb, **q2l, **qb, **q, **qm,
                    *scale, *expMLbase, expMLclosing, temp, prm_MLb, prmt1, prmt, *tmp,
                    Qmax;
  double            max_real;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;
  sc_int            sc_int_f;
  add_QI5           add_QI5_f;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  prml    = aux_arrays->prml;
  prm_l   = aux_arrays->prm_l;
  prm_l1  = aux_arrays->prm_l1;

  n             = vc->length;
  winSize       = vc->window_size;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype_local;
  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  rtype         = &(md->rtype[0]);
  expMLclosing  = pf_params->expMLclosing;
  scale         = vc->exp_matrices->scale;
  expMLbase     = vc->exp_matrices->expMLbase;
  hc            = vc->hc;
  sc            = vc->sc;

  pR  = vc->exp_matrices->pR;
  QI5 = vc->exp_matrices->QI5;
  qmb = vc->exp_matrices->qmb;
  q2l = vc->exp_matrices->q2l;
  q   = vc->exp_matrices->q_local;
  qb  = vc->exp_matrices->qb_local;
  qm  = vc->exp_matrices->qm_local;

  Qmax = 0;

  /* assign helper functions */
  if (sc)
    sc_int_f = &sc_contribution;
  else
    sc_int_f = &sc_dummy;

  if (options & VRNA_PROBS_WINDOW_UP)
    add_QI5_f = &add_QI5_contribution;
  else
    add_QI5_f = &add_QI5_dummy;

  /* start recursion */

  /*
   * i=j-winSize;
   * initialize multiloopfs
   */
  for (k = j - winSize; k <= MIN2(n, j); k++) {
    prml[k]   = 0;
    prm_l[k]  = 0;
    /*        prm_l1[k]=0;  others stay */
  }
  k         = j - winSize;
  prm_l1[k] = 0;

  for (l = k + 1; l <= MIN2(n, k + winSize - 1); l++) {
    int a;
    pR[k][l]  = 0; /* set zero at start */
    type      = vrna_get_ptype_window(k, l + k, ptype);

    if (qb[k][l] == 0)
      continue;

    /* Exterior loop cases */
    if (hc->matrix_local[k][l - k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
      for (a = MAX2(1, l - winSize + 2); a < MIN2(k, n - winSize + 2); a++)
        pR[k][l] += q[a][k - 1] *
                    q[l + 1][a + winSize - 1] /
                    q[a][a + winSize - 1];

      if (l - k + 1 == winSize) {
        pR[k][l] += 1. / q[k][l];
      } else {
        if (k + winSize - 1 <= n)    /* k outermost */
          pR[k][l] += q[l + 1][k + winSize - 1] /
                      q[k][k + winSize - 1];

        if (l - winSize + 1 >= 1) /* l outermost */
          pR[k][l] += q[l - winSize + 1][k - 1] /
                      q[l - winSize + 1][l];
      }

      pR[k][l] *= vrna_exp_E_exterior_stem(type,
                                      (k > 1) ? S1[k - 1] : -1,
                                      (l < n) ? S1[l + 1] : -1,
                                      pf_params);

      if ((sc) &&
          (sc->exp_f))
        pR[k][l] *= sc->exp_f(k, l, k, l, VRNA_DECOMP_EXT_STEM, sc->data);
    }

    if (hc->matrix_local[k][l - k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
      FLT_OR_DBL ppp;

      type_2  = rtype[vrna_get_ptype_window(k, l + k, ptype)];
      ppp     = 0.;
      start_i = k - MAXLOOP - 1;

      if (start_i < l - winSize + 1)
        start_i = l - winSize + 1;

      if (start_i < 1)
        start_i = 1;

      unsigned int  u1 = 0;
      short         sk1, sl1, si1;

      sk1 = S1[k - 1];
      sl1 = S1[l + 1];
      for (i = k - 1; i >= start_i; i--, u1++) {
        int max_m = i + winSize - 1;

        if (hc->up_int[i + 1] < u1)
          break;

        si1 = S1[i + 1];

        if (max_m > l + MAXLOOP - u1 + 1)
          max_m = l + MAXLOOP - u1 + 1;

        if (max_m > n)
          max_m = n;

        for (m = l + 1; m <= max_m; m++) {
          unsigned int u2 = m - l - 1;

          if (hc->up_int[l + 1] < u2)
            break;

          if (hc->matrix_local[i][m - i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
            type = vrna_get_ptype_window(i, m + i, ptype);
            if (pR[i][m] > 0) {
              temp = pR[i][m] *
                     vrna_exp_E_internal(u1,
                                   u2,
                                   type,
                                   type_2,
                                   si1,
                                   S1[m - 1],
                                   sk1,
                                   sl1,
                                   pf_params) *
                     sc_int_f(vc, i, m, k, l) *
                     scale[u1 + u2 + 2];

              add_QI5_f(QI5, i, k - i - 1, temp, qb[k][l]);
              add_QI5_f(QI5, l, m - l - 1, temp, qb[k][l]);

              ppp += temp;
            }
          }
        }
      }

      pR[k][l] += ppp;
    }
  }

  /* 3. bonding k,l as substem of multi-loop enclosed by i,m */
  prm_MLb = 0.;
  if (k > 1) {
    /* sonst nix! */
    for (l = MIN2(n - 1, k + winSize - 2); l >= k + 1; l--) {
      FLT_OR_DBL ppp;

      /* opposite direction */
      m     = l + 1;
      prmt  = prmt1 = 0.0;

      for (i = MAX2(1, l - winSize + 2); i < k - 1 /* turn */; i++) {
        if (hc->matrix_local[i][m - i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
          tt  = rtype[vrna_get_ptype_window(i, m + i, ptype)];
          ppp = pR[i][m] *
                vrna_exp_E_multibranch_stem(tt, S1[m - 1], S1[i + 1], pf_params) *
                qm[i + 1][k - 1];

          if (sc) {
            if (sc->exp_energy_bp_local)
              ppp *= sc->exp_energy_bp_local[i][m - i];

            if (sc->exp_f)
              ppp *= sc->exp_f(i, m, i + 1, m - 1, VRNA_DECOMP_PAIR_ML, sc->data);
          }

          prmt += ppp;
        }
      }
      prmt *= expMLclosing;

      prml[m] = prmt;

      if (hc->matrix_local[k - 1][m - k + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        tt    = rtype[vrna_get_ptype_window(k - 1, m + k - 1, ptype)];
        prmt1 = pR[k - 1][m] *
                expMLclosing *
                vrna_exp_E_multibranch_stem(tt,
                             S1[l],
                             S1[k],
                             pf_params);


        if (sc) {
          if (sc->exp_energy_bp_local)
            prmt1 *= sc->exp_energy_bp_local[k - 1][m - k + 1];

          if (sc->exp_f)
            prmt1 *= sc->exp_f(k - 1, m, k, m - 1, VRNA_DECOMP_PAIR_ML, sc->data);
        }
      }

      /* k-1 is unpaired */
      if (hc->up_ml[k - 1]) {
        ppp = prm_l1[m] * expMLbase[1];

        if (sc) {
          if (sc->exp_energy_up)
            ppp *= sc->exp_energy_up[k - 1][1];

          if (sc->exp_f)
            ppp *= sc->exp_f(k - 1, l, k, l, VRNA_DECOMP_ML_ML, sc->data);
        }

        prm_l[m] = ppp + prmt1;
      } else {
        /* skip configuration where k-1 is unpaired */
        prm_l[m] = prmt1;
      }

      /* m is unpaired */
      if (hc->up_ml[m]) {
        ppp = prm_MLb * expMLbase[1];

        if (sc) {
          if (sc->exp_energy_up)
            ppp *= sc->exp_energy_up[m][1];

        }

        prm_MLb = ppp + prml[m];
      } else {
        prm_MLb = prml[m];
      }

      /*
       * same as:    prm_MLb = 0;
       * for (i=n; i>k; i--)  prm_MLb += prml[i]*expMLbase[k-i-1];
       */
      prml[m] = prml[m] + prm_l[m];

      if (qb[k][l] == 0.)
        continue;

      if (hc->matrix_local[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
        tt = vrna_get_ptype_window(k, l + k, ptype);

        if (options & VRNA_PROBS_WINDOW_UP) {
          double dang;
          /* coefficient for computations of unpairedarrays */
          dang = qb[k][l] *
                 vrna_exp_E_multibranch_stem(tt,
                              (k > 1) ? S1[k - 1] : -1,
                              (l < n) ? S1[l + 1] : -1,
                              pf_params) *
                 scale[2];

          if ((sc) &&
              (sc->exp_f)) {
            dang *= sc->exp_f(k, l, k, l, VRNA_DECOMP_ML_STEM, sc->data);
          }

          for (m = MIN2(k + winSize - 2, n); m >= l + 2; m--) {
            qmb[l][m - l - 1] += prml[m] * dang;
            q2l[l][m - l - 1] += (prml[m] - prm_l[m]) * dang;
          }
        }

        temp = prm_MLb;

        for (m = MIN2(k + winSize - 2, n); m >= l + 2; m--)
          temp += prml[m] * qm[l + 1][m - 1];


        temp *= vrna_exp_E_multibranch_stem(tt,
                             (k > 1) ? S1[k - 1] : -1,
                             (l < n) ? S1[l + 1] : -1,
                             pf_params) *
                scale[2];

        if ((sc) &&
            (sc->exp_f)) {
          temp *= sc->exp_f(k, l, k, l, VRNA_DECOMP_ML_STEM, sc->data);
        }

        pR[k][l] += temp;
      }

      if (pR[k][l] > Qmax) {
        Qmax = pR[k][l];
        if (Qmax > max_real / 10.)
          vrna_log_warning("P close to overflow: %d %d %g %g\n",
                               i, m, pR[k][l], qb[k][l]);
      }

      if (pR[k][l] >= max_real) {
        (*ov)++;
        pR[k][l] = FLT_MAX;
      }
    } /* end for (l=..) */
  }

  tmp                 = prm_l1;
  aux_arrays->prm_l1  = prm_l;
  aux_arrays->prm_l   = tmp;
}


PRIVATE void
probability_correction(vrna_fold_compound_t *vc,
                       int                  i)
{
  int         j, howoften, pairdist, n, winSize;
  FLT_OR_DBL  **qb, **pR;

  n         = vc->length;
  winSize   = vc->window_size;
  howoften  = 0; /* how many samples do we have for this pair */

  qb  = vc->exp_matrices->qb_local;
  pR  = vc->exp_matrices->pR;

  for (j = i; j < MIN2(i + winSize, n + 1); j++) {
    pairdist = (j - i + 1);
    /* 4cases */
    howoften  = MIN2(winSize - pairdist + 1, i);  /* pairdist,start */
    howoften  = MIN2(howoften, n - j + 1);        /* end */
    howoften  = MIN2(howoften, n - winSize + 1);  /* windowsize */
    pR[i][j]  *= qb[i][j] / howoften;
  }
  return;
}


PRIVATE void
make_ptypes(vrna_fold_compound_t  *vc,
            int                   i)
{
  /* make new entries in ptype array */
  char        **ptype;
  const short *S;
  int         j, type, pairSize, n;
  vrna_md_t   *md;

  ptype     = vc->ptype_local;
  md        = &(vc->exp_params->model_details);
  pairSize  = md->max_bp_span;
  S         = vc->sequence_encoding2;
  n         = vc->length;

  for (j = i; j <= MIN2(i + pairSize, n); j++) {
    type        = md->pair[S[i]][S[j]];
    ptype[i][j] = (char)type;
  }
  return;
}


#if 0
PRIVATE vrna_ep_t *
get_deppp(vrna_fold_compound_t  *vc,
          vrna_ep_t             *pl,
          int                   start)
{
  /* compute dependent pair probabilities */
  int               i, j, count = 0;
  double            tmp;
  vrna_ep_t         *temp;
  char              **ptype;
  short             *S1;
  FLT_OR_DBL        **qb, *scale;
  int               *rtype, turn, pairsize, length;

  vrna_exp_param_t  *pf_params;

  S1        = vc->sequence_encoding;
  pf_params = vc->exp_params;
  ptype     = vc->ptype_local;
  qb        = vc->exp_matrices->qb_local;
  scale     = vc->exp_matrices->scale;
  rtype     = &(pf_params->model_details.rtype[0]);
  turn      = pf_params->model_details.min_loop_size;
  pairsize  = pf_params->model_details.max_bp_span;
  length    = vc->length;

  temp = (vrna_ep_t *)vrna_alloc(pairsize * sizeof(vrna_ep_t)); /* holds temporary deppp */
  for (j = start + turn; j < MIN2(start + pairsize, length); j++) {
    if ((qb[start][j] * qb[start - 1][(j + 1)]) > 10e-200) {
      int type    = ptype[start - 1][j + 1];
      int type_2  = rtype[(unsigned char)ptype[start][j]];
      tmp = qb[start][j] / qb[start - 1][(j + 1)] * vrna_exp_E_internal(0,
                                                                  0,
                                                                  type,
                                                                  type_2,
                                                                  S1[start],
                                                                  S1[j],
                                                                  S1[start - 1],
                                                                  S1[j + 1],
                                                                  pf_params) * scale[2];
      temp[count].i   = start;
      temp[count].j   = j;
      temp[count++].p = tmp;
    }
  }
  /* write it to list of deppps */
  for (i = 0; pl[i].i != 0; i++);
  pl = (vrna_ep_t *)vrna_realloc(pl, (i + count + 1) * sizeof(vrna_ep_t));
  for (j = 0; j < count; j++) {
    pl[i + j].i = temp[j].i;
    pl[i + j].j = temp[j].j;
    pl[i + j].p = temp[j].p;
  }
  pl[i + count].i = 0;
  pl[i + count].j = 0;
  pl[i + count].p = 0;
  free(temp);
  return pl;
}


#endif

PRIVATE FLT_OR_DBL *
compute_stack_probabilities(vrna_fold_compound_t  *vc,
                            int                   start)
{
  /* compute dependent pair probabilities */
  char              **ptype;
  short             *S1;
  int               j, max_j, *rtype, pairsize, length, type, type_2;
  FLT_OR_DBL        **qb, *scale, *probs;
  double            tmp;
  vrna_exp_param_t  *pf_params;
  vrna_sc_t         *sc;

  length    = vc->length;
  S1        = vc->sequence_encoding;
  pf_params = vc->exp_params;
  ptype     = vc->ptype_local;
  qb        = vc->exp_matrices->qb_local;
  scale     = vc->exp_matrices->scale;
  rtype     = &(pf_params->model_details.rtype[0]);
  pairsize  = pf_params->model_details.max_bp_span;
  sc        = vc->sc;

  max_j = MIN2(start + pairsize, length) - 1;

  probs = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (max_j - start + 1));

  for (j = start + 1; j <= max_j; j++) {
    if ((qb[start][j] * qb[start - 1][(j + 1)]) > 10e-200) {
      type    = vrna_get_ptype_window(start - 1, j + 1 + start - 1, ptype);
      type_2  = rtype[vrna_get_ptype_window(start, j + start, ptype)];
      tmp     = qb[start][j] /
                qb[start - 1][(j + 1)] *
                vrna_exp_E_internal(0,
                              0,
                              type,
                              type_2,
                              S1[start],
                              S1[j],
                              S1[start - 1],
                              S1[j + 1],
                              pf_params) *
                scale[2];

      if (sc) {
        if (sc->exp_energy_stack) {
          tmp *= sc->exp_energy_stack[start] *
                 sc->exp_energy_stack[j] *
                 sc->exp_energy_stack[start - 1] *
                 sc->exp_energy_stack[j + 1];
        }

        if (sc->exp_f)
          tmp *= sc->exp_f(start - 1, j + 1, start, j, VRNA_DECOMP_PAIR_IL, sc->data);
      }

      probs[j - start - 1] = tmp;
    }
  }
  return probs;
}


/*
 * Here: Space for questions...
 */
PRIVATE void
compute_pU(vrna_fold_compound_t       *vc,
           int                        k,
           int                        ulength,
           helper_arrays              *aux_arrays,
           vrna_probs_window_f cb,
           void                       *data,
           unsigned int               options)
{
  /*
   *  here, we try to add a function computing all unpaired probabilities starting at some i,
   *  going down to $unpaired, to be unpaired, i.e. a list with entries from 1 to unpaired for
   *  every i, with the probability of a stretch of length x, starting at i-x+1, to be unpaired
   */
  char              **ptype;
  short             *S1;
  int               startu, i5, j3, len, obp, *rtype, turn, winSize, n, leftmost,
                    rightmost, tt;
  FLT_OR_DBL        expMLclosing, *expMLbase, **q, **qm, **qm2, *scale, **pR, **QI5,
                    **q2l, **qmb;
  double            qqq, temp, *QBE, *QBI, *QBM, *QBH, **pU, **pUO, **pUH, **pUI, **pUM;
  vrna_exp_param_t  *pf_params;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;

  n             = vc->length;
  winSize       = vc->window_size;
  S1            = vc->sequence_encoding;
  pf_params     = vc->exp_params;
  ptype         = vc->ptype_local;
  rtype         = &(pf_params->model_details.rtype[0]);
  scale         = vc->exp_matrices->scale;
  q             = vc->exp_matrices->q_local;
  qm            = vc->exp_matrices->qm_local;
  qm2           = vc->exp_matrices->qm2_local;
  expMLbase     = vc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;
  pR            = vc->exp_matrices->pR;
  QI5           = vc->exp_matrices->QI5;
  q2l           = vc->exp_matrices->q2l;
  qmb           = vc->exp_matrices->qmb;
  turn          = pf_params->model_details.min_loop_size;
  hc            = vc->hc;
  sc            = vc->sc;

  pU  = aux_arrays->pU;
  pUO = aux_arrays->pUO;
  pUH = aux_arrays->pUH;
  pUI = aux_arrays->pUI;
  pUM = aux_arrays->pUM;

  QBE = (double *)vrna_alloc((MAX2(ulength, MAXLOOP) + 2) * sizeof(double));
  QBM = (double *)vrna_alloc((MAX2(ulength, MAXLOOP) + 2) * sizeof(double));
  QBI = (double *)vrna_alloc((MAX2(ulength, MAXLOOP) + 2) * sizeof(double));
  QBH = (double *)vrna_alloc((MAX2(ulength, MAXLOOP) + 2) * sizeof(double));

  /*
   * first, we will
   * for k<=ulength, pU[k][k]=0, because no bp can enclose it
   */

  /* compute pu[k+ulength][ulength] */
  for (i5 = MAX2(k + ulength - winSize + 1, 1); i5 <= k; i5++) {
    for (j3 = k + ulength + 1; j3 <= MIN2(n, i5 + winSize - 1); j3++) {
      /* Multiloops */
      if (hc->matrix_local[i5][j3 - i5] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        tt    = rtype[vrna_get_ptype_window(i5, j3 + i5, ptype)];
        temp  = 0.;
        /*
         * (.. >-----|..........)
         * i5  j     j+ulength  j3
         */
        /* (..{}{}-----|......) */
        if ((hc->up_ml[k + 1] >= j3 - k - 1) && (i5 < k)) {
          qqq = qm2[i5 + 1][k] * expMLbase[j3 - k - 1];

          if (sc) {
            if (sc->exp_energy_up)
              qqq *= sc->exp_energy_up[k + 1][j3 - k - 1];

            if (sc->exp_f)
              qqq *= sc->exp_f(i5 + 1, j3 - 1, i5 + 1, k, VRNA_DECOMP_ML_ML, sc->data);
          }

          temp += qqq;
        }

        /* (..|-----|{}{}) */
        if ((hc->up_ml[i5 + 1] >= k + ulength - i5) && (j3 - 1 > k + ulength)) {
          qqq = qm2[k + ulength + 1][j3 - 1] *
                expMLbase[k + ulength - i5];

          if (sc) {
            if (sc->exp_energy_up)
              qqq *= sc->exp_energy_up[i5 + 1][k + ulength - i5];

            if (sc->exp_f)
              qqq *= sc->exp_f(i5 + 1, j3 - 1, k + ulength + 1, j3 - 1, VRNA_DECOMP_ML_ML, sc->data);
          }

          temp += qqq;
        }

        /* ({}|-----|{}) */
        if ((hc->up_ml[k + 1] >= ulength) && (i5 < k) && (j3 - 1 > k + ulength)) {
          qqq = qm[i5 + 1][k] *
                qm[k + ulength + 1][j3 - 1] *
                expMLbase[ulength];

          if (sc) {
            if (sc->exp_energy_up)
              qqq *= sc->exp_energy_up[k + 1][ulength];

            if (sc->exp_f)
              qqq *= sc->exp_f(i5 + 1, j3 - 1, k, k + ulength + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          }

          temp += qqq;
        }

        /* add dangles, multloopclosing etc. */
        qqq = vrna_exp_E_multibranch_stem(tt,
                           S1[j3 - 1],
                           S1[i5 + 1],
                           pf_params) *
              scale[2] *
              expMLclosing;

        if (sc) {
          if (sc->exp_energy_bp_local)
            qqq *= sc->exp_energy_bp_local[i5][j3 - i5];

          if (sc->exp_f)
            qqq *= sc->exp_f(i5, j3, i5 + 1, j3 - 1, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        temp *= qqq;

        pU[k + ulength][ulength] += temp * pR[i5][j3];

        if (options & VRNA_PROBS_WINDOW_UP_SPLIT)
          pUM[k + ulength][ulength] += temp * pR[i5][j3];
      }

      /* add hairpins */
      if (hc->matrix_local[i5][j3 - i5] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
        temp                      = vrna_exp_eval_hairpin(vc, i5, j3, VRNA_EVAL_LOOP_DEFAULT);
        pU[k + ulength][ulength]  += temp * pR[i5][j3];

        if (options & VRNA_PROBS_WINDOW_UP_SPLIT)
          pUH[k + ulength][ulength] += temp * pR[i5][j3];
      }
    }
  }

  /* Add Interior loop contribution to QBE (and QBI) */
  temp = 0.;
  for (len = winSize; len > MAX2(ulength, MAXLOOP); len--)
    temp += QI5[k][len];

  for (; len > 0; len--) {
    temp      += QI5[k][len];
    QBI[len]  += temp;
    QBE[len]  += temp;
  }

  /* Add Hairpin loop contribution to QBE (and QBH) */
  temp = 0.;
  for (obp = MIN2(n, k + winSize - 1); obp > k + ulength; obp--)
    temp += pR[k][obp] *
            vrna_exp_eval_hairpin(vc, k, obp, VRNA_EVAL_LOOP_DEFAULT);

  for (obp = MIN2(n, MIN2(k + winSize - 1, k + ulength)); obp > k + 1; obp--) {
    temp += pR[k][obp] *
            vrna_exp_eval_hairpin(vc, k, obp, VRNA_EVAL_LOOP_DEFAULT);
    QBH[obp - k - 1]  += temp;
    QBE[obp - k - 1]  += temp;
  }

  /*
   * Add up Multiloopterms  qmb[l][m]+=prml[m]*dang;
   * q2l[l][m]+=(prml[m]-prm_l[m])*dang;
   */

  temp = 0.;

  /* add (()()____) type cont. to I3 */
  if (sc && sc->exp_energy_up) {
    for (len = winSize; len >= ulength; len--)
      if (hc->up_ml[k + 1] >= len) {
        temp += q2l[k][len] *
                expMLbase[len] *
                sc->exp_energy_up[k + 1][len];
      }

    for (; len > 0; len--) {
      if (hc->up_ml[k + 1] >= len) {
        temp += q2l[k][len] *
                expMLbase[len] *
                sc->exp_energy_up[k + 1][len];
      }

      QBM[len]  += temp;
      QBE[len]  += temp;
    }
  } else {
    for (len = winSize; len >= ulength; len--)
      if (hc->up_ml[k + 1] >= len)
        temp += q2l[k][len] *
                expMLbase[len];

    for (; len > 0; len--) {
      if (hc->up_ml[k + 1] >= len)
        temp += q2l[k][len] *
                expMLbase[len];

      QBM[len]  += temp;
      QBE[len]  += temp;
    }
  }

  /* add (()___()) */
  for (len = 1; len < ulength; len++) {
    if (hc->up_ml[k + 1] >= len) {
      for (obp = k + len + turn; obp <= MIN2(n, k + winSize - 1); obp++) {
        temp = qmb[k][obp - k - 1] *
               qm[k + len + 1 /*2*/][obp - 1] *
               expMLbase[len];

        if (sc)
          if (sc->exp_energy_up)
            temp *= sc->exp_energy_up[k + 1][len];

        QBM[len]  += temp;
        QBE[len]  += temp;
      }
    }
  }

  /* add (___()()) */
  for (len = 1; len < ulength; len++) {
    if (hc->up_ml[k + 1] >= len) {
      for (obp = k + len + turn + turn; obp <= MIN2(n, k + winSize - 1); obp++) {
        if (hc->matrix_local[k][obp - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
          tt    = rtype[vrna_get_ptype_window(k, obp + k, ptype)];
          temp  = vrna_exp_E_multibranch_stem(tt, S1[obp - 1], S1[k + 1], pf_params) *
                  scale[2] *
                  expMLbase[len] *
                  expMLclosing *
                  pR[k][obp] *
                  qm2[k + len + 1][obp - 1]; /* k:obp */

          if (sc) {
            if (sc->exp_energy_up)
              temp *= sc->exp_energy_up[k + 1][len];

            if (sc->exp_energy_bp)
              temp *= sc->exp_energy_bp_local[k][obp - k];

            if (sc->exp_f)
              temp *= sc->exp_f(k, obp, k + len + 1, obp - 1, VRNA_DECOMP_PAIR_ML, sc->data);
          }

          QBM[len]  += temp;
          QBE[len]  += temp;
        }
      }
    }
  }

  /*
   * After computing all these contributions in QBE[len], that k is paired
   * and the unpaired stretch is AT LEAST len long, we start to add that to
   * the old unpaired thingies;
   */
  for (len = 1; len <= MIN2(MAX2(ulength, MAXLOOP), n - k); len++)
    pU[k + len][len] += pU[k + len][len + 1] + QBE[len];

  if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
    for (len = 1; len <= MIN2(MAX2(ulength, MAXLOOP), n - k); len++) {
      pUH[k + len][len] += pUH[k + len][len + 1] + QBH[len];
      pUM[k + len][len] += pUM[k + len][len + 1] + QBM[len];
      pUI[k + len][len] += pUI[k + len][len + 1] + QBI[len];
    }
    /* open chain */
    if ((ulength >= winSize) &&
        (k >= ulength) &&
        (hc->up_ext[k - winSize + 1] >= winSize))
      pUO[k][winSize] = scale[winSize] / q[k - winSize + 1][k];
  }

  /* open chain */
  if ((ulength >= winSize) &&
      (k >= ulength) &&
      (hc->up_ext[k - winSize + 1] >= winSize)) {
    temp = scale[winSize] / q[k - winSize + 1][k];

    if (sc) {
      if (sc->exp_energy_up)
        temp *= sc->exp_energy_up[k][winSize];

      if (sc->exp_f)
        temp *= sc->exp_f(k - winSize + 1, k, k - winSize + 1, k, VRNA_DECOMP_EXT_UP, sc->data);
    }

    pU[k][winSize] = temp;
  }

  /*
   * now the not enclosed by any base pair terms for whatever it is we do not need anymore...
   * ... which should be e.g; k, again
   */
  for (startu = MIN2(ulength, k); startu > 0; startu--) {
    temp = 0.;
    /* check whether soft constraint unpaired contributions available */
    if (sc && sc->exp_energy_up) {
      if (hc->up_ext[k - startu + 1] >= startu) {
        for (i5 = MAX2(1, k - winSize + 2); i5 <= MIN2(k - startu, n - winSize + 1); i5++)
          temp += q[i5][k - startu] *
                  q[k + 1][i5 + winSize - 1] *
                  scale[startu] *
                  sc->exp_energy_up[k - startu + 1][startu] /
                  q[i5][i5 + winSize - 1];

        /* the 2 Cases where the borders are on the edge of the interval */
        if ((k >= winSize) && (startu + 1 <= winSize)) {
          temp += q[k - winSize + 1][k - startu] *
                  scale[startu] *
                  sc->exp_energy_up[k - startu + 1][startu] /
                  q[k - winSize + 1][k];
        }

        if ((k <= n - winSize + startu) && (k - startu >= 0) && (k < n) &&
            (startu + 1 <= winSize)) {
          temp += q[k + 1][k - startu + winSize] *
                  scale[startu] *
                  sc->exp_energy_up[k - startu + 1][startu] /
                  q[k - startu + 1][k - startu + winSize];
        }
      }
    } else {
      if (hc->up_ext[k - startu + 1] >= startu) {
        for (i5 = MAX2(1, k - winSize + 2); i5 <= MIN2(k - startu, n - winSize + 1); i5++)
          temp += q[i5][k - startu] *
                  q[k + 1][i5 + winSize - 1] *
                  scale[startu] /
                  q[i5][i5 + winSize - 1];

        /* the 2 Cases where the borders are on the edge of the interval */
        if ((k >= winSize) && (startu + 1 <= winSize))
          temp += q[k - winSize + 1][k - startu] *
                  scale[startu] /
                  q[k - winSize + 1][k];

        if ((k <= n - winSize + startu) && (k - startu >= 0) && (k < n) && (startu + 1 <= winSize))
          temp += q[k + 1][k - startu + winSize] *
                  scale[startu] /
                  q[k - startu + 1][k - startu + winSize];
      }
    }

    /* Divide by number of possible windows */
    leftmost  = MAX2(1, k - winSize + 1);
    rightmost = MIN2(n - winSize + 1, k - startu + 1);

    pU[k][startu] += temp;
    pU[k][startu] /= (rightmost - leftmost + 1);

    if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
      pUO[k][startu] += temp;

      /* Do we want to make a distinction between those? */
      pUO[k][startu]  /= (rightmost - leftmost + 1);
      pUH[k][startu]  /= (rightmost - leftmost + 1);
      pUI[k][startu]  /= (rightmost - leftmost + 1);
      pUM[k][startu]  /= (rightmost - leftmost + 1);
    }
  }
  free(QBE);
  free(QBI);
  free(QBH);
  free(QBM);

  /* call return callback */
  return_pU(MIN2(ulength, k), k, ulength, aux_arrays, cb, data, options);

  return;
}


PRIVATE void
print_bpp_callback(FLT_OR_DBL *pr,
                   int        size,
                   int        k,
                   void       *data)
{
  int         j;
  FILE        *fp     = ((default_cb_data *)data)->fp_bpp;
  FLT_OR_DBL  cutoff  = ((default_cb_data *)data)->bpp_cutoff;

  for (j = k + 1; j <= size; j++) {
    if (pr[j] < cutoff)
      continue;

    fprintf(fp, "%d  %d  %g\n", k, j, pr[j]);
  }
}


PRIVATE void
store_bpp_callback(FLT_OR_DBL *pr,
                   int        size,
                   int        k,
                   void       *data)
{
  int           j;
  vrna_ep_t     *pl         = ((default_cb_data *)data)->bpp;
  unsigned int  pl_size     = ((default_cb_data *)data)->bpp_size;
  unsigned int  pl_max_size = ((default_cb_data *)data)->bpp_max_size;
  FLT_OR_DBL    cutoff      = ((default_cb_data *)data)->bpp_cutoff;

  if (pl_max_size == 0) {
    /* init if necessary */
    pl_max_size = 100;
    pl          = (vrna_ep_t *)vrna_realloc(pl, sizeof(vrna_ep_t) * pl_max_size);
  }

  for (j = k + 1; j <= size; j++) {
    if (pr[j] < cutoff)
      continue;

    /* resize vrna_ep_t memory if necessary */
    if (pl_size >= pl_max_size - 1) {
      pl_max_size *= 1.5;
      pl          = (vrna_ep_t *)vrna_realloc(pl, sizeof(vrna_ep_t) * pl_max_size);
    }

    pl[pl_size].i     = k;
    pl[pl_size].j     = j;
    pl[pl_size].type  = VRNA_PLIST_TYPE_BASEPAIR;
    pl[pl_size++].p   = pr[j];
  }

  /* mark end of vrna_ep_t */
  pl[pl_size].i     = 0;
  pl[pl_size].j     = 0;
  pl[pl_size].type  = VRNA_PLIST_TYPE_BASEPAIR;
  pl[pl_size].p     = 0.;

  /* update data */
  ((default_cb_data *)data)->bpp          = pl;
  ((default_cb_data *)data)->bpp_size     = pl_size;
  ((default_cb_data *)data)->bpp_max_size = pl_max_size;
}


#if 0
PRIVATE void
store_stack_prob_callback(FLT_OR_DBL  *pr,
                          int         size,
                          int         k,
                          void        *data)
{
  int           j;
  vrna_ep_t     *pl         = ((default_cb_data *)data)->stack_prob;
  unsigned int  pl_size     = ((default_cb_data *)data)->stack_prob_size;
  unsigned int  pl_max_size = ((default_cb_data *)data)->stack_prob_max_size;
  FLT_OR_DBL    cutoff      = ((default_cb_data *)data)->bpp_cutoff;

  if (pl_max_size == 0) {
    /* init if necessary */
    pl_max_size = 100;
    pl          = (vrna_ep_t *)vrna_realloc(pl, sizeof(vrna_ep_t) * pl_max_size);
  }

  for (j = k + 1; j <= size; j++) {
    if (pr[j] < cutoff)
      continue;

    /* resize vrna_ep_t memory if necessary */
    if (pl_size >= pl_max_size - 1) {
      pl_max_size *= 1.5;
      pl          = (vrna_ep_t *)vrna_realloc(pl, sizeof(vrna_ep_t) * pl_max_size);
    }

    pl[pl_size].i     = k;
    pl[pl_size].j     = j;
    pl[pl_size].type  = VRNA_PLIST_TYPE_BASEPAIR;
    pl[pl_size++].p   = pr[j];
  }

  /* mark end of vrna_ep_t */
  pl[pl_size].i     = 0;
  pl[pl_size].j     = 0;
  pl[pl_size].type  = VRNA_PLIST_TYPE_BASEPAIR;
  pl[pl_size].p     = 0.;

  /* update data */
  ((default_cb_data *)data)->stack_prob           = pl;
  ((default_cb_data *)data)->stack_prob_size      = pl_size;
  ((default_cb_data *)data)->stack_prob_max_size  = pl_max_size;
}


#endif


PRIVATE void
print_pU_callback(double        *pU,
                  int           size,
                  int           k,
                  int           ulength VRNA_UNUSED,
                  unsigned int  type,
                  void          *data)
{
  if (type & VRNA_PROBS_WINDOW_UP) {
    int   i;
    FILE  *fp = ((default_cb_data *)data)->fp_pU;

    fprintf(fp, "%d\t", k);

    for (i = 1; i < size; i++)
      fprintf(fp, "%.7g\t", pU[i]);
    fprintf(fp, "%.7g", pU[size]);

    if ((type & VRNA_ANY_LOOP) == VRNA_ANY_LOOP)
      fprintf(fp, "\n");
    else if (type & VRNA_EXT_LOOP)
      fprintf(fp, "\tE\n");
    else if (type & VRNA_HP_LOOP)
      fprintf(fp, "\tH\n");
    else if (type & VRNA_INT_LOOP)
      fprintf(fp, "\tI\n");
    else if (type & VRNA_MB_LOOP)
      fprintf(fp, "\tM\n");
    else
      vrna_log_warning("unknown loop type");
  }
}


PRIVATE void
store_pU_callback(double        *pU,
                  int           size,
                  int           k,
                  int           ulength,
                  unsigned int  type,
                  void          *data)
{
  int     i;
  double  **pU_storage = ((default_cb_data *)data)->pU;

  if ((type & VRNA_PROBS_WINDOW_UP) && ((type & VRNA_ANY_LOOP) == VRNA_ANY_LOOP)) {
    pU_storage[k] = (double *)vrna_alloc(sizeof(double) * (ulength + 1));
    for (i = 1; i <= size; i++)
      pU_storage[k][i] = pU[i];
  }
}


PRIVATE void
backward_compat_callback(FLT_OR_DBL   *pr,
                         int          pr_size,
                         int          i,
                         int          max,
                         unsigned int type,
                         void         *data)
{
  default_cb_data *d = (default_cb_data *)data;

  if (type & VRNA_PROBS_WINDOW_BPP) {
    if (d->bpp_print)
      print_bpp_callback(pr, pr_size, i, data);
    else
      store_bpp_callback(pr, pr_size, i, data);
  } else if (type & VRNA_PROBS_WINDOW_UP) {
    if (d->up_print)
      print_pU_callback(pr, pr_size, i, max, type, data);
    else
      store_pU_callback(pr, pr_size, i, max, type, data);
  }
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */

PRIVATE void
putoutpU_prob_old(double            **pU,
                  int               length,
                  int               ulength,
                  FILE              *fp,
                  int               energies,
                  vrna_exp_param_t  *parameters);


PRIVATE void
putoutpU_prob_bin_old(double            **pU,
                      int               length,
                      int               ulength,
                      FILE              *fp,
                      int               energies,
                      vrna_exp_param_t  *parameters);


PRIVATE vrna_ep_t *
wrap_pf_foldLP(char             *sequence,
               int              winSize,
               int              pairSize,
               float            cutoffb,
               double           **pU,
               vrna_ep_t        **dpp2,
               FILE             *pUfp,
               FILE             *spup,
               vrna_exp_param_t *parameters)
{
  int                   ulength, r;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;
  default_cb_data       data;

  vc      = NULL;
  ulength = 0;

  /*
   *  if present, extract model details from provided parameters variable,
   *  to properly initialize the fold compound. Otherwise use default
   *  settings taken from deprecated global variables
   */
  if (parameters)
    vrna_md_copy(&md, &(parameters->model_details));
  else
    set_model_details(&md);

  md.compute_bpp  = 1;        /* turn on base pair probability computations */
  md.window_size  = winSize;  /* set size of sliding window */
  md.max_bp_span  = pairSize; /* set maximum base pair span */

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT | VRNA_OPTION_WINDOW);

  /*
   *  if present, attach a copy of the parameters structure instead of the
   *  default parameters but take care of re-setting it to (initialized)
   *  model details
   */
  free(vc->exp_params);
  if (parameters) {
    vrna_md_copy(&(parameters->model_details), &(vc->params->model_details));
    vc->exp_params = vrna_exp_params_copy(parameters);
  } else {
    vc->exp_params = vrna_exp_params(&(vc->params->model_details));
  }

  /* propagate global pf_scale into vc->exp_params */
  vc->exp_params->pf_scale = pf_scale;

  if (backward_compat_compound && backward_compat)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;
  iindx                     = backward_compat_compound->iindx; /* for backward compatibility and Perl wrapper */

  if (pU)
    ulength = (int)pU[0][0] + 0.49;

  data.fp_pU                = pUfp;
  data.pU                   = pU;
  data.bpp_cutoff           = (FLT_OR_DBL)cutoffb;
  data.fp_bpp               = spup;
  data.bpp                  = NULL;
  data.bpp_max_size         = 0;
  data.bpp_size             = 0;
  data.stack_prob           = NULL;
  data.stack_prob_max_size  = 0;
  data.stack_prob_size      = 0;

  data.bpp_print  = (spup) ? 1 : 0;
  data.up_print   = (pUfp) ? 1 : 0;

  unsigned int options = VRNA_PROBS_WINDOW_BPP; /* always compute base pair probabilities */

  if (dpp2 && (*dpp2))
    options |= VRNA_PROBS_WINDOW_STACKP;

  if (ulength > 0)
    options |= VRNA_PROBS_WINDOW_UP;

  r = vrna_probs_window(vc, ulength, options, &backward_compat_callback, (void *)&data);

  if (!r)
    return NULL;

  if (dpp2 && (*dpp2)) {
    data.stack_prob = (vrna_ep_t *)vrna_realloc(data.stack_prob,
                                                sizeof(vrna_ep_t) *
                                                (data.stack_prob_size + 1));
    data.stack_prob[data.stack_prob_size].i     = 0;
    data.stack_prob[data.stack_prob_size].j     = 0;
    data.stack_prob[data.stack_prob_size].type  = VRNA_PLIST_TYPE_BASEPAIR;
    data.stack_prob[data.stack_prob_size].p     = 0;
    free(*dpp2); /* free already occupied memory */
    *dpp2 = data.stack_prob;
  }

  if (!spup) {
    data.bpp =
      (vrna_ep_t *)vrna_realloc(data.bpp, sizeof(vrna_ep_t) * (data.bpp_size + 1));
    data.bpp[data.bpp_size].i     = 0;
    data.bpp[data.bpp_size].j     = 0;
    data.bpp[data.bpp_size].type  = VRNA_PLIST_TYPE_BASEPAIR;
    data.bpp[data.bpp_size].p     = 0;
    return data.bpp;
  } else {
    return NULL;
  }
}


PUBLIC void
init_pf_foldLP(int length VRNA_UNUSED)
{
  ;/* DO NOTHING */
}


PUBLIC void
update_pf_paramsLP(int length VRNA_UNUSED)
{
  if (backward_compat_compound && backward_compat) {
    vrna_md_t md;
    set_model_details(&md);
    vrna_exp_params_reset(backward_compat_compound, &md);

    /* compatibility with RNAup, may be removed sometime */
    pf_scale = backward_compat_compound->exp_params->pf_scale;
  }
}


PUBLIC void
update_pf_paramsLP_par(int              length VRNA_UNUSED,
                       vrna_exp_param_t *parameters)
{
  if (backward_compat_compound && backward_compat) {
    vrna_md_t md;
    if (parameters) {
      vrna_exp_params_subst(backward_compat_compound, parameters);
    } else {
      set_model_details(&md);
      vrna_exp_params_reset(backward_compat_compound, &md);
    }

    /* compatibility with RNAup, may be removed sometime */
    pf_scale = backward_compat_compound->exp_params->pf_scale;
  }
}


PUBLIC vrna_ep_t *
pfl_fold(char       *sequence,
         int        winSize,
         int        pairSize,
         float      cutoffb,
         double     **pU,
         vrna_ep_t  **dpp2,
         FILE       *pUfp,
         FILE       *spup)
{
  return wrap_pf_foldLP(sequence, winSize, pairSize, cutoffb, pU, dpp2, pUfp, spup, NULL);
}


PUBLIC vrna_ep_t *
pfl_fold_par(char             *sequence,
             int              winSize,
             int              pairSize,
             float            cutoffb,
             double           **pU,
             vrna_ep_t        **dpp2,
             FILE             *pUfp,
             FILE             *spup,
             vrna_exp_param_t *parameters)
{
  return wrap_pf_foldLP(sequence, winSize, pairSize, cutoffb, pU, dpp2, pUfp, spup, parameters);
}


PUBLIC void
putoutpU_prob(double  **pU,
              int     length,
              int     ulength,
              FILE    *fp,
              int     energies)
{
  if (backward_compat_compound && backward_compat)
    putoutpU_prob_old(pU, length, ulength, fp, energies, backward_compat_compound->exp_params);
  else
    vrna_log_warning("putoutpU_prob: Not doing anything! First, run pfl_fold()!");
}


PUBLIC void
putoutpU_prob_par(double            **pU,
                  int               length,
                  int               ulength,
                  FILE              *fp,
                  int               energies,
                  vrna_exp_param_t  *parameters)
{
  if ((pU) && (fp) && (parameters))
    putoutpU_prob_old(pU, length, ulength, fp, energies, parameters);
}


PRIVATE void
putoutpU_prob_old(double            **pU,
                  int               length,
                  int               ulength,
                  FILE              *fp,
                  int               energies,
                  vrna_exp_param_t  *parameters)
{
  /* put out unpaireds */
  int     i, k;
  double  temp, kT = parameters->kT / 1000.0;

  if (energies)
    fprintf(fp, "#opening energies\n #i$\tl=");
  else
    fprintf(fp, "#unpaired probabilities\n #i$\tl=");

  for (i = 1; i <= ulength; i++)
    fprintf(fp, "%d\t", i);
  fprintf(fp, "\n");

  for (k = 1; k <= length; k++) {
    fprintf(fp, "%d\t", k);
    for (i = 1; i <= ulength; i++) {
      if (i > k) {
        fprintf(fp, "NA\t");
        continue;
      }

      if (energies)
        temp = -log(pU[k][i]) * kT;
      else
        temp = pU[k][i];

      fprintf(fp, "%.7g\t", temp);
    }
    fprintf(fp, "\n");
    free(pU[k]);
  }
  fflush(fp);
}


PUBLIC void
putoutpU_prob_bin(double  **pU,
                  int     length,
                  int     ulength,
                  FILE    *fp,
                  int     energies)
{
  if (backward_compat_compound && backward_compat)
    putoutpU_prob_bin_old(pU, length, ulength, fp, energies, backward_compat_compound->exp_params);
  else
    vrna_log_warning("putoutpU_prob_bin: Not doing anything! First, run pfl_fold()!");
}


PUBLIC void
putoutpU_prob_bin_par(double            **pU,
                      int               length,
                      int               ulength,
                      FILE              *fp,
                      int               energies,
                      vrna_exp_param_t  *parameters)
{
  if ((pU) && (fp) && (parameters))
    putoutpU_prob_bin_old(pU, length, ulength, fp, energies, parameters);
}


PRIVATE void
putoutpU_prob_bin_old(double            **pU,
                      int               length,
                      int               ulength,
                      FILE              *fp,
                      int               energies VRNA_UNUSED,
                      vrna_exp_param_t  *parameters)
{
  /* put out unpaireds */
  int     i, k, *p;
  double  kT = parameters->kT / 1000.0;

  p = (int *)vrna_alloc(sizeof(int) * 1);
  /* write first line */
  p[0] = ulength; /* u length */
  fwrite(p, sizeof(int), 1, fp);
  p[0] = length;  /* seq length */
  fwrite(p, sizeof(int), 1, fp);
  for (k = 3; k <= (length + 20); k++) {
    /* all the other lines are set to 1000000 because we are at ulength=0 */
    p[0] = 1000000;
    fwrite(p, sizeof(int), 1, fp);
  }
  /* data */
  for (i = 1; i <= ulength; i++) {
    for (k = 1; k <= 11; k++) {
      /* write first ten entries to 1000000 */
      p[0] = 1000000;
      fwrite(p, sizeof(int), 1, fp);
    }
    for (k = 1; k <= length; k++) {
      /* write data now */
      if (i > k) {
        p[0] = 1000000;         /* check if u > pos */
        fwrite(p, sizeof(int), 1, fp);
        continue;
      } else {
        p[0] = (int)rint(100 * (-log(pU[k][i]) * kT));
        fwrite(p, sizeof(int), 1, fp);
      }
    }
    for (k = 1; k <= 9; k++) {
      /* finish by writing the last 10 entries */
      p[0] = 1000000;
      fwrite(p, sizeof(int), 1, fp);
    }
  }
  /* free pU array; */
  for (k = 1; k <= length; k++)
    free(pU[k]);
  free(p);
  fflush(fp);
}


#endif
