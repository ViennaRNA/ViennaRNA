/*
 *                Equilibrium probabilities of structures and structure ensembles
 *
 *                Ivo L Hofacker + Ronny Lorenz
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/eval/structures.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/partfunc/global.h"
#include "ViennaRNA/probabilities/structures.h"

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
PRIVATE double
wrap_mean_bp_distance(FLT_OR_DBL  *p,
                      int         length,
                      int         *index);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC double
vrna_pr_structure(vrna_fold_compound_t  *fc,
                  const char            *structure)
{
  if ((fc) &&
      (fc->exp_params) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->q)) {
    unsigned int      n;
    double            e, kT, Q, dG, p;
    vrna_exp_param_t  *params = fc->exp_params;
    n = fc->length;

    if (fc->params->model_details.dangles % 2) {
      /* only compute probabilities with dangles = 2 || 0 */
      int dang_bak = fc->params->model_details.dangles;
      fc->params->model_details.dangles = 2;
      e                                 = (double)vrna_eval_structure(fc, structure);
      fc->params->model_details.dangles = dang_bak;
    } else {
      e = (double)vrna_eval_structure(fc, structure);
    }

    kT  = params->kT / 1000.;
    Q   = params->model_details.circ ? fc->exp_matrices->qo : fc->exp_matrices->q[fc->iindx[1] - n];

    dG = (-log(Q) - n * log(params->pf_scale)) * kT;

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      /* add covariance term */
      e -= vrna_eval_covar_structure(fc, structure);

      /* divide ensemble free energy by number of sequences */
      dG /= fc->n_seq;
    }

    p = exp((dG - e) / kT);

    return p;
  }

  return -1.;
}


PUBLIC double
vrna_pr_energy(vrna_fold_compound_t *fc,
               double               e)
{
  if ((fc) &&
      (fc->exp_params) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->q)) {
    unsigned int      n;
    double            kT, Q, dG, p;
    vrna_exp_param_t  *params = fc->exp_params;
    n = fc->length;

    kT  = params->kT / 1000.;
    Q   = params->model_details.circ ? fc->exp_matrices->qo : fc->exp_matrices->q[fc->iindx[1] - n];

    dG = (-log(Q) - n * log(params->pf_scale)) * kT;

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      dG /= fc->n_seq;

    p = exp((dG - e) / kT);

    return p;
  }

  return -1.;
}


PUBLIC double
vrna_mean_bp_distance_pr(int        length,
                         FLT_OR_DBL *p)
{
  int     *index = vrna_idx_row_wise((unsigned int)length);
  double  d;

  if (!p) {
    vrna_log_warning("vrna_mean_bp_distance_pr: "
                         "p == NULL. "
                         "You need to supply a valid probability matrix");
    return (double)INF / 100.;
  }

  d = wrap_mean_bp_distance(p, length, index);

  free(index);
  return d;
}


PUBLIC double
vrna_mean_bp_distance(vrna_fold_compound_t *vc)
{
  if (!vc) {
    vrna_log_warning("vrna_mean_bp_distance: run vrna_pf_fold first!");
  } else if (!vc->exp_matrices) {
    vrna_log_warning("vrna_mean_bp_distance: exp_matrices == NULL!");
  } else if (!vc->exp_matrices->probs) {
    vrna_log_warning("vrna_mean_bp_distance: probs==NULL!");
  } else {
    return wrap_mean_bp_distance(vc->exp_matrices->probs,
                                 vc->length,
                                 vc->iindx);
  }

  return (double)INF / 100.;
}


PUBLIC double
vrna_ensemble_defect_pt(vrna_fold_compound_t  *fc,
                        const short           *pt)
{
  unsigned int  i, j, n;
  int           ii;
  double        ed = -1.;

  if ((fc) &&
      (pt) &&
      ((unsigned int)pt[0] == fc->length) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->probs)) {
    n = fc->length;

    FLT_OR_DBL  *probs  = fc->exp_matrices->probs;
    int         *idx    = fc->iindx;
    ed = 0.;

    for (i = 1; i <= n; i++) {
      ii = idx[i];
      double pi;

      /* compute probability to be paired */
      for (pi = 0., j = 1; j < i; j++)
        pi += probs[idx[j] - i];

      for (j = i + 1; j <= n; j++)
        pi += probs[ii - j];

      if (pt[i] == 0)
        ed += pi;
      else if ((unsigned int)pt[i] > i)
        ed += 1 - probs[ii - pt[i]];
      else
        ed += 1 - probs[idx[pt[i]] - i];
    }

    ed /= (double)n;
  }

  return ed;
}


PUBLIC double
vrna_ensemble_defect(vrna_fold_compound_t *fc,
                     const char           *structure)
{
  double  ed  = -1.;
  short   *pt = vrna_ptable(structure);

  ed = vrna_ensemble_defect_pt(fc, pt);

  free(pt);

  return ed;
}


PUBLIC double *
vrna_positional_entropy(vrna_fold_compound_t *fc)
{
  double *pos_ent = NULL;

  if ((fc) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->probs)) {
    unsigned int  i, j, n;
    int           *my_iindx, ii;
    FLT_OR_DBL    *probs;
    double        log2, a, p, *pp;

    log2      = log(2.);
    n         = fc->length;
    my_iindx  = fc->iindx;
    probs     = fc->exp_matrices->probs;
    pos_ent   = (double *)vrna_alloc(sizeof(double) * (n + 1));
    pp        = (double *)vrna_alloc(sizeof(double) * (n + 1));

    pos_ent[0] = (double)n;

    for (i = 1; i <= n; i++) {
      ii = my_iindx[i];
      for (j = i + 1; j <= n; j++) {
        p           = (double)probs[ii - j];
        a           = (p > 0.) ? p * log(p) : 0.;
        pos_ent[i]  += a;
        pos_ent[j]  += a;
        pp[i]       += p;
        pp[j]       += p;
      }
    }

    for (i = 1; i <= n; i++) {
      pos_ent[i]  += (pp[i] < 1.) ? (1 - pp[i]) * log(1 - pp[i]) : 0;
      pos_ent[i]  /= -log2;
    }

    free(pp);
  }

  return pos_ent;
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */

PRIVATE double
wrap_mean_bp_distance(FLT_OR_DBL  *p,
                      int         length,
                      int         *index)
{
  int     i, j;
  double  d = 0.;

  /* compute the mean base pair distance in the thermodynamic ensemble */
  /* <d> = \sum_{a,b} p_a p_b d(S_a,S_b)
   * this can be computed from the pair probs p_ij as
   * <d> = \sum_{ij} p_{ij}(1-p_{ij}) */

  for (i = 1; i <= length; i++)
    for (j = i + 1; j <= length; j++)
      d += p[index[i] - j] * (1 - p[index[i] - j]);

  return 2 * d;
}
