/*
 *                partiton function and base pair probabilities
 *                for RNA secvondary structures
 *                of a set of aligned sequences
 *
 *                Ivo L Hofacker
 *                Vienna RNA package
 */

/**
*** \file alipfold.c
**/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MIN */
#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/plotting/probabilities.h"
#include "ViennaRNA/params/ribosum.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/eval/structures.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/backtrack/global.h"
#include "ViennaRNA/partfunc/global.h"
#include "ViennaRNA/structures/problist.h"
#include "ViennaRNA/alifold.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 #################################
 # PUBLIC GLOBAL VARIABLES       #
 #################################
 */

/*
 #################################
 # PRIVATE GLOBAL VARIABLES      #
 #################################
 */

/* some backward compatibility stuff */
PRIVATE vrna_fold_compound_t  *backward_compat_compound = NULL;
PRIVATE int                   backward_compat           = 0;
PRIVATE unsigned short        **backward_compat_a2s     = NULL;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat, backward_compat_a2s)

#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE float
wrap_alipf_fold(const char        **sequences,
                char              *structure,
                plist             **pl,
                vrna_exp_param_t  *parameters,
                int               calculate_bppm,
                int               is_constrained,
                int               is_circular);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
/*-----------------------------------------------------------------*/
PRIVATE float
wrap_alipf_fold(const char        **sequences,
                char              *structure,
                plist             **pl,
                vrna_exp_param_t  *parameters,
                int               calculate_bppm,
                int               is_constrained,
                int               is_circular)
{
  unsigned int          i, n_seq;
  float                 free_energy;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  if (sequences == NULL)
    return 0.;

  for (n_seq = 0; sequences[n_seq]; n_seq++);  /* count the sequences */

  vc = NULL;

  /*
   *  if present, extract model details from provided parameters variable,
   *  to properly initialize the fold compound. Otherwise use default
   *  settings taken from deprecated global variables
   */
  if (parameters)
    vrna_md_copy(&md, &(parameters->model_details));
  else
    set_model_details(&md);

  /* set circular and backtracing options */
  md.circ         = is_circular;
  md.compute_bpp  = calculate_bppm;

  vc = vrna_fold_compound_comparative(sequences, &md, VRNA_OPTION_DEFAULT);

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
    vc->exp_params = vrna_exp_params_comparative(n_seq, &(vc->params->model_details));
  }

  /* propagate global pf_scale into vc->exp_params */
  vc->exp_params->pf_scale = pf_scale;

  if (is_constrained && structure) {
    unsigned int constraint_options = 0;
    constraint_options |= VRNA_CONSTRAINT_DB
                          | VRNA_CONSTRAINT_DB_PIPE
                          | VRNA_CONSTRAINT_DB_DOT
                          | VRNA_CONSTRAINT_DB_X
                          | VRNA_CONSTRAINT_DB_ANG_BRACK
                          | VRNA_CONSTRAINT_DB_RND_BRACK;

    vrna_constraints_add(vc, (const char *)structure, constraint_options);
  }

  if (backward_compat && backward_compat_compound) {
    for (n_seq = 0; n_seq < backward_compat_compound->n_seq; n_seq++)
      free(backward_compat_a2s[n_seq]);
    free(backward_compat_a2s);
    vrna_fold_compound_free(backward_compat_compound);
  }

  backward_compat_compound  = vc;
  iindx                     = backward_compat_compound->iindx;

  /* create alignment-column to sequence position mapping compatibility array */
  backward_compat_a2s = (unsigned short **)vrna_alloc(sizeof(unsigned short *) * (vc->n_seq + 1));
  for (n_seq = 0; n_seq < vc->n_seq; n_seq++) {
    backward_compat_a2s[n_seq] =
      (unsigned short *)vrna_alloc(sizeof(unsigned short) * (vc->length + 2));
    for (i = 1; i <= vc->length; i++)
      backward_compat_a2s[n_seq][i] = (unsigned short)vc->a2s[n_seq][i];
  }

  backward_compat = 1;

  free_energy = vrna_pf(vc, structure);

  /* fill plist */
  if (pl && calculate_bppm)
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);

  return free_energy;
}


/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

PUBLIC float
alipf_fold(const char **sequences,
           char       *structure,
           plist      **pl)
{
  return wrap_alipf_fold(sequences, structure, pl, NULL, do_backtrack, fold_constrained, 0);
}


PUBLIC float
alipf_circ_fold(const char  **sequences,
                char        *structure,
                plist       **pl)
{
  return wrap_alipf_fold(sequences, structure, pl, NULL, do_backtrack, fold_constrained, 1);
}


PUBLIC float
alipf_fold_par(const char       **sequences,
               char             *structure,
               plist            **pl,
               vrna_exp_param_t *parameters,
               int              calculate_bppm,
               int              is_constrained,
               int              is_circular)
{
  return wrap_alipf_fold(sequences,
                         structure,
                         pl,
                         parameters,
                         calculate_bppm,
                         is_constrained,
                         is_circular);
}


PUBLIC FLT_OR_DBL *
alipf_export_bppm(void)
{
  if (backward_compat_compound)
    if (backward_compat_compound->exp_matrices)
      if (backward_compat_compound->exp_matrices->probs)
        return backward_compat_compound->exp_matrices->probs;

  return NULL;
}


PUBLIC FLT_OR_DBL *
export_ali_bppm(void)
{
  if (backward_compat_compound)
    if (backward_compat_compound->exp_matrices)
      if (backward_compat_compound->exp_matrices->probs)
        return backward_compat_compound->exp_matrices->probs;

  return NULL;
}


/*brauch ma nurnoch pscores!*/
PUBLIC char *
alipbacktrack(double *prob)
{
  if (backward_compat_compound) {
    if (backward_compat_compound->exp_matrices) {
      vrna_exp_param_t  *params = backward_compat_compound->exp_params;
      int               n       = backward_compat_compound->length;
      int               n_seq   = backward_compat_compound->n_seq;
      int               *idx    = backward_compat_compound->iindx;
      double            Q       = (double)backward_compat_compound->exp_matrices->q[idx[1] - n];
      char              *s      = vrna_pbacktrack(backward_compat_compound);
      double            e       = (double)vrna_eval_structure(backward_compat_compound, s);
      e -= (double)vrna_eval_covar_structure(backward_compat_compound, s);
      double            fe = (-log(Q) - n * log(params->pf_scale)) * params->kT / (1000.0 * n_seq);
      *prob = exp((fe - e) / params->kT);
      return s;
    }
  }

  return NULL;
}


/*
 * -------------------------------------------------------------------------
 * make arrays used for alipf_fold available to other routines
 */
PUBLIC int
get_alipf_arrays(short          ***S_p,
                 short          ***S5_p,
                 short          ***S3_p,
                 unsigned short ***a2s_p,
                 char           ***Ss_p,
                 FLT_OR_DBL     **qb_p,
                 FLT_OR_DBL     **qm_p,
                 FLT_OR_DBL     **q1k_p,
                 FLT_OR_DBL     **qln_p,
                 short          **pscore_p)
{
  if (backward_compat_compound) {
    if (backward_compat_compound->exp_matrices) {
      if (backward_compat_compound->exp_matrices->qb) {
        *S_p      = backward_compat_compound->S;
        *S5_p     = backward_compat_compound->S5;
        *S3_p     = backward_compat_compound->S3;
        *Ss_p     = backward_compat_compound->Ss;
        *qb_p     = backward_compat_compound->exp_matrices->qb;
        *qm_p     = backward_compat_compound->exp_matrices->qm;
        *q1k_p    = backward_compat_compound->exp_matrices->q1k;
        *qln_p    = backward_compat_compound->exp_matrices->qln;
        *pscore_p = backward_compat_compound->pscore_pf_compat;
        *a2s_p    = backward_compat_a2s;
        return 1;
      }
    }
  }

  return 0;
}


PUBLIC void
free_alipf_arrays(void)
{
  if (backward_compat_compound && backward_compat) {
    vrna_fold_compound_free(backward_compat_compound);
    backward_compat_compound  = NULL;
    backward_compat           = 0;
    iindx                     = NULL;
  }
}


#endif
