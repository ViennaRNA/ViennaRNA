/*
                  partiton function and base pair probabilities
                  for RNA secvondary structures
                  of a set of aligned sequences

                  Ivo L Hofacker
                  Vienna RNA package
*/

/**
*** \file alipfold.c
**/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MIN */
#include <limits.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/ribo.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/mfe.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/structure_utils.h"
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

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE float     wrap_alipf_fold(const char **sequences,
                                  char *structure,
                                  plist **pl,
                                  vrna_exp_param_t *parameters,
                                  int calculate_bppm,
                                  int is_constrained,
                                  int is_circular);



/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC float
vrna_pf_alifold(const char **strings,
                char *structure,
                vrna_plist_t **pl){

  float                 free_energy;
  double                mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  /* no need to backtrack MFE structure */
  md.backtrack = 0;

  if(!pl){ /* no need for pair probability computations if we do not store them somewhere */
    md.compute_bpp = 0;
  }

  vc  = vrna_fold_compound_comparative(strings, &md, VRNA_OPTION_DEFAULT);
  mfe = (double)vrna_pf(vc, structure);
  vrna_exp_params_rescale(vc, &mfe);
  free_energy = vrna_pf(vc, structure);

  /* fill plist */
  if(pl){
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);
  }

  vrna_fold_compound_free(vc);

  return free_energy;
}

PUBLIC float
vrna_pf_circalifold(const char **sequences,
                    char *structure,
                    vrna_plist_t **pl){

  float                 free_energy;
  double                mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  md.circ = 1;

  /* no need to backtrack MFE structure */
  md.backtrack = 0;

  if(!pl){ /* no need for pair probability computations if we do not store them somewhere */
    md.compute_bpp = 0;
  }

  vc  = vrna_fold_compound_comparative(sequences, &md, VRNA_OPTION_DEFAULT);
  mfe = (double)vrna_mfe(vc, structure);
  vrna_exp_params_rescale(vc, &mfe);
  free_energy = vrna_pf(vc, structure);

  /* fill plist */
  if(pl){
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);
  }

  vrna_fold_compound_free(vc);

  return free_energy;
}



/*-----------------------------------------------------------------*/
PRIVATE float
wrap_alipf_fold(const char **sequences,
                char *structure,
                plist **pl,
                vrna_exp_param_t *parameters,
                int calculate_bppm,
                int is_constrained,
                int is_circular){

  int                 n_seq;
  float               free_energy;
  vrna_fold_compound_t  *vc;
  vrna_exp_param_t    *exp_params;

  if(sequences == NULL) return 0.;

  for(n_seq=0;sequences[n_seq];n_seq++); /* count the sequences */
  
  vc                  = NULL;

  /* we need vrna_exp_param_t datastructure to correctly init default hard constraints */
  if(parameters)
    exp_params = vrna_exp_params_copy(parameters);
  else{
    vrna_md_t md;
    set_model_details(&md); /* get global default parameters */
    exp_params = vrna_exp_params_comparative(n_seq, &md);
  }
  exp_params->model_details.circ        = is_circular;
  exp_params->model_details.compute_bpp = calculate_bppm;

  vc = vrna_fold_compound_comparative(sequences, &(exp_params->model_details), VRNA_OPTION_PF);

  if(parameters){ /* replace exp_params if necessary */
    free(vc->exp_params);
    vc->exp_params = exp_params;
  } else {
    free(exp_params);
  }
  vc->exp_params->pf_scale = pf_scale;

  if(is_constrained && structure){
    unsigned int constraint_options = 0;
    constraint_options |= VRNA_CONSTRAINT_DB
                          | VRNA_CONSTRAINT_DB_PIPE
                          | VRNA_CONSTRAINT_DB_DOT
                          | VRNA_CONSTRAINT_DB_X
                          | VRNA_CONSTRAINT_DB_ANG_BRACK
                          | VRNA_CONSTRAINT_DB_RND_BRACK;

    vrna_constraints_add(vc, (const char *)structure, constraint_options);
  }

  if(backward_compat_compound && backward_compat_compound)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  iindx                     = backward_compat_compound->iindx;
  backward_compat           = 1;

  free_energy = vrna_pf(vc, structure);
  
  /* fill plist */
  if(pl && calculate_bppm){
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);
  }

  return free_energy;
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC float
alipf_fold( const char **sequences,
                  char *structure,
                  plist **pl){

  return wrap_alipf_fold(sequences, structure, pl, NULL, do_backtrack, fold_constrained, 0);
}

PUBLIC float
alipf_circ_fold(const char **sequences,
                      char *structure,
                      plist **pl){

  return wrap_alipf_fold(sequences, structure, pl, NULL, do_backtrack, fold_constrained, 1);
}

PUBLIC float
alipf_fold_par( const char **sequences,
                char *structure,
                plist **pl,
                vrna_exp_param_t *parameters,
                int calculate_bppm,
                int is_constrained,
                int is_circular){

  return wrap_alipf_fold(sequences, structure, pl, parameters, calculate_bppm, is_constrained, is_circular);
}

PUBLIC FLT_OR_DBL *
alipf_export_bppm(void){

  if(backward_compat_compound)
    if(backward_compat_compound->exp_matrices)
      if(backward_compat_compound->exp_matrices->probs)
        return backward_compat_compound->exp_matrices->probs;

  return NULL;
}

PUBLIC FLT_OR_DBL *
export_ali_bppm(void){

  if(backward_compat_compound)
    if(backward_compat_compound->exp_matrices)
      if(backward_compat_compound->exp_matrices->probs)
        return backward_compat_compound->exp_matrices->probs;

  return NULL;
}

/*brauch ma nurnoch pscores!*/
PUBLIC char *
alipbacktrack(double *prob){

  if(backward_compat_compound)
    if(backward_compat_compound->exp_matrices){
      vrna_exp_param_t *params = backward_compat_compound->exp_params;
      int n     = backward_compat_compound->length;
      int n_seq = backward_compat_compound->n_seq;
      int *idx  = backward_compat_compound->iindx;
      double Q  = (double)backward_compat_compound->exp_matrices->q[idx[1]-n];
      char *s   = vrna_pbacktrack(backward_compat_compound);
      double e  = (double)vrna_eval_structure(backward_compat_compound, s);
      e        -= (double)vrna_eval_covar_structure(backward_compat_compound, s);
      double fe = (-log(Q)-n*log(params->pf_scale))*params->kT/(1000.0 * n_seq);
      *prob     = exp((fe - e)/params->kT);
      return s;
    }
  return NULL;
}

/*-------------------------------------------------------------------------*/
/* make arrays used for alipf_fold available to other routines */
PUBLIC int
get_alipf_arrays( short ***S_p,
                  short ***S5_p,
                  short ***S3_p,
                  unsigned short ***a2s_p,
                  char ***Ss_p,
                  FLT_OR_DBL **qb_p,
                  FLT_OR_DBL **qm_p,
                  FLT_OR_DBL **q1k_p,
                  FLT_OR_DBL **qln_p,
                  short **pscore_p) {

  if(backward_compat_compound){
    if(backward_compat_compound->exp_matrices)
      if(backward_compat_compound->exp_matrices->qb){
        *S_p      = backward_compat_compound->S;
        *S5_p     = backward_compat_compound->S5;
        *S3_p     = backward_compat_compound->S3;
        *a2s_p    = backward_compat_compound->a2s;
        *Ss_p     = backward_compat_compound->Ss;
        *qb_p     = backward_compat_compound->exp_matrices->qb;
        *qm_p     = backward_compat_compound->exp_matrices->qm;
        *q1k_p    = backward_compat_compound->exp_matrices->q1k;
        *qln_p    = backward_compat_compound->exp_matrices->qln;
        *pscore_p = backward_compat_compound->pscore_pf_compat;
        return 1;
      }
  }
  return 0;
}

PUBLIC void
free_alipf_arrays(void){

  if(backward_compat_compound && backward_compat){
    vrna_fold_compound_free(backward_compat_compound);
    backward_compat_compound  = NULL;
    backward_compat           = 0;
    iindx                     = NULL;
  }
}
