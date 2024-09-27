/*
 *                partiton function for single RNA secondary structures
 *
 *                Simplified interfaces and backward compatibility
 *                wrappers
 *
 *                Ivo L Hofacker + Ronny Lorenz
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

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
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/probabilities/basepairs.h"
#include "ViennaRNA/probabilities/structures.h"
#include "ViennaRNA/structures/centroid.h"
#include "ViennaRNA/sampling/basic.h"
#include "ViennaRNA/partfunc/global.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */
PUBLIC int st_back = 0;

/*
 #################################
 # PRIVATE VARIABLES             #
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

PRIVATE float
wrap_pf_fold(const char       *sequence,
             char             *structure,
             vrna_exp_param_t *parameters,
             int              calculate_bppm,
             int              is_constrained,
             int              is_circular);


PRIVATE double
wrap_mean_bp_distance(FLT_OR_DBL  *p,
                      int         length,
                      int         *index,
                      int         turn);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE double
wrap_mean_bp_distance(FLT_OR_DBL  *p,
                      int         length,
                      int         *index,
                      int         turn)
{
  int     i, j;
  double  d = 0.;

  /* compute the mean base pair distance in the thermodynamic ensemble */
  /* <d> = \sum_{a,b} p_a p_b d(S_a,S_b)
   * this can be computed from the pair probs p_ij as
   * <d> = \sum_{ij} p_{ij}(1-p_{ij}) */

  for (i = 1; i <= length; i++)
    for (j = i + turn + 1; j <= length; j++)
      d += p[index[i] - j] * (1 - p[index[i] - j]);

  return 2 * d;
}


PRIVATE float
wrap_pf_fold(const char       *sequence,
             char             *structure,
             vrna_exp_param_t *parameters,
             int              calculate_bppm,
             int              is_constrained,
             int              is_circular)
{
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vc = NULL;

  /* we need vrna_exp_param_t datastructure to correctly init default hard constraints */
  if (parameters)
    md = parameters->model_details;
  else
    set_model_details(&md); /* get global default parameters */

  md.circ         = is_circular;
  md.compute_bpp  = calculate_bppm;

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT);

  /* prepare exp_params and set global pf_scale */
  vc->exp_params            = vrna_exp_params(&(vc->params->model_details));
  vc->exp_params->pf_scale  = pf_scale;

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

  if (backward_compat_compound && backward_compat)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;
  iindx                     = backward_compat_compound->iindx;

  return vrna_pf(vc, structure);
}


PUBLIC vrna_ep_t *
stackProb(double cutoff)
{
  if (!(backward_compat_compound && backward_compat)) {
    vrna_log_warning("stackProb: "
                         "run pf_fold() first!");
    return NULL;
  } else if (!backward_compat_compound->exp_matrices->probs) {
    vrna_log_warning("stackProb: "
                         "probs == NULL!");
    return NULL;
  }

  return vrna_stack_prob(backward_compat_compound, cutoff);
}


PUBLIC char *
centroid(int    length,
         double *dist)
{
  if (pr == NULL) {
    vrna_log_warning("centroid: "
                         "pr == NULL. You need to call pf_fold() before centroid()");
    return NULL;
  }

  return vrna_centroid_from_probs(length, dist, pr);
}


PUBLIC double
mean_bp_dist(int length)
{
  /* compute the mean base pair distance in the thermodynamic ensemble */
  /* <d> = \sum_{a,b} p_a p_b d(S_a,S_b)
   * this can be computed from the pair probs p_ij as
   * <d> = \sum_{ij} p_{ij}(1-p_{ij}) */

  int     i, j, *my_iindx;
  double  d = 0;

  if (pr == NULL) {
    vrna_log_warning("mean_bp_dist: "
                         "pr == NULL. You need to call pf_fold() before mean_bp_dist()");
    return d;
  }

  my_iindx = vrna_idx_row_wise(length);

  for (i = 1; i <= length; i++)
    for (j = i + TURN + 1; j <= length; j++)
      d += pr[my_iindx[i] - j] * (1 - pr[my_iindx[i] - j]);

  free(my_iindx);
  return 2 * d;
}


/* get the free energy of a subsequence from the q[] array */
PUBLIC double
get_subseq_F(int  i,
             int  j)
{
  if (backward_compat_compound) {
    if (backward_compat_compound->exp_matrices) {
      if (backward_compat_compound->exp_matrices->q) {
        int               *my_iindx   = backward_compat_compound->iindx;
        vrna_exp_param_t  *pf_params  = backward_compat_compound->exp_params;
        FLT_OR_DBL        *q          = backward_compat_compound->exp_matrices->q;
        return (-log(q[my_iindx[i] - j]) - (j - i + 1) * log(pf_params->pf_scale)) * pf_params->kT /
               1000.0;
      }
    }
  }

  vrna_log_warning("get_subseq_F: "
                       "call pf_fold() to fill q[] array before calling get_subseq_F()");

  return 0.; /* we will never get to this point */
}


/*----------------------------------------------------------------------*/
PUBLIC double
expHairpinEnergy(int        u,
                 int        type,
                 short      si1,
                 short      sj1,
                 const char *string)
{
  /* compute Boltzmann weight of a hairpin loop, multiply by scale[u+2] */

  vrna_exp_param_t  *pf_params = backward_compat_compound->exp_params;

  double            q, kT;

  kT = pf_params->kT;   /* kT in cal/mol  */
  if (u <= 30)
    q = pf_params->exphairpin[u];
  else
    q = pf_params->exphairpin[30] * exp(-(pf_params->lxc * log(u / 30.)) * 10. / kT);

  if ((tetra_loop) && (u == 4)) {
    char tl[7] = {
      0
    }, *ts;
    strncpy(tl, string, 6);
    if ((ts = strstr(pf_params->Tetraloops, tl)))
      return pf_params->exptetra[(ts - pf_params->Tetraloops) / 7];
  }

  if ((tetra_loop) && (u == 6)) {
    char tl[9] = {
      0
    }, *ts;
    strncpy(tl, string, 6);
    if ((ts = strstr(pf_params->Hexaloops, tl)))
      return pf_params->exphex[(ts - pf_params->Hexaloops) / 9];
  }

  if (u == 3) {
    char tl[6] = {
      0
    }, *ts;
    strncpy(tl, string, 5);
    if ((ts = strstr(pf_params->Triloops, tl)))
      return pf_params->exptri[(ts - pf_params->Triloops) / 6];

    if (type > 2)
      q *= pf_params->expTermAU;
  } else {
    /* no mismatches for tri-loops */
    q *= pf_params->expmismatchH[type][si1][sj1];
  }

  return q;
}


PUBLIC double
expLoopEnergy(int   u1,
              int   u2,
              int   type,
              int   type2,
              short si1,
              short sj1,
              short sp1,
              short sq1)
{
  /* compute Boltzmann weight of internal loop,
   * multiply by scale[u1+u2+2] for scaling */
  double            z           = 0;
  int               no_close    = 0;
  vrna_exp_param_t  *pf_params  = backward_compat_compound->exp_params;


  if ((no_closingGU) && ((type2 == 3) || (type2 == 4) || (type == 2) || (type == 4)))
    no_close = 1;

  if ((u1 == 0) && (u2 == 0)) {
    /* stack */
    z = pf_params->expstack[type][type2];
  } else if (no_close == 0) {
    if ((u1 == 0) || (u2 == 0)) {
      /* bulge */
      int u;
      u = (u1 == 0) ? u2 : u1;
      z = pf_params->expbulge[u];
      if (u2 + u1 == 1) {
        z *= pf_params->expstack[type][type2];
      } else {
        if (type > 2)
          z *= pf_params->expTermAU;

        if (type2 > 2)
          z *= pf_params->expTermAU;
      }
    } else {
      /* internal loop */
      if (u1 + u2 == 2) {
        /* size 2 is special */
        z = pf_params->expint11[type][type2][si1][sj1];
      } else if ((u1 == 1) && (u2 == 2)) {
        z = pf_params->expint21[type][type2][si1][sq1][sj1];
      } else if ((u1 == 2) && (u2 == 1)) {
        z = pf_params->expint21[type2][type][sq1][si1][sp1];
      } else if ((u1 == 2) && (u2 == 2)) {
        z = pf_params->expint22[type][type2][si1][sp1][sq1][sj1];
      } else if (((u1 == 2) && (u2 == 3)) || ((u1 == 3) && (u2 == 2))) {
        /*2-3 is special*/
        z = pf_params->expinternal[5] *
            pf_params->expmismatch23I[type][si1][sj1] *
            pf_params->expmismatch23I[type2][sq1][sp1];
        z *= pf_params->expninio[2][1];
      } else if ((u1 == 1) || (u2 == 1)) {
        /*1-n is special*/
        z = pf_params->expinternal[u1 + u2] *
            pf_params->expmismatch1nI[type][si1][sj1] *
            pf_params->expmismatch1nI[type2][sq1][sp1];
        z *= pf_params->expninio[2][abs(u1 - u2)];
      } else {
        z = pf_params->expinternal[u1 + u2] *
            pf_params->expmismatchI[type][si1][sj1] *
            pf_params->expmismatchI[type2][sq1][sp1];
        z *= pf_params->expninio[2][abs(u1 - u2)];
      }
    }
  }

  return z;
}


PUBLIC void
init_pf_circ_fold(int length VRNA_UNUSED)
{
  ;/* DO NOTHING */
}


PUBLIC void
init_pf_fold(int length VRNA_UNUSED)
{
  ;/* DO NOTHING */
}


/**
*** Allocate memory for all matrices and other stuff
**/
PUBLIC void
free_pf_arrays(void)
{
  if (backward_compat_compound && backward_compat) {
    vrna_fold_compound_free(backward_compat_compound);
    backward_compat_compound  = NULL;
    backward_compat           = 0;
    iindx                     = NULL;
  }
}


PUBLIC FLT_OR_DBL *
export_bppm(void)
{
  if (backward_compat_compound)
    if (backward_compat_compound->exp_matrices)
      if (backward_compat_compound->exp_matrices->probs)
        return backward_compat_compound->exp_matrices->probs;

  return NULL;
}


/*
 * -------------------------------------------------------------------------
 * make arrays used for pf_fold available to other routines
 */
PUBLIC int
get_pf_arrays(short       **S_p,
              short       **S1_p,
              char        **ptype_p,
              FLT_OR_DBL  **qb_p,
              FLT_OR_DBL  **qm_p,
              FLT_OR_DBL  **q1k_p,
              FLT_OR_DBL  **qln_p)
{
  if (backward_compat_compound) {
    if (backward_compat_compound->exp_matrices) {
      if (backward_compat_compound->exp_matrices->qb) {
        *S_p      = backward_compat_compound->sequence_encoding2;
        *S1_p     = backward_compat_compound->sequence_encoding;
        *ptype_p  = backward_compat_compound->ptype_pf_compat;
        *qb_p     = backward_compat_compound->exp_matrices->qb;
        *qm_p     = backward_compat_compound->exp_matrices->qm;
        *q1k_p    = backward_compat_compound->exp_matrices->q1k;
        *qln_p    = backward_compat_compound->exp_matrices->qln;
        return 1;
      }
    }
  }

  return 0;
}


/*-----------------------------------------------------------------*/
PUBLIC float
pf_fold(const char  *sequence,
        char        *structure)
{
  return wrap_pf_fold(sequence, structure, NULL, do_backtrack, fold_constrained, 0);
}


PUBLIC float
pf_circ_fold(const char *sequence,
             char       *structure)
{
  return wrap_pf_fold(sequence, structure, NULL, do_backtrack, fold_constrained, 1);
}


PUBLIC float
pf_fold_par(const char        *sequence,
            char              *structure,
            vrna_exp_param_t  *parameters,
            int               calculate_bppm,
            int               is_constrained,
            int               is_circular)
{
  return wrap_pf_fold(sequence, structure, parameters, calculate_bppm, is_constrained, is_circular);
}


PUBLIC char *
pbacktrack(char *seq)
{
  int n = (int)strlen(seq);

  return vrna_pbacktrack5(backward_compat_compound, n);
}


PUBLIC char *
pbacktrack5(char  *seq VRNA_UNUSED,
            int   length)
{
  /* the seq parameter must no differ to the one stored globally anyway, so we just ignore it */
  return vrna_pbacktrack5(backward_compat_compound, length);
}


PUBLIC char *
pbacktrack_circ(char *seq VRNA_UNUSED)
{
  char      *structure;
  vrna_md_t *md;

  structure = NULL;

  if (backward_compat_compound) {
    md = &(backward_compat_compound->exp_params->model_details);
    if (md->circ && backward_compat_compound->exp_matrices->qm2)
      structure = vrna_pbacktrack(backward_compat_compound);
  }

  return structure;
}


PUBLIC void
update_pf_params(int length VRNA_UNUSED)
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
update_pf_params_par(int              length VRNA_UNUSED,
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


PUBLIC char *
get_centroid_struct_gquad_pr(int    length VRNA_UNUSED,
                             double *dist)
{
  return vrna_centroid(backward_compat_compound, dist);
}


PUBLIC void
assign_plist_gquad_from_pr(vrna_ep_t  **pl,
                           int        length VRNA_UNUSED, /* ignored */
                           double     cut_off)
{
  if (!backward_compat_compound)
    *pl = NULL;
  else if (!backward_compat_compound->exp_matrices->probs)
    *pl = NULL;
  else
    *pl = vrna_plist_from_probs(backward_compat_compound, cut_off);
}


PUBLIC double
mean_bp_distance(int length VRNA_UNUSED)
{
  if (backward_compat_compound)
    if (backward_compat_compound->exp_matrices)
      if (backward_compat_compound->exp_matrices->probs)
        return vrna_mean_bp_distance(backward_compat_compound);

  vrna_log_warning("mean_bp_distance: "
                       "you need to call vrna_pf_fold first");

  return 0.; /* we will never get to this point */
}


PUBLIC double
mean_bp_distance_pr(int         length,
                    FLT_OR_DBL  *p)
{
  double  d       = 0;
  int     *index  = vrna_idx_row_wise((unsigned int)length);

  if (p == NULL) {
    vrna_log_warning("mean_bp_distance_pr: "
                         "p == NULL. You need to supply a valid probability matrix for mean_bp_distance_pr()");
    return d;
  }

  d = wrap_mean_bp_distance(p, length, index, TURN);

  free(index);
  return d;
}


#endif
