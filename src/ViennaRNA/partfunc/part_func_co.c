/*
 *                partiton function for RNA secondary structures
 *
 *                Ivo L Hofacker
 *                Stephan Bernhart
 *                Ronny Lorenz
 *                ViennaRNA package
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
#include "ViennaRNA/structures/problist.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/plotting/probabilities.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/partfunc/global.h"
#include "ViennaRNA/part_func_co.h"

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

int     mirnatog      = 0;
double  F_monomer[2]  = {
  0, 0
};                             /* free energies of the two monomers */

#endif

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

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

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PRIVATE vrna_dimer_pf_t
wrap_co_pf_fold(char              *sequence,
                char              *structure,
                vrna_exp_param_t  *parameters,
                int               calculate_bppm,
                int               is_constrained);


#endif

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_dimer_pf_t
vrna_pf_co_fold(const char  *seq,
                char        *structure,
                vrna_ep_t   **pl)
{
  double                mfe;
  vrna_dimer_pf_t       X;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  /* no need to backtrack MFE structure */
  md.backtrack = 0;
  if (pl)
    md.compute_bpp = 1;
  else
    md.compute_bpp = 0;

  vc = vrna_fold_compound(seq, &md, VRNA_OPTION_DEFAULT);

  mfe = (double)vrna_mfe_dimer(vc, NULL);
  vrna_exp_params_rescale(vc, &mfe);

  X = vrna_pf_dimer(vc, structure);

  if (pl)
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);

  vrna_fold_compound_free(vc);

  return X;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 *****************************************
 * BEGIN backward compatibility wrappers *
 *****************************************
 */
PRIVATE vrna_dimer_pf_t
wrap_co_pf_fold(char              *sequence,
                char              *structure,
                vrna_exp_param_t  *parameters,
                int               calculate_bppm,
                int               is_constrained)
{
  int                   length;
  char                  *seq;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vc      = NULL;
  length  = strlen(sequence);

  seq = (char *)vrna_alloc(sizeof(char) * (length + 2));
  if (cut_point > -1) {
    int i;
    for (i = 0; i < cut_point - 1; i++)
      seq[i] = sequence[i];
    seq[i] = '&';
    for (; i < (int)length; i++)
      seq[i + 1] = sequence[i];
  } else {
    /* this ensures the allocation of all cofold matrices via vrna_fold_compound_t */
    free(seq);
    seq = strdup(sequence);
  }

  /*
   *  if present, extract model details from provided parameters variable,
   *  to properly initialize the fold compound. Otherwise use default
   *  settings taken from deprecated global variables
   */
  if (parameters)
    vrna_md_copy(&md, &(parameters->model_details));
  else
    set_model_details(&md);

  /* set min_loop_size and backtracing options */
  md.compute_bpp    = calculate_bppm;
  md.min_loop_size  = 0;

  vc = vrna_fold_compound(seq, &md, VRNA_OPTION_DEFAULT);

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

  if (backward_compat_compound)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;
  iindx                     = backward_compat_compound->iindx;

  free(seq);
  return vrna_pf_dimer(vc, structure);
}


/*
 *****************************************
 * END backward compatibility wrappers   *
 *****************************************
 */


/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

PUBLIC vrna_dimer_pf_t
co_pf_fold(char *sequence,
           char *structure)
{
  return wrap_co_pf_fold(sequence, structure, NULL, do_backtrack, fold_constrained);
}


PUBLIC vrna_dimer_pf_t
co_pf_fold_par(char             *sequence,
               char             *structure,
               vrna_exp_param_t *parameters,
               int              calculate_bppm,
               int              is_constrained)
{
  return wrap_co_pf_fold(sequence, structure, parameters, calculate_bppm, is_constrained);
}


PUBLIC vrna_ep_t *
get_plist(vrna_ep_t *pl,
          int       length,
          double    cut_off)
{
  int i, j, n, count, *my_iindx;

  my_iindx = backward_compat_compound->iindx;
  /* get pair probibilities out of pr array */
  count = 0;
  n     = 2;
  for (i = 1; i < length; i++) {
    for (j = i + 1; j <= length; j++) {
      if (pr[my_iindx[i] - j] < cut_off)
        continue;

      if (count == n * length - 1) {
        n   *= 2;
        pl  = (vrna_ep_t *)vrna_realloc(pl, n * length * sizeof(vrna_ep_t));
      }

      pl[count].i   = i;
      pl[count].j   = j;
      pl[count++].p = pr[my_iindx[i] - j];
      /* printf("gpl: %2d %2d %.9f\n",i,j,pr[my_iindx[i]-j]); */
    }
  }
  pl[count].i   = 0;
  pl[count].j   = 0; /* ->?? */
  pl[count++].p = 0.;
  pl            = (vrna_ep_t *)vrna_realloc(pl, (count) * sizeof(vrna_ep_t));
  return pl;
}


PUBLIC void
compute_probabilities(double    FAB,
                      double    FA,
                      double    FB,
                      vrna_ep_t *prAB,
                      vrna_ep_t *prA,
                      vrna_ep_t *prB,
                      int       Alength)
{
  if (backward_compat_compound && backward_compat) {
    vrna_pf_dimer_probs(FAB,
                        FA,
                        FB,
                        prAB,
                        (const vrna_ep_t *)prA,
                        (const vrna_ep_t *)prB,
                        Alength,
                        (const vrna_exp_param_t *)backward_compat_compound->exp_params);
  }
}


PUBLIC vrna_dimer_conc_t *
get_concentrations(double FcAB,
                   double FcAA,
                   double FcBB,
                   double FEA,
                   double FEB,
                   double *startconc)
{
  return vrna_pf_dimer_concentrations(FcAB,
                                      FcAA,
                                      FcBB,
                                      FEA,
                                      FEB,
                                      (const double *)startconc,
                                      (const vrna_exp_param_t *)backward_compat_compound->exp_params);
}


PUBLIC void
init_co_pf_fold(int length VRNA_UNUSED)
{
  ;/* DO NOTHING */
}


PUBLIC void
free_co_pf_arrays(void)
{
  if (backward_compat_compound && backward_compat) {
    vrna_fold_compound_free(backward_compat_compound);
    backward_compat_compound  = NULL;
    backward_compat           = 0;
  }
}


PUBLIC FLT_OR_DBL *
export_co_bppm(void)
{
  if (backward_compat_compound)
    return backward_compat_compound->exp_matrices->probs;
  else
    return NULL;
}


PUBLIC void
update_co_pf_params(int length VRNA_UNUSED)
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
update_co_pf_params_par(int               length VRNA_UNUSED,
                        vrna_exp_param_t  *parameters)
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


#endif
