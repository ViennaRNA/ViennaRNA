/*
 *                minimum free energy
 *                RNA secondary structure prediction
 *
 *                c Ivo Hofacker, Chrisoph Flamm
 *                original implementation by
 *                Walter Fontana
 *
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/structures/pairtable.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/backtrack/global.h"
#include "ViennaRNA/subopt/zuker.h"
#include "ViennaRNA/cofold.h"

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#ifdef _OPENMP
#include <omp.h>
#endif

#endif

#define MAXSECTORS        500     /* dimension for a backtrack array */

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

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/* some backward compatibility stuff */
PRIVATE int                   backward_compat           = 0;
PRIVATE vrna_fold_compound_t  *backward_compat_compound = NULL;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/* wrappers for old API compatibility */
PRIVATE void
wrap_array_export(int   **f5_p,
                  int   **c_p,
                  int   **fML_p,
                  int   **fM1_p,
                  int   **fc_p,
                  int   **indx_p,
                  char  **ptype_p);


PRIVATE float
wrap_cofold(const char    *string,
            char          *structure,
            vrna_param_t  *parameters,
            int           is_constrained);


PRIVATE SOLUTION *
wrap_zukersubopt(const char   *string,
                 vrna_param_t *parameters);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */


/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */
PRIVATE void
wrap_array_export(int   **f5_p,
                  int   **c_p,
                  int   **fML_p,
                  int   **fM1_p,
                  int   **fc_p,
                  int   **indx_p,
                  char  **ptype_p)
{
  /* make the DP arrays available to routines such as subopt() */
  if (backward_compat_compound) {
    *f5_p     = backward_compat_compound->matrices->f5;
    *c_p      = backward_compat_compound->matrices->c;
    *fML_p    = backward_compat_compound->matrices->fML;
    *fM1_p    = backward_compat_compound->matrices->fM1;
    *fc_p     = NULL; /* this breaks compatibility, but there is no way around it */
    *indx_p   = backward_compat_compound->jindx;
    *ptype_p  = backward_compat_compound->ptype;
  }
}


/*--------------------------------------------------------------------------*/

PRIVATE float
wrap_cofold(const char    *string,
            char          *structure,
            vrna_param_t  *parameters,
            int           is_constrained)
{
  unsigned int          length;
  char                  *seq;
  vrna_fold_compound_t  *vc;
  vrna_param_t          *P;
  float                 mfe;

  vc      = NULL;
  length  = strlen(string);

#ifdef _OPENMP
  /* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  /* we need the parameter structure for hard constraints */
  if (parameters) {
    P = vrna_params_copy(parameters);
  } else {
    vrna_md_t md;
    set_model_details(&md);
    md.temperature  = temperature;
    P               = vrna_params(&md);
  }

  P->model_details.min_loop_size = 0;  /* set min loop length to 0 */

  /* dirty hack to reinsert the '&' according to the global variable 'cut_point' */
  seq = vrna_cut_point_insert(string, cut_point);

  /* get compound structure */
  vc = vrna_fold_compound(seq, &(P->model_details), 0);

  if (parameters) {
    /* replace params if necessary */
    free(vc->params);
    vc->params = P;
  } else {
    free(P);
  }

  /* handle hard constraints in pseudo dot-bracket format if passed via simple interface */
  if (is_constrained && structure) {
    unsigned int constraint_options = 0;
    constraint_options |= VRNA_CONSTRAINT_DB
                          | VRNA_CONSTRAINT_DB_PIPE
                          | VRNA_CONSTRAINT_DB_DOT
                          | VRNA_CONSTRAINT_DB_X
                          | VRNA_CONSTRAINT_DB_ANG_BRACK
                          | VRNA_CONSTRAINT_DB_RND_BRACK
                          | VRNA_CONSTRAINT_DB_INTRAMOL
                          | VRNA_CONSTRAINT_DB_INTERMOL;

    vrna_constraints_add(vc, (const char *)structure, constraint_options);
  }

  if (backward_compat_compound)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;

  /* cleanup */
  free(seq);

  /* compute mfe without backtracing */
  mfe = vrna_mfe(vc, NULL);

  /* now we backtrace in a backward compatible way */
  if (structure && vc->params->model_details.backtrack) {
    char            *s;
    sect            bt_stack[MAXSECTORS];
    vrna_bp_stack_t *bp;

    bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + length / 2))); /* add a guess of how many G's may be involved in a G quadruplex */

    vrna_backtrack_from_intervals(vc, bp, bt_stack, 0);

    s = vrna_db_from_bp_stack(bp, length);
    strncpy(structure, s, length + 1);
    free(s);

    if (base_pair)
      free(base_pair);

    base_pair = bp;
  }

  return mfe;
}


PRIVATE SOLUTION *
wrap_zukersubopt(const char   *string,
                 vrna_param_t *parameters)
{
  vrna_fold_compound_t  *vc;
  vrna_param_t          *P;

  vc = NULL;

#ifdef _OPENMP
  /* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  /* we need the parameter structure for hard constraints */
  if (parameters) {
    P = vrna_params_copy(parameters);
  } else {
    vrna_md_t md;
    set_model_details(&md);
    md.temperature  = temperature;
    P               = vrna_params(&md);
  }

  /* get compound structure */
  vc = vrna_fold_compound(string, &(P->model_details), VRNA_OPTION_DEFAULT);

  if (parameters) {
    /* replace params if necessary */
    free(vc->params);
    vc->params = P;
  } else {
    free(P);
  }

  if (backward_compat_compound)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;

  return vrna_subopt_zuker(vc);
}


PUBLIC void
initialize_cofold(int length VRNA_UNUSED)
{
  ;/* DO NOTHING */
}


PUBLIC void
free_co_arrays(void)
{
  if (backward_compat_compound && backward_compat) {
    vrna_fold_compound_free(backward_compat_compound);
    backward_compat_compound  = NULL;
    backward_compat           = 0;
  }
}


/*--------------------------------------------------------------------------*/

PUBLIC void
export_cofold_arrays_gq(int   **f5_p,
                        int   **c_p,
                        int   **fML_p,
                        int   **fM1_p,
                        int   **fc_p,
                        int   **ggg_p,
                        int   **indx_p,
                        char  **ptype_p)
{
  /* make the DP arrays available to routines such as subopt() */
  wrap_array_export(f5_p, c_p, fML_p, fM1_p, fc_p, indx_p, ptype_p);
  if (backward_compat_compound) {
    *ggg_p = NULL; /* This will break backward compatibility */
  }
}


PUBLIC void
export_cofold_arrays(int  **f5_p,
                     int  **c_p,
                     int  **fML_p,
                     int  **fM1_p,
                     int  **fc_p,
                     int  **indx_p,
                     char **ptype_p)
{
  wrap_array_export(f5_p, c_p, fML_p, fM1_p, fc_p, indx_p, ptype_p);
}


PUBLIC float
cofold(const char *string,
       char       *structure)
{
  return wrap_cofold(string, structure, NULL, fold_constrained);
}


PUBLIC float
cofold_par(const char   *string,
           char         *structure,
           vrna_param_t *parameters,
           int          is_constrained)
{
  return wrap_cofold(string, structure, parameters, is_constrained);
}


PUBLIC SOLUTION *
zukersubopt(const char *string)
{
  return wrap_zukersubopt(string, NULL);
}


PUBLIC SOLUTION *
zukersubopt_par(const char    *string,
                vrna_param_t  *parameters)
{
  return wrap_zukersubopt(string, parameters);
}


PUBLIC void
update_cofold_params(void)
{
  vrna_fold_compound_t *v;

  if (backward_compat_compound && backward_compat) {
    vrna_md_t md;
    v = backward_compat_compound;

    if (v->params)
      free(v->params);

    set_model_details(&md);
    v->params = vrna_params(&md);
  }
}


PUBLIC void
update_cofold_params_par(vrna_param_t *parameters)
{
  vrna_fold_compound_t *v;

  if (backward_compat_compound && backward_compat) {
    v = backward_compat_compound;

    if (v->params)
      free(v->params);

    if (parameters) {
      v->params = vrna_params_copy(parameters);
    } else {
      vrna_md_t md;
      set_model_details(&md);
      md.temperature  = temperature;
      v->params       = vrna_params(&md);
    }
  }
}


#endif
