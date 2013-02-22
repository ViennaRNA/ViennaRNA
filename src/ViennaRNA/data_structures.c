/** \file data_structures.c **/

/*
                  Data structure creation/destruction

                  This file contains everything which is necessary to
                  obtain and destroy datastructures used in the folding
                  recurrences throughout the VienneRNA paclage

                  c Ronny Lorenx

                  Vienna RNA package
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"

#ifdef _OPENMP
#include <omp.h>
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
PRIVATE char *tokenize(char *line, int *cut_point);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC mfe_matrices  *
get_mfe_matrices_alloc( unsigned int n,
                        unsigned int alloc_vector){

  if(n >= (unsigned int)sqrt((double)INT_MAX))
    nrerror("get_mfe_matrices_alloc@data_structures.c: sequence length exceeds addressable range");

  mfe_matrices  *vars   = (mfe_matrices *)space(sizeof(mfe_matrices));

  vars->allocated       = 0;
  vars->f5              = NULL;
  vars->f3              = NULL;
  vars->fc              = NULL;
  vars->c               = NULL;
  vars->fML             = NULL;
  vars->fM1             = NULL;
  vars->fM2             = NULL;
  vars->FcH             = INF;
  vars->FcI             = INF;
  vars->FcM             = INF;
  vars->Fc              = INF;

  if(alloc_vector){
    vars->allocated = alloc_vector;
    unsigned int size     = ((n + 1) * (n + 2)) >> 1;
    unsigned int lin_size = n + 2;

    if(alloc_vector & ALLOC_F5)
      vars->f5  = (int *) space(sizeof(int) * lin_size);

    if(alloc_vector & ALLOC_F3)
      vars->f3  = (int *) space(sizeof(int) * lin_size);

    if(alloc_vector & ALLOC_FC)
      vars->fc  = (int *) space(sizeof(int) * lin_size);

    if(alloc_vector & ALLOC_C)
      vars->c      = (int *) space(sizeof(int) * size);

    if(alloc_vector & ALLOC_FML)
      vars->fML    = (int *) space(sizeof(int) * size);

    if(alloc_vector & ALLOC_FM1)
      vars->fM1    = (int *) space(sizeof(int) * size);

    if(alloc_vector & ALLOC_FM2)
      vars->fM2    = (int *) space(sizeof(int) * lin_size);

  }

  return vars;
}

PUBLIC void
destroy_mfe_matrices(mfe_matrices *self){

  if(self){
    if(self->allocated){
      if(self->allocated & ALLOC_F5)
        free(self->f5);
      if(self->allocated & ALLOC_F3)
        free(self->f3);
      if(self->allocated & ALLOC_FC)
        free(self->fc);
      if(self->allocated & ALLOC_C)
        free(self->c);
      if(self->allocated & ALLOC_FML)
        free(self->fML);
      if(self->allocated & ALLOC_FM1)
        free(self->fM1);
      if(self->allocated & ALLOC_FM2)
        free(self->fM2);
    }
    free(self);
  }
}

PUBLIC vrna_fold_compound *
get_fold_compound_mfe(const char *sequence, paramT *P){

  return get_fold_compound_mfe_constrained(sequence, NULL, NULL, P);
}

PUBLIC vrna_fold_compound *
get_fold_compound_mfe_constrained(const char *sequence,
                                  hard_constraintT *hc,
                                  soft_constraintT *sc,
                                  paramT *P){

  vrna_fold_compound  *vc;
  paramT              *params;
  unsigned int        alloc_vector, length;
  int                 cut_point;
  char                *seq;

  cut_point = -1;

  /* sanity check */
  length = (sequence) ? strlen(sequence) : 0;
  if(length == 0)
    nrerror("get_fold_compound_mfe_constraint@data_structures.c: sequence length must be greater 0");

  seq = tokenize(sequence, &cut_point); /* splice out the '&' if concatenated sequences and reset cut_point */


    /* prepare the parameters datastructure */
  if(P){
    params = get_parameter_copy(P);
  } else { /* this fallback relies on global parameters and thus is not threadsafe */
    model_detailsT md;
    set_model_details(&md);
    params = get_scaled_parameters(temperature, md);
  }

  /* prepare the allocation vector for the folding matrices */
  if(cut_point > -1){
    alloc_vector = ALLOC_MFE_HYBRID;
    if(params->model_details.uniq_ML)
      alloc_vector |= ALLOC_MFE_HYBRID_UNIQ_ML;
  } else {
    alloc_vector = ALLOC_MFE_DEFAULT;
    if(params->model_details.circ)
      alloc_vector |= ALLOC_MFE_CIRC;
    if(params->model_details.uniq_ML)
      alloc_vector |= ALLOC_MFE_UNIQ_ML;
  }

  /* start making the fold compound */
  vc                      = (vrna_fold_compound *)space(sizeof(vrna_fold_compound));

  vc->params              = params;
  vc->cut_point           = cut_point;
  vc->sequence            = strdup(sequence);
  vc->length              = strlen(sequence);
  vc->sequence_encoding   = get_sequence_encoding(sequence, 1, &(params->model_details));
  vc->sequence_encoding2  = get_sequence_encoding(sequence, 0, &(params->model_details));

  vc->ptype               = get_ptypes(vc->sequence_encoding2, &(P->model_details), 0);
  vc->exp_params          = NULL;

  vc->matrices            = get_mfe_matrices_alloc(vc->length, alloc_vector);

  vc->hc                  = hc ? hc : get_hard_constraints(sequence, NULL, &(vc->params->model_details), TURN, (unsigned int)0);
  vc->sc                  = sc ? sc : NULL;

  vc->iindx               = NULL;
  vc->jindx               = get_indx(vc->length);

  return vc;
}

PUBLIC  void
destroy_fold_compound(vrna_fold_compound *vc){

  if(vc){
    if(vc->matrices)
      destroy_mfe_matrices(vc->matrices);
    if(vc->sequence)
      free(vc->sequence);
    if(vc->sequence_encoding)
      free(vc->sequence_encoding);
    if(vc->sequence_encoding2)
      free(vc->sequence_encoding2);
    if(vc->ptype)
      free(vc->ptype);
    if(vc->iindx)
      free(vc->iindx);
    if(vc->jindx)
      free(vc->jindx);
    if(vc->params)
      free(vc->params);
    if(vc->exp_params)
      free(vc->exp_params);
    if(vc->hc)
      destroy_hard_constraints(vc->hc);
    if(vc->sc)
      destroy_soft_constraints(vc->sc);

    free(vc);
  }
}


PRIVATE char *tokenize(char *line, int *cut_point){
  char *pos, *copy;
  int cut = -1;
  copy = NULL;
  if(line){
    copy = (char *) space(strlen(line)+1);
    (void) sscanf(line, "%s", copy);
    pos = strchr(copy, '&');
    if (pos) {
      cut = (int) (pos-copy)+1;
      if (cut >= strlen(copy)) cut = -1;
      if (strchr(pos+1, '&')) nrerror("more than one cut-point in input");
      for (;*pos;pos++) *pos = *(pos+1); /* splice out the & */
    }
    if (cut > -1) {
      if (*cut_point==-1) *cut_point = cut;
      else if (*cut_point != cut) {
        fprintf(stderr,"cut_point = %d cut = %d\n", *cut_point, cut);
        nrerror("Sequence and Structure have different cut points.");
      }
    }
  }
  return copy;
}
