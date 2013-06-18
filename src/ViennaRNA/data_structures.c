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
#include "ViennaRNA/gquad.h"
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

PUBLIC mfe_matricesT  *
get_mfe_matrices_alloc( unsigned int n,
                        unsigned int alloc_vector){

  if(n >= (unsigned int)sqrt((double)INT_MAX))
    nrerror("get_mfe_matrices_alloc@data_structures.c: sequence length exceeds addressable range");

  mfe_matricesT *vars   = (mfe_matricesT *)space(sizeof(mfe_matricesT));

  vars->allocated       = 0;
  vars->length          = 0;
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
  vars->ggg             = NULL;

  if(alloc_vector){
    vars->allocated = alloc_vector;
    vars->length    = n;
    unsigned int size     = ((n + 1) * (n + 2)) / 2;
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
destroy_mfe_matrices(mfe_matricesT *self){

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
      if(self->ggg);
        free(self->ggg);
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
  char                *seq, *seq2;

  cut_point = -1;

  /* sanity check */
  length = (sequence) ? strlen(sequence) : 0;
  if(length == 0)
    nrerror("get_fold_compound_mfe_constraint@data_structures.c: sequence length must be greater 0");

  seq2 = strdup(sequence);
  seq = tokenize(seq2, &cut_point); /* splice out the '&' if concatenated sequences and reset cut_point */

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
  vc->cutpoint            = cut_point;
  vc->sequence            = seq;
  vc->length              = strlen(seq);
  vc->sequence_encoding   = get_sequence_encoding(seq, 1, &(params->model_details));
  vc->sequence_encoding2  = get_sequence_encoding(seq, 0, &(params->model_details));

  vc->ptype               = get_ptypes(vc->sequence_encoding2, &(params->model_details), 0);
  vc->exp_params          = NULL;

  vc->matrices            = get_mfe_matrices_alloc(vc->length, alloc_vector);

  /* get gquadruplex matrix if needed */
  if(vc->params->model_details.gquad)
    vc->matrices->ggg = get_gquad_matrix(vc->sequence_encoding2, vc->params);

  vc->hc                  = hc ? hc : get_hard_constraints(seq, NULL, &(vc->params->model_details), TURN, (unsigned int)0);
  vc->sc                  = sc;

  vc->iindx               = NULL;
  vc->jindx               = get_indx(vc->length);

  free(seq2);
  return vc;
}

PUBLIC  void
destroy_fold_compound(vrna_fold_compound *vc){

  if(vc){
    if(vc->matrices)
      destroy_mfe_matrices(vc->matrices);
    if(vc->exp_matrices)
      destroy_pf_matrices(vc->exp_matrices);
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

PUBLIC pf_matricesT  *
get_pf_matrices_alloc(unsigned int n,
                      unsigned int alloc_vector){

  if(n >= (unsigned int)sqrt((double)INT_MAX))
    nrerror("get_pf_matrices_alloc@data_structures.c: sequence length exceeds addressable range");

  pf_matricesT  *vars   = (pf_matricesT *)space(sizeof(pf_matricesT));

  vars->allocated       = 0;
  vars->length          = 0;
  vars->q               = NULL;
  vars->qb              = NULL;
  vars->qm              = NULL;
  vars->qm1             = NULL;
  vars->qm2             = NULL;
  vars->qho             = 0.;
  vars->qio             = 0.;
  vars->qmo             = 0.;
  vars->qo              = 0.;
  vars->G               = NULL;
  vars->q1k             = NULL;
  vars->qln             = NULL;

  vars->scale           = NULL;
  vars->expMLbase       = NULL;

  if(alloc_vector){
    vars->allocated = alloc_vector;
    vars->length    = n;
    unsigned int size     = ((n + 1) * (n + 2)) / 2;
    unsigned int lin_size = n + 2;

    if(alloc_vector & ALLOC_F)
      vars->q     = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * size);

    if(alloc_vector & ALLOC_C)
      vars->qb    = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * size);

    if(alloc_vector & ALLOC_FML)
      vars->qm    = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * size);

    if(alloc_vector & ALLOC_FM1)
      vars->qm1   = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * size);

    if(alloc_vector & ALLOC_FM2)
      vars->qm2   = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * lin_size);

    if(alloc_vector & ALLOC_PROBS)
      vars->probs = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * size);

    if(alloc_vector & ALLOC_AUX){
      vars->q1k   = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * lin_size);
      vars->qln   = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * lin_size);
    }

    /*  always alloc the helper arrays for unpaired nucleotides in multi-
        branch loops and scaling
    */
    vars->scale     = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * lin_size);
    vars->expMLbase = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * lin_size);
  }

  return vars;
}

PUBLIC void
destroy_pf_matrices(pf_matricesT *self){

  if(self){
    if(self->allocated){
      if(self->allocated & ALLOC_F)
        free(self->q);
      if(self->allocated & ALLOC_C)
        free(self->qb);
      if(self->allocated & ALLOC_FML)
        free(self->qm);
      if(self->allocated & ALLOC_FM1)
        free(self->qm1);
      if(self->allocated & ALLOC_FM2)
        free(self->qm2);
      if(self->allocated & ALLOC_PROBS)
        free(self->probs);
      if(self->G)
        free(self->G);
      if(self->allocated & ALLOC_AUX){
        free(self->q1k);
        free(self->qln);
      }
      if(self->scale)
        free(self->scale);
      if(self->expMLbase)
        free(self->expMLbase);
    }
    free(self);
  }
}

PUBLIC vrna_fold_compound *
get_fold_compound_pf_constrained( const char *sequence,
                                  hard_constraintT *hc,
                                  soft_constraintT *sc,
                                  pf_paramT *P){

  vrna_fold_compound  *vc;
  pf_paramT           *params;
  unsigned int        alloc_vector, length;
  int                 cut_point, i;
  char                *seq, *seq2;
  double              scaling_factor;

  cut_point = -1;

  /* sanity check */
  length = (sequence) ? strlen(sequence) : 0;
  if(length == 0)
    nrerror("get_fold_compound_mfe_constraint@data_structures.c: sequence length must be greater 0");

  seq2 = strdup(sequence);
  seq = tokenize(seq2, &cut_point); /* splice out the '&' if concatenated sequences and reset cut_point */

    /* prepare the parameters datastructure */
  if(P){
    params = get_boltzmann_factor_copy(P);
  } else { /* this fallback relies on global parameters and thus is not threadsafe */
    model_detailsT md;
    set_model_details(&md);
    params = get_boltzmann_factors(temperature, 1.0, md, pf_scale);
  }

  /* prepare the allocation vector for the folding matrices */
  if(cut_point > -1){
    alloc_vector = ALLOC_PF_HYBRID;
    if(params->model_details.uniq_ML)
      alloc_vector |= ALLOC_PF_HYBRID_UNIQ_ML;
  } else {
    alloc_vector = ALLOC_PF_DEFAULT;
    if(params->model_details.circ)
      alloc_vector |= ALLOC_PF_CIRC;
    if(params->model_details.uniq_ML)
      alloc_vector |= ALLOC_PF_UNIQ_ML;
  }

  /* start making the fold compound */
  vc                      = (vrna_fold_compound *)space(sizeof(vrna_fold_compound));

  vc->exp_params          = params;
  vc->cutpoint            = cut_point;
  vc->sequence            = seq;
  vc->length              = strlen(seq);
  vc->sequence_encoding   = get_sequence_encoding(seq, 1, &(params->model_details));
  vc->sequence_encoding2  = get_sequence_encoding(seq, 0, &(params->model_details));

  vc->ptype               = get_ptypes(vc->sequence_encoding2, &(params->model_details), 1);
  vc->params              = NULL;

  vc->exp_matrices        = get_pf_matrices_alloc(vc->length, alloc_vector);

  /* get gquadruplex matrix if needed */
  if(params->model_details.gquad)
    vc->exp_matrices->G   = get_gquad_pf_matrix(vc->sequence_encoding2, vc->exp_matrices->scale, params);

  /* fill additional helper arrays for scaling etc. */
  scaling_factor = params->pf_scale;
  if (scaling_factor == -1) { /* mean energy for random sequences: 184.3*length cal */
    scaling_factor = exp(-(-185+(params->temperature-37.)*7.27)/params->kT);
    if (scaling_factor<1) scaling_factor=1;
    params->pf_scale = scaling_factor;
    pf_scale = params->pf_scale; /* compatibility with RNAup, may be removed sometime */
  }
  vc->exp_matrices->scale[0] = 1.;
  vc->exp_matrices->scale[1] = 1./scaling_factor;
  vc->exp_matrices->expMLbase[0] = 1;
  vc->exp_matrices->expMLbase[1] = params->expMLbase/scaling_factor;
  for (i=2; i<=length; i++) {
    vc->exp_matrices->scale[i] = vc->exp_matrices->scale[i/2]*vc->exp_matrices->scale[i-(i/2)];
    vc->exp_matrices->expMLbase[i] = pow(params->expMLbase, (double)i) * vc->exp_matrices->scale[i];
  }

  /* get gquadruplex matrix if needed */
  if(params->model_details.gquad)
    vc->exp_matrices->G   = get_gquad_pf_matrix(vc->sequence_encoding2, vc->exp_matrices->scale, params);



  vc->hc                  = hc ? hc : get_hard_constraints(seq, NULL, &(params->model_details), TURN, VRNA_CONSTRAINT_IINDX);
  vc->sc                  = sc;

  vc->iindx               = get_iindx(vc->length);
  iindx                   = get_iindx(vc->length); /* for backward compatibility and Perl wrapper (they need the global iindx) */
  vc->jindx               = get_indx(vc->length);

  free(seq2);

  return vc;
}

PUBLIC vrna_alifold_compound *
get_alifold_compound_mfe(const char **sequences, paramT *P){

  return get_alifold_compound_mfe_constrained(sequences, NULL, NULL, P);
}

PUBLIC vrna_alifold_compound *
get_alifold_compound_mfe_constrained( const char **sequences,
                                      hard_constraintT **hc,
                                      soft_constraintT **sc,
                                      paramT *P){

  vrna_alifold_compound *vc;
  paramT                *params;
  unsigned int          alloc_vector, length;

  /* sanity check */
  length = (sequences) ? strlen(sequences[0]) : 0;
  if(length == 0)
    nrerror("get_alifold_compound_mfe_constraint@data_structures.c: sequence length must be greater 0");

    /* prepare the parameters datastructure */
  if(P){
    params = get_parameter_copy(P);
  } else { /* this fallback relies on global parameters and thus is not threadsafe */
    model_detailsT md;
    set_model_details(&md);
    params = get_scaled_parameters(temperature, md);
  }

  /* prepare the allocation vector for the folding matrices */
  alloc_vector = ALLOC_MFE_DEFAULT;
  if(params->model_details.circ)
    alloc_vector |= ALLOC_MFE_CIRC;
  if(params->model_details.uniq_ML)
    alloc_vector |= ALLOC_MFE_UNIQ_ML;

  /* start making the alfold compound */
  vc                      = (vrna_alifold_compound *)space(sizeof(vrna_alifold_compound));

  vc->params              = params;
  vc->sequences           = sequences;
  vc->length              = strlen(sequences[0]);

  /* reserve more memory */
  vc->pscore              = (int *) space(sizeof(int)*((size*(size+1))/2+2));


  vc->exp_params          = NULL;
  vc->matrices            = get_mfe_matrices_alloc(vc->length, alloc_vector);

  /* get gquadruplex matrix if needed */
  if(vc->params->model_details.gquad){
    vc->cons_seq  = consensus(sequences);
    vc->S_cons    = encode_sequence(vc->cons_seq, 0);
    vc->matrices->ggg = get_gquad_ali_matrix(vc->S_cons, vc->S, vc->n_seq,  vc->params);
  }

  vc->hc                  = hc ? hc : get_hard_constraints(seq, NULL, &(vc->params->model_details), TURN, (unsigned int)0);
  vc->sc                  = sc;

  vc->iindx               = NULL;
  vc->jindx               = get_indx(vc->length);

  return vc;
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
