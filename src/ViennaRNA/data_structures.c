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
#include "ViennaRNA/aln_util.h"
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
PRIVATE void vrna_add_pf_matrices( vrna_fold_compound *vc, unsigned int alloc_vector);
PRIVATE void vrna_add_mfe_matrices(vrna_fold_compound *vc, unsigned int alloc_vector);
PRIVATE void set_fold_compound(vrna_fold_compound *vc, model_detailsT *md_p, unsigned int options);

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

    if(alloc_vector & ALLOC_HYBRID)
      vars->fc  = (int *) space(sizeof(int) * lin_size);

    if(alloc_vector & ALLOC_C)
      vars->c      = (int *) space(sizeof(int) * size);

    if(alloc_vector & ALLOC_FML)
      vars->fML    = (int *) space(sizeof(int) * size);

    if(alloc_vector & ALLOC_UNIQ)
      vars->fM1    = (int *) space(sizeof(int) * size);

    if(alloc_vector & ALLOC_CIRC)
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
      if(self->allocated & ALLOC_HYBRID)
        free(self->fc);
      if(self->allocated & ALLOC_C)
        free(self->c);
      if(self->allocated & ALLOC_FML)
        free(self->fML);
      if(self->allocated & ALLOC_UNIQ)
        free(self->fM1);
      if(self->allocated & ALLOC_CIRC)
        free(self->fM2);
      free(self->ggg);
    }
    free(self);
  }
}

PUBLIC  void
destroy_fold_compound(vrna_fold_compound *vc){

  int s;

  if(vc){

    /* first destroy common attributes */
    destroy_mfe_matrices(vc->matrices);
    destroy_pf_matrices(vc->exp_matrices);
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

    /* now distinguish the vc type */
    if(vc->type == VRNA_VC_TYPE_SINGLE){
      if(vc->sequence)
        free(vc->sequence);
      if(vc->sequence_encoding)
        free(vc->sequence_encoding);
      if(vc->sequence_encoding2)
        free(vc->sequence_encoding2);
      if(vc->ptype)
        free(vc->ptype);
      if(vc->ptype_pf_compat)
        free(vc->ptype_pf_compat);
      if(vc->sc)
        vrna_sc_destroy(vc->sc);

    } else if (vc->type == VRNA_VC_TYPE_ALIGNMENT){
      for(s=0;s<vc->n_seq;s++){
        free(vc->sequences[s]);
        free(vc->S[s]);
        free(vc->S5[s]);
        free(vc->S3[s]);
        free(vc->Ss[s]);
        free(vc->a2s[s]);
      }
      free(vc->sequences);
      free(vc->cons_seq);
      free(vc->S_cons);
      free(vc->S);
      free(vc->S5);
      free(vc->S3);
      free(vc->Ss);
      free(vc->a2s);
      if(vc->scs){
        for(s=0;s<vc->n_seq;s++)
          free(vc->scs[s]);
        free(vc->scs);
      }
    }

    free(vc);
  }
}


PUBLIC vrna_fold_compound*
vrna_get_fold_compound( const char *sequence,
                        model_detailsT *md_p,
                        unsigned int options){

  int length;
  vrna_fold_compound *vc;

  if(sequence == NULL) return NULL;

  /* sanity check */
  length = strlen(sequence);
  if(length == 0)
    nrerror("vrna_get_fold_compound: sequence length must be greater 0");

  vc            = (vrna_fold_compound *)space(sizeof(vrna_fold_compound));
  vc->type      = VRNA_VC_TYPE_SINGLE;
  vc->length    = length;
  vc->sequence  = strdup(sequence);
  set_fold_compound(vc, md_p, options);

  return vc;
}

PUBLIC vrna_fold_compound*
vrna_get_fold_compound_ali( const char **sequences,
                            model_detailsT *md_p,
                            unsigned int options){

  int s, n_seq, length;
  vrna_fold_compound *vc;
  
  if(sequences == NULL) return NULL;

  for(s=0;sequences[s];s++); /* count the sequences */

  n_seq = s;

  length = strlen(sequences[0]);
  /* sanity check */
  if(length == 0)
    nrerror("vrna_get_fold_compound_ali: sequence length must be greater 0");
  for(s = 0; s < n_seq; s++)
    if(strlen(sequences[s]) != length)
      nrerror("vrna_get_fold_compound_ali: uneqal sequence lengths in alignment");

  vc            = (vrna_fold_compound *)space(sizeof(vrna_fold_compound));
  vc->type      = VRNA_VC_TYPE_ALIGNMENT;

  vc->n_seq     = n_seq;
  vc->length    = length;
  vc->sequences = (char **)space(sizeof(char *) * (vc->n_seq + 1));
  for(s = 0; sequences[s]; s++)
    vc->sequences[s] = strdup(sequences[s]);

  set_fold_compound(vc, md_p, options);

  return vc;
}

PRIVATE void
set_fold_compound(vrna_fold_compound *vc,
                  model_detailsT *md_p,
                  unsigned int options){


  char *sequence, **sequences;
  model_detailsT      md;
  unsigned int        alloc_vector, length, s;
  int                 cp;                     /* cut point for cofold */
  char                *seq, *seq2;

  sequence    = NULL;
  sequences   = NULL;

  /* default values */
  cp            = -1;
  alloc_vector  = ALLOC_NOTHING;

  /* get a copy of the model details */
  if(md_p)
    md = *md_p;
  else /* this fallback relies on global parameters and thus is not threadsafe */
    set_model_details(&md);

  length    = vc->length;

  if(vc->type == VRNA_VC_TYPE_SINGLE){
    sequence  = vc->sequence;

    seq2 = strdup(sequence);
    seq = tokenize(seq2, &cp); /*  splice out the '&' if concatenated sequences and
                                          reset cp... this should also be safe for
                                          single sequences */
    vc->cutpoint            = cp;
    free(vc->sequence);
    vc->sequence            = seq;
    vc->sequence_encoding   = get_sequence_encoding(seq, 1, &md);
    vc->sequence_encoding2  = get_sequence_encoding(seq, 0, &md);
    vc->ptype               = vrna_get_ptypes(vc->sequence_encoding2, &md);
    vc->ptype_pf_compat     = (options & VRNA_OPTION_PF) ? get_ptypes(vc->sequence_encoding2, &md, 1) : NULL;
    vc->sc                  = NULL;
    free(seq2);
  }
  else if(vc->type == VRNA_VC_TYPE_ALIGNMENT){
    sequences     = vc->sequences;

    vc->cons_seq  = consensus((const char **)sequences);
    vc->S_cons    = get_sequence_encoding(vc->cons_seq, 0, &md);

    vc->pscore    = (int *) space(sizeof(int)*((length*(length+1))/2+2));

    vc->oldAliEn  = md.oldAliEn;

    vc->S   = (short **)          space((vc->n_seq+1) * sizeof(short *));
    vc->S5  = (short **)          space((vc->n_seq+1) * sizeof(short *));
    vc->S3  = (short **)          space((vc->n_seq+1) * sizeof(short *));
    vc->a2s = (unsigned short **) space((vc->n_seq+1) * sizeof(unsigned short *));
    vc->Ss  = (char **)           space((vc->n_seq+1) * sizeof(char *));

    for (s = 0; s < vc->n_seq; s++) {
#if 0
      get_sequence_encoding_gapped( vc->sequences[s],
                                    &(vc->S[s]),
                                    &(vc->S5[s]),
                                    &(vc->S3[s]),
                                    &(vc->Ss[s]),
                                    &(vc->a2s[s]),
                                    &md);
#else
      vc->S5[s]  = (short *)         space((length + 2) * sizeof(short));
      vc->S3[s]  = (short *)         space((length + 2) * sizeof(short));
      vc->a2s[s] = (unsigned short *)space((length + 2) * sizeof(unsigned short));
      vc->Ss[s]  = (char *)          space((length + 2) * sizeof(char));
      vc->S[s]   = (short *)         space((length + 2) * sizeof(short));
      encode_ali_sequence(vc->sequences[s], vc->S[s], vc->S5[s], vc->S3[s], vc->Ss[s], vc->a2s[s], md.circ);
    }
#endif
    vc->S5[vc->n_seq]  = NULL;
    vc->S3[vc->n_seq]  = NULL;
    vc->a2s[vc->n_seq] = NULL;
    vc->Ss[vc->n_seq]  = NULL;
    vc->S[vc->n_seq]   = NULL;

    vc->scs       = NULL;
  }

  /* prepare the allocation vector for the DP matrices */
  if(options & VRNA_OPTION_HYBRID){
    md.min_loop_size = 0;
    alloc_vector |= ALLOC_HYBRID;
    if(cp < 0)
      warn_user("vrna_get_fold_compound@data_structures.c: hybrid DP matrices requested but single sequence provided");
  }

  if(options & VRNA_OPTION_MFE)
    alloc_vector |= ALLOC_MFE_DEFAULT;

  if(options & VRNA_OPTION_PF)
    alloc_vector |= (md.compute_bpp) ? ALLOC_PF_DEFAULT : ALLOC_PF_WO_PROBS;

  if(md.uniq_ML)
    alloc_vector |= ALLOC_UNIQ;

  if(md.circ)
    alloc_vector |= ALLOC_CIRC;


  if(options & VRNA_OPTION_MFE){
    vc->params    = vrna_get_energy_contributions(md);
    vrna_add_mfe_matrices(vc, alloc_vector);
  } else {
    vc->params    = NULL;
    vc->matrices  = NULL;
  }


  if(options & VRNA_OPTION_PF){
    vc->exp_params    = (vc->type == VRNA_VC_TYPE_SINGLE) ? \
                            vrna_get_boltzmann_factors(md) : \
                            vrna_get_boltzmann_factors_ali(vc->n_seq, md);
    vrna_add_pf_matrices(vc, alloc_vector);
  } else {
    vc->exp_params    = NULL;
    vc->exp_matrices  = NULL;
  }

  vrna_hc_add(vc, NULL, (unsigned int)0); /* add hard constraints according to canonical base pairs */


  vc->iindx               = (options & VRNA_OPTION_PF) ? get_iindx(vc->length) : NULL;
  vc->jindx               = get_indx(vc->length);

}


PRIVATE void
vrna_add_pf_matrices( vrna_fold_compound *vc,
                      unsigned int alloc_vector){

  if(vc){
    vc->exp_matrices  = get_pf_matrices_alloc(vc->length, alloc_vector);
    if(vc->exp_params->model_details.gquad)
      vc->exp_matrices->G = get_gquad_pf_matrix(vc->sequence_encoding2, vc->exp_matrices->scale, vc->exp_params);

    vrna_update_pf_params(vc, NULL);
  }
}

PRIVATE void
vrna_add_mfe_matrices(vrna_fold_compound *vc,
                      unsigned int alloc_vector){

  if(vc){
    vc->matrices = get_mfe_matrices_alloc(vc->length, alloc_vector);
    if(vc->params->model_details.gquad){
      if(vc->type == VRNA_VC_TYPE_SINGLE)
        vc->matrices->ggg = get_gquad_matrix(vc->sequence_encoding2, vc->params);
      else if(vc->type == VRNA_VC_TYPE_ALIGNMENT);
        vc->matrices->ggg = get_gquad_ali_matrix(vc->S_cons, vc->S, vc->n_seq,  vc->params);
    }

//    vrna_update_pf_params(vc, NULL);
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

    if(alloc_vector & ALLOC_UNIQ)
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
      if(self->allocated & ALLOC_UNIQ)
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

PUBLIC vrna_alifold_compound *
get_alifold_compound_mfe(const char **sequences, paramT *P){

  return get_alifold_compound_mfe_constrained(sequences, NULL, NULL, P);
}

PUBLIC vrna_alifold_compound *
get_alifold_compound_mfe_constrained( const char **sequences,
                                      hard_constraintT **hc,
                                      soft_constraintT **sc,
                                      paramT *P){

  return vrna_alifold_get_compund_constraints(sequences, hc, sc, NULL, P, NULL, VRNA_CONSTRAINT_SOFT_MFE);
}

PUBLIC vrna_alifold_compound *
vrna_alifold_get_compund_constraints( const char **sequences,
                                      hard_constraintT **hc,
                                      soft_constraintT **sc,
                                      model_detailsT *model_details,
                                      paramT *P,
                                      pf_paramT *pf,
                                      unsigned int options){

  vrna_alifold_compound *vc;
  paramT                *params;
  pf_paramT             *exp_params;
  model_detailsT        md;
  unsigned int          alloc_vector, length, s;

  alloc_vector  = 0;
  params        = NULL;
  exp_params    = NULL;

  /* sanity check */
  length = (sequences) ? strlen(sequences[0]) : 0;
  if(length == 0)
    nrerror("get_alifold_compound_mfe_constraint@data_structures.c: sequence length must be greater 0");

  /* start making the alfold compound */
  vc = (vrna_alifold_compound *)space(sizeof(vrna_alifold_compound));

  for(s=0;sequences[s];s++);

  vc->n_seq           = s;
  vc->sequences       = (char **)space(sizeof(char *) * (vc->n_seq + 1));
  for(s = 0; sequences[s]; s++)
    vc->sequences[s]  = strdup(sequences[s]);

  vc->cons_seq      = consensus(sequences);
  vc->S_cons        = get_sequence_encoding(vc->cons_seq, 0, &md);

  vc->length          = strlen(sequences[0]);
  vc->pscore          = (int *) space(sizeof(int)*((length*(length+1))/2+2));

  /* copy the user-provided model details or just create them from the default options */
  if(model_details){
    md = *model_details;
  } else {
    set_model_details(&md);
  }

  vc->oldAliEn            = oldAliEn;

  /* prepare the parameters datastructure */

  /* Minimum free energy specific things following */
  if(options & VRNA_CONSTRAINT_SOFT_MFE){
    if(P){
      params = get_parameter_copy(P);
    } else { /* this fallback relies on global parameters and thus is not threadsafe */
      params = get_scaled_parameters(temperature, md);
    }
    vc->params  = params;

    alloc_vector = ALLOC_MFE_DEFAULT;
    if(md.circ)
      alloc_vector |= ALLOC_CIRC;
    if(md.uniq_ML)
      alloc_vector |= ALLOC_UNIQ;
//    vc->matrices  = get_mfe_matrices_alloc(length, alloc_vector);
    vc->matrices = NULL;

    /* get gquadruplex matrix if needed */
    if(md.gquad){
      vc->matrices->ggg = get_gquad_ali_matrix(vc->S_cons, vc->S, vc->n_seq,  vc->params);
    }

    vc->jindx               = get_indx(vc->length);
  }


  /* Partition function specific things following */
  if(options & VRNA_CONSTRAINT_SOFT_PF){
    if(pf){
      exp_params = get_boltzmann_factor_copy(pf);
    } else { /* this fallback relies on global parameters and thus is not threadsafe */
      exp_params = get_boltzmann_factors_ali(s, temperature, 1.0, md, pf_scale);
    }
    vc->exp_params  = exp_params;

    alloc_vector = ALLOC_PF_DEFAULT;
    if(md.circ)
      alloc_vector |= ALLOC_CIRC;
    if(md.uniq_ML)
      alloc_vector |= ALLOC_UNIQ;
    vc->exp_matrices  = get_pf_matrices_alloc(length, alloc_vector);
    vc->exp_matrices = NULL;

    /* get gquadruplex matrix if needed */
    if(md.gquad){
      vc->exp_matrices->G = NULL;
    }

    vc->iindx = get_iindx(vc->length);
    if(!vc->jindx)
      vc->jindx = get_indx(vc->length);

  }


  vc->S   = (short **)          space((vc->n_seq+1) * sizeof(short *));
  vc->S5  = (short **)          space((vc->n_seq+1) * sizeof(short *));
  vc->S3  = (short **)          space((vc->n_seq+1) * sizeof(short *));
  vc->a2s = (unsigned short **) space((vc->n_seq+1) * sizeof(unsigned short *));
  vc->Ss  = (char **)           space((vc->n_seq+1) * sizeof(char *));

  for (s=0; s<vc->n_seq; s++) {
    if(strlen(vc->sequences[s]) != vc->length) nrerror("uneqal seqence lengths");

    get_sequence_encoding_gapped( vc->sequences[s],
                                  &(vc->S[s]),
                                  &(vc->S5[s]),
                                  &(vc->S3[s]),
                                  &(vc->Ss[s]),
                                  &(vc->a2s[s]),
                                  &md);
  }
  vc->S5[vc->n_seq]  = NULL;
  vc->S3[vc->n_seq]  = NULL;
  vc->a2s[vc->n_seq] = NULL;
  vc->Ss[vc->n_seq]  = NULL;
  vc->S[vc->n_seq]   = NULL;





#if 0
  if(hc)
    for(s = 0; sequences[s]; s++)
      vc->hc[s] = hc[s] ? hc[s] : get_hard_constraints(sequences[s], NULL, &(vc->params->model_details), TURN, (unsigned int)0);
  else
    for(s=0; sequences[s]; s++)
      vc->hc[s] = get_hard_constraints(sequences[s], NULL, &(vc->params->model_details), TURN, (unsigned int)0);
#endif

  vc->sc  = sc;
  vc->hc  = hc;

  return vc;
}

PUBLIC void
destroy_alifold_compound(vrna_alifold_compound *vc){
  int s;

  if(vc){
    destroy_mfe_matrices(vc->matrices);
    destroy_pf_matrices(vc->exp_matrices);

    for(s = 0;s < vc->n_seq; s++)
      free(vc->sequences[s]);
    free(vc->sequences);

    for (s=0; s < vc->n_seq; s++){
      free(vc->S[s]);
      free(vc->S5[s]);
      free(vc->S3[s]);
      free(vc->a2s[s]);
      free(vc->Ss[s]);
    }
    free(vc->S);
    free(vc->S5);
    free(vc->S3);
    free(vc->a2s);
    free(vc->Ss);

    free(vc->cons_seq);
    free(vc->S_cons);

    if(vc->pscore)
      free(vc->pscore);
    if(vc->iindx)
      free(vc->iindx);
    if(vc->jindx)
      free(vc->jindx);
    if(vc->params)
      free(vc->params);
    if(vc->exp_params)
      free(vc->exp_params);
    if(vc->hc)
      for (s=0; s < vc->n_seq; s++)
        destroy_hard_constraints(vc->hc[s]);
    if(vc->sc)
      for (s=0; s < vc->n_seq; s++)
        vrna_sc_destroy(vc->sc[s]);


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
