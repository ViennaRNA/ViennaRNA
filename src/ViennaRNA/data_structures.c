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
#include "ViennaRNA/ribo.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/part_func.h"

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
PRIVATE void make_pscores(vrna_fold_compound *vc);

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
      free(vc->pscore);
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

  if(vc->type == VRNA_VC_TYPE_SINGLE){
    sequence  = vc->sequence;

    seq2 = strdup(sequence);
    seq = tokenize(seq2, &cp); /*  splice out the '&' if concatenated sequences and
                                          reset cp... this should also be safe for
                                          single sequences */
    vc->cutpoint            = cp;
    free(vc->sequence);
    vc->sequence            = seq;
    vc->length              = length = strlen(seq);
    vc->sequence_encoding   = get_sequence_encoding(seq, 1, &md);
    vc->sequence_encoding2  = get_sequence_encoding(seq, 0, &md);
    vc->ptype               = vrna_get_ptypes(vc->sequence_encoding2, &md);
    vc->ptype_pf_compat     = (options & VRNA_OPTION_PF) ? get_ptypes(vc->sequence_encoding2, &md, 1) : NULL;
    vc->sc                  = NULL;
    free(seq2);
  }
  else if(vc->type == VRNA_VC_TYPE_ALIGNMENT){
    sequences     = vc->sequences;

    vc->length    = length = vc->length;

    vc->cons_seq  = consensus((const char **)sequences);
    vc->S_cons    = get_sequence_encoding(vc->cons_seq, 0, &md);

    vc->pscore    = (int *) space(sizeof(int)*((length*(length+1))/2+2));

    if(md.circ){
      warn_user("Not implemented yet: "
                "Falling back to oldAliEn model for alifold due to circular RNA");
      md.oldAliEn = 1;
    }
    oldAliEn = vc->oldAliEn  = md.oldAliEn;

    vc->S   = (short **)          space((vc->n_seq+1) * sizeof(short *));
    vc->S5  = (short **)          space((vc->n_seq+1) * sizeof(short *));
    vc->S3  = (short **)          space((vc->n_seq+1) * sizeof(short *));
    vc->a2s = (unsigned short **) space((vc->n_seq+1) * sizeof(unsigned short *));
    vc->Ss  = (char **)           space((vc->n_seq+1) * sizeof(char *));

    for (s = 0; s < vc->n_seq; s++) {
#if 1
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
#endif
    }
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

  if(md.circ){
    md.uniq_ML = 1;
    alloc_vector |= ALLOC_CIRC | ALLOC_UNIQ;
  }

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

  if(vc->type == VRNA_VC_TYPE_ALIGNMENT)
    make_pscores(vc);
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

    if(alloc_vector & ALLOC_CIRC)
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
    free(self->q);
    free(self->qb);
    free(self->qm);
    free(self->qm1);
    free(self->qm2);
    free(self->probs);
    free(self->G);
    free(self->q1k);
    free(self->qln);
    free(self->scale);
    free(self->expMLbase);
    free(self);
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

PRIVATE void
make_pscores(vrna_fold_compound *vc){

  /* calculate co-variance bonus for each pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */

#define NONE -10000 /* score for forbidden pairs */

  char *structure = NULL;
  int i,j,k,l,s, max_span;
  float **dm;
  int olddm[7][7]={{0,0,0,0,0,0,0}, /* hamming distance between pairs */
                  {0,0,2,2,1,2,2} /* CG */,
                  {0,2,0,1,2,2,2} /* GC */,
                  {0,2,1,0,2,1,2} /* GU */,
                  {0,1,2,2,0,2,1} /* UG */,
                  {0,2,2,1,2,0,2} /* AU */,
                  {0,2,2,2,1,2,0} /* UA */};

  short           **S         = vc->S;
  char            **AS        = vc->sequences;
  int             n_seq       = vc->n_seq;
  model_detailsT  *md         = (vc->params) ? &(vc->params->model_details) : &(vc->exp_params->model_details);
  int             *pscore     = vc->pscore;     /* precomputed array of pair types */             
  int             *indx       = vc->jindx;                                             
  int             n           = vc->length;                                            

  if (md->ribo) {
    if (RibosumFile !=NULL) dm=readribosum(RibosumFile);
    else dm=get_ribosum((const char **)AS, n_seq, n);
  }
  else { /*use usual matrix*/
    dm=(float **)space(7*sizeof(float*));
    for (i=0; i<7;i++) {
      dm[i]=(float *)space(7*sizeof(float));
      for (j=0; j<7; j++)
        dm[i][j] = (float) olddm[i][j];
    }
  }

  max_span = md->max_bp_span;
  if((max_span < TURN+2) || (max_span > n))
    max_span = n;
  for (i=1; i<n; i++) {
    for (j=i+1; (j<i+TURN+1) && (j<=n); j++)
      pscore[indx[j]+i] = NONE;
    for (j=i+TURN+1; j<=n; j++) {
      int pfreq[8]={0,0,0,0,0,0,0,0};
      double score;
      for (s=0; s<n_seq; s++) {
        int type;
        if (S[s][i]==0 && S[s][j]==0) type = 7; /* gap-gap  */
        else {
          if ((AS[s][i] == '~')||(AS[s][j] == '~')) type = 7;
          else type = md->pair[S[s][i]][S[s][j]];
        }
        pfreq[type]++;
      }
      if (pfreq[0]*2+pfreq[7]>n_seq) { pscore[indx[j]+i] = NONE; continue;}
      for (k=1,score=0; k<=6; k++) /* ignore pairtype 7 (gap-gap) */
        for (l=k; l<=6; l++)
          score += pfreq[k]*pfreq[l]*dm[k][l];
      /* counter examples score -1, gap-gap scores -0.25   */
      pscore[indx[j]+i] = md->cv_fact *
        ((UNIT*score)/n_seq - md->nc_fact*UNIT*(pfreq[0] + pfreq[7]*0.25));

      if((j - i + 1) > max_span){
        pscore[indx[j]+i] = NONE;
      }
    }
  }

  if (md->noLP) /* remove unwanted pairs */
    for (k=1; k<n-TURN-1; k++)
      for (l=1; l<=2; l++) {
        int type,ntype=0,otype=0;
        i=k; j = i+TURN+l;
        type = pscore[indx[j]+i];
        while ((i>=1)&&(j<=n)) {
          if ((i>1)&&(j<n)) ntype = pscore[indx[j+1]+i-1];
          if ((otype<md->cv_fact*MINPSCORE)&&(ntype<md->cv_fact*MINPSCORE))  /* too many counterexamples */
            pscore[indx[j]+i] = NONE; /* i.j can only form isolated pairs */
          otype =  type;
          type  = ntype;
          i--; j++;
        }
      }


  if (fold_constrained&&(structure!=NULL)) {
    int psij, hx, hx2, *stack, *stack2;
    stack = (int *) space(sizeof(int)*(n+1));
    stack2 = (int *) space(sizeof(int)*(n+1));

    for(hx=hx2=0, j=1; j<=n; j++) {
      switch (structure[j-1]) {
      case 'x': /* can't pair */
        for (l=1; l<j-TURN; l++) pscore[indx[j]+l] = NONE;
        for (l=j+TURN+1; l<=n; l++) pscore[indx[l]+j] = NONE;
        break;
      case '(':
        stack[hx++]=j;
        /* fallthrough */
      case '[':
        stack2[hx2++]=j;
        /* fallthrough */
      case '<': /* pairs upstream */
        for (l=1; l<j-TURN; l++) pscore[indx[j]+l] = NONE;
        break;
      case ']':
        if (hx2<=0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced brackets in constraints");
        }
        i = stack2[--hx2];
        pscore[indx[j]+i]=NONE;
        break;
      case ')':
        if (hx<=0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced brackets in constraints");
        }
        i = stack[--hx];
        psij = pscore[indx[j]+i]; /* store for later */
        for (k=j; k<=n; k++)
          for (l=i; l<=j; l++)
            pscore[indx[k]+l] = NONE;
        for (l=i; l<=j; l++)
          for (k=1; k<=i; k++)
            pscore[indx[l]+k] = NONE;
        for (k=i+1; k<j; k++)
          pscore[indx[k]+i] = pscore[indx[j]+k] = NONE;
        pscore[indx[j]+i] = (psij>0) ? psij : 0;
        /* fallthrough */
      case '>': /* pairs downstream */
        for (l=j+TURN+1; l<=n; l++) pscore[indx[l]+j] = NONE;
        break;
      }
    }
    if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in constraint string");
    }
    free(stack); free(stack2);
  }
  /*free dm */
  for (i=0; i<7;i++) {
    free(dm[i]);
  }
  free(dm);
}
