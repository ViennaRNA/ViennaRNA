/* Last changed Time-stamp: <2009-02-24 14:37:05 ivo> */
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

#include <config.h>
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
#include "ViennaRNA/part_func.h"
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
PRIVATE vrna_fold_compound  *backward_compat_compound = NULL;
PRIVATE int                 backward_compat           = 0;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE void      alipf_linear(vrna_fold_compound *vc, char *structure);
PRIVATE void      alipf_create_bppm(vrna_fold_compound *vc, char *structure, struct plist **pl);
PRIVATE void      backtrack(vrna_fold_compound *vc, char *pstruc, int i, int j, double *prob);
PRIVATE void      backtrack_qm1(vrna_fold_compound *vc, char *pstruc, int i,int j, double *prob);
PRIVATE float     wrap_alipf_fold(const char **sequences,
                                  char *structure,
                                  plist **pl,
                                  pf_paramT *parameters,
                                  int calculate_bppm,
                                  int is_constrained,
                                  int is_circular);
PRIVATE void      wrap_alipf_circ(vrna_fold_compound *vc,
                                  char *structure);



/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/


/*-----------------------------------------------------------------*/
PRIVATE float
wrap_alipf_fold(const char **sequences,
                char *structure,
                plist **pl,
                pf_paramT *parameters,
                int calculate_bppm,
                int is_constrained,
                int is_circular){

  int                 n_seq;
  vrna_fold_compound  *vc;
  pf_paramT           *exp_params;

  if(sequences == NULL) return 0.;

  for(n_seq=0;sequences[n_seq];n_seq++); /* count the sequences */
  
  vc                  = NULL;

  /* we need pf_paramT datastructure to correctly init default hard constraints */
  if(parameters)
    exp_params = get_boltzmann_factor_copy(parameters);
  else{
    model_detailsT md;
    set_model_details(&md); /* get global default parameters */
    exp_params = vrna_get_boltzmann_factors_ali(n_seq, md);
  }
  exp_params->model_details.circ        = is_circular;
  exp_params->model_details.compute_bpp = calculate_bppm;

  vc = vrna_get_fold_compound_ali(sequences, &(exp_params->model_details), VRNA_OPTION_PF);

  if(parameters){ /* replace exp_params if necessary */
    free(vc->exp_params);
    vc->exp_params = exp_params;
  } else {
    free(exp_params);
  }

  if(is_constrained && structure){
    unsigned int constraint_options = 0;
    constraint_options |= VRNA_CONSTRAINT_DB
                          | VRNA_CONSTRAINT_PIPE
                          | VRNA_CONSTRAINT_DOT
                          | VRNA_CONSTRAINT_X
                          | VRNA_CONSTRAINT_ANG_BRACK
                          | VRNA_CONSTRAINT_RND_BRACK;

    vrna_hc_add(vc, (const char *)structure, constraint_options);
  }

  if(backward_compat_compound && backward_compat_compound)
    destroy_fold_compound(backward_compat_compound);

  backward_compat_compound  = vc;
  iindx                     = backward_compat_compound->iindx;
  backward_compat           = 1;

  return vrna_ali_pf_fold(vc, structure, pl);
}

PUBLIC float
vrna_ali_pf_fold( vrna_fold_compound *vc,
                  char *structure,
                  plist **pl){

  FLT_OR_DBL  Q;
  float       free_energy;

  int             n         = vc->length;
  int             n_seq     = vc->n_seq;
  pf_paramT       *params   = vc->exp_params;
  model_detailsT  *md       = &(params->model_details);
  pf_matricesT    *matrices = vc->exp_matrices;

  alipf_linear(vc, structure);

  /* calculate post processing step for circular  */
  /* RNAs                                         */
  if(md->circ)
    wrap_alipf_circ(vc, structure);

  if (md->backtrack_type=='C')      Q = matrices->qb[vc->iindx[1]-n];
  else if (md->backtrack_type=='M') Q = matrices->qm[vc->iindx[1]-n];
  else Q = (md->circ) ? matrices->qo : matrices->q[vc->iindx[1]-n];

  /* ensemble free energy in Kcal/mol */
  if (Q<=FLT_MIN) fprintf(stderr, "pf_scale too large\n");
  free_energy = (-log(Q)-n*log(params->pf_scale))*params->kT/(1000.0 * n_seq);
  /* in case we abort because of floating point errors */
  if (n>1600) fprintf(stderr, "free energy = %8.2f\n", free_energy);

  /* backtracking to construct binding probabilities of pairs*/
  if(md->compute_bpp){
    alipf_create_bppm(vc, structure, pl);
    /*
    *  Backward compatibility:
    *  This block may be removed if deprecated functions
    *  relying on the global variable "pr" vanish from within the package!
    */
    pr = matrices->probs;
  }

  return free_energy;
}



PRIVATE void
alipf_linear( vrna_fold_compound *vc,
              char *structure){

  int         s, i,j,k,l, ij, u, u1, d, ii, *type, type_2, tt;
  FLT_OR_DBL  temp;
  FLT_OR_DBL  qbt1, *tmp;
  FLT_OR_DBL  *qqm = NULL, *qqm1 = NULL, *qq = NULL, *qq1 = NULL;
  double      kTn;

  int             n_seq         = vc->n_seq;
  int             n             = vc->length;


  short             **S           = vc->S;                                                                   
  short             **S5          = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/            
  short             **S3          = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/            
  char              **Ss          = vc->Ss;                                                                   
  unsigned short    **a2s         = vc->a2s;                                                                   
  pf_paramT         *pf_params    = vc->exp_params;
  pf_matricesT      *matrices     = vc->exp_matrices;
  model_detailsT    *md           = &(pf_params->model_details);
  hard_constraintT  *hc           = vc->hc;
  soft_constraintT  **sc          = vc->scs;
  int               *my_iindx     = vc->iindx;
  int               *jindx        = vc->jindx;
  FLT_OR_DBL        *q            = matrices->q;
  FLT_OR_DBL        *qb           = matrices->qb;
  FLT_OR_DBL        *qm           = matrices->qm;
  FLT_OR_DBL        *qm1          = matrices->qm1;
  int               *pscore       = vc->pscore;     /* precomputed array of pair types */                      
  int               *rtype        = &(md->rtype[0]);
  int               circular      = md->circ;
  FLT_OR_DBL        *scale        = matrices->scale;
  FLT_OR_DBL        *expMLbase    = matrices->expMLbase;
  FLT_OR_DBL        expMLclosing  = pf_params->expMLclosing;

  kTn   = pf_params->kT/10.;   /* kT in cal/mol  */
  type  = (int *)space(sizeof(int) * n_seq);

  int max_bpspan = (md->max_bp_span > 0) ? md->max_bp_span : n;

  /* allocate memory for helper arrays */
  qq        = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));
  qq1       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));
  qqm       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));
  qqm1      = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));


  /* array initialization ; qb,qm,q
     qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  for (d=0; d<=TURN; d++)
    for (i=1; i<=n-d; i++) {
      j=i+d;
      ij = my_iindx[i]-j;
      q[ij]=1.0*scale[d+1];
      qb[ij]=qm[ij]=0.0;
    }

  for (i=1; i<=n; i++)
    qq[i]=qq1[i]=qqm[i]=qqm1[i]=0;

  for (j=TURN+2;j<=n; j++) {
    for (i=j-TURN-1; i>=1; i--) {
      int ij, psc;
      /* construction of partition function for segment i,j */
      /* calculate pf given that i and j pair: qb(i,j)      */
       ij = my_iindx[i]-j;

      for (s=0; s<n_seq; s++) {
        type[s] = md->pair[S[s][i]][S[s][j]];
        if (type[s]==0) type[s]=7;
      }
      psc = pscore[jindx[j]+i];
      if (psc>=md->cv_fact*MINPSCORE && ((j-i) < max_bpspan)) {   /* otherwise ignore this pair */

        /* hairpin contribution */
        for (qbt1=1,s=0; s<n_seq; s++) {
          u = a2s[s][j-1]-a2s[s][i];
          if (a2s[s][i]<1) continue;
          char loopseq[10];
          if (u<7){
            strncpy(loopseq, Ss[s]+a2s[s][i]-1, 10);
          }
          qbt1 *= exp_E_Hairpin(u, type[s], S3[s][i], S5[s][j], loopseq, pf_params);

        }
        qbt1 *= scale[j-i+1];

        /* interior loops with interior pair k,l */
        for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++){

          for (l=MAX2(k+TURN+1,j-1-MAXLOOP+k-i-1); l<=j-1; l++){
            double qloop=1.;
            if (qb[my_iindx[k]-l]==0) {qloop=0; continue;}
            for (s=0; s<n_seq; s++) {
              u1 = a2s[s][k-1]-a2s[s][i]/*??*/;
              type_2 = md->pair[S[s][l]][S[s][k]]; if (type_2 == 0) type_2 = 7;
              qloop *= exp_E_IntLoop( u1, a2s[s][j-1]-a2s[s][l],
                                  type[s], type_2, S3[s][i],
                                  S5[s][j], S5[s][k], S3[s][l],
                                  pf_params
                                );

              if(sc)
                if((i + 1 == k) && (j - 1 == l))
                  if(sc[s])
                    if(sc[s]->exp_en_stack)
                      if((i + 1 == k) && (j - 1 == l)){
                        qloop *=    sc[s]->exp_en_stack[i]
                                  * sc[s]->exp_en_stack[k]
                                  * sc[s]->exp_en_stack[l]
                                  * sc[s]->exp_en_stack[j];
                      }

            }
            
            qbt1 += qb[my_iindx[k]-l] * qloop * scale[k-i+j-l];
          }
        }

        /* multi-loop loop contribution */
        ii = my_iindx[i+1]; /* ii-k=[i+1,k-1] */
        temp = 0.0;
        for (k=i+2; k<=j-1; k++) temp += qm[ii-(k-1)]*qqm1[k];
        for (s=0; s<n_seq; s++) {
          tt = rtype[type[s]];
          temp *= exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params)* expMLclosing;
        }
        temp *= scale[2] ;
        qbt1 += temp;
        qb[ij] = qbt1;
        qb[ij] *= exp(psc/kTn);
      } /* end if (type!=0) */
      else qb[ij] = 0.0;
      /* construction of qqm matrix containing final stem
         contributions to multiple loop partition function
         from segment i,j */
      qqm[i] = qqm1[i]*expMLbase[1];  /* expMLbase[1]^n_seq */
      for (qbt1=1, s=0; s<n_seq; s++) {
        qbt1 *= exp_E_MLstem(type[s], (i>1) || circular ? S5[s][i] : -1, (j<n) || circular ? S3[s][j] : -1, pf_params);
      }
      qqm[i] += qb[ij]*qbt1;
      if(qm1) qm1[jindx[j]+i] = qqm[i]; /* for circ folding and stochBT */

      /* construction of qm matrix containing multiple loop
         partition function contributions from segment i,j */
      temp = 0.0;
      ii = my_iindx[i];  /* ii-k=[i,k-1] */
      for (k=i+1; k<=j; k++)
        temp += (qm[ii-(k-1)]+expMLbase[k-i])*qqm[k];
      qm[ij] = (temp + qqm[i]);

      /* auxiliary matrix qq for cubic order q calculation below */
      qbt1 = qb[ij];
      if (qbt1>0)
        for (s=0; s<n_seq; s++) {
          qbt1 *= exp_E_ExtLoop(type[s], (i>1) || circular ? S5[s][i] : -1, (j<n) || circular ? S3[s][j] : -1, pf_params);
        }
      qq[i] = qq1[i]*scale[1] + qbt1;

      /* construction of partition function for segment i,j */
      temp = 1.0*scale[1+j-i] + qq[i];
      for (k=i; k<=j-1; k++) temp += q[ii-k]*qq[k+1];
      q[ij] = temp;

#ifndef LARGE_PF
      if (temp>Qmax) {
        Qmax = temp;
        if (Qmax>FLT_MAX/10.)
          fprintf(stderr, "%d %d %g\n", i,j,temp);
      }
      if (temp>FLT_MAX) {
        PRIVATE char msg[128];
        sprintf(msg, "overflow in pf_fold while calculating q[%d,%d]\n"
          "use larger pf_scale", i,j);
        nrerror(msg);
      }
#endif
    }
    tmp = qq1;  qq1 =qq;  qq =tmp;
    tmp = qqm1; qqm1=qqm; qqm=tmp;
  }

  /* clean up */
  free(type);
  free(qq);
  free(qq1);
  free(qqm);
  free(qqm1);

}

PRIVATE void
alipf_create_bppm(vrna_fold_compound *vc,
                  char *structure,
                  plist **pl){

  int s;
  int i,j,k,l, ij, kl, ii, ll, tt, *type, ov=0;
  FLT_OR_DBL temp, prm_MLb;
  FLT_OR_DBL prmt,prmt1;
  FLT_OR_DBL qbt1, *tmp, tmp2, tmp3;

  char            **sequences   = vc->sequences;
  int             n_seq         = vc->n_seq;
  int             n             = vc->length;


  short             **S           = vc->S;                                                                   
  short             **S5          = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/            
  short             **S3          = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/            
  unsigned short    **a2s         = vc->a2s;                                                                   
  pf_paramT         *pf_params    = vc->exp_params;
  pf_matricesT      *matrices     = vc->exp_matrices;
  model_detailsT    *md           = &(pf_params->model_details);
  hard_constraintT  *hc           = vc->hc;
  soft_constraintT  **sc          = vc->scs;
  int               *my_iindx     = vc->iindx;
  int               *jindx        = vc->jindx;
  FLT_OR_DBL        *q            = matrices->q;
  FLT_OR_DBL        *qb           = matrices->qb;
  FLT_OR_DBL        *qm           = matrices->qm;
  FLT_OR_DBL        *qm1          = matrices->qm1;
  FLT_OR_DBL        qo            = matrices->qo;
  int               *pscore       = vc->pscore;     /* precomputed array of pair types */                      
  int               *rtype        = &(md->rtype[0]);
  int               circular      = md->circ;
  FLT_OR_DBL        *scale        = matrices->scale;
  FLT_OR_DBL        *expMLbase    = matrices->expMLbase;
  FLT_OR_DBL        expMLclosing  = pf_params->expMLclosing;
  FLT_OR_DBL        *probs        = matrices->probs;

  double kTn;

  FLT_OR_DBL *q1k    = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+1));
  FLT_OR_DBL *qln    = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));
  FLT_OR_DBL *prm_l  = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));
  FLT_OR_DBL *prm_l1 = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));
  FLT_OR_DBL *prml   = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));

    type  = (int *)space(sizeof(int) * n_seq);

  kTn = pf_params->kT/10.;   /* kT in cal/mol  */

  for (i=0; i<=n; i++)
    prm_l[i]=prm_l1[i]=prml[i]=0;

  /* backtracking to construct binding probabilities of pairs*/
  for (k=1; k<=n; k++) {
    q1k[k] = q[my_iindx[1] - k];
    qln[k] = q[my_iindx[k] -n];
  }
  q1k[0] = 1.0;
  qln[n+1] = 1.0;

  /* 1. exterior pair i,j and initialization of pr array */
  if(circular){
    for (i=1; i<=n; i++) {
      for (j=i; j<=MIN2(i+TURN,n); j++) probs[my_iindx[i]-j] = 0;
      for (j=i+TURN+1; j<=n; j++) {
        ij = my_iindx[i]-j;
        if (qb[ij]>0.) {
          probs[ij] =  exp(pscore[jindx[j]+i]/kTn)/qo;

          /* get pair types  */
          for (s=0; s<n_seq; s++) {
            type[s] = md->pair[S[s][j]][S[s][i]];
            if (type[s]==0) type[s]=7;
          }

          /* 1.1. Exterior Hairpin Contribution */
          int u = i + n - j -1;
          for (qbt1=1.,s=0; s<n_seq; s++) {

            char loopseq[10];
            if (u<7){
              strcpy(loopseq , sequences[s]+j-1);
              strncat(loopseq, sequences[s], i);
            }
            qbt1 *= exp_E_Hairpin(u, type[s], S[s][j+1], S[s][(i>1) ? i-1 : n], loopseq, pf_params);
          }
          tmp2 = qbt1 * scale[u];

          /* 1.2. Exterior Interior Loop Contribution */
          /* recycling of k and l... */
          /* 1.2.1. first we calc exterior loop energy with constraint, that i,j  */
          /* delimtis the "left" part of the interior loop                        */
          /* (j,i) is "outer pair"                                                */
          for(k=1; k < i-TURN-1; k++){
            /* so first, lets calc the length of loop between j and k */
            int ln1, lstart;
            ln1 = k + n - j - 1;
            if(ln1>MAXLOOP) break;
            lstart = ln1+i-1-MAXLOOP;
            if(lstart<k+TURN+1) lstart = k + TURN + 1;
            for(l=lstart; l < i; l++){
              int ln2,ln2a,ln1a, type_2;
              ln2 = i - l - 1;
                if(ln1+ln2>MAXLOOP) continue;

              double qloop=1.;
              if (qb[my_iindx[k]-l]==0.){ qloop=0.; continue;}

              for (s=0; s<n_seq; s++){
                ln2a= a2s[s][i-1];
                ln2a-=a2s[s][l];
                ln1a= a2s[s][n]-a2s[s][j];
                ln1a+=a2s[s][k-1];
                type_2 = md->pair[S[s][l]][S[s][k]];
                if (type_2 == 0) type_2 = 7;
                qloop *= exp_E_IntLoop(ln1a, ln2a, type[s], type_2,
                            S[s][j+1],
                            S[s][i-1],
                            S[s][(k>1) ? k-1 : n],
                            S[s][l+1], pf_params);
              }
              tmp2 += qb[my_iindx[k] - l] * qloop * scale[ln1+ln2];
            }
          }

          /* 1.2.2. second we calc exterior loop energy with constraint, that i,j  */
          /* delimtis the "right" part of the interior loop                        */
          /* (l,k) is "outer pair"                                                */
          for(k=j+1; k < n-TURN; k++){
            /* so first, lets calc the length of loop between l and i */
            int ln1, lstart;
            ln1 = k - j - 1;
            if((ln1 + i - 1)>MAXLOOP) break;
            lstart = ln1+i-1+n-MAXLOOP;
            if(lstart<k+TURN+1) lstart = k + TURN + 1;
            for(l=lstart; l <= n; l++){
              int ln2, type_2;
              ln2 = i - 1 + n - l;
              if(ln1+ln2>MAXLOOP) continue;
              double qloop=1.;
              if (qb[my_iindx[k]-l]==0.){ qloop=0.; continue;}

              for (s=0; s<n_seq; s++){
                    ln1 = a2s[s][k] - a2s[s][j+1];
                    ln2 = a2s[s][i-1] + a2s[s][n] - a2s[s][l];
                    type_2 = md->pair[S[s][l]][S[s][k]];
                    if (type_2 == 0) type_2 = 7;
                    qloop *= exp_E_IntLoop(ln2, ln1, type_2, type[s],
                            S3[s][l],
                            S5[s][k],
                            S5[s][i],
                            S3[s][j], pf_params);
              }
              tmp2 += qb[my_iindx[k] - l] * qloop * scale[(k-j-1)+(i-1+n-l)];
            }
          }

          /* 1.3 Exterior multiloop decomposition */
          /* 1.3.1 Middle part                    */
          if((i>TURN+2) && (j<n-TURN-1)){

            for (tmp3=1, s=0; s<n_seq; s++){
              tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);
            }
            tmp2 += qm[my_iindx[1]-i+1] * qm[my_iindx[j+1]-n] * tmp3 * pow(expMLclosing,n_seq);
          }
          /* 1.3.2 Left part    */
          for(k=TURN+2; k < i-TURN-2; k++){

            for (tmp3=1, s=0; s<n_seq; s++){
              tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);
            }
            tmp2 += qm[my_iindx[1]-k] * qm1[jindx[i-1]+k+1] * tmp3 * expMLbase[n-j] * pow(expMLclosing,n_seq);
          }
          /* 1.3.3 Right part    */
          for(k=j+TURN+2; k < n-TURN-1;k++){

            for (tmp3=1, s=0; s<n_seq; s++){
              tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);
            }
            tmp2 += qm[my_iindx[j+1]-k] * qm1[jindx[n]+k+1] * tmp3 * expMLbase[i-1] * pow(expMLclosing,n_seq);
          }
          probs[ij] *= tmp2;
        }
        else probs[ij] = 0;
      }  /* end for j=..*/
    }  /* end or i=...  */
  } /* end if(circular)  */
  else{
    for (i=1; i<=n; i++) {
      for (j=i; j<=MIN2(i+TURN,n); j++) probs[my_iindx[i]-j] = 0;
      for (j=i+TURN+1; j<=n; j++) {
        ij = my_iindx[i]-j;
        if (qb[ij]>0.){
          probs[ij] = q1k[i-1]*qln[j+1]/q1k[n] * exp(pscore[jindx[j]+i]/kTn);
          for (s=0; s<n_seq; s++) {
            int typ;
            typ = md->pair[S[s][i]][S[s][j]]; if (typ==0) typ=7;
            probs[ij] *= exp_E_ExtLoop(typ, (i>1) ? S5[s][i] : -1, (j<n) ? S3[s][j] : -1, pf_params);
          }
        } else
          probs[ij] = 0;
      }
    }
  } /* end if(!circular)  */
  for (l=n; l>TURN+1; l--) {

    /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
    for (k=1; k<l-TURN; k++) {
      double pp = 0;
      kl = my_iindx[k]-l;
      if (qb[kl]==0) continue;
      for (s=0; s<n_seq; s++) {
            type[s] = md->pair[S[s][l]][S[s][k]];
            if (type[s]==0) type[s]=7;
      }

      for (i=MAX2(1,k-MAXLOOP-1); i<=k-1; i++)
        for (j=l+1; j<=MIN2(l+ MAXLOOP -k+i+2,n); j++) {
          ij = my_iindx[i] - j;
          if ((probs[ij]>0.)) {
            double qloop=1;
            for (s=0; s<n_seq; s++) {
              int typ;
              typ = md->pair[S[s][i]][S[s][j]]; if (typ==0) typ=7;
              qloop *=  exp_E_IntLoop(a2s[s][k-1]-a2s[s][i], a2s[s][j-1]-a2s[s][l], typ, type[s], S3[s][i], S5[s][j], S5[s][k], S3[s][l], pf_params);

              if(sc)
                if((i + 1 == k) && (j - 1 == l))
                  if(sc[s])
                    if(sc[s]->exp_en_stack)
                      qloop *=    sc[s]->exp_en_stack[i]
                                * sc[s]->exp_en_stack[k]
                                * sc[s]->exp_en_stack[l]
                                * sc[s]->exp_en_stack[j];
            }
            pp += probs[ij]*qloop*scale[k-i + j-l];
          }
        }
      probs[kl] += pp * exp(pscore[jindx[l]+k]/kTn);
    }
    /* 3. bonding k,l as substem of multi-loop enclosed by i,j */
    prm_MLb = 0.;
    if (l<n) for (k=2; k<l-TURN; k++) {
      i = k-1;
      prmt = prmt1 = 0.0;

      ii = my_iindx[i];     /* ii-j=[i,j]     */
      ll = my_iindx[l+1];   /* ll-j=[l+1,j-1] */
      prmt1 = probs[ii-(l+1)];
      for (s=0; s<n_seq; s++) {
        tt = md->pair[S[s][l+1]][S[s][i]]; if (tt==0) tt=7;
        prmt1 *= exp_E_MLstem(tt, S5[s][l+1], S3[s][i], pf_params) * expMLclosing;
      }

      for (j=l+2; j<=n; j++) {
        double pp=1;
        if (probs[ii-j]==0) continue;
        for (s=0; s<n_seq; s++) {
          tt=md->pair[S[s][j]][S[s][i]]; if (tt==0) tt=7;
          pp *=  exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params)*expMLclosing;
        }
        prmt +=  probs[ii-j]*pp*qm[ll-(j-1)];
      }
      kl = my_iindx[k]-l;

      prml[ i] = prmt;
      prm_l[i] = prm_l1[i]*expMLbase[1]+prmt1; /* expMLbase[1]^n_seq */

      prm_MLb = prm_MLb*expMLbase[1] + prml[i];
      /* same as:    prm_MLb = 0;
         for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */

      prml[i] = prml[ i] + prm_l[i];

      if (qb[kl] == 0.) continue;

      temp = prm_MLb;

      for (i=1;i<=k-2; i++)
        temp += prml[i]*qm[my_iindx[i+1] - (k-1)];

      for (s=0; s<n_seq; s++) {
        tt=md->pair[S[s][k]][S[s][l]]; if (tt==0) tt=7;
        temp *= exp_E_MLstem(tt, S5[s][k], S3[s][l], pf_params);
      }
      probs[kl] += temp * scale[2] * exp(pscore[jindx[l]+k]/kTn);


#ifndef LARGE_PF
      if (probs[kl]>Qmax) {
        Qmax = probs[kl];
        if (Qmax>FLT_MAX/10.)
          fprintf(stderr, "%d %d %g %g\n", i,j,probs[kl],qb[kl]);
      }
      if (probs[kl]>FLT_MAX) {
        ov++;
        probs[kl]=FLT_MAX;
      }
#endif
    } /* end for (k=2..) */
    tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;

  }  /* end for (l=..)   */

  for (i=1; i<=n; i++)
    for (j=i+TURN+1; j<=n; j++) {
      ij = my_iindx[i]-j;
      probs[ij] *= qb[ij] *exp(-pscore[jindx[j]+i]/kTn);
    }

  /* did we get an adress where to save a pair-list? */
  if (pl != NULL)
    *pl = vrna_get_plist_from_pr(vc, /*cut_off:*/ 1e-6);

  if (structure!=NULL)
    bppm_to_structure(structure, probs, n);

  if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
        "you might try a smaller pf_scale than %g\n",
        ov, pf_params->pf_scale);

  free(type);
}

/*---------------------------------------------------------------------------*/
PRIVATE int
compare_pair_info(const void *pi1,
                  const void *pi2){

  pair_info *p1, *p2;
  int  i, nc1, nc2;
  p1 = (pair_info *)pi1;  p2 = (pair_info *)pi2;
  for (nc1=nc2=0, i=1; i<=6; i++) {
    if (p1->bp[i]>0) nc1++;
    if (p2->bp[i]>0) nc2++;
  }
  /* sort mostly by probability, add
     epsilon * comp_mutations/(non-compatible+1) to break ties */
  return (p1->p + 0.01*nc1/(p1->bp[0]+1.)) <
         (p2->p + 0.01*nc2/(p2->bp[0]+1.)) ? 1 : -1;
}

PUBLIC pair_info *
vrna_ali_get_pair_info( vrna_fold_compound *vc,
                        const char *structure,
                        double threshold){

  int i,j, num_p=0, max_p = 64;
  pair_info *pi;
  double *duck, p;
  short *ptable = NULL;

  short **S = vc->S;
  char **AS = vc->sequences;
  int n_seq = vc->n_seq;
  int n     = vc->length;
  int         *my_iindx = vc->iindx;
  FLT_OR_DBL  *probs    = vc->exp_matrices->probs;
  model_detailsT  *md = &(vc->exp_params->model_details);

  max_p = 64; pi = space(max_p*sizeof(pair_info));
  duck =  (double *) space((n+1)*sizeof(double));
  if(structure)
    ptable = vrna_pt_get(structure);

  for (i=1; i<n; i++)
    for (j=i+TURN+1; j<=n; j++) {
      if ((p=probs[my_iindx[i]-j])>=threshold) {
        duck[i] -=  p * log(p);
        duck[j] -=  p * log(p);

        int type, s;
        pi[num_p].i   = i;
        pi[num_p].j   = j;
        pi[num_p].p   = p;
        pi[num_p].ent = duck[i]+duck[j]-p*log(p);

        for (type=0; type<8; type++) pi[num_p].bp[type]=0;
        for (s=0; s<n_seq; s++) {
          type = md->pair[S[s][i]][S[s][j]];
          if(S[s][i]==0 && S[s][j]==0) type = 7; /* gap-gap  */
          if ((AS[s][i-1] == '-')||(AS[s][j-1] == '-')) type = 7;
          if ((AS[s][i-1] == '~')||(AS[s][j-1] == '~')) type = 7;
          pi[num_p].bp[type]++;
        }
        if(ptable)
          pi[num_p].comp = (ptable[i] == j) ? 1:0;

        num_p++;
        if (num_p>=max_p) {
          max_p *= 2;
          pi = xrealloc(pi, max_p * sizeof(pair_info));
        }
      }
    }
  free(duck);
  pi = xrealloc(pi, (num_p+1)*sizeof(pair_info));
  pi[num_p].i=0;
  qsort(pi, num_p, sizeof(pair_info), compare_pair_info );

  free(ptable);
  return pi;
}

/* calculate partition function for circular case   */
/* NOTE: this is the postprocessing step ONLY        */
/* You have to call alipf_linear first to calculate  */
/* circular case!!!                                  */

PRIVATE void
wrap_alipf_circ(vrna_fold_compound *vc,
                char *structure){

  int u, p, q, k, l, s, *type;
  FLT_OR_DBL qbt1, qot, qo, qho, qio, qmo;

  char              **sequences = vc->sequences;
  int               n_seq       = vc->n_seq;
  int               n           = vc->length;
  short             **S         = vc->S;                                                                   
  short             **S5        = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/            
  short             **S3        = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/            
  unsigned short    **a2s       = vc->a2s;                                                                   
  pf_paramT         *pf_params  = vc->exp_params;
  pf_matricesT      *matrices   = vc->exp_matrices;
  model_detailsT    *md         = &(pf_params->model_details);
  int               *my_iindx   = vc->iindx;
  int               *jindx      = vc->jindx;
  hard_constraintT  *hc         = vc->hc;
  soft_constraintT  **sc        = vc->scs;
  FLT_OR_DBL        *qb         = matrices->qb;
  FLT_OR_DBL        *qm         = matrices->qm;
  FLT_OR_DBL        *qm1        = matrices->qm1;
  FLT_OR_DBL        *qm2        = matrices->qm2;
  int               *pscore     = vc->pscore;     /* precomputed array of pair types */             
  FLT_OR_DBL        *scale      = matrices->scale;
  FLT_OR_DBL        expMLclosing  = pf_params->expMLclosing;

  type  = (int *)space(sizeof(int) * n_seq);

  qo = qho = qio = qmo = 0.;
  /* calculate the qm2 matrix  */
  for(k=1; k<n-TURN; k++){
    qot = 0.;
    for (u=k+TURN+1; u<n-TURN-1; u++)
      qot += qm1[jindx[u]+k]*qm1[jindx[n]+(u+1)];
    qm2[k] = qot;
  }

  for(p=1;p<n;p++){
    for(q=p+TURN+1;q<=n;q++){
      int psc;
      u = n-q + p-1;
      if (u<TURN) continue;

      psc = pscore[jindx[q]+p];

      if(psc<md->cv_fact*MINPSCORE) continue;

      /* 1. exterior hairpin contribution  */
      /* Note, that we do not scale Hairpin Energy by u+2 but by u cause the scale  */
      /* for the closing pair was already done in the forward recursion              */
      for (qbt1=1,s=0; s<n_seq; s++) {
        char loopseq[10];
            type[s] = md->pair[S[s][q]][S[s][p]];
            if (type[s]==0) type[s]=7;

            if (u<9){
              strcpy(loopseq , sequences[s]+q-1);
              strncat(loopseq, sequences[s], p);
            }
        qbt1 *= exp_E_Hairpin(u, type[s], S[s][q+1], S[s][(p>1) ? p-1 : n], loopseq, pf_params);
      }
      qho += qb[my_iindx[p]-q] * qbt1 * scale[u];

      /* 2. exterior interior loop contribution*/

      for(k=q+1; k < n; k++){
        int ln1, lstart;
        ln1 = k - q - 1;
        if(ln1+p-1>MAXLOOP) break;
        lstart = ln1+p-1+n-MAXLOOP;
        if(lstart<k+TURN+1) lstart = k + TURN + 1;
        for(l=lstart;l <= n; l++){
          int ln2, type_2;

              ln2 = (p - 1) + (n - l);
              if((ln1+ln2) > MAXLOOP) continue;
              double qloop=1.;
              if (qb[my_iindx[k]-l]==0.){ qloop=0.; continue;}

              for (s=0; s<n_seq; s++){
                int ln1a=a2s[s][k-1]-a2s[s][q];
                int ln2a=a2s[s][n]-a2s[s][l]+a2s[s][p-1];
                type_2 = md->pair[S[s][l]][S[s][k]];
                if (type_2 == 0) type_2 = 7;
                qloop *= exp_E_IntLoop(ln2a, ln1a, type_2, type[s], S3[s][l], S5[s][k], S5[s][p], S3[s][q], pf_params);
              }
              qio += qb[my_iindx[p]-q] * qb[my_iindx[k]-l] * qloop * scale[ln1+ln2];
            }
      } /* end of kl double loop */
    }
  } /* end of pq double loop */

  /* 3. exterior multiloop contribution  */
  for(k=TURN+2; k<n-2*TURN-3; k++)
    qmo += qm[my_iindx[1]-k] * qm2[k+1] * pow(expMLclosing,n_seq);

  /* add additional pf of 1.0 to take open chain into account */
  qo = qho + qio + qmo + 1.0*scale[n];

  matrices->qo    = qo;
  matrices->qho   = qho;
  matrices->qio   = qio;
  matrices->qmo   = qmo;

  free(type);
}


PUBLIC char *
vrna_ali_pbacktrack(vrna_fold_compound *vc,
                    double *prob){

  double r, gr, qt;
  int k,i,j, start,s;
  double probs=1;
  char  *pstruc = NULL;

  int               n_seq       = vc->n_seq;
  int               n           = vc->length;
  short             **S         = vc->S;                                                                   
  short             **S5        = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/            
  short             **S3        = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/            
  pf_paramT         *pf_params  = vc->exp_params;
  pf_matricesT      *matrices   = vc->exp_matrices;
  model_detailsT    *md         = &(pf_params->model_details);
  int               *my_iindx   = vc->iindx;
  hard_constraintT  *hc         = vc->hc;
  soft_constraintT  **sc        = vc->scs;
  FLT_OR_DBL        *q          = matrices->q;
  FLT_OR_DBL        *qb         = matrices->qb;

  FLT_OR_DBL *q1k   = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+1));
  FLT_OR_DBL *qln   = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));

  FLT_OR_DBL        *scale        = matrices->scale;

  pstruc = space((n+1)*sizeof(char));

  for (i=0; i<n; i++)
    pstruc[i] = '.';

  for (k=1; k<=n; k++) {
    q1k[k] = q[my_iindx[1] - k];
    qln[k] = q[my_iindx[k] - n];
  }
  q1k[0] = 1.0;
  qln[n+1] = 1.0;

  start = 1;
  while (start<n) {
  /* find i position of first pair */
    probs=1.;
    for (i=start; i<n; i++) {
      gr = urn() * qln[i];
      if (gr > qln[i+1]*scale[1]) {
        *prob=*prob*probs*(1-qln[i+1]*scale[1]/qln[i]);
        break; /* i is paired */
      }
      probs*=qln[i+1]*scale[1]/qln[i];
    }
    if (i>=n) {
      *prob=*prob*probs;
      break; /* no more pairs */
    }
    /* now find the pairing partner j */
    r = urn() * (qln[i] - qln[i+1]*scale[1]);
    for (qt=0, j=i+1; j<=n; j++) {
      int xtype;
      /*  type = ptype[my_iindx[i]-j];
          if (type) {*/
      double qkl;
      if (qb[my_iindx[i]-j]>0) {
        qkl = qb[my_iindx[i]-j]*qln[j+1];  /*if psc too small qb=0!*/
        for (s=0; s< n_seq; s++) {
          xtype=md->pair[S[s][i]][S[s][j]];
          if (xtype==0) xtype=7;
          qkl *= exp_E_ExtLoop(xtype, (i>1) ? S5[s][i] : -1, (j<n) ? S3[s][j] : -1, pf_params);
        }
        qt += qkl; /*?*exp(pscore[jindx[j]+i]/kTn)*/
        if (qt > r) {
          *prob=*prob*(qkl/(qln[i] - qln[i+1]*scale[1]));/*probs*=qkl;*/
          break; /* j is paired */
        }
      }
    }
    if (j==n+1) nrerror("backtracking failed in ext loop");
    start = j+1;
    backtrack(vc, pstruc, i, j, prob); /*?*/
  }

  return pstruc;
}


PRIVATE void
backtrack(vrna_fold_compound *vc,
          char *pstruc,
          int i,
          int j,
          double *prob){

  int               n_seq       = vc->n_seq;
  short             **S         = vc->S;                                                                   
  short             **S5        = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/            
  short             **S3        = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/            
  char              **Ss        = vc->Ss;                                                                   
  unsigned short    **a2s       = vc->a2s;                                                                   
  pf_paramT         *pf_params  = vc->exp_params;
  pf_matricesT      *matrices   = vc->exp_matrices;
  model_detailsT    *md         = &(pf_params->model_details);
  int               *my_iindx   = vc->iindx;
  int               *jindx      = vc->jindx;
  hard_constraintT  *hc         = vc->hc;
  soft_constraintT  **sc        = vc->scs;
  FLT_OR_DBL        *qb         = matrices->qb;
  FLT_OR_DBL        *qm         = matrices->qm;
  FLT_OR_DBL        *qm1        = matrices->qm1;
  int               *pscore     = vc->pscore;     /* precomputed array of pair types */             

  FLT_OR_DBL        *scale        = matrices->scale;
  FLT_OR_DBL        *expMLbase    = matrices->expMLbase;

  /*backtrack given i,j basepair!*/
  double kTn = pf_params->kT/10.;
  int *type = (int *)space(sizeof(int) * n_seq);

  do {
    double r, qbt1, max_k, min_l;
    int k, l, u, u1,s;
    pstruc[i-1] = '('; pstruc[j-1] = ')';
    for (s=0; s<n_seq; s++) {
      type[s] = md->pair[S[s][i]][S[s][j]];
      if (type[s]==0) type[s]=7;
    }
    r = urn() * (qb[my_iindx[i]-j]/exp(pscore[jindx[j]+i]/kTn)); /*?*exp(pscore[jindx[j]+i]/kTn)*/

    qbt1=1.;
    for (s=0; s<n_seq; s++){
      u = a2s[s][j-1]-a2s[s][i];
      if (a2s[s][i]<1) continue;
      char loopseq[10];
      if(u < 7){
        strncpy(loopseq, Ss[s]+a2s[s][i]-1, 10);
      }
      qbt1 *= exp_E_Hairpin(u, type[s], S3[s][i], S5[s][j], loopseq, pf_params);
    }
    qbt1 *= scale[j-i+1];

    if (qbt1>r) {
      *prob=*prob*qbt1/(qb[my_iindx[i]-j]/exp(pscore[jindx[j]+i]/kTn));/*probs*=qbt1;*/
      free(type);
      return; /* found the hairpin we're done */
    }

    
    max_k = MIN2(i+MAXLOOP+1,j-TURN-2);
    l = MAX2(i+TURN+2,j-MAXLOOP-1);
    for (k=i+1; k<=max_k; k++){
      min_l = MAX2(k+TURN+1,j-1-MAXLOOP+k-i-1);

      for (l=min_l; l<j; l++){
        double qloop=1;
        int type_2;
        if (qb[my_iindx[k]-l]==0) {qloop=0; continue;}
        for (s=0; s<n_seq; s++) {
          u1 = a2s[s][k-1]-a2s[s][i]/*??*/;
          type_2 = md->pair[S[s][l]][S[s][k]]; if (type_2 == 0) type_2 = 7;
          qloop *= exp_E_IntLoop(u1, a2s[s][j-1]-a2s[s][l], type[s], type_2,
                                 S3[s][i], S5[s][j],S5[s][k], S3[s][l], pf_params);
          if(sc)
            if((i + 1 == k) && (j - 1 == l))
              if(sc[s])
                if(sc[s]->exp_en_stack)
                  qloop *=    sc[s]->exp_en_stack[i]
                            * sc[s]->exp_en_stack[k]
                            * sc[s]->exp_en_stack[l]
                            * sc[s]->exp_en_stack[j];

        }
        qbt1 += qb[my_iindx[k]-l] * qloop * scale[k-i+j-l];

        if (qbt1 > r) {
         *prob=*prob*qb[my_iindx[k]-l] * qloop * scale[k-i+j-l]/(qb[my_iindx[i]-j]/exp(pscore[jindx[j]+i]/kTn));
         /*
          prob*=qb[my_iindx[k]-l] * qloop * scale[k-i+j-l];
         */
          break;
        }
      }
      if (qbt1 > r) break;
    }
    if (l<j) {
      i=k; j=l;
    }
    else {
       *prob=*prob*(1-qbt1/(qb[my_iindx[i]-j]/exp(pscore[jindx[j]+i]/kTn)));
      break;
    }
  } while (1);

  /* backtrack in multi-loop */
  {
    double r, qt;
    int k, ii, jj;
    double qttemp=0;;
    i++; j--;
    /* find the first split index */
    ii = my_iindx[i]; /* ii-j=[i,j] */
    jj = jindx[j]; /* jj+i=[j,i] */
    for (qt=0., k=i+1; k<j; k++) qttemp += qm[ii-(k-1)]*qm1[jj+k];
    r = urn() * qttemp;
    for (qt=0., k=i+1; k<j; k++) {
      qt += qm[ii-(k-1)]*qm1[jj+k];
      if (qt>=r){
        *prob=*prob*qm[ii-(k-1)]*qm1[jj+k]/qttemp;/*qttemp;*/
        /*        prob*=qm[ii-(k-1)]*qm1[jj+k];*/
        break;
      }
    }
    if (k>=j) nrerror("backtrack failed, can't find split index ");

    backtrack_qm1(vc, pstruc, k, j, prob);

    j = k-1;
    while (j>i) {
      /* now backtrack  [i ... j] in qm[] */
      jj = jindx[j];/*habides??*/
      ii = my_iindx[i];
      r = urn() * qm[ii - j];
      qt = qm1[jj+i]; k=i;
      if (qt<r)
        for (k=i+1; k<=j; k++) {
          qt += (qm[ii-(k-1)]+expMLbase[k-i]/*n_seq??*/)*qm1[jj+k];
          if (qt >= r) {
            *prob=*prob*(qm[ii-(k-1)]+expMLbase[k-i])*qm1[jj+k]/qm[ii - j];/*???*/
            /*            probs*=qt;*/
            break;
          }
        }
      else {
        *prob=*prob*qt/qm[ii - j];/*??*/
      }
      if (k>j) nrerror("backtrack failed in qm");

      backtrack_qm1(vc, pstruc, k, j, prob);

      if (k<i+TURN) break; /* no more pairs */
      r = urn() * (qm[ii-(k-1)] + expMLbase[k-i]);
      if (expMLbase[k-i] >= r) {
        break; /* no more pairs */
        *prob=*prob*expMLbase[k-i]/(qm[ii-(k-1)] + expMLbase[k-i]);
      }
      j = k-1;
      /* whatishere?? */
    }
  }
  free(type);
}

PRIVATE void
backtrack_qm1(vrna_fold_compound *vc,
              char *pstruc,
              int i,
              int j,
              double *prob){

  int               n_seq       = vc->n_seq;
  short             **S         = vc->S;                                                                   
  short             **S5        = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/            
  short             **S3        = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/            
  pf_paramT         *pf_params  = vc->exp_params;
  pf_matricesT      *matrices   = vc->exp_matrices;
  model_detailsT    *md         = &(pf_params->model_details);
  int               *my_iindx   = vc->iindx;
  int               *jindx      = vc->jindx;
  hard_constraintT  *hc         = vc->hc;
  soft_constraintT  **sc        = vc->scs;
  FLT_OR_DBL        *qb         = matrices->qb;
  FLT_OR_DBL        *qm1        = matrices->qm1;
  FLT_OR_DBL        *expMLbase    = matrices->expMLbase;

  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  int ii, l, xtype,s;
  double qt, r, tempz;
  r = urn() * qm1[jindx[j]+i];
  ii = my_iindx[i];
  for (qt=0., l=i+TURN+1; l<=j; l++) {
    if (qb[ii-l]==0) continue;
    tempz=1.;
    for (s=0; s<n_seq; s++) {
      xtype = md->pair[S[s][i]][S[s][l]];
      if (xtype==0) xtype=7;
      tempz*=exp_E_MLstem(xtype, S5[s][i], S3[s][l], pf_params);
    }
    qt +=  qb[ii-l]*tempz*expMLbase[j-l];
    if (qt>=r) {
      *prob=*prob*qb[ii-l]*tempz*expMLbase[j-l]/qm1[jindx[j]+i];
      /* probs*=qb[ii-l]*tempz*expMLbase[j-l];*/
      break;
    }
  }
  if (l>j) nrerror("backtrack failed in qm1");

  backtrack(vc, pstruc, i, l, prob);
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
                pf_paramT *parameters,
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

  return vrna_ali_pbacktrack(backward_compat_compound, prob);
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
                  int **pscore_p) {

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
        *pscore_p = backward_compat_compound->pscore;
        return 1;
      }
  }
  return 0;
}

PUBLIC void
free_alipf_arrays(void){

  if(backward_compat_compound && backward_compat){
    destroy_fold_compound(backward_compat_compound);
    backward_compat_compound = NULL;
    iindx = NULL;
  }
}
