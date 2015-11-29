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
PRIVATE int                 backward_compat           = 0;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE void      alipf_linear(vrna_fold_compound_t *vc, char *structure);
PRIVATE void      alipf_create_bppm(vrna_fold_compound_t *vc, char *structure, plist **pl);
PRIVATE float     wrap_alipf_fold(const char **sequences,
                                  char *structure,
                                  plist **pl,
                                  vrna_exp_param_t *parameters,
                                  int calculate_bppm,
                                  int is_constrained,
                                  int is_circular);
PRIVATE void      wrap_alipf_circ(vrna_fold_compound_t *vc,
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
                vrna_exp_param_t *parameters,
                int calculate_bppm,
                int is_constrained,
                int is_circular){

  int                 n_seq;
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

  return vrna_pf_comparative(vc, structure, pl);
}

PUBLIC float
vrna_pf_comparative( vrna_fold_compound_t *vc,
                  char *structure,
                  plist **pl){

  int         i;
  FLT_OR_DBL  Q;
  float       free_energy;

  int               n         = vc->length;
  int               n_seq     = vc->n_seq;
  vrna_exp_param_t  *params   = vc->exp_params;
  vrna_md_t         *md       = &(params->model_details);
  vrna_mx_pf_t      *matrices = vc->exp_matrices;

  /* call user-defined recursion status callback function */
  if(vc->stat_cb)
    vc->stat_cb(vc, VRNA_STATUS_PF_PRE);

  alipf_linear(vc, structure);

  /* calculate post processing step for circular  */
  /* RNAs                                         */
  if(md->circ)
    wrap_alipf_circ(vc, structure);

  /* call user-defined recursion status callback function */
  if(vc->stat_cb)
    vc->stat_cb(vc, VRNA_STATUS_PF_POST);

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
alipf_linear( vrna_fold_compound_t *vc,
              char *structure){

  int         s, i,j,k,l, ij, jij, u, u1, u2, d, ii, *type, type_2, tt;
  FLT_OR_DBL  temp, temp2;
  FLT_OR_DBL  qbt1, *tmp;
  FLT_OR_DBL  *qqm = NULL, *qqm1 = NULL, *qq = NULL, *qq1 = NULL;
  double      kTn;

  int               n_seq             = vc->n_seq;
  int               n                 = vc->length;
  short             **S               = vc->S;                                                               
  short             **S5              = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/        
  short             **S3              = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/        
  char              **Ss              = vc->Ss;                                                               
  unsigned short    **a2s             = vc->a2s;                                                               
  vrna_exp_param_t  *pf_params        = vc->exp_params;
  vrna_mx_pf_t      *matrices         = vc->exp_matrices;
  vrna_md_t         *md               = &(pf_params->model_details);
  vrna_hc_t         *hc               = vc->hc;
  vrna_sc_t         **sc              = vc->scs;
  int               *my_iindx         = vc->iindx;
  int               *jindx            = vc->jindx;
  FLT_OR_DBL        *q                = matrices->q;
  FLT_OR_DBL        *qb               = matrices->qb;
  FLT_OR_DBL        *qm               = matrices->qm;
  FLT_OR_DBL        *qm1              = matrices->qm1;
  int               *pscore           = vc->pscore;     /* precomputed array of pair types */                  
  int               *rtype            = &(md->rtype[0]);
  int               circular          = md->circ;
  FLT_OR_DBL        *scale            = matrices->scale;
  FLT_OR_DBL        *expMLbase        = matrices->expMLbase;
  FLT_OR_DBL        expMLclosing      = pf_params->expMLclosing;
  char              *hard_constraints = hc->matrix;

  kTn   = pf_params->kT/10.;   /* kT in cal/mol  */
  type  = (int *)vrna_alloc(sizeof(int) * n_seq);

  int max_bpspan = (md->max_bp_span > 0) ? md->max_bp_span : n;

  /* allocate memory for helper arrays */
  qq        = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
  qq1       = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
  qqm       = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
  qqm1      = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));


  /* array initialization ; qb,qm,q
     qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  for (d=0; d<=TURN; d++)
    for (i=1; i<=n-d; i++) {
      j=i+d;
      ij = my_iindx[i]-j;
      if(hc->up_ext[i] > d){
        q[ij]=1.0*scale[d+1];
        if(sc)
          for(s = 0; s < n_seq; s++)
            if(sc[s]){
              int u = d + 1 /* a2s[s][j] - a2s[s][i] + 1 */;
              if(sc[s]->exp_energy_up)
                q[ij] *= sc[s]->exp_energy_up[a2s[s][i]][u];
            }
      }
      qb[ij]=qm[ij]=0.0;
    }

  for (i=1; i<=n; i++)
    qq[i]=qq1[i]=qqm[i]=qqm1[i]=0;

  for (j=TURN+2;j<=n; j++) {
    for (i=j-TURN-1; i>=1; i--) {
      int psc;
      /* construction of partition function for segment i,j */
      /* calculate pf given that i and j pair: qb(i,j)      */
      ij  = my_iindx[i] - j;
      jij = jindx[j] + i;

      for (s=0; s<n_seq; s++) {
        type[s] = md->pair[S[s][i]][S[s][j]];
        if (type[s]==0) type[s]=7;
      }

      psc = pscore[jij];
      qbt1 = 0.;

      if(hard_constraints[jij]){
        /* process hairpin loop(s) */
        qbt1 += vrna_exp_E_hp_loop(vc, i, j);

        /* interior loops with interior pair k,l */
        if(hard_constraints[jij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
          for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++){
            if(hc->up_int[i+1] < k - i - 1)
              break;

            for (l=MAX2(k+TURN+1,j-1-MAXLOOP+k-i-1); l<=j-1; l++){
              double qloop=1.;

              if(!(hard_constraints[jindx[l] + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC))
                continue;
              if(hc->up_int[l+1] < j - l - 1)
                continue;
              if(qb[my_iindx[k]-l] == 0){
                qloop=0;
                continue;
              }

              for (s=0; s<n_seq; s++) {
                u1 = a2s[s][k-1] - a2s[s][i];
                u2 = a2s[s][j-1] - a2s[s][l];
                type_2 = md->pair[S[s][l]][S[s][k]];
                if (type_2 == 0) type_2 = 7;
                qloop *= exp_E_IntLoop(u1, u2,
                                    type[s], type_2, S3[s][i],
                                    S5[s][j], S5[s][k], S3[s][l],
                                    pf_params
                                  );
              }
              if(sc){
                for(s = 0; s < n_seq; s++){
                  if(sc[s]){
                    u1 = a2s[s][k-1] - a2s[s][i];
                    u2 = a2s[s][j-1] - a2s[s][l];
/*
                    u1 = k - i - 1;
                    u2 = j - l - 1;
*/
                    if(sc[s]->exp_energy_bp)
                      qloop *=    sc[s]->exp_energy_bp[jij];

                    if(sc[s]->exp_energy_up)
                      qloop *=    sc[s]->exp_energy_up[a2s[s][i]+1][u1]
                                * sc[s]->exp_energy_up[a2s[s][l]+1][u2];

                    if(sc[s]->exp_energy_stack)
                      if(u1 + u2 == 0){
                        qloop *=    sc[s]->exp_energy_stack[i]
                                  * sc[s]->exp_energy_stack[k]
                                  * sc[s]->exp_energy_stack[l]
                                  * sc[s]->exp_energy_stack[j];
                      }
                  }
                }
              }

              qbt1 += qb[my_iindx[k]-l] * qloop * scale[k-i+j-l];
            }
          }
        }

        if(hard_constraints[jij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
          /* multi-loop loop contribution */
          ii = my_iindx[i+1]; /* ii-k=[i+1,k-1] */
          temp = 0.0;
          for (k=i+2; k<=j-1; k++)
            temp += qm[ii-(k-1)] * qqm1[k];
          for (s=0; s<n_seq; s++) {
            tt = rtype[type[s]];
            temp *= exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params) * expMLclosing;
          }
          if(sc)
            for(s = 0; s < n_seq; s++){
              if(sc[s]){
                if(sc[s]->exp_energy_bp)
                  temp *= sc[s]->exp_energy_bp[jij];
              }
            }
          temp *= scale[2] ;
          qbt1 += temp;
        }
        qb[ij] = qbt1;
        qb[ij] *= exp(psc/kTn);
      } /* end if (type!=0) */
      else qb[ij] = 0.0;

      /* construction of qqm matrix containing final stem
         contributions to multiple loop partition function
         from segment i,j */
      qqm[i] = 0.;
      if(hc->up_ml[i]){
        temp = qqm1[i] * expMLbase[1]; /* expMLbase[1]^n_seq */
        if(sc)
          for(s = 0; s < n_seq; s++){
            if(sc[s]){
              if(sc[s]->exp_energy_up)
                temp *= sc[s]->exp_energy_up[a2s[s][i]][1];
            }
          }
        qqm[i] += temp;
      }
      if(hard_constraints[jij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
        for (qbt1=1, s=0; s<n_seq; s++) {
          qbt1 *= exp_E_MLstem(type[s], (i>1) || circular ? S5[s][i] : -1, (j<n) || circular ? S3[s][j] : -1, pf_params);
        }
        qqm[i] += qb[ij]*qbt1;
      }

      if(qm1)
        qm1[jij] = qqm[i]; /* for circ folding and stochBT */

      /* construction of qm matrix containing multiple loop
         partition function contributions from segment i,j */
      temp = 0.0;
      ii = my_iindx[i];  /* ii-k=[i,k-1] */
      for (k=i+1; k<=j; k++)
        temp += qm[ii-(k-1)] * qqm[k];

      for (k=i+1; k<=j; k++){
        if(hc->up_ml[i] < k - i)
          break;
        temp2 = expMLbase[k-i] * qqm[k];
        if(sc)
          for(s = 0; s < n_seq; s++){
            if(sc[s]){
              if(sc[s]->exp_energy_up)
                temp2 *= sc[s]->exp_energy_up[a2s[s][i]][a2s[s][k] - a2s[s][i]];
            }
          }
        temp += temp2;
      }
      qm[ij] = (temp + qqm[i]);

      /* auxiliary matrix qq for cubic order q calculation below */
      qbt1 = 0.;
      if((qb[ij] > 0) && (hard_constraints[jij] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)){
        qbt1 = qb[ij];
        for (s=0; s<n_seq; s++) {
          qbt1 *= exp_E_ExtLoop(type[s], (i>1) || circular ? S5[s][i] : -1, (j<n) || circular ? S3[s][j] : -1, pf_params);
        }
      }
      if(hc->up_ext[i]){
        temp = qq1[i]*scale[1];
        if(sc)
          for(s = 0; s < n_seq; s++){
            if(sc[s]){
              if(sc[s]->exp_energy_up)
                temp *= sc[s]->exp_energy_up[a2s[s][i]][1];
            }
          }
        qbt1 += temp;
      }
      qq[i] = qbt1;

      /* construction of partition function for segment i,j */
      temp = qq[i];
      if(hc->up_ext[i] >= j - i + 1){
        temp2 = 1.0 * scale[j - i + 1];
        if(sc)
          for(s = 0; s < n_seq; s++){
            if(sc[s]){
              if(sc[s]->exp_energy_up)
                temp2 *= sc[s]->exp_energy_up[a2s[s][i]][a2s[s][j] - a2s[s][i] + 1];
            }
          }
        temp += temp2;
      }

      for(k = i; k <= j - 1; k++)
        temp += q[ii - k] * qq[k + 1];

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
        vrna_message_error(msg);
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
alipf_create_bppm(vrna_fold_compound_t *vc,
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
  char              **Ss          = vc->Ss;
  unsigned short    **a2s         = vc->a2s;                                                                   
  vrna_exp_param_t  *pf_params    = vc->exp_params;
  vrna_mx_pf_t      *matrices     = vc->exp_matrices;
  vrna_md_t         *md           = &(pf_params->model_details);
  vrna_hc_t         *hc           = vc->hc;
  vrna_sc_t         **sc          = vc->scs;
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
  char              *hard_constraints = hc->matrix;

  double kTn, pp;

  FLT_OR_DBL *prm_l   = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
  FLT_OR_DBL *prm_l1  = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
  FLT_OR_DBL *prml    = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
  type                = (int *)vrna_alloc(sizeof(int) * n_seq);

  if((matrices->q1k == NULL) || (matrices->qln == NULL)){
    free(matrices->q1k);
    matrices->q1k = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+1));
    free(matrices->qln);
    matrices->qln = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
  }

  FLT_OR_DBL *q1k    = matrices->q1k;
  FLT_OR_DBL *qln    = matrices->qln;

  for (k=1; k<=n; k++) {
    q1k[k] = q[my_iindx[1] - k];
    qln[k] = q[my_iindx[k] - n];
  }
  q1k[0] = 1.0;
  qln[n+1] = 1.0;


  kTn = pf_params->kT/10.;   /* kT in cal/mol  */

  for (i=0; i<=n; i++)
    prm_l[i]=prm_l1[i]=prml[i]=0;

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

          tmp2 = 0.;

          /* 1.1. Exterior Hairpin Contribution */
          if(hard_constraints[jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP){
            int u = i + n - j -1;
            if(hc->up_hp[j + 1] >= u){
              for (qbt1=1.,s=0; s<n_seq; s++) {
                int u1 = a2s[s][i] - 1 + a2s[s][n] - a2s[s][j];

                char loopseq[10];
                if (u1<9){
                  strcpy(loopseq , Ss[s] + a2s[s][j] - 1);
                  strncat(loopseq, Ss[s], a2s[s][i]);
                }
                qbt1 *= exp_E_Hairpin(u1, type[s], S3[s][j], S5[s][i], loopseq, pf_params);
              }
              if(sc)
                for(s = 0; s < n_seq; s++){
                  if(sc[s]){
                    if(sc[s]->exp_energy_up)
                      qbt1 *=   ((i > 1) ? sc[s]->exp_energy_up[a2s[s][j]+1][a2s[s][n] - a2s[s][j]] : 1.)
                              * ((j < n) ? sc[s]->exp_energy_up[a2s[s][1]][a2s[s][i] - a2s[s][1]] : 1.);
                   }
                }
            }
            tmp2 = qbt1 * scale[u];
          }

          /* 1.2. Exterior Interior Loop Contribution */
          /* recycling of k and l... */
          if(hard_constraints[jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){

            /* 1.2.1. first we calc exterior loop energy with constraint, that i,j  */
            /* delimtis the "right" part of the interior loop                       */
            /* (l,k) is "outer pair"                                                */
            for(k=1; k < i-TURN-1; k++){
              /* so first, lets calc the length of loop between j and k */
              int ln1, lstart;
              ln1 = k + n - j - 1;
              if(ln1>MAXLOOP)
                break;
              if(hc->up_int[j+1] < ln1)
                break;

              lstart = ln1+i-1-MAXLOOP;
              if(lstart<k+TURN+1) lstart = k + TURN + 1;
              for(l=lstart; l < i; l++){
                int ln2,ln2a,ln1a, type_2;
                ln2 = i - l - 1;
                if(ln1+ln2>MAXLOOP)
                  continue;
                if(hc->up_int[l+1] < ln2)
                  continue;
                if(!(hard_constraints[jindx[l] + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP))
                  continue;
                
                double qloop=1.;
                if(qb[my_iindx[k]-l]==0.){
                  qloop=0.;
                  continue;
                }

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
                if(sc)
                  for(s = 0; s < n_seq; s++){
                    if(sc[s]){
                      ln2a= a2s[s][i-1];
                      ln2a-=a2s[s][l];
                      ln1a= a2s[s][n]-a2s[s][j];
                      ln1a+=a2s[s][k-1];

                      if(sc[s]->exp_energy_up)
                        qloop *=    sc[s]->exp_energy_up[a2s[s][l]+1][ln2a]
                                  * ((j < n) ? sc[s]->exp_energy_up[a2s[s][j]+1][a2s[s][n] - a2s[s][j]] : 1.)
                                  * ((k > 1) ? sc[s]->exp_energy_up[1][a2s[s][k]-1] : 1.);

                      if((ln1a + ln2a == 0) && sc[s]->exp_energy_stack)
                        qloop *=    sc[s]->exp_energy_stack[a2s[s][k]]
                                  * sc[s]->exp_energy_stack[a2s[s][l]]
                                  * sc[s]->exp_energy_stack[a2s[s][i]]
                                  * sc[s]->exp_energy_stack[a2s[s][j]];
                    }
                  }
                tmp2 += qb[my_iindx[k] - l] * qloop * scale[ln1+ln2];
              }
            }

            /* 1.2.2. second we calc exterior loop energy with constraint, that i,j */
            /* delimtis the "left" part of the interior loop                        */
            /* (j,i) is "outer pair"                                                */
            for(k=j+1; k < n-TURN; k++){
              /* so first, lets calc the length of loop between l and i */
              int ln1, lstart;
              ln1 = k - j - 1;
              if((ln1 + i - 1)>MAXLOOP)
                break;
              if(hc->up_int[j+1] < ln1)
                break;

              lstart = ln1+i-1+n-MAXLOOP;
              if(lstart<k+TURN+1) lstart = k + TURN + 1;
              for(l=lstart; l <= n; l++){
                int ln2, type_2;
                ln2 = i - 1 + n - l;
                if(ln1+ln2>MAXLOOP)
                  continue;
                if(hc->up_int[l+1] < ln2)
                  continue;
                if(!(hard_constraints[jindx[l] + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP))
                  continue;

                double qloop=1.;
                if(qb[my_iindx[k]-l]==0.){
                  qloop=0.;
                  continue;
                }

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
                if(sc)
                  for(s = 0; s < n_seq; s++){
                    if(sc[s]){
                      ln1 = a2s[s][k] - a2s[s][j+1];
                      ln2 = a2s[s][i-1] + a2s[s][n] - a2s[s][l];

                      if(sc[s]->exp_energy_up)
                        qloop *=    sc[s]->exp_energy_up[a2s[s][j]+1][ln1]
                                  * ((l < n) ? sc[s]->exp_energy_up[a2s[s][l]+1][a2s[s][n] - a2s[s][l]] : 1.)
                                  * ((i > 1) ? sc[s]->exp_energy_up[1][a2s[s][i]-1] : 1.);

                      if((ln1 + ln2 == 0) && sc[s]->exp_energy_stack)
                        qloop *=    sc[s]->exp_energy_stack[a2s[s][k]]
                                  * sc[s]->exp_energy_stack[a2s[s][l]]
                                  * sc[s]->exp_energy_stack[a2s[s][i]]
                                  * sc[s]->exp_energy_stack[a2s[s][j]];
                    }
                  }
                tmp2 += qb[my_iindx[k] - l] * qloop * scale[ln1+ln2];
              }
            }
          }
          /* 1.3 Exterior multiloop decomposition */
          if(hard_constraints[jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
            /* 1.3.1 Middle part                    */
            if((i>TURN+2) && (j<n-TURN-1)){

              for (tmp3=1, s=0; s<n_seq; s++){
                tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);
              }
              tmp2 += qm[my_iindx[1]-i+1] * qm[my_iindx[j+1]-n] * tmp3 * pow(expMLclosing,n_seq);
            }
            /* 1.3.2 Left part    */
            for(k=TURN+2; k < i-TURN-2; k++){
              if(hc->up_ml[j+1] < n-j)
                break;

              for (tmp3=1, s=0; s<n_seq; s++){
                tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);
              }

              if(sc)
                for(s = 0; s < n_seq; s++){
                  if(sc[s]){
                    if(sc[s]->exp_energy_bp)
                      tmp3 *= sc[s]->exp_energy_bp[jindx[j] + i];

                    if(sc[s]->exp_energy_up)
                      tmp3 *= sc[s]->exp_energy_up[a2s[s][j]+1][a2s[s][n]-a2s[s][j]];
                  }
                }

              tmp2 += qm[my_iindx[1]-k] * qm1[jindx[i-1]+k+1] * tmp3 * expMLbase[n-j] * pow(expMLclosing,n_seq);
            }
            /* 1.3.3 Right part    */
            for(k=j+TURN+2; k < n-TURN-1;k++){
              if(hc->up_ml[1] < i-1)
                break;

              for (tmp3=1, s=0; s<n_seq; s++){
                tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);
              }

              if(sc)
                for(s = 0; s < n_seq; s++){
                  if(sc[s]){
                    if(sc[s]->exp_energy_bp)
                      tmp3 *= sc[s]->exp_energy_bp[jindx[j] + i];

                    if(sc[s]->exp_energy_up)
                      tmp3 *= sc[s]->exp_energy_up[a2s[s][1]][a2s[s][i]-a2s[s][1]];
                  }
                }

              tmp2 += qm[my_iindx[j+1]-k] * qm1[jindx[n]+k+1] * tmp3 * expMLbase[i-1] * pow(expMLclosing,n_seq);
            }
          }
          probs[ij] *= tmp2;
        }
        else probs[ij] = 0;
      }  /* end for j=..*/
    }  /* end or i=...  */
  } /* end if(circular)  */
  else{
    for (i=1; i<=n; i++) {
      for (j=i; j<=MIN2(i+TURN,n); j++)
        probs[my_iindx[i]-j] = 0;

      for (j=i+TURN+1; j<=n; j++) {
        ij = my_iindx[i]-j;
        if ((qb[ij] > 0.) && (hard_constraints[jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)){
          probs[ij] = q1k[i-1] * qln[j+1]/q1k[n] * exp(pscore[jindx[j]+i]/kTn);
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
      pp = 0.;
      kl = my_iindx[k]-l;
      if (qb[kl] == 0.) continue;
      if(!(hard_constraints[jindx[l] + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) continue;

      for (s=0; s<n_seq; s++) {
        type[s] = md->pair[S[s][l]][S[s][k]];
        if (type[s]==0) type[s]=7;
      }

      for (i=MAX2(1,k-MAXLOOP-1); i<=k-1; i++){
        if(hc->up_int[i+1] < k - i - 1)
          continue;

        for (j=l+1; j<=MIN2(l+ MAXLOOP -k+i+2,n); j++) {
          double qloop=1;
          ij = my_iindx[i] - j;

          if(probs[ij] == 0.) continue;
          if(!(hard_constraints[jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)) continue;
          if(hc->up_int[l+1] < j - l - 1) break;

          for (s=0; s<n_seq; s++) {
            int typ, u1, u2;
            u1 = a2s[s][k-1] - a2s[s][i];
            u2 = a2s[s][j-1] - a2s[s][l];
            typ = md->pair[S[s][i]][S[s][j]]; if (typ==0) typ=7;
            qloop *=  exp_E_IntLoop(u1, u2, typ, type[s], S3[s][i], S5[s][j], S5[s][k], S3[s][l], pf_params);
          }

          if(sc){
            for(s = 0; s < n_seq; s++){
              if(sc[s]){
                int u1, u2;
                u1 = a2s[s][k-1] - a2s[s][i];
                u2 = a2s[s][j-1] - a2s[s][l];
/*
                u1 = k - i - 1;
                u2 = j - l - 1;
*/
                if(sc[s]->exp_energy_bp)
                  qloop *= sc[s]->exp_energy_bp[jindx[j] + i];

                if(sc[s]->exp_energy_up)
                  qloop *=    sc[s]->exp_energy_up[a2s[s][i]+1][u1]
                              * sc[s]->exp_energy_up[a2s[s][l]+1][u2];

                if(sc[s]->exp_energy_stack)
                  if(u1 + u2 == 0)
                    qloop *=    sc[s]->exp_energy_stack[i]
                              * sc[s]->exp_energy_stack[k]
                              * sc[s]->exp_energy_stack[l]
                              * sc[s]->exp_energy_stack[j];

              }
            }
          }
          pp += probs[ij]*qloop*scale[k-i + j-l];
        }
      }
      probs[kl] += pp * exp(pscore[jindx[l]+k]/kTn);
    }
    /* 3. bonding k,l as substem of multi-loop enclosed by i,j */
    prm_MLb = 0.;
    if (l<n)
      for (k=2; k<l-TURN; k++) {
      i = k-1;
      prmt = prmt1 = 0.;

      if(1 /* hard_constraints[jindx[l] + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC */){
        ii = my_iindx[i];     /* ii-j=[i,j]     */
        ll = my_iindx[l+1];   /* ll-j=[l+1,j-1] */
        if(hard_constraints[jindx[l+1] + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
          prmt1 = probs[ii-(l+1)];
          for (s=0; s<n_seq; s++) {
            tt = md->pair[S[s][l+1]][S[s][i]]; if (tt==0) tt=7;
            prmt1 *= exp_E_MLstem(tt, S5[s][l+1], S3[s][i], pf_params) * expMLclosing;
          }

          if(sc)
            for(s = 0; s < n_seq; s++){
              if(sc[s]){
                if(sc[s]->exp_energy_bp)
                  prmt1 *= sc[s]->exp_energy_bp[jindx[l+1] + i];
              }
            }
        }

        for (j=l+2; j<=n; j++){
          pp = 1.;
          if(probs[ii-j]==0) continue;
          if(!(hard_constraints[jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP)) continue;

          for (s=0; s<n_seq; s++) {
            tt = md->pair[S[s][j]][S[s][i]]; if (tt==0) tt=7;
            pp *=  exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params) * expMLclosing;
          }

          if(sc)
            for(s = 0; s < n_seq; s++){
              if(sc[s]){
                if(sc[s]->exp_energy_bp)
                  pp *= sc[s]->exp_energy_bp[jindx[j] + i];
              }
            }

          prmt +=  probs[ii-j] * pp * qm[ll-(j-1)];
        }
        kl = my_iindx[k]-l;

        prml[ i] = prmt;

        pp = 0.;
        if(hc->up_ml[l+1]){
          pp = prm_l1[i] * expMLbase[1];
          if(sc)
            for(s = 0; s < n_seq; s++){
              if(sc[s]){
                if(sc[s]->exp_energy_up)
                  pp *= sc[s]->exp_energy_up[a2s[s][l+1]][1];
              }
            }
        }
        prm_l[i] = pp + prmt1; /* expMLbase[1]^n_seq */

        pp = 0.;
        if(hc->up_ml[i]){
          pp = prm_MLb * expMLbase[1];
          if(sc)
            for(s = 0; s < n_seq; s++){
              if(sc[s]){
                if(sc[s]->exp_energy_up)
                  pp *= sc[s]->exp_energy_up[a2s[s][i]][1];
              }
            }
        }
        prm_MLb = pp + prml[i];

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
      } else { /* (k,l) not allowed to be substem of multiloop closed by (i,j) */
        prml[i] = prm_l[i] = prm_l1[i] = 0.;
      }

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
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);

  if (structure!=NULL){
    char *s = vrna_db_from_probs(probs, (unsigned int)n);
    memcpy(structure, s, n);
    structure[n] = '\0';
    free(s);
  }

  if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
        "you might try a smaller pf_scale than %g\n",
        ov, pf_params->pf_scale);

  free(type);
  free(prm_l);
  free(prm_l1);
  free(prml);
}

/* calculate partition function for circular case   */
/* NOTE: this is the postprocessing step ONLY        */
/* You have to call alipf_linear first to calculate  */
/* circular case!!!                                  */

PRIVATE void
wrap_alipf_circ(vrna_fold_compound_t *vc,
                char *structure){

  int u, p, q, pq, k, l, s, *type;
  FLT_OR_DBL qbt1, qot, qo, qho, qio, qmo;

  char              **sequences = vc->sequences;
  int               n_seq       = vc->n_seq;
  int               n           = vc->length;
  short             **S         = vc->S;                                                                   
  short             **S5        = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/            
  short             **S3        = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/            
  char              **Ss        = vc->Ss;                                                                   
  unsigned short    **a2s       = vc->a2s;                                                                   
  vrna_exp_param_t  *pf_params  = vc->exp_params;
  vrna_mx_pf_t      *matrices   = vc->exp_matrices;
  vrna_md_t         *md         = &(pf_params->model_details);
  int               *my_iindx   = vc->iindx;
  int               *jindx      = vc->jindx;
  vrna_hc_t         *hc         = vc->hc;
  vrna_sc_t         **sc        = vc->scs;
  FLT_OR_DBL        *qb         = matrices->qb;
  FLT_OR_DBL        *qm         = matrices->qm;
  FLT_OR_DBL        *qm1        = matrices->qm1;
  FLT_OR_DBL        *qm2        = matrices->qm2;
  int               *pscore     = vc->pscore;     /* precomputed array of pair types */             
  FLT_OR_DBL        *scale      = matrices->scale;
  FLT_OR_DBL        expMLclosing      = pf_params->expMLclosing;
  char              *hard_constraints = hc->matrix;
  int               *rtype            = &(md->rtype[0]);

  type  = (int *)vrna_alloc(sizeof(int) * n_seq);

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
      pq  = jindx[q] + p;
      psc = pscore[pq];

      if(!hard_constraints[pq]) continue;

      for(s = 0; s < n_seq; s++){
        type[s] = md->pair[S[s][p]][S[s][q]];
        if (type[s]==0) type[s]=7;
      }

      /* 1. exterior hairpin contribution  */
      /* Note, that we do not scale Hairpin Energy by u+2 but by u cause the scale  */
      /* for the closing pair was already done in the forward recursion              */
      if(hard_constraints[pq] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP){
        if(hc->up_hp[q+1] > u){
          for (qbt1=1,s=0; s<n_seq; s++) {
            int rt;
            char loopseq[10];
            u   = a2s[s][n] - a2s[s][q] + a2s[s][p] - 1;
            rt  = rtype[type[s]];

            if (u<9){
              strcpy(loopseq , Ss[s] + a2s[s][q] - 1);
              strncat(loopseq, Ss[s], a2s[s][p]);
            }
            qbt1 *= exp_E_Hairpin(u, rt, S3[s][q], S5[s][p], loopseq, pf_params);
          }
          if(sc)
            for(s = 0; s < n_seq; s++){
              if(sc[s]){
                if(sc[s]->exp_energy_up){
                  qbt1 *=   ((p > 1) ? sc[s]->exp_energy_up[1][a2s[s][p]-1] : 1.)
                          * ((q < n) ? sc[s]->exp_energy_up[a2s[s][q]+1][a2s[s][n] - a2s[s][q]] : 1.);
                }
              }
            }
          qho += qb[my_iindx[p]-q] * qbt1 * scale[u];
        }
      }
      /* 2. exterior interior loop contribution*/

      if(hard_constraints[pq] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
        for(k=q+1; k < n; k++){
          int ln1, lstart;
          ln1 = k - q - 1;
          if(ln1 + p - 1 > MAXLOOP)
            break;
          if(hc->up_int[q+1] < ln1)
            break;

          lstart = ln1+p-1+n-MAXLOOP;
          if(lstart<k+TURN+1) lstart = k + TURN + 1;
          for(l=lstart;l <= n; l++){
            int ln2, type_2;

            ln2 = (p - 1) + (n - l);

            if(!(hard_constraints[jindx[l]+k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP))
              continue;
            if((ln1+ln2) > MAXLOOP)
              continue;
            if(hc->up_int[l+1] < ln2)
              continue;

            double qloop=1.;
            if (qb[my_iindx[k]-l]==0.){ qloop=0.; continue;}

            for (s=0; s<n_seq; s++){
              int ln1a = a2s[s][k] - 1 - a2s[s][q];
              int ln2a = a2s[s][n] - a2s[s][l] + a2s[s][p] - 1;
              int rt = rtype[type[s]];
              type_2 = md->pair[S[s][l]][S[s][k]];
              if (type_2 == 0) type_2 = 7;
              qloop *= exp_E_IntLoop(ln1a, ln2a, rt, type_2, S3[s][q], S5[s][p], S5[s][k], S3[s][l], pf_params);
            }
            if(sc)
              for(s = 0; s < n_seq; s++){
                int ln1a = a2s[s][k] - 1 - a2s[s][q];
                int ln2a = a2s[s][n] - a2s[s][l] + a2s[s][p] - 1;
                if(sc[s]){
                  if((ln1a+ln2a == 0) && (sc[s]->exp_energy_stack))
                    qloop *=    sc[s]->exp_energy_stack[a2s[s][p]]
                              * sc[s]->exp_energy_stack[a2s[s][q]]
                              * sc[s]->exp_energy_stack[a2s[s][k]]
                              * sc[s]->exp_energy_stack[a2s[s][l]];

                  if(sc[s]->exp_energy_up)
                    qloop *=    sc[s]->exp_energy_up[a2s[s][q] + 1][ln1a]
                              * ((l < n) ? sc[s]->exp_energy_up[a2s[s][l]+1][a2s[s][n] - a2s[s][l]] : 1.)
                              * ((p > 1) ? sc[s]->exp_energy_up[1][a2s[s][p]-1] : 1.);
                }
              }

            qio += qb[my_iindx[p]-q] * qb[my_iindx[k]-l] * qloop * scale[ln1+ln2];
          }
        } /* end of kl double loop */
      }
    }
  } /* end of pq double loop */

  /* 3. exterior multiloop contribution  */
  for(k=TURN+2; k<n-2*TURN-3; k++)
    qmo += qm[my_iindx[1]-k] * qm2[k+1] * pow(expMLclosing,n_seq);

  /* add additional pf of 1.0 to take open chain into account */
  qo = qho + qio + qmo;
  if(hc->up_ext[1] >= n)
     qo += 1.0 * scale[n];

  matrices->qo    = qo;
  matrices->qho   = qho;
  matrices->qio   = qio;
  matrices->qmo   = qmo;

  free(type);
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

  return vrna_pbacktrack_comparative(backward_compat_compound, prob);
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
    vrna_fold_compound_free(backward_compat_compound);
    backward_compat_compound  = NULL;
    backward_compat           = 0;
    iindx                     = NULL;
  }
}
