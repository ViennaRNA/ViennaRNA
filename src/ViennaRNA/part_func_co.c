/* Last changed Time-stamp: <2007-05-09 16:11:21 ivo> */
/*
                  partiton function for RNA secondary structures

                  Ivo L Hofacker
                  Stephan Bernhart
                  Vienna RNA package
*/
/*
  $Log: part_func_co.c,v $
  Revision 1.10  2007/05/10 17:27:01  ivo
  make sure the relative error eps is positive in newton iteration

  Revision 1.9  2006/05/10 15:12:11  ivo
  some compiler choked on  double semicolon after declaration

  Revision 1.8  2006/04/05 12:52:31  ivo
  Fix performance bug (O(n^4) loop)

  Revision 1.7  2006/01/19 11:30:04  ivo
  compute_probabilities should only look at one dimer at a time

  Revision 1.6  2006/01/18 12:55:40  ivo
  major cleanup of berni code
  fix bugs related to confusing which free energy is returned by co_pf_fold()

  Revision 1.5  2006/01/16 11:32:25  ivo
  small bug in multiloop pair probs

  Revision 1.4  2006/01/05 18:13:40  ivo
  update

  Revision 1.3  2006/01/04 15:14:29  ivo
  fix bug in concentration calculations

  Revision 1.2  2004/12/23 12:14:41  berni
  *** empty log message ***

  Revision 1.1  2004/12/22 10:46:17  berni

  Partition function Cofolding 0.9, Computation of concentrations.

  Revision 1.16  2003/08/04 09:14:09  ivo
  finish up stochastic backtracking

  Revision 1.15  2002/03/19 16:51:12  ivo
  more on stochastic backtracking (still incomplete)

  Revision 1.13  2001/11/16 17:30:04  ivo
  add stochastic backtracking (still incomplete)
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include <limits.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/part_func_co.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define ISOLATED  256.0
#undef TURN
#define TURN 0
#define ON_SAME_STRAND(I,J,C)  (((I)>=(C))||((J)<(C)))

/* #define SAME_STRAND(I,J) (((J)<cut_point)||((I)>=cut_point2)||(((I)>=cut_point)&&((J)<cut_point2)))
 */

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
int     mirnatog      = 0;
double  F_monomer[2]  = {0,0}; /* free energies of the two monomers */

/*
#################################
# PRIVATE VARIABLES             #
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
PRIVATE void    pf_co(vrna_fold_compound *vc);
PRIVATE void    pf_co_bppm(vrna_fold_compound *vc, char *structure);
PRIVATE double  *Newton_Conc(double ZAB, double ZAA, double ZBB, double concA, double concB,double* ConcVec);
PRIVATE void    get_arrays(unsigned int length);
PRIVATE void    backtrack(vrna_fold_compound *vc, int i, int j, char *pstruc);
PRIVATE cofoldF wrap_co_pf_fold(char *sequence,
                                char *structure,
                                pf_paramT *parameters,
                                int calculate_bppm,
                                int is_constrained);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

/*-----------------------------------------------------------------*/
PRIVATE cofoldF
wrap_co_pf_fold(char *sequence,
                char *structure,
                pf_paramT *parameters,
                int calculate_bppm,
                int is_constrained){

  int                 length;
  char                *seq;
  vrna_fold_compound  *vc;
  model_detailsT      md;

  vc      = NULL;
  length  = strlen(sequence);

  /* we need pf_paramT datastructure to correctly init default hard constraints */
  if(parameters)
    md = parameters->model_details;
  else{
    set_model_details(&md); /* get global default parameters */
  }
  md.compute_bpp    = calculate_bppm;
  md.min_loop_size  = 0;

  seq = (char *)space(sizeof(char) * (length + 2));
  if(cut_point > -1){
    int i;
    for(i = 0; i < cut_point-1; i++)
      seq[i] = sequence[i];
    seq[i] = '&';
    for(;i<(int)length;i++)
      seq[i+1] = sequence[i];
  } else { /* this ensures the allocation of all cofold matrices via vrna_get_fold_compound */
    free(seq);
    seq = strdup(sequence);
  }

  vc = vrna_get_fold_compound(seq, &md, VRNA_OPTION_PF | VRNA_OPTION_HYBRID);

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

  if(backward_compat_compound)
    vrna_free_fold_compound(backward_compat_compound);

  backward_compat_compound = vc;
  backward_compat           = 1;
  iindx = backward_compat_compound->iindx;

  free(seq);
  return vrna_co_pf_fold(vc, structure);
}

PUBLIC cofoldF
vrna_co_pf_fold(vrna_fold_compound *vc,
                char *structure){
                
  int             n;
  FLT_OR_DBL      Q;
  cofoldF         X;
  double          free_energy;
  char            *sequence;
  model_detailsT  *md;
  pf_paramT       *params;
  pf_matricesT    *matrices;

  params    = vc->exp_params;
  n         = vc->length;
  md        = &(params->model_details);
  matrices  = vc->exp_matrices;
  sequence  = vc->sequence;

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

#ifdef SUN4
  nonstandard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(1);
#endif
#endif

  pf_co(vc);

  if (md->backtrack_type=='C')
    Q = matrices->qb[vc->iindx[1]-n];
  else if (md->backtrack_type=='M')
    Q = matrices->qm[vc->iindx[1]-n];
  else Q = matrices->q[vc->iindx[1]-n];

  /* ensemble free energy in Kcal/mol */
  if (Q<=FLT_MIN)
    fprintf(stderr, "pf_scale too large\n");
  free_energy = (-log(Q)-n*log(params->pf_scale))*params->kT/1000.0;
  /* in case we abort because of floating point errors */
  if (n>1600) fprintf(stderr, "free energy = %8.2f\n", free_energy);
  /*probability of molecules being bound together*/

  /*Computation of "real" Partition function*/
  /*Need that for concentrations*/
  if (vc->cutpoint > 0){
    double kT, QAB, QToT, Qzero;

    kT = params->kT/1000.0;
    Qzero = matrices->q[vc->iindx[1] - n];
    QAB = (matrices->q[vc->iindx[1] - n]- matrices->q[vc->iindx[1] - (vc->cutpoint - 1)] * matrices->q[vc->iindx[vc->cutpoint] - n]) * params->expDuplexInit;
    /*correction for symmetry*/
    if((n - (vc->cutpoint - 1) * 2) == 0){
      if((strncmp(sequence, sequence + vc->cutpoint - 1, vc->cutpoint - 1)) == 0){
        QAB/=2;
      }
    }

    QToT    = matrices->q[vc->iindx[1] - (vc->cutpoint - 1)] * matrices->q[vc->iindx[vc->cutpoint] - n] + QAB;
    X.FAB   = -kT * (log(QToT) + n * log(params->pf_scale));
    X.F0AB  = -kT * (log(Qzero)+ n * log(params->pf_scale));
    X.FcAB  = (QAB>1e-17) ? -kT * (log(QAB) + n * log(params->pf_scale)) : 999;
    X.FA    = -kT * (log(matrices->q[vc->iindx[1] - (vc->cutpoint - 1)]) + (vc->cutpoint - 1) * log(params->pf_scale));
    X.FB    = -kT * (log(matrices->q[vc->iindx[vc->cutpoint] - n]) + (n - vc->cutpoint + 1) * log(params->pf_scale));

    /* printf("QAB=%.9f\tQtot=%.9f\n",QAB/scale[n],QToT/scale[n]);*/
  }
  else {
    X.FA    = X.FB = X.FAB = X.F0AB = free_energy;
    X.FcAB  = 0;
  }

  /* backtracking to construct binding probabilities of pairs*/
  if(md->compute_bpp){
    pf_co_bppm(vc, structure);
    /*
    *  Backward compatibility:
    *  This block may be removed if deprecated functions
    *  relying on the global variable "pr" vanish from within the package!
    */
    pr = vc->exp_matrices->probs;
    /*
    {
      if(pr) free(pr);
      pr = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * ((n+1)*(n+2)/2));
      memcpy(pr, probs, sizeof(FLT_OR_DBL) * ((n+1)*(n+2)/2));
    }
    */
  }

#ifdef SUN4
  standard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(0);
#endif
#endif

  return X;
}

/* forward recursion of pf cofolding */
PRIVATE void
pf_co(vrna_fold_compound *vc){

  int               n, i,j,k,l, ij, kl, u,u1,u2,ii, type, type_2, tt, cp, turn, maxk, minl;
  FLT_OR_DBL        *qqm = NULL, *qqm1 = NULL, *qq = NULL, *qq1 = NULL;
  FLT_OR_DBL        temp, q_temp, Qmax=0;
  FLT_OR_DBL        qbt1, *tmp;
  FLT_OR_DBL        *q, *qb, *qm, *qm1;
  FLT_OR_DBL        *scale;
  FLT_OR_DBL        *expMLbase;
  short             *S1;
  int               *my_iindx, *jindx;
  char              *ptype, *sequence;
  model_detailsT    *md;
  hard_constraintT  *hc;
  soft_constraintT  *sc;
  FLT_OR_DBL        expMLclosing;
  int               noGUclosure;
  double            max_real;
  int               *rtype;
  pf_paramT         *pf_params;
  pf_matricesT      *matrices;
  int               hc_decompose;
  char              *hard_constraints;
  int               *hc_up_ext;
  int               *hc_up_hp;
  int               *hc_up_int;
  int               *hc_up_ml;

  sequence          = vc->sequence;
  S1                = vc->sequence_encoding;
  n                 = vc->length;
  cp                = vc->cutpoint;
  my_iindx          = vc->iindx;
  jindx             = vc->jindx;
  ptype             = vc->ptype;
  pf_params         = vc->exp_params;
  md                = &(pf_params->model_details);
  rtype             = &(md->rtype[0]);
  hc                = vc->hc;
  sc                = vc->sc;
  expMLclosing      = pf_params->expMLclosing;
  noGUclosure       = md->noGUclosure;
  matrices          = vc->exp_matrices;
  turn              = md->min_loop_size;

  q                 = matrices->q;
  qb                = matrices->qb;
  qm                = matrices->qm;
  qm1               = matrices->qm1;
  scale             = matrices->scale;
  expMLbase         = matrices->expMLbase;

  hard_constraints  = hc->matrix;
  hc_up_ext         = hc->up_ext;
  hc_up_hp          = hc->up_hp;
  hc_up_int         = hc->up_int;
  hc_up_ml          = hc->up_ml;

  max_real          = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  /* allocate memory for helper arrays */
  qq        = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));
  qq1       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));
  qqm       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));
  qqm1      = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));

  /*array initialization ; qb,qm,q
    qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  /* for (d=0; d<=TURN; d++) */
  for (i=1; i<=n/*-d*/; i++) {
      ij      = my_iindx[i]-i;
      if(hc_up_ext[i]){
        q[ij] = scale[1];

        if(sc){
          if(sc->boltzmann_factors)
            q[ij] *= sc->boltzmann_factors[i][1];
        }
      } else {
        q[ij] = 0.;
      }

      qb[ij]  = qm[ij] = 0.0;
    }

  for (i=0; i<=n; i++)
    qq[i] = qq1[i] = qqm[i] = qqm1[i] = 0;

  for (j = turn + 2; j <= n; j++) {
    for (i = j - turn - 1; i >= 1; i--) {
      /* construction of partition function of segment i,j*/
       /*firstly that given i bound to j : qb(i,j) */
      u             = j - i - 1;
      ij            = my_iindx[i] - j;
      type          = ptype[jindx[j] + i];
      hc_decompose  = hard_constraints[jindx[j] + i];
      qbt1          = 0;
      q_temp        = 0.;

      if(hc_decompose){
        /*hairpin contribution*/
        if(hc_decompose & VRNA_HC_CONTEXT_HP_LOOP){
          if(hc_up_hp[i+1] >= u){
            if (ON_SAME_STRAND(i,j,cp)){
              if(((type==3)||(type==4)) && noGUclosure)
                qbt1 = 0;
              else{
                qbt1 =  exp_E_Hairpin(u, type, S1[i+1], S1[j-1], sequence+i-1, pf_params)
                        * scale[u+2];
                if(sc){
                  if(sc->boltzmann_factors)
                    qbt1 *= sc->boltzmann_factors[i+1][u];

                  if(sc->exp_en_basepair)
                    qbt1 *= sc->exp_en_basepair[ij];

                  if(sc->exp_f)
                    qbt1 *= sc->exp_f(i, j, n, n, VRNA_DECOMP_PAIR_HP, sc->data);
                }
              }
            }
          }
        }

        /* interior loops with interior pair k,l */
        if(hc_decompose & VRNA_HC_CONTEXT_INT_LOOP){
          maxk = i + MAXLOOP + 1;
          maxk = MIN2(maxk, j - turn - 2);
          maxk = MIN2(maxk, i + 1 + hc_up_int[i+1]);

          for (k = i + 1; k <= maxk; k++) {
            u1    = k-i-1;
            minl  = MAX2(k + turn + 1, j - 1 - MAXLOOP + u1);
            kl    = my_iindx[k] - j + 1;
            for (u2 = 0, l=j-1; l>=minl; l--, kl++, u2++){
              if(hc_up_int[l+1] < u2) break;
              if(hard_constraints[jindx[l] + k] & VRNA_HC_CONTEXT_INT_LOOP_ENC){

                if ((ON_SAME_STRAND(i,k,cp))&&(ON_SAME_STRAND(l,j,cp))){
                  type_2  =   ptype[jindx[l] + k];
                  type_2  =   rtype[type_2];
                  q_temp  +=  qb[my_iindx[k]-l]
                              * exp_E_IntLoop(u1, j-l-1, type, type_2, S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params)
                              * scale[u1+j-l+1];

                  if(sc){
                    if(sc->boltzmann_factors)
                      q_temp *= sc->boltzmann_factors[i+1][u1]
                                * sc->boltzmann_factors[l+1][u2];

                    if(sc->exp_en_basepair)
                      q_temp *= sc->exp_en_basepair[ij];

                    if(sc->exp_f)
                      q_temp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);

                    if(sc->exp_en_stack)
                      if((i+1 == k) && (j-1 == l)){
                        q_temp *=   sc->exp_en_stack[i]
                                  * sc->exp_en_stack[k]
                                  * sc->exp_en_stack[l]
                                  * sc->exp_en_stack[j];
                      }
                  }

                  qbt1 += q_temp;
                }
              }
            }
          }
        }

        /*multiple stem loop contribution*/
        if(hc_decompose & VRNA_HC_CONTEXT_MB_LOOP){
          ii = my_iindx[i+1]; /* ii-k=[i+1,k-1] */
          temp = 0.0;
          kl = my_iindx[i+1]-(i+1);
          if (ON_SAME_STRAND(i,i+1,cp) && ON_SAME_STRAND(j-1,j,cp)) {
            for (k = i + 2; k <= j - 1; k++,kl--){
              if (ON_SAME_STRAND(k-1,k,cp)){
                q_temp = qm[kl] * qqm1[k];

                if(sc){
                  if(sc->exp_en_basepair)
                    q_temp *= sc->exp_en_basepair[ij];

                  if(sc->exp_f)
                    q_temp *= sc->exp_f(i, j, k, n, VRNA_DECOMP_PAIR_ML, sc->data);
                }

                temp += q_temp;
              }
            }
            tt    =   rtype[type];
            qbt1  +=  temp
                      * expMLclosing
                      * exp_E_MLstem(tt, S1[j-1], S1[i+1], pf_params)
                      * scale[2];
          }
        }

        /*qc contribution*/
        temp=0.0;
        if (!ON_SAME_STRAND(i,j,cp)){
          tt = rtype[type];
          temp=q[my_iindx[i+1]-(cp-1)]*q[my_iindx[cp]-(j-1)];
          if ((j==cp)&&(i==cp-1)) temp=scale[2];
          else if (i==cp-1) temp=q[my_iindx[cp]-(j-1)]*scale[1];
          else if (j==cp) temp=q[my_iindx[i+1]-(cp-1)]*scale[1];
          if (j>cp) temp*=scale[1];
          if (i<cp-1) temp*=scale[1];
          temp  *=  exp_E_ExtLoop(tt, ON_SAME_STRAND(j-1,j,cp) ? S1[j-1] : -1, ON_SAME_STRAND(i,i+1,cp) ? S1[i+1] : -1, pf_params);
          qbt1  +=  temp;
        }
        qb[ij] = qbt1;
      } else  /* end if allowed to be paired */
        qb[ij] = 0.0;

      /* construction of qqm matrix containing final stem
         contributions to multiple loop partition function
         from segment i,j */
      qqm[i] = 0.;

      if(hc_up_ml[j]){
        if (ON_SAME_STRAND(j-1,j,cp)) {
          q_temp  =  qqm1[i] * expMLbase[1];

          if(sc){
            if(sc->boltzmann_factors)
              q_temp *= sc->boltzmann_factors[j][1];

            if(sc->exp_f)
              q_temp *= sc->exp_f(i, j, j-1, n, VRNA_DECOMP_ML_UP_3, sc->data);
          }

          qqm[i] = q_temp;
        }
      }

      if(hc_decompose & VRNA_HC_CONTEXT_MB_LOOP_ENC){
        if(ON_SAME_STRAND(i-1,i,cp) && ON_SAME_STRAND(j,j+1,cp)){
          qbt1    =   qb[ij];
          qbt1    *=  exp_E_MLstem(type, (i>1) ? S1[i-1] : -1, (j<n) ? S1[j+1] : -1, pf_params);
          qqm[i]  +=  qbt1;
        }
      }

      if (qm1) qm1[jindx[j]+i] = qqm[i]; /* for stochastic backtracking */


      /*construction of qm matrix containing multiple loop
        partition function contributions from segment i,j */
      temp  = 0.0;
      kl = my_iindx[i] - j + 1; /* ii-k=[i,k-1] */
      for (k=j; k>i; k--, kl++){
        if (ON_SAME_STRAND(k-1,k,cp)){
          q_temp = qm[kl] * qqm[k];

          if(sc){
            if(sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, n, VRNA_DECOMP_ML_ML_ML, sc->data);
          }

          temp += q_temp;
        }
      }

      maxk = MIN2(i+hc_up_ml[i], j);
      ii = 1; /* length of unpaired stretch */
      for (k=i+1; k<=maxk; k++, ii++){
        if (ON_SAME_STRAND(i,k,cp)){
          q_temp = expMLbase[ii] * qqm[k];

          if(sc){
            if(sc->boltzmann_factors)
              q_temp *= sc->boltzmann_factors[i][ii];

            if(sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, n, VRNA_DECOMP_ML_UP_5, sc->data);
          }

          temp += q_temp;
        }
      }

      qm[ij] = (temp + qqm[i]);

      /*auxiliary matrix qq for cubic order q calculation below */
      qbt1 = 0.;

      if(hc_decompose & VRNA_HC_CONTEXT_EXT_LOOP){
        qbt1 =  qb[ij]
                * exp_E_ExtLoop(type, ((i>1)&&(ON_SAME_STRAND(i-1,i,cp))) ? S1[i-1] : -1, ((j<n)&&(ON_SAME_STRAND(j,j+1,cp))) ? S1[j+1] : -1, pf_params);
      }

      if(hc_up_ext[j]){
        q_temp = qq1[i] * scale[1];

        if(sc){
          if(sc->boltzmann_factors)
            q_temp *= sc->boltzmann_factors[j][1];

          if(sc->exp_f)
            q_temp *= sc->exp_f(i, j, j-1, n, VRNA_DECOMP_EXT_UP_3, sc->data);
        }

        qbt1 += q_temp;
      }
      qq[i] = qbt1;

       /*construction of partition function for segment i,j */
      temp = qq[i];

      /* the whole stretch [i,j] is unpaired */
      if(hc_up_ext[i] >= (j-i+1)){
        q_temp = 1.0 * scale[j-i+1];

        if(sc){
          if(sc->boltzmann_factors)
            q_temp *= sc->boltzmann_factors[i][j-i+1];

          if(sc->exp_f)
            q_temp *= sc->exp_f(i, j, n, n, VRNA_DECOMP_EXT_UP, sc->data);
        }

        temp += q_temp;
      }

      kl = my_iindx[i] - i;
      for (k=i; k<j; k++, kl--){
        q_temp = q[kl] * qq[k+1];

        if(sc){
          if(sc->exp_f)
            q_temp *= sc->exp_f(i, j, k, n, VRNA_DECOMP_EXT_EXT, sc->data);
        }

        temp += q_temp;
      }
      q[ij] = temp;

      if (temp>Qmax) {
        Qmax = temp;
        if (Qmax>max_real/10.)
          fprintf(stderr, "Q close to overflow: %d %d %g\n", i,j,temp);
      }
      if (temp>=max_real) {
        PRIVATE char msg[128];
        snprintf(msg, 127, "overflow in co_pf_fold while calculating q[%d,%d]\n"
                "use larger pf_scale", i,j);
        nrerror(msg);
      }
    }
    tmp = qq1;  qq1 =qq;  qq =tmp;
    tmp = qqm1; qqm1=qqm; qqm=tmp;
  }

  /* clean up */
  free(qq);
  free(qq1);
  free(qqm);
  free(qqm1);

}

/* backward recursion of pf cofolding */
PRIVATE void
pf_co_bppm(vrna_fold_compound *vc, char *structure){

  int               n, i,j,k,l, ij, kl, ii, ll, lj, u1, u2, type, type_2, tt, turn, ov=0, *my_iindx, *jindx, cp;
  FLT_OR_DBL        temp, Qmax=0, prm_MLb, tmp2, ppp;
  FLT_OR_DBL        prmt,prmt1, *expMLbase;
  FLT_OR_DBL        *tmp;
  FLT_OR_DBL        expMLclosing, *probs, *q1k, *qln, *q, *qb, *qm, *scale;
  double            max_real;
  pf_paramT         *pf_params;
  model_detailsT    *md;
  short             *S,*S1;
  char              *ptype;
  hard_constraintT  *hc;
  soft_constraintT  *sc;
  pf_matricesT      *matrices;
  char              *sequence;
  char              *hard_constraints;
  int               *hc_up_ext;
  int               *hc_up_hp;
  int               *hc_up_int;
  int               *hc_up_ml;

  sequence          = vc->sequence;
  n                 = vc->length;
  cp                = vc->cutpoint;
  pf_params         = vc->exp_params;
  md                = &(pf_params->model_details);
  expMLclosing      = pf_params->expMLclosing;
  S                 = vc->sequence_encoding2;
  S1                = vc->sequence_encoding;
  jindx             = vc->jindx;
  my_iindx          = vc->iindx;
  ptype             = vc->ptype;
  turn              = md->min_loop_size;

  matrices          = vc->exp_matrices;
  probs             = matrices->probs;
  scale             = matrices->scale;
  q1k               = matrices->q1k;
  qln               = matrices->qln;
  q                 = matrices->q;
  qb                = matrices->qb;
  qm                = matrices->qm;
  expMLbase         = matrices->expMLbase;

  hard_constraints  = hc->matrix;
  hc_up_ext         = hc->up_ext;
  hc_up_hp          = hc->up_hp;
  hc_up_int         = hc->up_int;
  hc_up_ml          = hc->up_ml;

  hc            = vc->hc;
  sc            = vc->sc;

  max_real      = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  /* backtracking to construct binding probabilities of pairs*/
  if ((S != NULL) && (S1 != NULL)) {
    FLT_OR_DBL   *Qlout, *Qrout;
    FLT_OR_DBL *prm_l  = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));
    FLT_OR_DBL *prm_l1 = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));
    FLT_OR_DBL *prml   = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(n+2));

    Qmax  = 0;
    Qrout = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * (n+2));
    Qlout = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * (cp+2));

    for (k=1; k<=n; k++) {
      q1k[k] = q[my_iindx[1] - k];
      qln[k] = q[my_iindx[k] -n];
    }
    q1k[0] = 1.0;
    qln[n+1] = 1.0;

    /* 1. exterior pair i,j and initialization of pr array */
    for (i=1; i<=n; i++) {
      for (j=i; j<=MIN2(i + turn, n); j++)
        probs[my_iindx[i]-j] = 0;

      for (j = i + turn + 1; j <= n; j++){
        ij = my_iindx[i]-j;
        if(hard_constraints[jindx[j] + i] & VRNA_HC_CONTEXT_EXT_LOOP){
          type  = ptype[jindx[j] + i];
          probs[ij] = q1k[i-1]*qln[j+1]/q1k[n];
          probs[ij] *= exp_E_ExtLoop(type, ((i>1)&&(ON_SAME_STRAND(i-1,i,cp))) ? S1[i-1] : -1, ((j<n)&&(ON_SAME_STRAND(j,j+1,cp))) ? S1[j+1] : -1, pf_params);
        } else
          probs[ij] = 0;
      }
    }

    for(l = n; l > turn + 1; l--){

      /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
      for(k = 1; k < l - turn; k++){
        kl      = my_iindx[k]-l;
        type_2  = ptype[jindx[l] + k];
        type_2  = rtype[type_2];

        if(qb[kl]==0.) continue;

        if(hard_constraints[jindx[l] + k] & VRNA_HC_CONTEXT_INT_LOOP_ENC){
          for(i = MAX2(1, k - MAXLOOP - 1); i <= k - 1; i++){
            u1 = k - i - 1;
            if(hc_up_int[i+1] < u1) continue;

            for(j = l + 1; j <= MIN2(l + MAXLOOP - k + i + 2, n); j++){
              u2 = j-l-1;
              if(hc_up_int[l+1] < u2) break;

              ij = my_iindx[i] - j;
              if(hard_constraints[jindx[j] + i] & VRNA_HC_CONTEXT_INT_LOOP){

                if ((ON_SAME_STRAND(i,k,cp)) && (ON_SAME_STRAND(l,j,cp))){
                  type  = ptype[jindx[j] + i];
                  if ((probs[ij]>0)) {
                    tmp2  = probs[ij]
                            * exp_E_IntLoop(k-i-1, j-l-1, type, type_2, S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params)
                            * scale[k-i+j-l];

                    if(sc){
                      if(sc->boltzmann_factors)
                        tmp2 *=   sc->boltzmann_factors[i+1][u1]
                                * sc->boltzmann_factors[l+1][u2];

                      if(sc->exp_en_basepair)
                        tmp2 *=   sc->exp_en_basepair[ij];

                      if(sc->exp_en_stack){
                        if((i+1 == k) && (j-1 == l)){
                          tmp2 *=   sc->exp_en_stack[i]
                                  * sc->exp_en_stack[k]
                                  * sc->exp_en_stack[l]
                                  * sc->exp_en_stack[j];
                        }
                      }

                      if(sc->exp_f)
                        tmp2 *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);
                    }

                    probs[kl] += tmp2;
                  }
                }
              }
            }
          }
        }
      }

      /* 3. bonding k,l as substem of multi-loop enclosed by i,j */
      prm_MLb = 0.;
      if((l < n) && (ON_SAME_STRAND(l, l + 1, cp)))
        for (k = 2; k < l - turn; k++) {
          kl    = my_iindx[k]-l;
          i     = k - 1;
          prmt  = prmt1 = 0.0;

          ii    = my_iindx[i];     /* ii-j=[i,j]     */
          ll    = my_iindx[l+1];   /* ll-j=[l+1,j] */
          tt    = ptype[jindx[l+1] + i];
          tt    = rtype[tt];
          if (ON_SAME_STRAND(i,k,cp)){
            if(hard_constraints[jindx[l+1] + i] & VRNA_HC_CONTEXT_MB_LOOP){
              prmt1 = probs[ii-(l+1)]
                      * expMLclosing
                      * exp_E_MLstem(tt, S1[l], S1[i+1], pf_params);

              if(sc){
                /* which decompositions are covered here? => (i, l+1) -> enclosing pair, (k,l) -> enclosed pair, */
                if(sc->exp_en_basepair)
                  prmt1 *= sc->exp_en_basepair[ii - (l+1)];

/*
                if(sc->exp_f)
                  prmt1 *= sc->exp_f(i, l+1, k, l, , sc->data);
*/
              }
            }
            ij = my_iindx[i] - (l+2);
            lj = my_iindx[l+1]-(l+1);

            for(j = l + 2; j <= n; j++, ij--, lj--){
              if(hard_constraints[jindx[j] + i] & VRNA_HC_CONTEXT_MB_LOOP){
                if (ON_SAME_STRAND(j-1,j,cp)){ /*??*/
                  tt    =   ptype[jindx[j] + i];
                  tt    =   rtype[tt];
                  /* which decomposition is covered here? =>
                    i + 1 = k < l < j:
                    (i,j)       -> enclosing pair
                    (k, l)      -> enclosed pair
                    (l+1, j-1)  -> multiloop part with at least one stem
                  */
                  ppp +=  probs[ii-j]
                          * exp_E_MLstem(tt, S1[j-1], S1[i+1], pf_params)
                          * qm[ll-(j-1)];

                  if(sc){
                    if(sc->exp_en_basepair)
                      ppp *= sc->exp_en_basepair[ij];
/*
                    if(sc->exp_f)
                      ppp *= sc->exp_f(i, j, l+1, j-1, , sc->data);
*/
                  }
                  prmt += ppp;
                }
              }
            }
          }
          prmt *= expMLclosing;

          tt        =   ptype[jindx[l] + k];

          prml[ i]  =   prmt;

          /* l+1 is unpaired */
          ppp = prm_l1[i] * expMLbase[1];
          if(sc){
            if(sc->boltzmann_factors)
              ppp *= sc->boltzmann_factors[l+1][1];

/*
            if(sc_exp_f)
              ppp *= sc->exp_f(, sc->data);
*/
          }
          prm_l[i]  =   ppp + prmt1;

          /* i is unpaired */
          ppp = prm_MLb*expMLbase[1];
          if(sc){
            if(sc->boltzmann_factors)
              ppp *= sc->boltzmann_factors[i][1];

/*
            if(sc->exp_f)
              ppp *= sc->exp_f(, sc->data);
*/
          }
          prm_MLb = ppp + prml[i];
          /* same as:    prm_MLb = 0;
             for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */

          prml[i]   =   prml[ i] + prm_l[i];

          if (qb[kl] == 0.) continue;

          if(hard_constraints[jindx[l] + k] & VRNA_HC_CONTEXT_MB_LOOP_ENC){
            temp = prm_MLb;

            for (i=1;i<=k-2; i++) {
              if ((ON_SAME_STRAND(i,i+1,cp)) && (ON_SAME_STRAND(k-1,k,cp))){
                temp += prml[i]*qm[my_iindx[i+1] - (k-1)];
              }
            }
            temp *= exp_E_MLstem( tt,
                                  ((k>1) && ON_SAME_STRAND(k-1,k,cp)) ? S1[k-1] : -1,
                                  ((l<n) && ON_SAME_STRAND(l,l+1,cp)) ? S1[l+1] : -1,
                                  pf_params) * scale[2];
            probs[kl] += temp;

            if (probs[kl]>Qmax) {
              Qmax = probs[kl];
              if (Qmax>max_real/10.)
                fprintf(stderr, "P close to overflow: %d %d %g %g\n",
                        i, j, probs[kl], qb[kl]);
            }
            if (probs[kl]>=max_real) {
              ov++;
              probs[kl]=FLT_MAX;
            }
          }
        } /* end for (k=..) multloop*/
      else  /* set prm_l to 0 to get prm_l1 to be 0 */
        for (i=0; i<=n; i++) prm_l[i]=0;

      tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;
      /*computation of .(..(...)..&..). type features?*/
      if (cp<=0) continue;  /* no .(..(...)..&..). type features*/
      if ((l==n)||(l<=2)) continue; /* no .(..(...)..&..). type features*/
      /*new version with O(n^3)??*/
      if (l>cp) {
        if (l<n) {
          int t,kt;
          for (t=n; t>l; t--) {
            for (k=1; k<cp; k++) {
              kt    = my_iindx[k]-t;
              type  = rtype[(unsigned char)ptype[jindx[t] + k]];
              temp  = probs[kt]
                      * exp_E_ExtLoop(type, S1[t-1], (ON_SAME_STRAND(k,k+1,cp)) ? S1[k+1] : -1, pf_params)
                      * scale[2];

              if (l+1<t)
                temp    *=  q[my_iindx[l+1]-(t-1)];

              if (ON_SAME_STRAND(k,k+1,cp))
                temp    *=  q[my_iindx[k+1]-(cp-1)];

              Qrout[l]  +=  temp;
            }
          }
        }
        for (k=l-1; k>=cp; k--) {
          if (qb[my_iindx[k]-l]) {
            kl        =   my_iindx[k]-l;
            type      =   ptype[jindx[l] + k];
            temp      =   Qrout[l];
            temp      *=  exp_E_ExtLoop(type, (k>cp) ? S1[k-1] : -1, (l < n) ? S1[l+1] : -1, pf_params);
            if (k>cp)
              temp    *=  q[my_iindx[cp]-(k-1)];
            probs[kl] +=  temp;
          }
        }
      }
      else if (l==cp ) {
        int t, sk,s;
        for (t=2; t<cp;t++) {
          for (s=1; s<t; s++) {
            for (k=cp; k<=n; k++) {
              sk=my_iindx[s]-k;
              if (qb[sk]) {
                type      =   rtype[(unsigned char)ptype[jindx[k] + s]];
                temp      =   probs[sk]
                              * exp_E_ExtLoop(type, (ON_SAME_STRAND(k-1,k,cp)) ? S1[k-1] : -1, S1[s+1], pf_params)
                              * scale[2];
                if (s+1<t)
                  temp    *=  q[my_iindx[s+1]-(t-1)];
                if (ON_SAME_STRAND(k-1,k,cp))
                  temp    *=  q[my_iindx[cp]-(k-1)];
                Qlout[t]  +=  temp;
              }
            }
          }
        }
      }
      else if (l<cp) {
        for (k=1; k<l; k++) {
          if (qb[my_iindx[k]-l]) {
            type    =   ptype[jindx[l] + k];
            temp    =   Qlout[k];
            temp    *=  exp_E_ExtLoop(type, (k>1) ? S1[k-1] : -1, (l<(cp-1)) ? S1[l+1] : -1, pf_params);
            if (l+1<cp)
              temp  *=  q[my_iindx[l+1]-(cp-1)];
            probs[my_iindx[k]-l]  +=  temp;
          }
        }
      }
    }  /* end for (l=..)   */
    free(Qlout);
    free(Qrout);
    for (i=1; i<=n; i++)
      for (j=i+TURN+1; j<=n; j++) {
        ij        =   my_iindx[i]-j;
        probs[ij] *=  qb[ij];
      }

    if (structure!=NULL)
      bppm_to_structure(structure, probs, n);

    /* clean up */
    free(prm_l);
    free(prm_l1);
    free(prml);

  }   /* end if (do_backtrack)*/

  if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
                    "you might try a smaller pf_scale than %g\n",
                    ov, pf_params->pf_scale);
}



PRIVATE void
wrap_update_pf_params(int length,
                      pf_paramT *parameters){

  pf_paramT *p = NULL;
  if(parameters == NULL){
    model_detailsT md;
    set_model_details(&md);
    p = get_boltzmann_factors(temperature, 1.0, md, pf_scale);
  }

  vrna_update_pf_params(backward_compat_compound, p);
  free(p);
}

/*----------------------------------------------------------------------*/
PUBLIC void
update_co_pf_params(int length){

/*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

  wrap_update_pf_params(length, NULL);
}

PUBLIC void
update_co_pf_params_par(int length,
                        pf_paramT *parameters){

/*---------------------------------------------------------------------------*/
  wrap_update_pf_params(length, parameters);
}

/*
  stochastic backtracking in pf_fold arrays
  returns random structure S with Boltzman probabilty
  p(S) = exp(-E(S)/kT)/Z
*/
PRIVATE void
backtrack_qm1(vrna_fold_compound *vc,
              int i,
              int j,
              char *pstruc){

  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  int           ii, l, type, *jindx, *my_iindx;
  double        qt, r;
  FLT_OR_DBL    *qm, *qm1, *qb, *expMLbase;
  short         *S1;
  char          *ptype;

  pf_paramT     *pf_params;
  pf_matricesT  *matrices;

  pf_params     = vc->exp_params;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;

  matrices      = vc->exp_matrices;
  qb            = matrices->qb;
  qm            = matrices->qm;
  qm1           = matrices->qm1;
  expMLbase     = matrices->expMLbase;

  jindx         = vc->jindx;
  my_iindx      = vc->iindx;

  r   = urn() * qm1[jindx[j]+i];
  ii  = my_iindx[i];
  for (qt=0., l=i+TURN+1; l<=j; l++) {
    type = ptype[jindx[l] + i];
    if (type)
      qt +=  qb[ii-l]*exp_E_MLstem(type, S1[i-1], S1[l+1], pf_params) * expMLbase[j-l];
    if (qt>=r) break;
  }
  if (l>j) nrerror("backtrack failed in qm1");
  backtrack(vc, i,l, pstruc);
}

PRIVATE void
backtrack(vrna_fold_compound *vc,
          int i,
          int j,
          char *pstruc){

  int           *jindx, *my_iindx;
  FLT_OR_DBL    *qm, *qm1, *qb, *expMLbase, *scale;
  pf_paramT     *pf_params;
  pf_matricesT  *matrices;
  short         *S1;
  char          *ptype, *sequence;
  int           noGUclosure;

  sequence      = vc->sequence;
  pf_params     = vc->exp_params;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;

  matrices      = vc->exp_matrices;
  qb            = matrices->qb;
  qm            = matrices->qm;
  qm1           = matrices->qm1;
  expMLbase     = matrices->expMLbase;
  scale         = matrices->scale;
  jindx         = vc->jindx;
  my_iindx      = vc->iindx;
  noGUclosure   = pf_params->model_details.noGUclosure;

  do {
    double r, qbt1;
    int k, l, type, u, u1;

    pstruc[i-1] = '('; pstruc[j-1] = ')';

    r     = urn() * qb[my_iindx[i]-j];
    type  = ptype[jindx[j] + i];
    u     = j - i - 1;
    /*hairpin contribution*/
    if (((type==3)||(type==4))&&noGUclosure) qbt1 = 0;
    else
      qbt1 = exp_E_Hairpin(u, type, S1[i+1], S1[j-1], sequence+i-1, pf_params)*scale[u+2];

    if (qbt1>r) return; /* found the hairpin we're done */

    for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++) {
      u1 = k-i-1;
      for (l=MAX2(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
        int type_2;
        type_2 = ptype[jindx[l] + k];
        if (type_2) {
          type_2  =   rtype[type_2];
          qbt1    +=  qb[my_iindx[k]-l] *
            exp_E_IntLoop(u1, j-l-1, type, type_2,
                          S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params)*scale[u1+j-l+1];
        }
        if (qbt1 > r) break;
      }
      if (qbt1 > r) break;
    }
    if (l<j) {
      i=k; j=l;
    }
    else break;
  } while (1);

  /* backtrack in multi-loop */
  {
    double r, qt;
    int k, ii, jj;

    i++; j--;
    /* find the first split index */
    ii = my_iindx[i]; /* ii-j=[i,j] */
    jj = jindx[j]; /* jj+i=[j,i] */
    for (qt=0., k=i+1; k<j; k++) qt += qm[ii-(k-1)]*qm1[jj+k];
    r = urn() * qt;
    for (qt=0., k=i+1; k<j; k++) {
      qt += qm[ii-(k-1)]*qm1[jj+k];
      if (qt>=r) break;
    }
    if (k>=j) nrerror("backtrack failed, can't find split index ");

    backtrack_qm1(vc, k, j, pstruc);

    j = k-1;
    while (j>i) {
      /* now backtrack  [i ... j] in qm[] */
      jj  = jindx[j];
      ii  = my_iindx[i];
      r   = urn() * qm[ii - j];
      qt  = qm1[jj+i]; k=i;
      if (qt<r)
        for (k=i+1; k<=j; k++) {
          qt += (qm[ii-(k-1)]+expMLbase[k-i])*qm1[jj+k];
          if (qt >= r) break;
        }
      if (k>j) nrerror("backtrack failed in qm");

      backtrack_qm1(vc, k,j, pstruc);

      if (k<i+TURN) break; /* no more pairs */
      r = urn() * (qm[ii-(k-1)] + expMLbase[k-i]);
      if (expMLbase[k-i] >= r) break; /* no more pairs */
      j = k-1;
    }
  }
}

PUBLIC void compute_probabilities(double FAB, double FA,double FB,
                                  struct plist *prAB,
                                  struct plist *prA, struct plist *prB,
                                  int Alength) {
  /*computes binding probabilities and dimer free energies*/
  int     i, j;
  double  pAB;
  double  mykT;
  struct  plist  *lp1, *lp2;
  int     offset;

  mykT = backward_compat_compound->exp_params->kT/1000.;

  /* pair probabilities in pr are relative to the null model (without DuplexInit) */

  /*Compute probabilities pAB, pAA, pBB*/

  pAB = 1. - exp((1/mykT)*(FAB-FA-FB));

  /* compute pair probabilities given that it is a dimer */
  /* AB dimer */
  offset  = 0;
  lp2     = prA;
  if (pAB>0)
    for (lp1=prAB; lp1->j>0; lp1++) {
      float pp=0;
      i = lp1->i;
      j = lp1->j;
      while (offset+lp2->i < i && lp2->i>0) lp2++;
      if (offset+lp2->i == i)
        while ((offset+lp2->j) < j  && (lp2->j>0)) lp2++;
      if (lp2->j == 0) {lp2=prB; offset=Alength;}/* jump to next list */
      if ((offset+lp2->i==i) && (offset+lp2->j ==j)) {
        pp = lp2->p;
        lp2++;
      }
      lp1->p=(lp1->p-(1-pAB)*pp)/pAB;
      if(lp1->p < 0.){
        warn_user("part_func_co: numeric instability detected, probability below zero!");
        lp1->p = 0.;
      }
    }
  return;
}

PRIVATE double *
Newton_Conc(double KAB,
            double KAA,
            double KBB,
            double concA,
            double concB,
            double* ConcVec){

  double  TOL, EPS, xn, yn, det, cA, cB;
  int     i;

  i       = 0;
  /*Newton iteration for computing concentrations*/
  cA      = concA;
  cB      = concB;
  TOL     = 1e-6; /*Tolerance for convergence*/
  ConcVec = (double*)space(5*sizeof(double)); /* holds concentrations */
  do {
    /* det = (4.0 * KAA * cA + KAB *cB + 1.0) * (4.0 * KBB * cB + KAB *cA + 1.0) - (KAB *cB) * (KAB *cA); */
    det = 1 + 16. *KAA*KBB*cA*cB + KAB*(cA+cB) + 4.*KAA*cA + 4.*KBB*cB + 4.*KAB*(KBB*cB*cB + KAA*cA*cA);
    /* xn  = ( (2.0 * KBB * cB*cB + KAB *cA *cB + cB - concB) * (KAB *cA) -
       (2.0 * KAA * cA*cA + KAB *cA *cB + cA - concA) * (4.0 * KBB * cB + KAB *cA + 1.0) ) /det; */
    xn  = ( (2.0 * KBB * cB*cB + cB - concB) * (KAB *cA) - KAB*cA*cB*(4. * KBB*cB + 1.) -
	    (2.0 * KAA * cA*cA + cA - concA) * (4.0 * KBB * cB + KAB *cA + 1.0) ) /det;
    /* yn  = ( (2.0 * KAA * cA*cA + KAB *cA *cB + cA - concA) * (KAB *cB) -
       (2.0 * KBB * cB*cB + KAB *cA *cB + cB - concB) * (4.0 * KAA * cA + KAB *cB + 1.0) ) /det; */
    yn  = ( (2.0 * KAA * cA*cA + cA - concA) * (KAB *cB) - KAB*cA*cB*(4. * KAA*cA + 1.) -
            (2.0 * KBB * cB*cB + cB - concB) * (4.0 * KAA * cA + KAB *cB + 1.0) ) /det;
    EPS = fabs(xn/cA) + fabs(yn/cB);
    cA += xn;
    cB += yn;
    i++;
    if (i>10000) {
      fprintf(stderr, "Newton did not converge after %d steps!!\n",i);
      break;
    }
  } while(EPS>TOL);

  ConcVec[0] = cA*cB*KAB ;/*AB concentration*/
  ConcVec[1] = cA*cA*KAA ;/*AA concentration*/
  ConcVec[2] = cB*cB*KBB ;/*BB concentration*/
  ConcVec[3] = cA;        /* A concentration*/
  ConcVec[4] = cB;        /* B concentration*/

  return ConcVec;
}

PUBLIC struct ConcEnt *
get_concentrations( double FcAB,
                    double FcAA,
                    double FcBB,
                    double FEA,
                    double FEB,
                    double *startconc){

  /*takes an array of start concentrations, computes equilibrium concentrations of dimers, monomers, returns array of concentrations in strucutre ConcEnt*/
  double          *ConcVec;
  int             i;
  struct  ConcEnt *Concentration;
  double          KAA, KAB, KBB, kT;

  kT            = backward_compat_compound->exp_params->kT/1000.;
  Concentration = (struct ConcEnt *)space(20*sizeof(struct ConcEnt));
 /* Compute equilibrium constants */
  /* again note the input free energies are not from the null model (without DuplexInit) */

  KAA = exp(( 2.0 * FEA - FcAA)/kT);
  KBB = exp(( 2.0 * FEB - FcBB)/kT);
  KAB = exp(( FEA + FEB - FcAB)/kT);
  /* printf("Kaa..%g %g %g\n", KAA, KBB, KAB); */
  for (i=0; ((startconc[i]!=0)||(startconc[i+1]!=0));i+=2) {
    ConcVec                 = Newton_Conc(KAB, KAA, KBB, startconc[i], startconc[i+1], ConcVec);
    Concentration[i/2].A0   = startconc[i];
    Concentration[i/2].B0   = startconc[i+1];
    Concentration[i/2].ABc  = ConcVec[0];
    Concentration[i/2].AAc  = ConcVec[1];
    Concentration[i/2].BBc  = ConcVec[2];
    Concentration[i/2].Ac   = ConcVec[3];
    Concentration[i/2].Bc   = ConcVec[4];

    if (!(((i+2)/2)%20))  {
      Concentration = (struct ConcEnt *)xrealloc(Concentration,((i+2)/2+20)*sizeof(struct ConcEnt));
    }
    free(ConcVec);
  }

  return Concentration;
}

PUBLIC FLT_OR_DBL *
export_co_bppm(void){

  if(backward_compat_compound)
    return backward_compat_compound->exp_matrices->probs;
  else
    return NULL;
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC cofoldF
co_pf_fold(char *sequence, char *structure){

  return wrap_co_pf_fold(sequence, structure, NULL, do_backtrack, fold_constrained);
}

PUBLIC cofoldF
co_pf_fold_par( char *sequence,
                char *structure,
                pf_paramT *parameters,
                int calculate_bppm,
                int is_constrained){

  return wrap_co_pf_fold(sequence, structure, parameters, calculate_bppm, is_constrained);
}


PUBLIC struct plist *
get_plist(struct plist *pl,
          int length,
          double cut_off){

  int i, j,n, count, *my_iindx;

  my_iindx = backward_compat_compound->iindx;
  /*get pair probibilities out of pr array*/
  count=0;
  n=2;
  for (i=1; i<length; i++) {
    for (j=i+1; j<=length; j++) {
      if (pr[my_iindx[i]-j]<cut_off) continue;
      if (count==n*length-1) {
        n*=2;
        pl=(struct plist *)xrealloc(pl,n*length*sizeof(struct plist));
      }
      pl[count].i=i;
      pl[count].j=j;
      pl[count++].p=pr[my_iindx[i]-j];
      /*      printf("gpl: %2d %2d %.9f\n",i,j,pr[my_iindx[i]-j]);*/
    }
  }
  pl[count].i=0;
  pl[count].j=0; /*->??*/
  pl[count++].p=0.;
  pl=(struct plist *)xrealloc(pl,(count)*sizeof(struct plist));
  return pl;
}

PUBLIC void
init_co_pf_fold(int length){

 /* DO NOTHING */
}

PUBLIC void
free_co_pf_arrays(void){

  if(backward_compat_compound && backward_compat){
    vrna_free_fold_compound(backward_compat_compound);
    backward_compat_compound  = NULL;
    backward_compat           = 0;
  }
}
