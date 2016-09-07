/*
                  partiton function for RNA secondary structures

                  Ivo L Hofacker

                  Vienna RNA package
*/
/*
  $Log: part_func_up.c,v $
  Revision 1.4  2008/07/04 14:27:36  ivo
  Modify output (again)

  Revision 1.3  2008/05/08 14:11:55  ivo
  minor output changes

  Revision 1.2  2007/12/13 10:19:54  ivo
  major RNAup update from Ulli

  Revision 1.1  2007/04/30 15:13:13  ivo
  merge RNAup into package

  Revision 1.11  2006/07/17 11:11:43  ulim
  removed all globals from fold_vars.h,c, cleaned code

  Revision 1.10  2006/07/12 09:19:29  ulim
  global variables w, incr3 and incr5 are now local

  Revision 1.9  2006/07/11 12:45:02  ulim
  remove redundancy in function pf_interact(...)

  Revision 1.8  2006/03/08 15:26:37  ulim
  modified -o[1|2], added meaningful default

  Revision 1.5  2006/01/23 11:27:04  ulim
  include file into new package version. cleaned it

  Revision 1.2  2005/07/29 15:13:37  ulim
  put the function, calculating the probability of an unpaired region in
  an RNA and the function calculating the prob. of interaction between 2 RNAs
  in a seperate file (pf_two.c)

  Revision 1.1  2005/07/26 13:27:12  ulim
  Initial revision

  Revision 1.2  2005/07/01 13:14:57  ulim
  fixed error in scaling, included new commandline options -incr5, -incr3 to
  allow a variable number of unpaired positions 5' and 3' of the site of
  interaction between the two RNAs

  Revision 1.1  2005/04/19 08:16:38  ulim
  Initial revision
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include <unistd.h>
#include "ViennaRNA/fold.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/part_func_up.h"
#include "ViennaRNA/duplex.h"


#define CO_TURN 0
#define ZERO(A) (fabs(A) < DBL_EPSILON)
#define EQUAL(A,B) (fabs((A)-(B)) < 1000*DBL_EPSILON)
#define ISOLATED  256.0
/* #define NUMERIC 1 */

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
PRIVATE short       *S=NULL, *S1=NULL, *SS=NULL, *SS2=NULL;
PRIVATE vrna_exp_param_t   *Pf = NULL;/* use this structure for all the exp-arrays*/
PRIVATE FLT_OR_DBL  *qb=NULL, *qm=NULL, *prpr=NULL; /* add arrays for pf_unpaired()*/
PRIVATE FLT_OR_DBL  *probs=NULL;
PRIVATE FLT_OR_DBL  *q1k=NULL, *qln=NULL;
PRIVATE double      *qqm2=NULL, *qq_1m2=NULL, *qqm=NULL, *qqm1=NULL;
PRIVATE FLT_OR_DBL  *scale=NULL, *expMLbase=NULL;
PRIVATE char        *ptype=NULL; /* precomputed array of pair types */
PRIVATE int         init_length;  /* length in last call to init_pf_fold()*/
PRIVATE double      init_temp; /* temperature in last call to scale_pf_params */
PRIVATE int         *my_iindx = NULL;
/* make iptypes array for intermolecular constrains (ipidx for indexing)*/


/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE pu_out      *get_u_vals(pu_contrib *p_c,
                                int **unpaired_values,
                                char *select_contrib);

PRIVATE int         plot_free_pu_out( pu_out* res,
                                      interact *pint,
                                      char *ofile,
                                      char *head);

PRIVATE void        scale_stru_pf_params(unsigned int length);

PRIVATE void        init_pf_two(int length);

PRIVATE void        scale_int(const char *s,
                              const char *sl,
                              double *sc_int);

PRIVATE void        encode_seq( const char *s1,
                                const char *s2);

PRIVATE constrain   *get_ptypes_up(char *S,
                                const char *structure);

PRIVATE void        get_up_arrays(unsigned int length);

PRIVATE void        free_up_arrays(void);

PRIVATE void        set_encoded_seq(const char *sequence,
                                    short **S,
                                    short **S1);

PRIVATE void        get_interact_arrays(unsigned int n1,
                                        unsigned int n2,
                                        pu_contrib *p_c,
                                        pu_contrib *p_c2,
                                        int w,
                                        int incr5,
                                        int incr3,
                                        double ***p_c_S,
                                        double ***p_c2_S);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC pu_contrib *get_pu_contrib_struct(unsigned int n, unsigned int w){
  unsigned int i;
  pu_contrib  *pu = (pu_contrib *)vrna_alloc(sizeof(pu_contrib));
  pu->length      = n;
  pu->w           = w;
  /* contributions to probability of being unpaired witihin a(n)
   H hairpin,
   I interior loop,
   M muliloop,
   E exterior loop*/
  /* pu_test->X[i][j] where i <= j and i [1...n], j = [1...w[ */
  pu->H           = (double **)vrna_alloc(sizeof(double *) * (n + 1));
  pu->I           = (double **)vrna_alloc(sizeof(double *) * (n + 1));
  pu->M           = (double **)vrna_alloc(sizeof(double *) * (n + 1));
  pu->E           = (double **)vrna_alloc(sizeof(double *) * (n + 1));
  for(i=0;i<=n;i++){
    pu->H[i]  = (double *)vrna_alloc(sizeof(double) * (w + 1));
    pu->I[i]  = (double *)vrna_alloc(sizeof(double) * (w + 1));
    pu->M[i]  = (double *)vrna_alloc(sizeof(double) * (w + 1));
    pu->E[i]  = (double *)vrna_alloc(sizeof(double) * (w + 1));
  }
  return pu;
}

PUBLIC  void
free_pu_contrib(pu_contrib *pu){

  free_pu_contrib_struct(pu);
}

PUBLIC  void
free_pu_contrib_struct(pu_contrib *pu){

  unsigned int i;
  if(pu != NULL){
    for(i=0;i<=pu->length;i++){
      free(pu->H[i]);
      free(pu->I[i]);
      free(pu->M[i]);
      free(pu->E[i]);
    }
    free(pu->H);
    free(pu->I);
    free(pu->M);
    free(pu->E);
    free(pu);
  }
}

/* you have to call pf_fold(sequence, structure); befor pf_unstru */
PUBLIC pu_contrib *pf_unstru(char *sequence, int w){
  int           n, i, j, v, k, l, o, p, ij, kl, po, u, u1, d, type, type_2, tt;
  unsigned int  size;
  double        temp, tqm2;
  double        qbt1, *tmp, sum_l, *sum_M;
  double        *store_H, *store_Io, **store_I2o; /* hairp., interior contribs */
  double        *store_M_qm_o,*store_M_mlbase;    /* multiloop contributions */
  pu_contrib    *pu_test;

  sum_l           = 0.0;
  temp            = 0;
  n               = (int) strlen(sequence);
  sum_M           = (double *)  vrna_alloc((n+1) * sizeof(double));
  pu_test         = get_pu_contrib_struct((unsigned)n, (unsigned)w);
  size            = ((n+1)*(n+2))>>1;

  get_up_arrays((unsigned) n);
  init_pf_two(n);

  /* init everything */
  for (d=0; d<=TURN; d++)
    for (i=1; i<=n-d; i++){
      j=i+d;
      ij = my_iindx[i]-j;
      if(d < w) {
        pu_test->H[i][d]=pu_test->I[i][d]=pu_test->M[i][d]=pu_test->E[i][d]=0.;
      }
    }


  for (i=0; i<size; i++)
    prpr[i]= probs[i];

  sum_M[0] = 0.;
  for (i=1; i<=n; i++){
    /* set auxillary arrays to 0, reuse qqm and qqm1, reuse qqm2 and qq_1m2*/
    sum_M[i] = qqm[i] = qqm1[i] = qqm2[i] = qq_1m2[i] = 0;
    for (j=i+TURN+1; j<=n; j++){
      ij = my_iindx[i]-j;
      /* i need the part_func of all structures outside bp[ij] */
      if(qb[ij] > 0.0) prpr[ij]= (probs[ij]/qb[ij]);
    }
  }

  /* alloc even more memory */
  store_I2o = (double **)vrna_alloc(sizeof(double *) * (n + 1)); /* for p,k */
  for(i=0;i<=n;i++)
    store_I2o[i] = (double *)vrna_alloc(sizeof(double) * (MAXLOOP + 2));

  /* expMLbase[i-p]*dangles_po */
  store_M_mlbase = (double *)vrna_alloc(sizeof(double) * (size + 1));

  /* 2. exterior bp (p,o) encloses unpaired region [i,i+w[*/
  for (o=TURN+2;o<=n; o++) {
    double sum_h;
    /*allocate space for arrays to store different contributions to H, I & M */
    store_H       = (double *)vrna_alloc(sizeof(double) * (o+2));
    /* unpaired between ]l,o[ */
    store_Io      = (double *)vrna_alloc(sizeof(double) * (o+2));
    /* qm[p+1,i-1]*dangles_po */
    store_M_qm_o  = (double *)vrna_alloc(sizeof(double) * (n+1));

    for (p=o-TURN-1; p>=1; p--) {
      /* construction of partition function of segment [p,o], given that
         an unpaired region [i,i+w[ exists within [p,o] */
      u = o-p-1;
      po = my_iindx[p]-o;
      type = ptype[po];
      if(type){

        /*hairpin contribution*/
        if (((type==3)||(type==4))&&no_closingGU)
          temp = 0.;
        else
          temp = prpr[po] * exp_E_Hairpin(u, type, S1[p+1], S1[o-1], sequence+p-1, Pf) * scale[u+2];
        /* all H contribs are collect for the longest unpaired region */
        store_H[p+1] = temp;

        /* interior loops with interior pair k,l and an unpaired region of
         length w between p and k || l and o*/
        for (k=p+1; k<=MIN2(p+MAXLOOP+1,o-TURN-2); k++) {
          u1    = k-p-1;
          sum_l = 0.;
          for (l=MAX2(k+TURN+1,o-1-MAXLOOP+u1); l<o; l++) {
            kl      = my_iindx[k]-l;
            type_2  = ptype[kl];
            if((l+1) < o) store_Io[l+1] += sum_l;

            temp=0.;
            if (type_2){
              type_2 = rtype[type_2];
              temp = prpr[po] * qb[kl] * exp_E_IntLoop(u1, o-l-1, type, type_2, S1[p+1], S1[o-1], S1[k-1], S1[l+1], Pf) *scale[u1+o-l+1];
              if((l+1) < o) store_Io[l+1] += temp; /* unpaired region between ]l,o[ */
              sum_l += temp;
            } /* end of if pair(k,l) */
          } /* end of l */
          /* unpaired in region ]p,k[  */
          for(i=p+1;i <= k-1;i++)
            store_I2o[i][MIN2(w-1,k-i-1)] += sum_l;
        } /* end of k */
      } /*end of if(type) test for bp (p,o) */

      /* multiple stem loop contribution
         calculate qm2[my_iindx[i]-j] in the course of the calculation
         of the multiple stem loop contribution:
         advantage: you save memory:
         instead of a (n+1)*n array for qqm2 you only need 2*n arrays
         disadvantage: you have to use two times the op-loop for the full
         multiloop contribution
         first op-loop: index o goes from 1...n and
                        index p from o-TURN-1 ... 1
         second op-loop: index o goes from n...1 and
                         index p from o+TURN+1 ... n !!
         HERE index o goes from 1...n and index p o-TURN-1 ... 1 ,
         we calculate the contributions to multiple stem loop
         where exp(i+w-1-p)*(qqm2 values between i+w and o-1)
         AND qm[iindex[p+1]-(i-1)]*exp(beta*w)*qm[iindex[i+w]-(o-1)]
         you have to recalculate of qqm matrix containing final stem
         contributions to multiple loop partition function
         from segment p,o */

      /* recalculate qqm[]
         qqm[p] := (contribution with exact one loop in region (p,o)*/
      qqm[p]  = qqm1[p] * expMLbase[1];
      if(type){
        qbt1    =   qb[po] * exp_E_MLstem(type, (p>1) ? S1[p-1] : -1, (o<n) ? S1[o+1] : -1, Pf);
        qqm[p]  +=  qbt1;
        /* reverse dangles for prpr[po]*... */
        temp    =   0.;
        tt      =   rtype[type];
        temp    =   prpr[po] * exp_E_MLstem(tt, S1[o-1], S1[p+1], Pf) * scale[2] * Pf->expMLclosing;
        for(i=p+1; i < o; i++) {
          int p1i = (p+1) < (i-1)  ? my_iindx[p+1]-(i-1)  : 0;
          /*unpaired region expMLbase[i-p] left of structured
            region qq_1m2[i+1]*/
          /* @expMLbase:  note distance of i-p == i-(p+1)+1 */
          store_M_mlbase[my_iindx[p+1]-i] += expMLbase[i-p] * temp * qq_1m2[i+1];
          /* structured region qm[p1i] left of unpaired region */
          /* contribition for unpaired region is added after the p-loop */
          store_M_qm_o[i] += qm[p1i] * temp;
        } /*end of for i ... */
      }

      for(tqm2 = 0., i=p+1; i < o; i++)
        tqm2  +=  qm[my_iindx[p]-i] * qqm[i+1];

      /* qqm2[p] contrib with at least 2 loops in region (p,o) */
      qqm2[p] = tqm2;
    } /* end for (p=..) */

    for(sum_h = 0., i=1; i < o; i++) {
      int max_v, vo;
      sum_h +=  store_H[i];
      max_v =   MIN2(w-1,o-i-1);
      for(v=max_v; v >= 0; v--){
        /* Hairpins */
        pu_test->H[i][v] += sum_h;/* store_H[i][v] + store_H[i][max_v]; */
        /* Interior loops: unpaired region between  ]l,o[ calculated here !*/
        /* unpaired region between ]p,k[ collected after after o-loop */
        if(v <= MIN2(max_v,MAXLOOP)) {
          pu_test->I[i][v] += store_Io[i]; /* ]l,o[ */
        }
        /* Multiloops:*/
        /* unpaired region [i,v] between structured regions ]p,i[ and ]v,o[. */
        /* store_M_qm_o[i] = part. funct over all structured regions ]p,i[ */
        vo = (i+v+1) <= (o-1) ? my_iindx[i+v+1]-(o-1): 0;
        pu_test->M[i][v] += store_M_qm_o[i]*expMLbase[v+1]*qm[vo];
      }
    }
    tmp = qqm1; qqm1=qqm; qqm=tmp;
    tmp = qqm2; qqm2=qq_1m2; qq_1m2=tmp;

    free(store_Io);
    free(store_H);
    free(store_M_qm_o);
  }/* end for (o=..) */

  for(i=1; i < n; i++) {
    int     max_v;
    double  sum_iv;
    sum_iv  = 0.;
    max_v   = MIN2(w-1,n-i);
    for(v=n; v >=0; v--) {
      if(v <= MIN2(max_v,MAXLOOP)) {
        /* all unpaired regions [i,v] between p and k in interior loops */
        /* notice v runs from max_v -> 0, sum_iv sums all int. l. contribs */
        /* for each x, v < x =< max_v, since they contribute to [i,v] */
        sum_iv            += store_I2o[i][v];
        pu_test->I[i][v]  += sum_iv;
      }
      /* all unpaired region [i,v] for a fixed v, given that */
      /* region ]v,o[ contains at least 2 structures qq_1m2[v+1]; */
      if(v >= i) {
        sum_M[v] += store_M_mlbase[my_iindx[i]-v];
        if(v-i<=max_v) {
          pu_test->M[i][v-i] += sum_M[v];
        }
      }
    }
  }

  for(i=0;i<=n;i++) {
    free(store_I2o[i]);
  }
  free(store_I2o);

  for (i=1; i<=n; i++) {
    /* set auxillary arrays to 0 */
    qqm[i] = qqm1[i] = qqm2[i] = qq_1m2[i] = 0;
  }

  /* 2. exterior bp (p,o) encloses unpaired region [i,j]
     HERE index o goes from n...1 and index p from o+TURN+1 ... n,
     that is, we add the one multiloop contribution that we
     could not calculate before  */

/* is free'ing plus allocating faster than looping over all entries an setting them to 0? */
#if 0
  free(store_M_mlbase);
  store_M_mlbase = (double *) vrna_alloc(sizeof(double) * (size + 1));
#else
  /* this should be the fastest way to set everything to 0 */
  memset(store_M_mlbase, 0, sizeof(double) * (size + 1));
#endif

  for (o=n-TURN-1;o>=1; o--) {
    for (p=o+TURN+1; p<=n; p++) {
      po    = my_iindx[o]-p;
      type  = ptype[po];
      /* recalculate of qqm matrix containing final stem
         contributions to multiple loop partition function
         from segment [o,p] */
      qqm[p] = qqm1[p] * expMLbase[1];
      if (type) {
        qbt1 = qb[po];
        qbt1 *= exp_E_MLstem(type, (o>1) ? S1[o-1] : -1, (p<n) ? S1[p+1] : -1, Pf);
        qqm[p] += qbt1;
        /* revers dangles for prpr[po]...  */
        temp=0.;
        tt=rtype[type];
        temp = prpr[po]*exp_E_MLstem(tt, S1[p-1], S1[o+1], Pf) * Pf->expMLclosing * scale[2];
      }
      tqm2=0.;
      for(i=o+1; i < p; i++) {
        tqm2+=qqm[i]*qm[my_iindx[i+1]-p];

        if(type !=0) {
          /* structured region qq_1m2[i-1] left of unpaired r. expMLbase[p-i]*/
          /* @expMLbase:  note distance of p-i == p+1-i+1 */
           store_M_mlbase[my_iindx[i]-p+1] +=  qq_1m2[i-1]*expMLbase[p-i]*temp;
        }
      }/*end of for i ....*/
      qqm2[p] = tqm2;
    }/* end for (p=..) */
    tmp = qqm1; qqm1=qqm; qqm=tmp;
    tmp = qqm2; qqm2=qq_1m2; qq_1m2=tmp;
  }/* end for (o=..) */
  /* now collect the missing multiloop contributions */
  for(i=0;i<=n;i++) { sum_M[i]=0.; }
  for(i=1; i<=n;i++) {
    int v_max = MIN2(w-1,n-i);
    for(v=n; v>=i; v--){
      sum_M[i]  += store_M_mlbase[my_iindx[i]-v];
      if ((v-i <= v_max) ) {
        pu_test->M[i][v-i] += sum_M[i];
      }
    }
  }

  /* 1. region [i,j] exterior to all loops */
  for (i=1; i<=n; i++) {
    for(j=i; j<MIN2(i+w,n+1);j++){
      ij=my_iindx[i]-j;
      temp=q1k[i-1]*1*scale[j-i+1]*qln[j+1]/q1k[n];
      pu_test->E[i][j-i]+=temp;

    }
  }

  free(sum_M);
  free(store_M_mlbase);
  free_up_arrays();
  return pu_test;
}


PRIVATE void  get_interact_arrays(unsigned int n1,
                                  unsigned int n2,
                                  pu_contrib *p_c,
                                  pu_contrib *p_c2,
                                  int w,
                                  int incr5,
                                  int incr3,
                                  double ***p_c_S,
                                  double ***p_c2_S){

  unsigned int i;
  int pc_size, j;
  *p_c_S = (double **)vrna_alloc(sizeof(double *)*(n1+1));

  for (i=1; i<=n1; i++){
    pc_size = MIN2((w + incr5 + incr3), (int)n1);
    (*p_c_S)[i] = (double *)vrna_alloc(sizeof(double) * (pc_size + 1));
    for (j=0; j < pc_size; j++)
      (*p_c_S)[i][j] = p_c->H[i][j] + p_c->I[i][j] + p_c->M[i][j] + p_c->E[i][j];
  }

  if(p_c2 != NULL){
    (*p_c2_S) = (double **)vrna_alloc(sizeof(double *) * (n2 + 1));
    for (i=1; i<=n2; i++){
      pc_size = MIN2(w, (int)n2);
      (*p_c2_S)[i]  = (double *)vrna_alloc(sizeof(double) * (pc_size + 2));
      for (j=0; j < pc_size; j++)
        (*p_c2_S)[i][j] = p_c2->H[i][j] + p_c2->I[i][j] + p_c2->M[i][j] + p_c2->E[i][j];
    }
  }
}

/*------------------------------------------------------------------------*/
/* s1 is the longer seq */
PUBLIC interact *pf_interact( const char *s1,
                              const char *s2,
                              pu_contrib *p_c,
                              pu_contrib *p_c2,
                              int w,
                              char *cstruc,
                              int incr3,
                              int incr5){

  int         i, j, k,l,n1,n2,add_i5,add_i3, pc_size;
  double      temp, Z, rev_d, E, Z2,**p_c_S, **p_c2_S, int_scale;
  FLT_OR_DBL  ****qint_4, **qint_ik;
  /* PRIVATE double **pint; array for pf_up() output */
  interact    *Int;
  double      G_min, G_is,Gi_min;
  int         gi,gj,gk,gl,ci,cj,ck,cl,prev_k,prev_l;
  FLT_OR_DBL  **int_ik;
  double      Z_int, temp_int, temppfs;
  double      const_scale, const_T;
  constrain   *cc = NULL;  /* constrains for cofolding */
  char        *Seq, *i_long,*i_short,*pos=NULL; /* short seq appended to long one */
  /* int ***pu_jl; */ /* positions of interaction in the short RNA */

  G_min = G_is = Gi_min = 100.0;
  gi = gj = gk = gl = ci = cj = ck = cl = 0;

  n1      = (int) strlen(s1);
  n2      = (int) strlen(s2);
  prev_k  = 1;
  prev_l  = n2;

  i_long  = (char *) vrna_alloc(sizeof(char)*(n1+1));
  i_short = (char *) vrna_alloc(sizeof(char)*(n2+1));
  Seq     = (char *) vrna_alloc(sizeof(char)*(n1+n2+2));

  strcpy(Seq,s1);
  strcat(Seq,s2);

  set_encoded_seq(s1, &S, &S1);
  set_encoded_seq(s2, &SS, &SS2);

  cc = get_ptypes_up(Seq,cstruc);

  get_interact_arrays(n1, n2, p_c, p_c2, w, incr5, incr3, &p_c_S, &p_c2_S);

  /*array for pf_up() output */
  Int = (interact *) vrna_alloc(sizeof(interact)*1);
  Int->Pi = (double *) vrna_alloc(sizeof(double)*(n1+2));
  Int->Gi = (double *) vrna_alloc(sizeof(double)*(n1+2));

  /* use a different scaling for pf_interact*/
  scale_int(s2, s1, &int_scale);

  /* set the global scale array and the global variable pf_scale to the
     values used to scale the interaction, keep their former values !! */
  temppfs = pf_scale;
  pf_scale = int_scale;

  /* in order to scale expLoopEnergy correctly call*/
  /* we also pass twice the seq-length to avoid bogus access to scale[] array */
  scale_stru_pf_params((unsigned) 2*n1);

  qint_ik = (FLT_OR_DBL **) vrna_alloc(sizeof(FLT_OR_DBL *) * (n1+1));
  for (i=1; i<=n1; i++) {
    qint_ik[i] = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * (n1+1));
  }
/* int_ik */
  int_ik = (FLT_OR_DBL **) vrna_alloc(sizeof(FLT_OR_DBL *) * (n1+1));
  for (i=1; i<=n1; i++) {
    int_ik[i] = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * (n1+1));
  }
  Z_int=0.;
  /*  Gint = ( -log(int_ik[gk][gi])-( ((int) w/2)*log(pf_scale)) )*((Pf->temperature+K0)*GASCONST/1000.0); */
  const_scale = ((int) w/2)*log(pf_scale);
  const_T = (Pf->kT/1000.0);
  encode_seq(s1, s2);
  /* static  short *S~S1, *S1~SS1, *SS~S2, *SS2; */
  for (i=0; i<=n1; i++) {
    Int->Pi[i]=Int->Gi[i]=0.;
  }
  E=0.;
  Z=0.;

  if ( fold_constrained && cstruc != NULL) {
    pos = strchr(cstruc,'|');
    if(pos) {
      ci=ck=cl=cj=0;
      /* long seq              & short seq
         .........||..|||||....&....||||...  w = maximal interaction length
                 ck       ci       cj  cl    */
      strncpy(i_long,cstruc,n1);
      i_long[n1] = '\0';
      strncpy(i_short,&cstruc[n1],n2);
      i_short[n2] ='\0';
      pos = strchr(i_long,'|');
      if(pos) ck = (int) (pos-i_long)+1; /* k */
      pos = strrchr(i_long,'|');
      if(pos) ci = (int) (pos-i_long)+1; /* i */
      pos = strrchr(i_short,'|');
      if(pos) cl = (int) (pos-i_short)+1; /* l */
      pos = strchr(i_short,'|');
      if(pos) cj = (int) (pos-i_short)+1; /* j */

      if(ck > 0 && ci > 0 && ci-ck+1 > w) {
        vrna_message_warning("distance between constrains in longer seq, %d, larger than -w = %d",ci-ck+1,w);
        vrna_message_error("pf_interact: could not satisfy all constraints");
      }
      if(cj > 0 && cl > 0 && cl-cj+1 > w) {
        vrna_message_warning("distance between constrains in shorter seq, %d, larger than -w = %d",cl-cj+1,w);
        vrna_message_error("pf_interact: could not satisfy all constraints");
      }
    }

  } else if ( fold_constrained && cstruc == NULL) {
    vrna_message_error("option -C selected, but no constrained structure given\n");
  }
  if(fold_constrained) pos = strchr(cstruc,'|');

  /*  qint_4[i][j][k][l] contribution that region (k-i) in seq1 (l=n1)
      is paired to region (l-j) in seq 2(l=n2) that is
      a region closed by bp k-l  and bp i-j */
  qint_4 = (FLT_OR_DBL ****) vrna_alloc(sizeof(FLT_OR_DBL ***) * (n1+1));

  /* qint_4[i][j][k][l] */
  for (i=1; i<=n1; i++) {
    int end_k;
    end_k = i-w;
    if(fold_constrained && pos && ci) end_k= MAX2(i-w, ci-w);
    /* '|' constrains for long sequence: index i from 1 to n1 (5' to 3')*/
    /* interaction has to include 3' most '|' constrain, ci */
    if(fold_constrained && pos && ci && i==1 && i<ci)
      i= ci-w+1 > 1 ? ci-w+1 : 1;
    /* interaction has to include 5' most '|' constrain, ck*/
    if(fold_constrained && pos && ck && i > ck+w-1) break;

    /* note: qint_4[i] will be freed before we allocate qint_4[i+1] */
    qint_4[i] = (FLT_OR_DBL ***) vrna_alloc(sizeof(FLT_OR_DBL **) * (n2+1));
    for (j=n2; j>0; j--) {
      qint_4[i][j] = (FLT_OR_DBL **) vrna_alloc(sizeof(FLT_OR_DBL*) * (w+1));
      for (k=0; k<=w; k++) {
        qint_4[i][j][k] = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * (w+1));
      }
    }

     prev_k=1;
    for (j=n2; j>0; j--) {
      int type, type2,end_l;
      end_l = j+w;
      if(fold_constrained && pos && ci) end_l= MIN2(cj+w,j+w);
      /* '|' constrains for short sequence: index j from n2 to 1 (3' to 5')*/
      /* interaction has to include 5' most '|' constrain, cj */
      if(fold_constrained && pos && cj && j==n2 && j>cj)
        j = cj+w-1 > n2 ? n2 : cj+w-1;
      /* interaction has to include 3' most '|' constrain, cl*/
      if(fold_constrained && pos && cl && j < cl-w+1) break;
      type = cc->ptype[cc->indx[i]-(n1+j)];
      qint_4[i][j][0][0] = type ? Pf->expDuplexInit : 0;

      if (!type) continue;
      qint_4[i][j][0][0] *= exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, (j<n2) ? SS2[j+1] : -1, Pf);

      rev_d = exp_E_ExtLoop(rtype[type], (j>1) ? SS2[j-1] : -1, (i<n1) ? S1[i+1] : -1, Pf);

      /* add inc5 and incr3 */
      if((i-incr5) > 0 ) add_i5=i-incr5;
      else add_i5=1;
      add_i3=incr3;
      pc_size = MIN2((w+incr3+incr5),n1);
      if(incr3 < pc_size) add_i3=incr3;
      else add_i3=pc_size-1;

      /* only one bp (no interior loop) */
      if(p_c2 == NULL) {/* consider only structure of longer seq. */
        qint_ik[i][i]+=qint_4[i][j][0][0]*rev_d*p_c_S[add_i5][add_i3]*scale[((int) w/2)];
        Z+=qint_4[i][j][0][0]*rev_d*p_c_S[add_i5][add_i3]*scale[((int) w/2)];
      } else {/* consider structures of both seqs. */
        qint_ik[i][i]+=qint_4[i][j][0][0]*rev_d*p_c_S[add_i5][add_i3]*p_c2_S[j][0]*scale[((int) w/2)];
        Z+=qint_4[i][j][0][0]*rev_d*p_c_S[add_i5][add_i3]*p_c2_S[j][0]*scale[((int) w/2)];
      }

/* int_ik */
      /* check deltaG_ges = deltaG_int + deltaG_unstr; */
      int_ik[i][i]+=qint_4[i][j][0][0]*rev_d*scale[((int) w/2)];
      Z_int+=qint_4[i][j][0][0]*rev_d*scale[((int) w/2)];
      temp_int=0.;

      temp=0.;
      prev_l = n2;
      for (k=i-1; k>end_k && k>0; k--) {
        if (fold_constrained && pos && cstruc[k-1] == '|' && k > prev_k)
          prev_k=k;
        for (l=j+1; l< end_l && l<=n2; l++) {
          int a,b,ia,ib,isw;
          double scalew, tt, intt;

          type2 = cc->ptype[cc->indx[k]-(n1+l)];
          /* '|' : l HAS TO be paired: not pair (k,x) where x>l allowed */
          if(fold_constrained && pos && cstruc[n1+l-1] == '|' && l < prev_l)
            prev_l=l; /*break*/
          if(fold_constrained && pos && (k<=ck || i>=ci) && !type2) continue;
          if(fold_constrained && pos && ((cstruc[k-1] == '|') || (cstruc[n1+l-1] == '|')) && !type2) break;

          if (!type2) continue;
          /* to save memory keep only qint_4[i-w...i][][][] in memory
             use indices qint_4[i][j][a={0,1,...,w-1}][b={0,1,...,w-1}] */
          a=i-k;/* k -> a from 1...w-1*/
          b=l-j;/* l -> b from 1...w-1 */

          /* scale everything to w/2 */
          isw = ((int) w/2);
          if ((a+b) < isw ){
            scalew = ( scale[isw - (a+b)] );
          } else if ( (a+b) > isw ) {
            scalew = 1/( scale[(a+b) - isw] );
          } else {
            scalew = 1;
          }

          if (i-k+l-j-2<=MAXLOOP) {
            if(k >= prev_k && l <= prev_l) { /* don't violate constrains */
              E = exp_E_IntLoop(i-k-1,l-j-1, type2, rtype[type],
                                S1[k+1], SS2[l-1], S1[i-1], SS2[j+1], Pf) *
                                scale[i-k+l-j]; /* add *scale[u1+u2+2] */

              qint_4[i][j][a][b] += ( qint_4[k][l][0][0]*E);

              /* use ia and ib to go from a....w-1 and from b....w-1  */
              ia=ib=1;
              while((a+ia)<w && i-(a+ia)>=1 && (b+ib)<w && (j+b+ib)<=n2) {
                int iaa,ibb;

                qint_4[i][j][a+ia][b+ib] += qint_4[k][l][ia][ib]*E;

                iaa=ia+1;
                while(a+iaa<w && i-(a+iaa)>=1) {
                  qint_4[i][j][a+iaa][b+ib] += qint_4[k][l][iaa][ib]*E;
                  ++iaa;
                }

                ibb=ib+1;
                while( (b+ibb)<w && (j+b+ibb)<=n2 ) {
                  qint_4[i][j][a+ia][b+ibb] += qint_4[k][l][ia][ibb]*E;
                  ++ibb;
                }
                ++ia;
                ++ib;
              }
            }
          }
          /* '|' constrain in long sequence */
          /* collect interactions starting before 5' most '|' constrain */
          if ( fold_constrained && pos && ci && i < ci) continue;
          /* collect interactions ending after 3' most '|' constrain*/
          if ( fold_constrained && pos && ck &&  k > ck) continue;
          /* '|' constrain in short sequence */
          /* collect interactions starting before 5' most '|' constrain */
          if ( fold_constrained && pos && cj && j > cj) continue;
          /* collect interactions ending after 3' most '|' constrain*/
          if ( fold_constrained && pos && cl && l < cl) continue;

          /* scale everything to w/2*/
          /* qint_ik[k][i] all interactions where k and i both are paired */
          /* substract incr5 from k */
          if(k-incr5 > 0) add_i5=k-incr5;
          else add_i5=1;
          /* add incr3 to i */
          pc_size = MIN2((w+incr3+incr5),n1);
          if(i-k+incr3 < pc_size) add_i3=i-k+incr3;
          else add_i3=pc_size-1;

          if(p_c2 == NULL) {/* consider only structure of longer seq. */
            tt = qint_4[i][j][a][b]*p_c_S[add_i5][add_i3]*scalew*rev_d;
          } else { /* consider structures of both seqs. */
            tt = qint_4[i][j][a][b]*p_c_S[add_i5][add_i3]*p_c2_S[j][b]*scalew*rev_d;
          }
          temp+= tt;
          qint_ik[k][i]+= tt;
          /* int_ik */
          /* check deltaG_ges = deltaG_int + deltaG_unstr; */
          intt = qint_4[i][j][a][b]*scalew*rev_d;
          temp_int += intt;
          int_ik[k][i]+= intt;
          G_is = (-log(tt)-const_scale)*(const_T);
          if (G_is < G_min || EQUAL(G_is,G_min)) {
            G_min = G_is;
            Gi_min =(-log(intt)-const_scale)*(const_T);
            gi=i;
            gj=j;
            gk=k;
            gl=l;
          }
        }
      }
      Z+=temp;
      /* int_ik */
      Z_int+=temp_int;
    }

    /* free qint_4 values not needed any more */
    if(i > w) {
      int bla;
      bla=i-w;
      if (fold_constrained && pos && ci && i-w < ci-w+1) continue;
      if (fold_constrained && pos && ci) bla = MAX2(ci-w+1,i-w);
      for (j=n2; j>0; j--) {
        for (k=0; k<=w; k++){
          free(qint_4[bla][j][k]);
        }
        free(qint_4[bla][j]);
      }
      free(qint_4[bla]);
      qint_4[bla] = NULL;
    }
  }


  Z2=0.0;
  for (i=1; i<=n1; i++) {
    for (k=i; k<=n1 && k<i+w; k++) {
      Z2+=qint_ik[i][k];
      for(l=i;l<=k;l++) {
        /* Int->Pi[l]: prob that position l is within a paired region */
        /* qint_ik[i][k] as well as Z are scaled to scale[((int) w/2) */
        Int->Pi[l]+=qint_ik[i][k]/Z;
        /* Int->Gi[l]: minimal delta G at position [l] */
        Int->Gi[l]=MIN2(Int->Gi[l],
                       ( -log(qint_ik[i][k])-( ((int) w/2)*log(pf_scale)) )*
                       (Pf->kT/1000.0) );
      }
    }
  }
  if(n1 > w){
    int start_i,end_i;
    start_i = n1-w+1;
    end_i=n1;
    if (fold_constrained && pos && ci) {
      /* a break in the k loop might result in unfreed values */
      start_i = ci-w+1 < n1-w+1 ? ci-w+1 : n1-w+1;
      start_i = start_i > 0 ? start_i : 1;
      /* start_i = ck; */
      end_i = ck+w-1 > n1 ? n1 : ck+w-1;
    }
    for (i=start_i; i<=end_i; i++) {
      if(qint_4[i] == NULL ) continue;
      for (j=n2; j>0; j--) {
        for (k=0; k<=w; k++) {
          free(qint_4[i][j][k]);
        }
        free(qint_4[i][j]);
      }
      free(qint_4[i]);
    }
    free(qint_4);
  } else {
    int start_i,end_i;
    start_i = 1;
    end_i=n1;
    if (fold_constrained && pos) {
      start_i = ci-w+1 > 0 ? ci-w+1 : 1;
      end_i = ck+w-1 > n1 ? n1 : ck+w-1;
    }

    for (i=start_i; i<=end_i; i++) {
      for (j=n2; j>0; j--) {
        for (k=0; k<=w; k++) {
          free(qint_4[i][j][k]);
        }
        free(qint_4[i][j]);
      }
      free(qint_4[i]);
    }
    free(qint_4);
  }
  if(fold_constrained && (gi==0 || gk==0 ||  gl==0 || gj==0)) {
    vrna_message_error("pf_interact: could not satisfy all constraints");
  }
  /* fill structure interact */
  Int->length = n1;
  Int->i = gi;
  Int->j = gj;
  Int->k = gk;
  Int->l = gl;
  Int->Gikjl = G_min;
  Int->Gikjl_wo = Gi_min;

  free(i_long);
  free(i_short);

  for (i=1; i<=n1; i++) {
    free(int_ik[i]);
  }
  free(int_ik);
  for (i=1; i<=n1; i++) {
    free(qint_ik[i]);
  }
  free(qint_ik);

  /* reset the global variables pf_scale and scale to their original values */
  pf_scale = temppfs;/* reset pf_scale */
  scale_stru_pf_params((unsigned) n1);/* reset the scale array */
  free_pf_arrays(); /* for arrays for pf_fold(...) */

  if(expMLbase != NULL) {
    free(expMLbase);
    expMLbase = NULL;
  }
  if(scale != NULL) {
    free(scale);
    scale = NULL;
  }
  for (i=1; i<=n1; i++) {
    free(p_c_S[i]);
  }
  free(p_c_S);
  if(p_c2 != NULL) {
    for (i=1; i<=n2; i++) {
      free(p_c2_S[i]);
    }
    free(p_c2_S);
  }
  free(Seq);
  free(cc->indx);
  free(cc->ptype);
  free(cc);
  return(Int);
}
/*------------------------------------------------------------------------*/
/* use an extra scale for pf_interact, here sl is the longer sequence */
PRIVATE void scale_int(const char *s, const char *sl, double *sc_int){
  int       n,nl;
  duplexT   mfe;
  double    kT;

  n         = strlen(s);
  nl        = strlen(sl);

  free(expMLbase);
  free(scale);

  expMLbase = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*((nl+1)*2));
  scale     = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*((nl+1)*2));

  /* use RNA duplex to get a realistic estimate for the best possible
     interaction energy between the short RNA s and its target sl */
  mfe = duplexfold(s,sl);

  kT = Pf->kT/1000.0;   /* in Kcal */

  /* sc_int is similar to pf_scale: i.e. one time the scale */
  *sc_int = exp(-(mfe.energy)/kT/n);

  /* free the structure returned by duplexfold */
  free(mfe.structure);
}

/*----------------------------------------------------------------------*/
/* init_pf_two(n) :gets the arrays, that you need, from part_func.c */
/* get_pf_arrays(&S, &S1, &ptype, &qb, &qm, &q1k, &qln);*/
/* init_pf_fold(), update_pf_params, encode_char(), make_ptypes() are called by pf_fold() */
PRIVATE void init_pf_two(int length){
#ifdef SUN4
  nonstandard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(1);
#endif
#endif
  make_pair_matrix();

  /* gets the arrays, that we need, from part_func.c */
  if(!get_pf_arrays(&S, &S1, &ptype, &qb, &qm, &q1k, &qln))
    vrna_message_error("init_pf_two: pf_fold() has to be called before calling pf_unstru()\n");
  /* get a pointer to the base pair probs */
  probs = export_bppm();

  scale_stru_pf_params((unsigned) length);

  init_length=length;
  if(init_temp != Pf->temperature)
    vrna_message_error("init_pf_two: inconsistency with temperature");
}

PRIVATE void  get_up_arrays(unsigned int length){
  unsigned int l1 = length + 1;
  unsigned int l2 = length + 2;
  prpr      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL)  * ((l1*l2)>>1));
  expMLbase = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL)  * l2);
  scale     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL)  * l2);
  qqm2      = (double *)    vrna_alloc(sizeof(double)      * l2);
  qq_1m2    = (double *)    vrna_alloc(sizeof(double)      * l2);
  qqm       = (double *)    vrna_alloc(sizeof(double)      * l2);
  qqm1      = (double *)    vrna_alloc(sizeof(double)      * l2);
  my_iindx  = vrna_idx_row_wise(length);
}

PRIVATE void  free_up_arrays(void){
  if(prpr       != NULL){ free(prpr);       prpr      = NULL;}
  if(expMLbase  != NULL){ free(expMLbase);  expMLbase = NULL;}
  if(scale      != NULL){ free(scale);      scale     = NULL;}
  if(qqm        != NULL){ free(qqm);        qqm       = NULL;}
  if(qqm1       != NULL){ free(qqm1);       qqm1      = NULL;}
  if(qqm2       != NULL){ free(qqm2);       qqm2      = NULL;}
  if(qq_1m2     != NULL){ free(qq_1m2);     qq_1m2    = NULL;}
  if(my_iindx   != NULL){ free(my_iindx);   my_iindx  = NULL;}
}

PUBLIC void free_interact(interact *pin) {
  if(S != NULL && pin != NULL){
    free(S);
    S=NULL;
  }
  if(S1 != NULL && pin != NULL){
    free(S1);
    S1=NULL;
  }
  if(pin != NULL){
    free(pin->Pi);
    free(pin->Gi);
    free(pin);
    pin=NULL;
  }
}
/*---------------------------------------------------------------------------*/

PRIVATE void encode_seq(const char *s1, const char *s2) {
  unsigned int i,l;

  l = strlen(s1);
  /* S and S1 are freed by free_pf_arrays(); ! */
  S = (short *) vrna_alloc(sizeof(short)*(l+1));
  S1= (short *) vrna_alloc(sizeof(short)*(l+1));
  /* S1 exists only for the special X K and I bases and energy_set!=0 */
  S[0] = l;
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    S[i]= (short) encode_char(toupper(s1[i-1]));
    S1[i] = alias[S[i]];   /* for mismatches of nostandard bases */
  }
  if(s2 != NULL) {
    l = strlen(s2);
    /* SS2 exists only for the special X K and I bases and energy_set!=0 */
    SS[0] = l;
    for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
      SS[i]= (short) encode_char(toupper(s2[i-1]));
      SS2[i] = alias[SS[i]];   /* for mismatches of nostandard bases */
    }
  }
}

/*-------------------------------------------------------------------------*/
 /* scale energy parameters and pre-calculate Boltzmann weights:
  most of this is done in structure Pf see params.c,h (function:
  get_scaled_pf_parameters(), only arrays scale and expMLbase are handled here*/
PRIVATE void scale_stru_pf_params(unsigned int length)
{
  unsigned int i;
  double  kT;


  /* Do this only at the first call for get_scaled_pf_parameters()
     and/or if temperature has changed*/
  if(init_temp != temperature) {
    if(Pf) free(Pf);
    vrna_md_t   md;
    set_model_details(&md);
    Pf=vrna_exp_params(&md);
  }

  init_temp = Pf->temperature;

  kT = Pf->kT;   /* kT in cal/mol  */

   /* scaling factors (to avoid overflows) */
  if (pf_scale == -1) { /* mean energy for random sequences: 184.3*length cal */
    pf_scale = exp(-(-185+(Pf->temperature-37.)*7.27)/kT);
    if (pf_scale<1) pf_scale=1;
  }
  Pf->pf_scale = pf_scale;
  scale[0] = 1.;
  scale[1] = 1./pf_scale;
  expMLbase[0] = 1;
  expMLbase[1] = Pf->expMLbase/pf_scale;
  for (i=2; i<=length+1; i++) {
    scale[i] = scale[i/2]*scale[i-(i/2)];
    expMLbase[i] = pow(Pf->expMLbase, (double)i) * scale[i];
  }
}
/*-------------------------------------------------------------------------*/
/* make a results structure containing all u-values & the header */
PUBLIC pu_out *get_u_vals(pu_contrib *p_c, int **unpaired_values, char *select_contrib) {
  int i, j, k, l, num_u_vals,count,contribs,size,w,len;
  int S,E,H,I,M;
  int off_S, off_E, off_H, off_I, off_M;
  /* double **p_cont,**p_cont_sh, dG_u; p_u AND its contributions */
  pu_out* u_results;

  len = p_c->length;

  /* number of different -u values */
  for (num_u_vals = 0, i = 1; i <= unpaired_values[0][0]; i++) {
    j = unpaired_values[i][0];
    do num_u_vals++; while(++j <= unpaired_values[i][1]);
  }
  /* check which contributions ([-c "SHIME"] ) are desired by the user,
     set the offset for each contribution */
  contribs = 0;
  S = E = H = I = M = 0;
  off_S = off_E = off_H = off_I = off_M = 0;
  if(strchr(select_contrib, 'S')) {
    S=1;
    off_S = contribs;
    ++contribs;
  }
  if(strchr(select_contrib, 'E')) {
    E=1;
    off_E = contribs;
    ++contribs;
  }
  if(strchr(select_contrib, 'H')) {
    H=1;
    off_H = contribs;
    ++contribs;
  }
  if(strchr(select_contrib, 'I')) {
    I=1;
    off_I = contribs;
    ++contribs;
  }
  if(strchr(select_contrib, 'M')) {
    M=1;
    off_M = contribs;
    ++contribs;
  }

  if(contribs > 5) {
    vrna_message_error("get_u_vals: error with contribs!");
  }
  /* allocate the results structure */
  u_results = (pu_out *) vrna_alloc(1*sizeof(pu_out));
  u_results->len = len; /* sequence length */
  /*num_u_vals differnet -u values, contribs [-c "SHIME"] */
  u_results->u_vals = num_u_vals;
  u_results->contribs = contribs;
  /* add 1 column for position within the sequence and
     add 1 column for the free energy of interaction values */
  /* header e.g. u3I (contribution for u3 interior loops */
  size = 1 + (num_u_vals*contribs) + 1;
  u_results->header = (char **) vrna_alloc((size+1)*sizeof(char*));
  for(i=0;i<(size+1);i++){
    u_results->header[i] = (char *) vrna_alloc(10*sizeof(char));
  }
  /* different free energies for all  -u and -c combinations */
  u_results->u_values = (double**) vrna_alloc((size+1) *sizeof(double*));
  for(i=0;i<(size+1);i++){
    /* position within the sequence  */
    u_results->u_values[i] = (double*) vrna_alloc((len+3)*sizeof(double));
  }
  /* write the position within the sequence in the u_results array
     at column zerro */
  sprintf(u_results->header[0],"pos");
  for(i=0;i<=len;i++){
    /* add the position*/
    u_results->u_values[0][i] = i;
  }
  /* go over the different -u values, u_vals[] listy of different -u values*/
  for (count = k = 1; k <= unpaired_values[0][0]; k++) {
    l = unpaired_values[k][0];
    do{
      int offset; /* offset for the respective -u value (depents on the number
                   of the -u value and on the numbers of contribs */

      offset = ((count - 1) * contribs) + 1; /* first colum is the position */
      /* set the current value of -u : here we call it w */
      w = l; /* set w to the actual -u value */
      if(w > len) break; /* corr caro */
      /* make the header - look which contribitions are wanted */
      if(S) sprintf(u_results->header[offset+off_S],"u%dS",w);
      if(E) sprintf(u_results->header[offset+off_E],"u%dE",w);
      if(H) sprintf(u_results->header[offset+off_H],"u%dH",w);
      if(I) sprintf(u_results->header[offset+off_I],"u%dI",w);
      if(M) sprintf(u_results->header[offset+off_M],"u%dM",w);

      if(p_c != NULL) {
        for (i=1; i<=len; i++) { /* for each position */
          /* w goes form j to i (intervall end at i) */
          for (j=i; j < MIN2((i+w),len+1); j++) { /* for each -u value < w
                                                this is not necessay ->
                                                calculate j from i and w
                                                : (j-i+1) == w */
            double blubb;
            /* if (j-i+1) == w we have the -u = w value wanted */
            if( (j-i+1) == w && i+w-1 <= len) {
              blubb = p_c->H[i][j-i]+p_c->I[i][j-i]+p_c->M[i][j-i]+p_c->E[i][j-i];

              /* printf("len %d  blubb %.3f \n",len, blubb); */
              if(S) u_results->u_values[offset+off_S][i+w-1]+=blubb;
              if(E) u_results->u_values[offset+off_E][i+w-1]+=p_c->E[i][j-i];
              if(H) u_results->u_values[offset+off_H][i+w-1]+=p_c->H[i][j-i];
              if(I) u_results->u_values[offset+off_I][i+w-1]+=p_c->I[i][j-i];
              if(M) u_results->u_values[offset+off_M][i+w-1]+=p_c->M[i][j-i];

            }
            if(i<w && (j-i+1) != w && i+w-1 > len &&  i+w-1 < len+3) {
              if(S) u_results->u_values[offset+off_S][i+w-1]=-1;
              if(E) u_results->u_values[offset+off_E][i+w-1]=-1;
              if(H) u_results->u_values[offset+off_H][i+w-1]=-1;
              if(I) u_results->u_values[offset+off_I][i+w-1]=-1;
              if(M) u_results->u_values[offset+off_M][i+w-1]=-1;
            }
          }
        }
      } else return(NULL); /* error */
      count++;
    } while(++l <= unpaired_values[k][1]);
  }
  return(u_results); /*success*/
}
/* plot the results structure */
/* when plotting the results for the target seq we add a header */
/* when plotting the results for the interaction partner u want no header,
   set s1 to NULL to avoid plotting the header */
/* currently we plot the free energies to a file: the probability of
   being unpaired for region [i,j], p_u[i,j], is related to the free
   energy to open region [i,j], dG_u[i,j] by:
   dG_u[i,j] = -log(p_u[i,j])*(temperature+K0)*GASCONST/1000.0; */
PUBLIC int plot_free_pu_out(pu_out* res, interact *pint, char *ofile, char *head) {
  int size,s,i,len;
  double dG_u;
  char nan[4], *time, dg[11];
  FILE *wastl;
  double  kT = Pf->kT;
  wastl = fopen(ofile,"a");
  if (wastl==NULL) {
    vrna_message_warning("p_cont: can't open %s for Up_plot", ofile);
    return(0);
  }
  sprintf(dg,"dG");

  /* printf("T=%.16f \n(temperature+K0)*GASCONST/1000.0 = %.16f\n",temperature,(temperature+K0)*GASCONST/1000.0); */

  /* write the header of the output file:  */
  /*  # timestamp commandlineaufruf   */
  /*  # length and name of first sequence (target) */
  /*  # first seq */
  /*  # length and name of second sequence (interaction partner) */
  /*  # second seq */
  /* the next line is the output for the target: colums
     position in target | dG_unpaired values for target | interaction energy */
  /*  # pos   u1S   u1H  dg */
  /*  values for target */
  /* if -b was choosen: the next lines are the dG_unpaired values for
     the interaction partner */
  /*  # pos   u1S   u1H  */
  /*  values for the interaction partner */

  /* print header, if nh is zerro */
  if(head){
    time = vrna_time_stamp();
    fprintf(wastl,"# %s\n", time);
    fprintf(wastl,"%s\n",head);
  }
  fprintf(wastl,"# ");
  /* }  else { fprintf(wastl," "); } close if before  */
  len  = res->len;
  size = res->u_vals * res->contribs;

  sprintf(nan,"NA");
  nan[2] = '\0';

  for(i=0;i<=len; i++) {
    for(s=0;s<=size+1;s++) { /* that is for different contribution */
      if ( i== 0 && s > size && pint != NULL)
        fprintf(wastl,"%8s  ",dg);
      if(i != 0) {
        if(s>0 && s<=size) {
          if(res->u_values[s][i] > 0.0) {
            dG_u = -log(res->u_values[s][i])*kT/1000.0;
            fprintf(wastl,"%8.3f  ",dG_u);
          } else { /* no p_u value was defined print nan*/
            fprintf(wastl,"%8s  ",nan);
          }

        } else if (s > size && pint != NULL) {
          fprintf(wastl,"%8.3f  ",pint->Gi[i]);
        } else if (s == 0) {
          fprintf(wastl,"%8.0f  ",res->u_values[s][i]);
        }
      } else {
        if(s>1) {
          fprintf(wastl,"%8s  ",res->header[s]);
        } else {
          fprintf(wastl,"%7s  ",res->header[s]);
        }
      }
    }
    fprintf(wastl,"\n");
  }
  fclose(wastl);
  /*free pu_out* res */
  if(res != NULL) {
    for(i=0;i<=(size+2);i++) {
      free(res->u_values[i]);
      free(res->header[i]);
    }
    free(res->u_values);
    free(res->header);
    free(res);
    res = NULL;
  }

  return(1); /* success */
}

PUBLIC int Up_plot(pu_contrib *p_c, pu_contrib *p_c_sh, interact *pint, char *ofile, int **unpaired_values, char *select_contrib, char *head, unsigned int mode) {
  pu_out *dada;
  int ret;
  /* check what case we have */

  /* upmode = 1 only one seq */
  /* if(p_c != NULL && pint == NULL) { */
  if(mode & RNA_UP_MODE_1){
    dada = get_u_vals(p_c,unpaired_values,select_contrib);
    ret = plot_free_pu_out(dada,NULL,ofile,head);

  /* upmode > 1 cofolding */
  /* } else if (p_c != NULL && pint != NULL) { */
  } else if(mode & RNA_UP_MODE_2) {
    dada = get_u_vals(p_c,unpaired_values,select_contrib);
    ret = plot_free_pu_out(dada,pint,ofile,head);

  /* upmode = 3  cofolding*/
  /* } else if (p_c == NULL && p_c_sh != NULL) { */
  }
  if(mode & RNA_UP_MODE_3) {
    dada  = get_u_vals(p_c,unpaired_values, select_contrib);
    ret   = plot_free_pu_out(dada, pint, ofile, head);

    dada = get_u_vals(p_c_sh, unpaired_values, select_contrib);
    ret = plot_free_pu_out(dada,NULL,ofile, NULL);
  }
  return(ret);
}

/*-------------------------------------------------------------------------*/
/* copy from part_func_co.c */
PRIVATE constrain *get_ptypes_up(char *Seq, const char *structure) {
  int n,i,j,k,l, length;
  constrain *con;
  short *s, *s1;

  length = strlen(Seq);
  make_pair_matrix();
  con = (constrain *) vrna_alloc(sizeof(constrain));
  con->indx = (int *) vrna_alloc(sizeof(int)*(length+1));
  for (i=1; i<=length; i++) {
    con->indx[i] = ((length+1-i)*(length-i))/2 +length+1;
  }
  con->ptype = (char *) vrna_alloc(sizeof(char)*((length+1)*(length+2)/2));

  set_encoded_seq((const char *)Seq, &s, &s1);

  n=s[0];
  for (k=1; k<=n-CO_TURN-1; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+CO_TURN+l; if (j>n) continue;
      type = pair[s[i]][s[j]];
      while ((i>=1)&&(j<=n)) {
        if ((i>1)&&(j<n)) ntype = pair[s[i-1]][s[j+1]];
        if (noLonelyPairs && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */
        con->ptype[con->indx[i]-j] = (char) type;
        otype =  type;
        type  = ntype;
        i--; j++;
      }
    }

  if (fold_constrained&&(structure!=NULL)) {
    int hx, *stack;
    char type;
    stack = (int *) vrna_alloc(sizeof(int)*(n+1));
    for(hx=0, j=1; j<=n; j++) {
      switch (structure[j-1]) {
      case 'x': /* can't pair */
        for (l=1; l<j-CO_TURN; l++) con->ptype[con->indx[l]-j] = 0;
        for (l=j+CO_TURN+1; l<=n; l++) con->ptype[con->indx[j]-l] = 0;
        break;
      case '(':
        stack[hx++]=j;
        /* fallthrough */
      case '<': /* pairs upstream */
        break;
      case ')':
        if (hx<=0) {
          vrna_message_error("1. unbalanced brackets in constraints\n%s", structure);
        }
        i = stack[--hx];
        type = con->ptype[con->indx[i]-j];
        /* don't allow pairs i<k<j<l */
        for (k=i; k<=j; k++)
          for (l=j; l<=n; l++) con->ptype[con->indx[k]-l] = 0;
        /* don't allow pairs k<i<l<j */
        for (k=1; k<=i; k++)
          for (l=i; l<=j; l++) con->ptype[con->indx[k]-l] = 0;
        con->ptype[con->indx[i]-j] = (type==0)?7:type;
      case '>': /* pairs downstream */
        break;
      }
    }
    if (hx!=0) {
      vrna_message_error("2. unbalanced brackets in constraint string\n%s", structure);
    }
    free(stack);
  }
  free(s);
  free(s1);
  return con;
}
PRIVATE  void  set_encoded_seq(const char *sequence, short **S, short **S1){
  unsigned int i,l;
  l = strlen(sequence);
  if(S!= NULL){
    *S  = (short *)vrna_alloc(sizeof(short) * (l + 2));
    for(i=1; i<=l; i++) /* make numerical encoding of sequence */
      (*S)[i]= (short) encode_char(toupper(sequence[i-1]));
    (*S)[l+1] = (*S)[1];
    (*S)[0]   = (short) l;
  }
  /* S1 exists only for the special X K and I bases and energy_set!=0 */
  if(S1 != NULL){
    *S1 = (short *)vrna_alloc(sizeof(short) * (l + 2));
    for(i=1; i<=l; i++) /* make numerical encoding of sequence */
      (*S1)[i]  = alias[(short) encode_char(toupper(sequence[i-1]))]; /* for mismatches of nostandard bases */
    /* for circular folding add first base at position n+1 and last base at position 0 in S1 */
    (*S1)[l+1]  = (*S1)[1];
    (*S1)[0]    = (*S1)[l];
  }
}
