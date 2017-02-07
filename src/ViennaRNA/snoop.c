/*                
           compute the duplex structure of two RNA strands,
                allowing only inter-strand base pairs.
         see cofold() for computing hybrid structures without
                             restriction.

                             Ivo Hofacker
                          Vienna RNA package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/snofold.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/snoop.h"
#include "ViennaRNA/PS_dot.h"
/* #include "ViennaRNA/fold.h" */
#include "ViennaRNA/duplex.h"
#include "ViennaRNA/loop_energies.h"


#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */



PRIVATE void  encode_seqs(const char *s1, const char *s2);
PRIVATE short *encode_seq(const char *seq);

PRIVATE void find_max_snoop(const char *s1, const char *s2, const int max, 
                            const int alignment_length, const int* position, 
                            const int delta, const int distance,  const int penalty, 
                            const int threshloop, const int threshLE, const int threshRE, 
                            const int threshDE, const int threshTE, const int threshSE, const int threshD,
                            const int half_stem, const int max_half_stem,  const int min_s2, 
                            const int max_s2, const int min_s1, const int  max_s1, const int min_d1, const int min_d2, const char* name, const int fullStemEnergy);

PRIVATE void find_max_snoop_XS(const char *s1, const char *s2, const int **access_s1, const int max, 
                               const int alignment_length, const int* position, const int *position_j,
                               const int delta, const int distance,  const int penalty, 
                               const int threshloop, const int threshLE, const int threshRE, 
                               const int threshDE, const int threshTE, const int threshSE, const int threshD,
                               const int half_stem, const int max_half_stem,  const int min_s2, 
                               const int max_s2, const int min_s1, const int  max_s1, const int min_d1, const int min_d2, const char *name, const int fullStemEnergy);





PRIVATE char * alisnoop_backtrack(int i, int j, const char ** s2, 
                                  int* Duplex_El, int* Duplex_Er, int* Loop_E, int *Loop_D, int *u, 
                                  int *pscd, int *psct, int *pscg,
                                  const int penalty, const int threshloop, 
                                  const int threshLE, const int threshRE, const int threshDE, const int threshD,
                                  const int half_stem, const int max_half_stem, 
                                  const int min_s2, const int max_s2, const int min_s1, 
                                  const int max_s1, const int min_d1, const int min_d2,
                                  const short **S1, const short **S2);

   
PRIVATE char * snoop_backtrack(int i, int j, const char* s2, int* Duplex_El, int* Duplex_Er, int* Loop_E, int *Loop_D, int *u, 
                               const int penalty, const int threshloop, const int threshLE, const int threshRE, const int threshDE, 
                               const int threshD,
                               const int half_stem, const int max_half_stem, 
                               const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2);

PRIVATE char * snoop_backtrack_XS(int i, int j, const char* s2, int* Duplex_El, int* Duplex_Er, int* Loop_E, int *Loop_D, int *u, 
                               const int penalty, const int threshloop, const int threshLE, const int threshRE, const int threshDE, 
                               const int threshD,
                               const int half_stem, const int max_half_stem, 
                               const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2);




PRIVATE int compare(const void *sub1, const void *sub2);
PRIVATE int covscore(const int *types, int n_seq);
PRIVATE short * aliencode_seq(const char *sequence);

PUBLIC  int snoop_subopt_sorted=0; /* from subopt.c, default 0 */


/*@unused@*/




#define MAXLOOP_L        3
#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define MAX2(A, B)      ((A) > (B) ? (A) : (B))
#define ASS                1
PRIVATE vrna_param_t *P = NULL;

PRIVATE int   **c = NULL;      /* energy array, given that i-j pair */
PRIVATE int   **r = NULL;
PRIVATE int   **lc = NULL;      /* energy array, given that i-j pair */
PRIVATE int   **lr = NULL;
PRIVATE int   **c_fill = NULL;
PRIVATE int   **r_fill = NULL;
PRIVATE int   **lpair = NULL;


PRIVATE short  *S1 = NULL, *SS1 = NULL, *S2 = NULL, *SS2 = NULL;
PRIVATE short *S1_fill = NULL, *SS1_fill = NULL, *S2_fill = NULL, *SS2_fill = NULL;
PRIVATE int   n1,n2;    /* sequence lengths */

extern int cut_point;

PRIVATE int delay_free=0;
/*--------------------------------------------------------------------------*/

snoopT alisnoopfold(const char **s1, const char **s2, 
                    const int penalty, const int threshloop, 
                    const int threshLE, const int threshRE, const int threshDE, const int threshD,
                    const int half_stem, const int max_half_stem, 
                    const int min_s2, const int max_s2, const int min_s1, 
                    const int max_s1, const int min_d1, const int min_d2) {
  
  int s,n_seq;
  int i, j, E, l1,Emin=INF, i_min=0, j_min=0;
  char *struc;
  snoopT mfe;
  int *indx;
  int *mLoop;
  int *cLoop;
  folden  **foldlist; folden **foldlist_XS;
  int Duplex_El, Duplex_Er,pscd,psct,pscg;
  int Loop_D;
  int u;
  int Loop_E;
  short **Sali1,**Sali2;
  int *type,*type2,*type3;
  vrna_md_t md;
  Duplex_El=0;Duplex_Er=0;Loop_E=0; Loop_D=0;pscd=0;psct=0;pscg=0;
  snoexport_fold_arrays(&indx, &mLoop, &cLoop,&foldlist, &foldlist_XS); 
  n1 = (int) strlen(s1[0]);
  n2 = (int) strlen(s2[0]);
  
  for (s=0; s1[s]!=NULL; s++);
  n_seq = s;
  for (s=0; s2[s]!=NULL; s++);
  if (n_seq != s) vrna_message_error("unequal number of sequences in aliduplexfold()\n");
  
  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    snoupdate_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  
  c = (int **) vrna_alloc(sizeof(int *) * (n1+1));
  r = (int **) vrna_alloc(sizeof(int *) * (n1+1));
  for (i=0; i<=n1; i++) {
          c[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        r[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        for(j=n2; j>-1; j--){
                c[i][j]=INF;
                r[i][j]=INF;
        }
  }
  Sali1 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  Sali2 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if ((int)strlen(s1[s]) != n1) vrna_message_error("uneqal seqence lengths");
    if ((int)strlen(s2[s]) != n2) vrna_message_error("uneqal seqence lengths");
    Sali1[s] = aliencode_seq(s1[s]);
    Sali2[s] = aliencode_seq(s2[s]);
  }
  type = (int *) vrna_alloc(n_seq*sizeof(int));
  type2 = (int *) vrna_alloc(n_seq*sizeof(int));
  type3 = (int *) vrna_alloc(n_seq*sizeof(int));
  /*   encode_seqs(s1, s2); */
  for (i=6; i<=n1-5; i++) {
    int U; U=0;
    for (s=0; s<n_seq; s++) {
      U+=Sali1[s][i-2];
    }
    U = (U==(n_seq)*4?1:0);
    for (j=n2-min_d2; j>min_d1; j--) {
      int type4, k,l,psc,psc2,psc3;
      for (s=0; s<n_seq; s++) {
        type[s] = pair[Sali1[s][i]][Sali2[s][j]];
      }
      psc = covscore(type, n_seq);
      for (s=0; s<n_seq; s++) if (type[s]==0) type[s]=7;
      c[i][j] = (psc>=MINPSCORE) ? (n_seq*P->DuplexInit) : INF;
      if (psc<MINPSCORE) continue;
      if(/*  pair[Sali1[i+1]][Sali2[j-1]] &&  */
         U && j < max_s1 && j > min_s1 &&  
         j > n2 - max_s2 - max_half_stem && 
         j < n2 -min_s2 -half_stem ) { /*constraint on s2 and i*/
        folden *temp;
        temp=foldlist[j+1];
        while(temp->next){
          int k = temp->k;
          for (s=0; s<n_seq; s++) {
            type2[s]= pair[Sali1[s][i-3]][Sali2[s][k+1]];
            type3[s]= pair[Sali1[s][i-4]][Sali2[s][k+1]];
          }
          psc2 = covscore(type2, n_seq);
          psc3 = covscore(type3, n_seq);
          if(psc2 > MINPSCORE){
            r[i][j]=MIN2(r[i][j],c[i-3][k+1]+temp->energy);
          }
          if(psc3 > MINPSCORE){
            r[i][j]=MIN2(r[i][j],c[i-4][k+1]+temp->energy);
          }
          temp=temp->next;
        }
      }
      /* dangle 5'SIDE relative to the mRNA  */
      for (s=0; s<n_seq; s++) {
        c[i][j] += E_ExtLoop(type[s], Sali1[s][i-1],Sali2[s][j+1],P);
      }
      for (k=i-1; k>0 && (i-k)<MAXLOOP_L; k--) {
        for (l=j+1; l<=n2 ; l++) {
          if (i-k+l-j>2*MAXLOOP_L-2) break;
          if (abs(i-k-l+j) >= ASS ) continue;
          for (E=s=0; s<n_seq; s++) { 
            type4 = pair[Sali1[s][k]][Sali2[s][l]];
            if (type4==0) type4=7;
            E += E_IntLoop(i-k-1, l-j-1, type4, rtype[type[s]],
                           Sali1[s][k+1], Sali2[s][l-1], Sali1[s][i-1], Sali2[s][j+1],P);
          }
          c[i][j] = MIN2(c[i][j], c[k][l] + E);
          r[i][j] = MIN2(r[i][j], r[k][l] + E);
        }
      }
      c[i][j]-=psc;
      r[i][j]-=psc;
      E = r[i][j]; 
      for (s=0; s<n_seq; s++) {
        E+= E_ExtLoop(rtype[type[s]], Sali2[s][j-1], Sali1[s][i+1], P);
        /**
        *** if (i<n1) E += P->dangle3[rtype[type[s]]][Sali1[s][i+1]];
        *** if (j>1)  E += P->dangle5[rtype[type[s]]][Sali2[s][j-1]];
        *** if (type[s]>2) E += P->TerminalAU;
        **/
      }
      if (E<Emin) {
        Emin=E; i_min=i; j_min=j;
      } 
    }
  }
  if(Emin > 0){
          printf("no target found under the constraints chosen\n");
        for (i=0; i<=n1; i++) {free(r[i]);free(c[i]);}
        free(c);
        free(r);
        for(s=0; s<n_seq;s++){
          free(Sali1[s]);
          free(Sali2[s]);
        }
        free(Sali1); free(Sali2);
        free(S2); free(SS1); free(SS2);free(type);free(type2);free(type3);
        mfe.energy=INF;
        mfe.structure=NULL;
        return mfe;
  }
  struc = alisnoop_backtrack(i_min, j_min,(const char**) s2, 
                             &Duplex_El, &Duplex_Er, &Loop_E, 
                             &Loop_D, &u, &pscd, &psct, &pscg,
                             penalty, threshloop, threshLE, 
                             threshRE,threshDE, threshD,
                             half_stem, max_half_stem, min_s2, 
                             max_s2, min_s1, max_s1, min_d1, 
                             min_d2,(const short**) Sali1,(const short**) Sali2);
  /* if (i_min<n1-5) i_min++; */
  /* if (j_min>6 ) j_min--; */
  l1 = strchr(struc, '&')-struc;
  mfe.i = i_min-5;
  mfe.j = j_min-5;
  mfe.u = u -5;
  mfe.Duplex_Er = (float) Duplex_Er/100;
  mfe.Duplex_El = (float) Duplex_El/100;
  mfe.Loop_D = (float) Loop_D/100;
  mfe.Loop_E = (float) Loop_E/100;
  mfe.energy = (float) Emin/100 ;
  /* mfe.fullStemEnergy = (float) fullStemEnergy/100; */
  mfe.pscd = pscd;
  mfe.psct = psct;
  mfe.structure = struc;
  for(s=0; s<n_seq;s++){
    free(Sali1[s]);free(Sali2[s]);
  }
  free(Sali1);free(Sali2);free(type);free(type2);free(type3);

  if (!delay_free) {
    for (i=0; i<=n1; i++) {free(r[i]);free(c[i]);}
    free(c);
    free(r);
    free(S2); free(SS1); free(SS2);
  }
  return mfe;
}

PUBLIC snoopT *alisnoop_subopt(const char **s1, const char **s2, int delta, int w, 
                              const int penalty, const int threshloop, 
                               const int threshLE, const int threshRE, const int threshDE, const int threshTE, const int threshSE, const int threshD,
                              const int distance, const int half_stem, const int max_half_stem,
                               const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2) {


  short **Sali1, **Sali2;
  /* printf("%d %d\n", min_s2, max_s2); */
  int i,j,s,n_seq, n1, n2, E, n_subopt=0, n_max;
  char *struc;
  snoopT mfe;
  snoopT *subopt;
  int thresh;
  int *type;
  int Duplex_El, Duplex_Er, Loop_E,pscd,psct,pscg;
  int Loop_D;
  Duplex_El=0; Duplex_Er=0; Loop_E=0;Loop_D=0;pscd=0;psct=0;pscg=0;
  int u;
  u=0;
  n_max=16;
  subopt = (snoopT *) vrna_alloc(n_max*sizeof(snoopT));
  delay_free=1;
  mfe = alisnoopfold(s1, s2, penalty, threshloop, threshLE, threshRE, threshDE,threshD,
                  half_stem, max_half_stem,
                     min_s2, max_s2, min_s1, max_s1, min_d1, min_d2);
  if(mfe.energy > 0){
          free(subopt);
        delay_free=0;
        return NULL;
  }
  thresh = MIN2((int) ((mfe.Duplex_Er + mfe.Duplex_El + mfe.Loop_E)*100+0.1 + 410) + delta, threshTE );
 /* subopt[n_subopt++]=mfe; */
  free(mfe.structure);
  n1 = (int)strlen(s1[0]);
  n2 = (int)strlen(s2[0]);
  for (s=0; s1[s]!=NULL; s++);
  n_seq = s;
  Sali1 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  Sali2 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if ((int)strlen(s1[s]) != n1) vrna_message_error("uneqal seqence lengths");
    if ((int)strlen(s2[s]) != n2) vrna_message_error("uneqal seqence lengths");
    Sali1[s] = aliencode_seq(s1[s]);
    Sali2[s] = aliencode_seq(s2[s]);
  }
  Sali1[n_seq]=NULL;  Sali2[n_seq]=NULL;
  type = (int *) vrna_alloc(n_seq*sizeof(int));
  for (i=n1; i>1; i--){
    for (j=1; j<=n2; j++) {
      int  ii,jj, Ed,psc,skip;
      for (s=0; s<n_seq; s++) {
        type[s] = pair[Sali2[s][j]][Sali1[s][i]];
      }
      psc = covscore(type, n_seq);
      for (s=0; s<n_seq; s++) if (type[s]==0) type[s]=7;
      if (psc<MINPSCORE) continue;
      E = Ed = r[i][j];
      for  (s=0; s<n_seq; s++) {
        /*         if (i<n1-5) Ed += P->dangle3[type[s]][Sali1[s][i+1]]; */
        /*       if (j>6)  Ed += P->dangle5[type[s]][Sali2[s][j-1]]; */
        if (type[s]>2) Ed += P->TerminalAU;
      }
      if (Ed>thresh) continue;
      /* too keep output small, remove hits that are dominated by a
         better one close (w) by. For simplicity we do test without
         adding dangles, which is slightly inaccurate. 
      */ 
      w=1;
      for (skip=0, ii=MAX2(i-w,1); (ii<=MIN2(i+w,n1)) && type; ii++) { 
        for (jj=MAX2(j-w,1); jj<=MIN2(j+w,n2); jj++)
          if (r[ii][jj]<E) {skip=1; break;}
      }
      if (skip){continue;}
      psct=0;
      pscg=0;
      struc = alisnoop_backtrack(i,j,s2, &Duplex_El, 
                                 &Duplex_Er, &Loop_E, &Loop_D, &u, &pscd, &psct,&pscg, 
                                 penalty, threshloop,threshLE,threshRE,threshDE, threshD,
                                 half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2,(const short int**) Sali1,(const int short **) Sali2);
              
      if (Duplex_Er > threshRE || Duplex_El > threshLE || Loop_D > threshD ||
         (Duplex_Er + Duplex_El) > threshDE || 
         (Duplex_Er + Duplex_El + Loop_E) > threshTE ||
         (Duplex_Er + Duplex_El + Loop_E + Loop_D + 410) > threshSE) {
                 /* printf(" Duplex_Er %d threshRE %d Duplex_El %d threshLE %d \n" */
                /*        " Duplex_Er + Duplex_El %d  threshDE %d \n" */
                /*        " Duplex_Er + Duplex_El + Loop_E %d  threshTE %d \n" */
                /*        " Duplex_Er + Duplex_El + Loop_E + Loop_D %d  threshSE %d \n",  */
                /*          Duplex_Er , threshRE , Duplex_El ,threshLE, */
                /*          Duplex_Er + Duplex_El, threshDE, */
                /*          Duplex_Er + Duplex_El+  Loop_E , threshTE, */
                /*          Duplex_Er + Duplex_El+  Loop_E + Loop_D, threshSE);  */
                 Duplex_Er=0; 
                Duplex_El=0;
                Loop_E = 0;
                Loop_D = 0;
                u=0,
                free(struc);
                continue;
        }

      if (n_subopt+1>=n_max) {
        n_max *= 2;
        subopt = (snoopT *) vrna_realloc(subopt, n_max*sizeof(snoopT));
      }
      
      subopt[n_subopt].i = i-5;
      subopt[n_subopt].j = j-5;
      subopt[n_subopt].u = u-5;
      subopt[n_subopt].Duplex_Er = Duplex_Er * 0.01;
      subopt[n_subopt].Duplex_El = Duplex_El * 0.01;
      subopt[n_subopt].Loop_E = Loop_E * 0.01;
      subopt[n_subopt].Loop_D = Loop_D * 0.01;
      subopt[n_subopt].energy = (Duplex_Er +Duplex_El + Loop_E + Loop_D + 410) * 0.01 ;
      subopt[n_subopt].pscd = pscd * 0.01;
      subopt[n_subopt].psct = -psct * 0.01;
      subopt[n_subopt++].structure = struc;
      
      /* i=u; */
      Duplex_Er=0; Duplex_El=0; Loop_E=0; Loop_D=0;u=0;pscd=0;psct=0;
    }
  }
  
  for (i=0; i<=n1; i++) {free(c[i]);free(r[i]);}
  free(c);free(r);
  for (s=0; s<n_seq; s++) {
    free(Sali1[s]); free(Sali2[s]);
  }
  free(Sali1); free(Sali2); free(type);
  
  if (snoop_subopt_sorted) 
    qsort(subopt, n_subopt, sizeof(snoopT), compare);
  subopt[n_subopt].i =0;
  subopt[n_subopt].j =0;
  subopt[n_subopt].structure = NULL;
  return subopt;
}







PRIVATE char *alisnoop_backtrack(int i, int j, const char ** snoseq, int *Duplex_El, 
                                 int *Duplex_Er, int *Loop_E, int *Loop_D, int *u, 
                                 int *pscd, int *psct, int *pscg,
                                 const int penalty, const int threshloop, const int threshLE, 
                                 const int threshRE, const int threshDE, const int threshD, const int half_stem, 
                                 const int max_half_stem, 
                                 const int min_s2, const int max_s2, const int min_s1, 
                                 const int max_s1, 
                                 const int min_d1, const int min_d2,const short **Sali1, const short **Sali2) {
  /* backtrack structure going backwards from i, and forwards from j 
     return structure in bracket notation with & as separator */
  int k, l, *type,*type2,*type3,type4, E, traced, i0, j0,s,n_seq,psc;
  int traced_r=0; /* flag for following backtrack in c or r */
  char *st1, *st2, *struc;
  char *struc_loop;
  n1 = (int) Sali1[0][0];
  n2 = (int) Sali2[0][0];
  
  for (s=0; Sali1[s]!=NULL; s++);
  n_seq = s;
  for (s=0; Sali2[s]!=NULL; s++);
  if (n_seq != s) vrna_message_error("unequal number of sequences in alibacktrack()\n");
  
  st1 = (char *) vrna_alloc(sizeof(char)*(n1+1));
  st2 = (char *) vrna_alloc(sizeof(char)*(n2+1));
  type = (int *) vrna_alloc(n_seq*sizeof(int));
  type2 = (int *) vrna_alloc(n_seq*sizeof(int));
  type3 = (int *) vrna_alloc(n_seq*sizeof(int));
  int *indx;
  int *mLoop;
  int *cLoop;
  folden **foldlist, **foldlist_XS;
  snoexport_fold_arrays(&indx, &mLoop, &cLoop,&foldlist, &foldlist_XS ); 
  i0=i; j0=j; /* MIN2(i+1,n1); j0=MAX2(j-1,1);!modified */
  for (s=0; s<n_seq; s++) {
    type[s] = pair[Sali1[s][i]][Sali2[s][j]];
    if(type[s]==0) type[s] = 7;
    *Duplex_Er += E_ExtLoop(rtype[type[s]], (j>1) ? Sali2[s][j-1] : -1, (i<n1) ? Sali1[s][i+1] : -1, P);
    /**
    *** if (i<n1)   *Duplex_Er += P->dangle3[rtype[type[s]]][Sali1[s][i+1]];
    *** if (j>1)    *Duplex_Er += P->dangle5[rtype[type[s]]][Sali2[s][j-1]];
    *** if (type[s]>2) *Duplex_Er += P->TerminalAU;
    **/
  }
  while (i>0 && j<=n2-min_d2 ) {
    if(!traced_r) {
      E = r[i][j]; traced=0;
      st1[i-1] = '<';
      st2[j-1] = '>'; 
      for (s=0; s<n_seq; s++) {
        type[s] = pair[Sali1[s][i]][Sali2[s][j]];
      }
      psc = covscore(type,n_seq);
      for (s=0; s<n_seq; s++) if (type[s]==0) type[s] = 7;
      E += psc;
      *pscd +=psc;
      for (k=i-1; k>0 && (i-k)<MAXLOOP_L; k--) {
        for (l=j+1; l<=n2 ; l++) {
          int LE;
          if (i-k+l-j>2*MAXLOOP_L-2) break;
          if (abs(i-k-l+j) >= ASS) continue;
          for (s=LE=0; s<n_seq; s++) {
            type4 = pair[Sali1[s][k]][Sali2[s][l]];
            if (type4==0) type4=7;
            LE += E_IntLoop(i-k-1, l-j-1, type4, rtype[type[s]], Sali1[s][k+1], Sali2[s][l-1], Sali1[s][i-1], Sali2[s][j+1],P);
          }
          if (E == r[k][l]+LE) {
            traced=1; 
            i=k; j=l;
            *Duplex_Er+=LE;
            break;
          }
        }
        if (traced) break;
      }
      if(!traced){
        int U=0;
        for (s=0; s<n_seq; s++) {
          U+=Sali1[s][i-2];
        }
        U = (U==(n_seq)*4?1:0);
        if(/*  pair[Sali1[i+1]][Sali2[j-1]] && */   /* only U's are allowed */
            U && j < max_s1 && j > min_s1 && 
            j > n2 - max_s2 - max_half_stem && 
            j < n2 -min_s2 -half_stem ) {
          int min_k, max_k;
          max_k = MIN2(n2-min_s2,j+max_half_stem+1);
          min_k = MAX2(j+half_stem+1, n2-max_s2);
          folden * temp;
          temp=foldlist[j+1];
          while(temp->next) {
            int psc2, psc3;
            int k = temp->k;
            for (s=0; s<n_seq; s++) {
              type2[s]= pair[Sali1[s][i-3]][Sali2[s][k+1]];
              type3[s]= pair[Sali1[s][i-4]][Sali2[s][k+1]];
            }
            psc2 = covscore(type2, n_seq);
            psc3 = covscore(type3, n_seq);
            if(psc2>MINPSCORE /*&& pair[Sali1[i-4]][Sali2[k+2]]*/    ){  /* introduce structure from RNAfold */
              if(E==c[i-3][k+1]+temp->energy){
                *Loop_E=temp->energy;
                st1[i-3]= '|';
                *u=i-2;
                int a,b;
                /* int fix_ij=indx[k-1+1]+j+1; */
                for(a=0; a< MISMATCH ;a++){
                  for(b=0; b< MISMATCH ; b++){
                    int ij=indx[k-1-a+1]+j+1+b;
                    if(cLoop[ij]==temp->energy) {
                      /* int bla; */
                      struc_loop=alisnobacktrack_fold_from_pair(snoseq, j+1+b, k-a-1+1,psct);
                    a=INF; b=INF;        
                    }
                  }
                }
                traced=1;
                traced_r=1;
                i=i-3;j=k+1;
                break;
              }
            }
             if (psc3>MINPSCORE  /*&& pair[Sali1[i-5]][Sali2[k+2]]*/){  /* introduce structure from RNAfold */
              if(E==c[i-4][k+1]+temp->energy){
                *Loop_E=temp->energy;
                st1[i-3]= '|';
                *u=i-2;
                int a,b;
                /* int fix_ij=indx[k-1+1]+j+1; */
                for(a=0; a< MISMATCH ;a++){
                  for(b=0; b< MISMATCH ; b++){
                    int ij=indx[k-1-a+1]+j+1+b;
                    if(cLoop[ij]==temp->energy) {
                      /* int bla; */
                      struc_loop=alisnobacktrack_fold_from_pair(snoseq, j+1+b, k-a-1+1,psct);
                      a=INF; b=INF;        
                    }
                  }
                }
                traced=1;
                traced_r=1;
                i=i-4;j=k+1;
                break;
              }
            } /* else if */
            temp=temp->next;
          } /* while temp-> next */
        } /* test on j  */
      }/* traced? */
    }/* traced_r? */
    else{
      E = c[i][j]; traced=0;
      st1[i-1] = '<';
      st2[j-1] = '>'; 
      for (s=0; s<n_seq; s++) {
        type[s] = pair[Sali1[s][i]][Sali2[s][j]];
      }
      psc = covscore(type,n_seq);
      for (s=0; s<n_seq; s++) if (type[s]==0) type[s] = 7;
      E += psc;
      *pscd+=psc;
      if (!type) vrna_message_error("backtrack failed in fold duplex c");
      for (k=i-1; (i-k)<MAXLOOP_L; k--) {
        for (l=j+1; l<=n2; l++) {
          int LE;
          if (i-k+l-j>2*MAXLOOP_L-2) break;
          if (abs(i-k-l+j) >= ASS) continue;
          for (s=LE=0; s<n_seq; s++) {
            type4 = pair[Sali1[s][k]][Sali2[s][l]];
            if (type4==0) type4=7;
            LE += E_IntLoop(i-k-1, l-j-1, type4, rtype[type[s]], Sali1[s][k+1], Sali2[s][l-1], Sali1[s][i-1], Sali2[s][j+1],P);
          }
          if (E == c[k][l]+LE) {
            traced=1; 
            i=k; j=l;
            *Duplex_El+=LE;
            break;
          }
        }
        if (traced) break;
      }
    }
    if (!traced) { 
      for (s=0; s<n_seq; s++) {
        int correction;
        correction = E_ExtLoop(type[s], (i>1) ? Sali1[s][i-1] : -1, (j<n2) ? Sali2[s][j+1] : -1, P);
        *Duplex_El += correction;
        E          -= correction;
        /**
        *** if (i>1)    {E -= P->dangle5[type[s]][Sali1[s][i-1]]; *Duplex_El +=P->dangle5[type[s]][Sali1[s][i-1]];}
        *** if (j<n2)   {E -= P->dangle3[type[s]][Sali2[s][j+1]]; *Duplex_El +=P->dangle3[type[s]][Sali2[s][j+1]];}
        *** if (type[s]>2) {E -= P->TerminalAU;                      *Duplex_El +=P->TerminalAU;}
        **/
      }
      if (E != n_seq * P->DuplexInit) {
        vrna_message_error("backtrack failed in fold duplex end");
      } else break;
    }
  }
/*   if (i>1)  i--;  */
/*   if (j<n2) j++;  */
  /* struc = (char *) vrna_alloc(i0-i+1+j-j0+1+2); */ /* declare final duplex structure */
  struc = (char *) vrna_alloc(i0-i+1+n2-1+1+2); /* declare final duplex structure */
  char * struc2;
  struc2 = (char *) vrna_alloc(n2+1);
  /* char * struct_const; */
  for (k=MAX2(i,1); k<=i0; k++) if (!st1[k-1]) st1[k-1] = '.';
  /* for (k=j0; k<=j; k++) if (!st2[k-1]) st2[k-1] = struc_loop[k-1];*/ /* '.'; normal */
  /*  char * struct_const; */
  /*  struct_const = (char *) vrna_alloc(sizeof(char)*(n2+1));   */
  for (k=1; k<=n2; k++) {
    if (!st2[k-1]) st2[k-1] = struc_loop[k-1];/* '.'; */
    struc2[k-1] = st2[k-1];/* '.'; */
    /*      if (k>=j0 && k<=j){ */
    /*              struct_const[k-1]='x'; */
    /*      } */
    /*      else{ */
    /*              if(k<j0) {struct_const[k-1]='<';} */
    /*              if(k>j) {struct_const[k-1]='>';} */
    /*      } */
  }

  /*   char duplexseq_1[j0+1]; */
  /* char duplexseq_2[n2-j+3]; */
  if(j<n2){
    char **duplexseq_1, **duplexseq_2;
    duplexseq_1 = (char**) vrna_alloc((n_seq+1) * sizeof(char*));
    duplexseq_2 = (char**) vrna_alloc((n_seq+1) * sizeof(char*));
    for(s=0; s<n_seq; s++){
      duplexseq_1[s] = (char*) vrna_alloc((j0)*sizeof(char)); /* modfied j0+1 */
      duplexseq_2[s] = (char*) vrna_alloc((n2-j+2)*sizeof(char)); /* modified j+3 */
      strncpy(duplexseq_1[s], snoseq[s], j0-1); /* modified j0 */
      strcpy(duplexseq_2[s], snoseq[s] + j); /* modified j-1 */
      duplexseq_1[s][j0-1]='\0'; /* modified j0 */
      duplexseq_2[s][n2-j+1]='\0';/* modified j+2 */
    }
    duplexseq_1[n_seq]=NULL;
    duplexseq_2[n_seq]=NULL;
    duplexT temp;
    temp=aliduplexfold((const char**)duplexseq_1, (const char**)duplexseq_2);
    *Loop_D =  MIN2(0,-410 + (int) 100 * temp.energy*n_seq);
    if(*Loop_D){
      int l1,ibegin, iend, jbegin, jend;
      l1=strchr(temp.structure, '&')-temp.structure;
      ibegin=temp.i-l1;
      iend  =temp.i-1;
      jbegin=temp.j;
      jend  =temp.j+(int)strlen(temp.structure)-l1-2-1;
      for(k=ibegin+1; k<=iend+1; k++){
        struc2[k-1]=temp.structure[k-ibegin-1];
      }
      for(k=jbegin+j; k<=jend+j; k++){
        struc2[k-1]=temp.structure[l1+k-j-jbegin+1];
      }
    }
    for(s=0; s<n_seq; s++){
      free(duplexseq_1[s]);
      free(duplexseq_2[s]);
    }
    free(duplexseq_1);free(duplexseq_2);
    free(temp.structure);
  }
  strcpy(struc, st1+MAX2(i-1,0)); strcat(struc, "&"); 
  /* strcat(struc, st2); */
  strncat(struc, struc2+5, (int)strlen(struc2)-10);
  free(struc2);
  free(struc_loop);
  free(st1); free(st2);
  free(type);free(type2);free(type3);
  /* free_arrays(); */
  return struc;
}







void Lsnoop_subopt(const char *s1, const char *s2, int delta, int w, 
                   const int penalty, const int threshloop, 
                   const int threshLE, const int threshRE, const int threshDE, const int threshTE,const int threshSE,const int threshD,
                   const int distance,
                   const int half_stem, const int max_half_stem,
                   const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2, const int alignment_length, const char* name, const int fullStemEnergy)
{

 int min_colonne=INF;
 int max_pos;
 int max;max=INF;
 /* int temp; */
 /* int nsubopt=10; */
 n1 = (int) strlen(s1);
 n2 = (int) strlen(s2);
 int *position;
 position = (int*) vrna_alloc((n1+3)*sizeof(int));


  /* int Eminj, Emin_l; */
  int i, j; /* l1, Emin=INF, i_min=0, j_min=0; */
  /* char *struc; */
  /* snoopT mfe; */
  int *indx;
  int *mLoop;
  int *cLoop;
  folden **foldlist, **foldlist_XS;
  int Duplex_El, Duplex_Er;
  int Loop_D;
  /* int u; */
  int Loop_E;
  vrna_md_t md;

  Duplex_El=0;Duplex_Er=0;Loop_E=0, Loop_D=0;
  snoexport_fold_arrays(&indx, &mLoop, &cLoop, &foldlist, &foldlist_XS); 
  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    snoupdate_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  
  lc = (int **) vrna_alloc(sizeof(int *) * (5));
  lr = (int **) vrna_alloc(sizeof(int *) * (5));
  for (i=0; i<5; i++) {
          lc[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        lr[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        for(j=n2; j>-1; j--){
                lc[i][j]=INF;
                lr[i][j]=INF;
        }
  }
  encode_seqs(s1, s2);
  for (i=1; i<=n1; i++) {
      int idx=i%5;
      int idx_1=(i-1)%5;
      int idx_2=(i-2)%5;
      int idx_3=(i-3)%5;
      int idx_4=(i-4)%5;
    for (j=n2-min_d2; j>min_d1; j--) {
      int type, type2, k;
      type = pair[S1[i]][S2[j]];
      lc[idx][j] = (type) ? P->DuplexInit + 2*penalty : INF;
      lr[idx][j] = INF;
      if(!type) continue;
      if( /*pair[S1[i+1]][S2[j-1]]   &&  check that we have a solid base stack after the mLoop */
          j < max_s1 && j > min_s1 &&  
          j > n2 - max_s2 - max_half_stem && 
          j < n2 -min_s2 -half_stem && S1[i-2]==4) { /*constraint on s2 and i*/
        int min_k, max_k;
        max_k = MIN2(n2-min_s2,j+max_half_stem+1);
        min_k = MAX2(j+half_stem+1, n2-max_s2);
        for(k=min_k; k <= max_k ; k++){         
          if(mLoop[indx[k-1]+j+1] < 0){
                }
          if(pair[S1[i-3]][S2[k]] /*genau zwei ungepaarte nucleotiden --NU--*/
             && mLoop[indx[k-1]+j+1] < threshloop){  
            lr[idx][j]=MIN2(lr[idx][j], lc[idx_3][k]+mLoop[indx[k-1]+j+1]);
          }
          else if(pair[S1[i-4]][S2[k]] &&  mLoop[indx[k-1]+j+1] < threshloop){/*--NUN--*/
            lr[idx][j]=MIN2(lr[idx][j], lc[idx_4][k]+mLoop[indx[k-1]+j+1]);
          }
              }
      }
      /* dangle 5'SIDE relative to the mRNA  */
      lc[idx][j] += E_ExtLoop(type, (i>1) ? SS1[i-1] : -1, (j<n2) ? SS2[j+1] : -1, P);
      /**
      *** if (i>1)    lc[idx][j] += P->dangle5[type][SS1[i-1]];
      *** if (j<n2)   lc[idx][j] += P->dangle3[type][SS2[j+1]];
      *** if (type>2) lc[idx][j] += P->TerminalAU;
      **/
      
      if(j<n2 && i>1){
        type2=pair[S1[i-1]][S2[j+1]];
              if(type2>0){
          lc[idx][j]=MIN2(lc[idx_1][j+1]+E_IntLoop(0,0,type2, rtype[type],SS1[i], SS2[j], SS1[i-1], SS2[j+1], P)+2*penalty, lc[idx][j]);
          lr[idx][j]=MIN2(lr[idx_1][j+1]+E_IntLoop(0,0,type2, rtype[type],SS1[i], SS2[j], SS1[i-1], SS2[j+1], P)+2*penalty, lr[idx][j]);
              }
      }
      if(j<n2-1 && i>2){
        type2=pair[S1[i-2]][S2[j+2]];
        if(type2>0 ){
          lc[idx][j]=MIN2(lc[idx_2][j+2]+E_IntLoop(1,1,type2, rtype[type],SS1[i-1], SS2[j+1], SS1[i-1], SS2[j+1], P)+4*penalty, lc[idx][j]);
          lr[idx][j]=MIN2(lr[idx_2][j+2]+E_IntLoop(1,1,type2, rtype[type],SS1[i-1], SS2[j+1], SS1[i-1], SS2[j+1], P)+4*penalty, lr[idx][j]);
        }
      }
      if(j<n2-2 && i>3){
        type2 = pair[S1[i-3]][S2[j+3]];
        if(type2>0){
          lc[idx][j]=MIN2(lc[idx_3][j+3]+E_IntLoop(2,2,type2, rtype[type],SS1[i-2], SS2[j+2], SS1[i-1], SS2[j+1], P)+6*penalty,lc[idx][j]);
          lr[idx][j]=MIN2(lr[idx_3][j+3]+E_IntLoop(2,2,type2, rtype[type],SS1[i-2], SS2[j+2], SS1[i-1], SS2[j+1], P)+6*penalty,lr[idx][j]);
              }
      }
      /**
      *** (type>2?P->TerminalAU:0)+(i<(n1)?P->dangle3[rtype[type]][SS1[i+1]]+penalty:0)+(j>1?P->dangle5[rtype[type]][SS2[j-1]]+penalty:0)
      **/
      min_colonne=MIN2(lr[idx][j]+E_ExtLoop(rtype[type], (j > 1) ? SS2[j-1] : -1, (i<n1) ? SS1[i+1] : -1, P), min_colonne);
      }
      position[i]=min_colonne;
      if(max>=min_colonne){
        max=min_colonne;
        max_pos=i;
      }
      min_colonne=INF;
 }
 
  free(S1); free(S2); free(SS1); free(SS2);
  if(max<threshTE){
        find_max_snoop(s1, s2, max, alignment_length, position, delta, 
                       distance, penalty, threshloop, threshLE, threshRE, threshDE, 
                       threshTE, threshSE, threshD, half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2,name, fullStemEnergy);
   }
  for (i=1; i<5; i++) {free(lc[i]);free(lr[i]);}
  free(lc[0]);free(lr[0]);
  free(lc);free(lr);
  free(position);
  
}  



void Lsnoop_subopt_list(const char *s1, const char *s2, int delta, int w, 
                        const int penalty, const int threshloop, 
                        const int threshLE, const int threshRE, const int threshDE, const int threshTE,const int threshSE,const int threshD,
                        const int distance,
                        const int half_stem, const int max_half_stem,
                        const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2, const int alignment_length,const char *name,const int fullStemEnergy)
{
 
 int min_colonne=INF;
 int max_pos;
 int max;max=INF;
 /* int temp; */
 /* int nsubopt=10; */
 n1 = (int) strlen(s1);
 n2 = (int) strlen(s2);
 int *position;
 position = (int*) vrna_alloc((n1+3)*sizeof(int));


 /* int Eminj, Emin_l; */
  int i, j;/*  l1, Emin=INF, i_min=0, j_min=0; */
  /* char *struc; */
  /* snoopT mfe; */
  int *indx;
  int *mLoop;
  int *cLoop;
  folden **foldlist, **foldlist_XS;
  int Duplex_El, Duplex_Er;
  int Loop_D;
  /* int u; */
  int Loop_E;
  vrna_md_t md;

  Duplex_El=0;Duplex_Er=0;Loop_E=0, Loop_D=0;
  snoexport_fold_arrays(&indx, &mLoop, &cLoop, &foldlist,  &foldlist_XS); 

  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    snoupdate_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  
  lpair = (int **) vrna_alloc(sizeof(int *) * (6));
  lc    = (int **) vrna_alloc(sizeof(int *) * (6));
  lr    = (int **) vrna_alloc(sizeof(int *) * (6));
  for (i=0; i<6; i++) {
          lc[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        lr[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        lpair[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        for(j=n2; j>-1; j--){
                lc[i][j]=INF;
                lr[i][j]=INF;
                lpair[i][j]=0;
        }
  }
  encode_seqs(s1, s2);
  int lim_maxj=n2-min_d2 ;
  int lim_minj=min_d1;
  int lim_maxi=n1;
  for (i=5; i<=lim_maxi; i++) {
    int idx=i%5;
    int idx_1=(i-1)%5;
    int idx_2=(i-2)%5;
    int idx_3=(i-3)%5;
    int idx_4=(i-4)%5;

    for (j=lim_maxj; j>lim_minj; j--) {
      int type, type2;/*  E, k,l; */
      type = pair[S1[i]][S2[j]];
      lpair[idx][j] = type;
      lc[idx][j] = (type) ? P->DuplexInit + 2*penalty : INF;
      lr[idx][j] = INF;
      if(!type) continue;
      if( /*pair[S1[i+1]][S2[j-1]] && Be sure it binds*/
          j < max_s1 && j > min_s1 &&  
          j > n2 - max_s2 - max_half_stem && 
          j < n2 -min_s2 -half_stem && S1[i-2]==4 ) { /*constraint on s2 and i*/
        int min_k, max_k;
        max_k = MIN2(n2-min_s2,j+max_half_stem+1);
        min_k = MAX2(j+half_stem+1, n2-max_s2);
        folden * temp;
        temp=foldlist[j+1];
        while(temp->next){
          int k = temp->k;
          /* if(k >= min_k-1 && k < max_k){ comment to recover normal behaviour */
          if(lpair[idx_3][k+1] /*&& lpair[idx_4][k+2]*/){
            lr[idx][j]=MIN2(lr[idx][j], lc[idx_3][k+1]+temp->energy);/*--NU--*/
          }
          /*else*/ if(lpair[idx_4][k+1]){/*--NUN--*/
            lr[idx][j]=MIN2(lr[idx][j], lc[idx_4][k+1]+temp->energy);
          }
            /*  } */
          temp=temp->next;
        }
      }
      /* dangle 5'SIDE relative to the mRNA  */
      lc[idx][j] += E_ExtLoop(type,  SS1[i-1] , SS2[j+1] , P);
      /**
      ***      lc[idx][j] += P->dangle5[type][SS1[i-1]];
      ***      lc[idx][j] += P->dangle3[type][SS2[j+1]];
      ***      if (type>2) lc[idx][j] += P->TerminalAU;
      **/
      /*       if(j<n2 && i>1){ */
      /* type2=pair[S1[i-1]][S2[j+1]]; */
        type2=lpair[idx_1][j+1];
        if(type2>0 ){
          lc[idx][j]=MIN2(lc[idx_1][j+1]+E_IntLoop(0,0,type2, rtype[type],SS1[i], SS2[j], SS1[i-1], SS2[j+1], P)+2*penalty, lc[idx][j]);
          lr[idx][j]=MIN2(lr[idx_1][j+1]+E_IntLoop(0,0,type2, rtype[type],SS1[i], SS2[j], SS1[i-1], SS2[j+1], P)+2*penalty, lr[idx][j]);
        }
        /* } */
        /*       if(j<n2-1 && i>2){ */
        /* type2=pair[S1[i-2]][S2[j+2]]; */
        type2=lpair[idx_2][j+2];
        if(type2>0 ){
          lc[idx][j]=MIN2(lc[idx_2][j+2]+E_IntLoop(1,1,type2, rtype[type],SS1[i-1], SS2[j+1], SS1[i-1], SS2[j+1], P), lc[idx][j]);
          lr[idx][j]=MIN2(lr[idx_2][j+2]+E_IntLoop(1,1,type2, rtype[type],SS1[i-1], SS2[j+1], SS1[i-1], SS2[j+1], P), lr[idx][j]);
          /*         } */
      }
        /*       if(j<n2-2 && i>3){ */
        /* type2 = pair[S1[i-3]][S2[j+3]]; */
        type2 =lpair[idx_3][j+3];
        if(type2>0 ){
          lc[idx][j]=MIN2(lc[idx_3][j+3]+E_IntLoop(2,2,type2, rtype[type],SS1[i-2], SS2[j+2], SS1[i-1], SS2[j+1], P)+6*penalty,lc[idx][j]);
          lr[idx][j]=MIN2(lr[idx_3][j+3]+E_IntLoop(2,2,type2, rtype[type],SS1[i-2], SS2[j+2], SS1[i-1], SS2[j+1], P)+6*penalty,lr[idx][j]);
          /*         } */
      }
      /* min_colonne=MIN2(lr[idx][j]+(type>2?P->TerminalAU:0)+P->dangle3[rtype[type]][SS1[i+1]]+P->dangle5[rtype[type]][SS2[j-1]], min_colonne); */
      int bla;
      bla=lr[idx][j]+E_ExtLoop(rtype[type], SS2[j-1] , SS1[i+1], P)+2*penalty;
      min_colonne=MIN2(bla, min_colonne);
    }
    position[i]=min_colonne;
    if(max>=min_colonne){
      max=min_colonne;
      max_pos=i;
      }
    min_colonne=INF;
 }
 
  free(S1); free(S2); free(SS1); free(SS2);
  if(max<threshTE){
      find_max_snoop(s1, s2, max, alignment_length, position, 
                     delta, distance, penalty, threshloop, 
                     threshLE, threshRE, threshDE, threshTE, threshSE, threshD,
                     half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2,name, fullStemEnergy);
   }
  for (i=1; i<6; i++) {free(lc[i]);free(lr[i]);free(lpair[i]);}
  free(lc[0]);free(lr[0]);free(lpair[0]);
  free(lc);free(lr);free(lpair);
  free(position);
}  


PRIVATE void find_max_snoop(const char *s1, const char *s2,const int max,  const int alignment_length, const int* position, const int delta, 
                            const int distance, const int penalty, const int threshloop,  const int threshLE, const int threshRE, 
                            const int threshDE, const int threshTE, const int threshSE, const int threshD, 
                            const int half_stem, const int max_half_stem, const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2, const char* name, const int fullStemEnergy)
{
  int count=0;
  int pos=n1+1;
  int threshold = MIN2(threshTE , max + delta );
  /* printf("threshTE %d max %d\n", threshTE, max); */
  /* #pragma omp parallel for */
  /*   for(pos=n1+1;pos>distance;pos--){ */
  while(pos-- > 5){
    int temp_min=0;
    if(position[pos]<(threshold)){
      int search_range;
      search_range=distance+1;
      while(--search_range){
        if(position[pos-search_range]<=position[pos-temp_min]){
          temp_min=search_range;
        }
      }
      pos-=temp_min;
      int begin=MAX2(6, pos-alignment_length+1);
      char *s3 = (char*) vrna_alloc(sizeof(char)*(pos-begin+3+12));
      strcpy(s3, "NNNNN");
      strncat(s3, (s1+begin-1), pos-begin+2);
      strcat(s3,"NNNNN\0");
      /* printf("%s s3\n", s3);  */
      snoopT test;
      test = snoopfold(s3, s2, penalty, threshloop, threshLE, threshRE, threshDE, threshD, half_stem, max_half_stem, min_s2, max_s2, min_s1, 
                       max_s1, min_d1, min_d2,fullStemEnergy);
      if(test.energy==INF){
        free(s3);
        continue;
      }
      if(test.Duplex_El > threshLE * 0.01 || test.Duplex_Er > threshRE * 0.01  ||
         test.Loop_D > threshD * 0.01 || (test.Duplex_Er + test.Duplex_El) > threshDE * 0.01 ||
         (test.Duplex_Er + test.Duplex_El + test.Loop_E + test.Loop_D + 410) > threshSE*0.01) {
        free(test.structure);free(s3);
        continue;
      }
      int l1;
      l1 = strchr(test.structure, '&')-test.structure;
      
      int shift=0;
      if(test.i > (int)strlen(s3)-10){
        test.i--;
        l1--; 
      }
      if(test.i-l1<0){
        l1--;
        shift++;
      }
      char *target_struct =  (char*) vrna_alloc(sizeof(char) * (strlen(test.structure)+1));
      strncpy(target_struct, test.structure+shift, l1);
      strncat(target_struct, test.structure + (strchr(test.structure, '&')-
                                               test.structure), (int)strlen(test.structure) - (strchr(test.structure, '&')-
                                                                                          test.structure));
      strcat(target_struct,"\0");
      char *target; 
      target = (char *) vrna_alloc(l1+1);
      strncpy(target, (s3+test.i+5-l1), l1);
      target[l1]='\0';
      char *s4;
      s4 = (char*) vrna_alloc(sizeof(char) *(strlen(s2)-9));
      strncpy(s4, s2+5, (int)strlen(s2)-10);
      s4[(int)strlen(s2)-10]='\0';
      printf("%s %3d,%-3d;%3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f + %5.2f + 4.1 ) (%5.2f) \n%s&%s\n", 
             target_struct,begin + test.i-5-l1,begin + test.i -6 , begin + test.u -6, 
             test.j+1, test.j + (int)(strrchr(test.structure,'>') - strchr(test.structure,'>'))+1 ,
             test.Loop_D + test.Duplex_El + test.Duplex_Er + test.Loop_E + 4.10, test.Duplex_El,
             test.Duplex_Er, test.Loop_E, test.Loop_D,test.fullStemEnergy, target,s4);
      if(name){
        char *temp_seq;
        char *temp_struc;
        char *psoutput;
        temp_seq = (char*) vrna_alloc(sizeof(char)*(l1+n2-9));
        temp_struc = (char*) vrna_alloc(sizeof(char)*(l1+n2-9));
        strcpy(temp_seq, target);
        strcat(temp_seq, s4);
        strncpy(temp_struc, target_struct, l1);
        strcat(temp_struc, target_struct+l1+1);
        temp_seq[n2+l1-10]='\0';
        temp_struc[n2+l1-10]='\0';
        cut_point = l1+1;
        psoutput = vrna_strdup_printf("sno_%d_u_%d_%s.ps",
                                      count,
                                      begin + test.u - 6,
                                      name);

        PS_rna_plot_snoop_a(temp_seq, temp_struc, psoutput, NULL, NULL);
        cut_point = -1;
        free(temp_seq);
        free(temp_struc);
        free(psoutput);
        count++;
        /* free(psoutput); */
      }
      free(s4);
      free(test.structure);
      free(target_struct);
      free(target);
      free(s3);
    }
  }
  
}








snoopT snoopfold(const char *s1, const char *s2, 
                 const int penalty, const int threshloop, const int threshLE, const int threshRE, const int threshDE,
                 const int threshD,
                 const int half_stem, const int max_half_stem, 
                 const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2, const int fullStemEnergy) {
  /* int Eminj, Emin_l; */
  int i, j, l1, Emin=INF, i_min=0, j_min=0;
  char *struc;
  snoopT mfe;
  int *indx;
  int *mLoop;
  int *cLoop;
  folden** foldlist, **foldlist_XS;
  int Duplex_El, Duplex_Er;
  int Loop_D;
  int u;
  int Loop_E;
  vrna_md_t md;
  Duplex_El=0;Duplex_Er=0;Loop_E=0, Loop_D=0;
  snoexport_fold_arrays(&indx, &mLoop, &cLoop,&foldlist, &foldlist_XS ); 
  n1 = (int) strlen(s1);
  n2 = (int) strlen(s2);
  
  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    snoupdate_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  
  c = (int **) vrna_alloc(sizeof(int *) * (n1+1));
  r = (int **) vrna_alloc(sizeof(int *) * (n1+1));
  for (i=0; i<=n1; i++) {
          c[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        r[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        for(j=n2; j>-1; j--){
                c[i][j]=INF;
                r[i][j]=INF;
        }
  }
  encode_seqs(s1, s2);
  for (i=6; i<=n1-5; i++) {
    for (j=n2-min_d2; j>min_d1; j--) {
      int type, type2, E, k,l;
      type = pair[S1[i]][S2[j]];
      c[i][j] = (type ) ? P->DuplexInit : INF;
      if(!type) continue;
      if(/*  pair[S1[i+1]][S2[j-1]] &&  */
         j < max_s1 && j > min_s1 &&  
         j > n2 - max_s2 - max_half_stem && 
         j < n2 -min_s2 -half_stem && S1[i-2]==4  ) { /*constraint on s2 and i*/
        int min_k, max_k;
        max_k = MIN2(n2-min_s2,j+max_half_stem);
        min_k = MAX2(j+half_stem, n2-max_s2);
        folden * temp;
        temp=foldlist[j+1];
        while(temp->next){
          int k = temp->k;
          /* if(k >= min_k-1 && k < max_k){ uncomment to recovernormal behaviour */
          if(pair[S1[i-3]][S2[k+1]] /*&& pair[S1[i-4]][S2[k+2]]*/ ){
            r[i][j]=MIN2(r[i][j], c[i-3][k+1]+temp->energy);
          }
          /*else*/ if(pair[S1[i-4]][S2[k+1]] /*&& pair[S1[i-5]][S2[k+2]]*/ ){
            r[i][j]=MIN2(r[i][j], c[i-4][k+1]+temp->energy);
          }
          /* } */
          temp=temp->next;
        }
      }
      /* dangle 5'SIDE relative to the mRNA  */
      /**
      *** c[i][j] += P->dangle5[type][SS1[i-1]];
      *** c[i][j] += P->dangle3[type][SS2[j+1]];
      *** if (type>2) c[i][j] += P->TerminalAU;
      **/
      c[i][j]+=E_ExtLoop(type, SS1[i-1] , SS2[j+1], P);
      for (k=i-1; k>0 && (i-k)<MAXLOOP_L; k--) {
        for (l=j+1; l<=n2 ; l++) {
          if (i-k+l-j>2*MAXLOOP_L-2) break;
          if (abs(i-k-l+j) >= ASS ) continue;
          type2 = pair[S1[k]][S2[l]];
          if (!type2) continue;
          E = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                        SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
          c[i][j] = MIN2(c[i][j], c[k][l]+E+(i-k+l-j)*penalty);
          r[i][j] = MIN2(r[i][j], r[k][l]+E+(i-k+l-j)*penalty);
        }
      }
      E = r[i][j]; 
      /**
      *** if (i<n1) E += P->dangle3[rtype[type]][SS1[i+1]];              
      *** if (j>1)  E += P->dangle5[rtype[type]][SS2[j-1]]; 
      *** f (type>2) E += P->TerminalAU;
      **/
      E+=E_ExtLoop(rtype[type], (j > 1) ? SS2[j-1] : -1, (i<n1) ? SS1[i+1] : -1, P);
      if (E<Emin) {
        Emin=E; i_min=i; j_min=j;
      } 
    }
  }
  if(Emin > 0){
          printf("no target found under the constraints chosen\n");
        for (i=0; i<=n1; i++) {free(r[i]);free(c[i]);}
        free(c);
        free(r);
        free(S1); free(S2); free(SS1); free(SS2);
        mfe.energy=INF;
        return mfe;
  }
  struc = snoop_backtrack(i_min, j_min,s2, &Duplex_El, &Duplex_Er, &Loop_E, &Loop_D, 
                          &u, penalty, threshloop, threshLE, threshRE,threshDE, threshD,
                            half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2);
/*   if (i_min<n1-5) i_min++; */
/*   if (j_min>1 ) j_min--; */
  l1 = strchr(struc, '&')-struc;
  mfe.i = i_min-5;
  mfe.j = j_min-5;
  mfe.u = u -5;
  mfe.Duplex_Er = (float) Duplex_Er/100;
  mfe.Duplex_El = (float) Duplex_El/100;
  mfe.Loop_D = (float) Loop_D/100;
  mfe.Loop_E = (float) Loop_E/100;
  mfe.energy = (float) Emin/100 ;
  mfe.fullStemEnergy = (float) fullStemEnergy/100;
  mfe.structure = struc;
  if (!delay_free) {
    for (i=0; i<=n1; i++) {free(r[i]);free(c[i]);}
    free(c);
    free(r);
    free(S1); free(S2); free(SS1); free(SS2);
  }
  return mfe;
}

PRIVATE int snoopfold_XS_fill(const char *s1, const char *s2, const int **access_s1,
                 const int penalty, const int threshloop, const int threshLE, const int threshRE, const int threshDE,
                 const int threshD,
                 const int half_stem, const int max_half_stem, 
                 const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2) {
  /* int Eminj, Emin_l; */
  int i, j, Emin=INF, i_min=0, j_min=0;
  /* char *struc; */
  /* snoopT mfe; */
  int *indx;
  int *mLoop;
  int *cLoop;
  folden** foldlist, **foldlist_XS;
  int Duplex_El, Duplex_Er;
  int Loop_D;
  /* int u; */
  int Loop_E;
  vrna_md_t   md;
  Duplex_El=0;Duplex_Er=0;Loop_E=0, Loop_D=0;
  snoexport_fold_arrays(&indx, &mLoop, &cLoop,&foldlist, &foldlist_XS ); 
  n1 = (int) strlen(s1);
  n2 = (int) strlen(s2);
  
  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    snoupdate_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  
  c_fill = (int **) vrna_alloc(sizeof(int *) * (n1+1));
  r_fill = (int **) vrna_alloc(sizeof(int *) * (n1+1));
  for (i=0; i<=n1; i++) {
          c_fill[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        r_fill[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        for(j=n2; j>-1; j--){
                c_fill[i][j]=INF;
                r_fill[i][j]=INF;
        }
  }
  encode_seqs(s1, s2);

  int di[5];
  di[0]=0;  
  for (i=6; i<=n1-5; i++) {
    di[1]=access_s1[5][i]   - access_s1[4][i-1];           
    di[2]=access_s1[5][i-1] - access_s1[4][i-2] + di[1];
    di[3]=access_s1[5][i-2] - access_s1[4][i-3] + di[2];
    di[4]=access_s1[5][i-3] - access_s1[4][i-4] + di[3];
    di[1]=MIN2(di[1],165);
    di[2]=MIN2(di[2],330);
    di[3]=MIN2(di[3],495);
    di[4]=MIN2(di[4],660);
    for (j=n2-min_d2; j>min_d1; j--) {
      int type, type2, E, k,l;
      type = pair[S1[i]][S2[j]];
      c_fill[i][j] = (type ) ? P->DuplexInit : INF;
      if(!type) continue;
      if(/*  pair[S1[i+1]][S2[j-1]] &&  */
         j < max_s1 && j > min_s1 &&  
         j > n2 - max_s2 - max_half_stem && 
         j < n2 -min_s2 -half_stem && S1[i-2]==4  ) { /*constraint on s2 and i*/
        int min_k, max_k;
        max_k = MIN2(n2-min_s2,j+max_half_stem);
        min_k = MAX2(j+half_stem, n2-max_s2);
        folden * temp;
        temp=foldlist[j+1];
        while(temp->next){
          int k = temp->k;
          /* if(k >= min_k-1 && k < max_k){ uncomment to recovernormal behaviour */
          if(pair[S1[i-3]][S2[k+1]] /*&& pair[S1[i-4]][S2[k+2]]*/ ){
            r_fill[i][j]=MIN2(r_fill[i][j], c_fill[i-3][k+1]+temp->energy+ di[3]);
          }
          /*else*/ if(pair[S1[i-4]][S2[k+1]] /*&& pair[S1[i-5]][S2[k+2]]*/ ){
            r_fill[i][j]=MIN2(r_fill[i][j], c_fill[i-4][k+1]+temp->energy + di[4]);
          }
          /* } */
          temp=temp->next;
        }
      }
      /* dangle 5'SIDE relative to the mRNA  */
      /**
      *** c_fill[i][j] += P->dangle5[type][SS1[i-1]];
      *** c_fill[i][j] += P->dangle3[type][SS2[j+1]];
      *** if (type>2) c_fill[i][j] += P->TerminalAU;
      **/
      c_fill[i][j]+= E_ExtLoop(type, SS1[i-1], SS2[j+1],P);
      for (k=i-1; k>0 && (i-k)<MAXLOOP_L; k--) {
        for (l=j+1; l<=n2 ; l++) {
          if (i-k+l-j>2*MAXLOOP_L-2) break;
          if (abs(i-k-l+j) >= ASS ) continue;
          type2 = pair[S1[k]][S2[l]];
          if (!type2) continue;
          E = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                        SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
          c_fill[i][j] = MIN2(c_fill[i][j], c_fill[k][l]+E+di[i-k]);
          r_fill[i][j] = MIN2(r_fill[i][j], r_fill[k][l]+E+di[i-k]);
        }
      }
      E = r_fill[i][j]; 
      /**
      ***  if (i<n1) E += P->dangle3[rtype[type]][SS1[i+1]];              
      ***  if (j>1)  E += P->dangle5[rtype[type]][SS2[j-1]]; 
      *** if (type>2) E += P->TerminalAU;
      **/
      E+= E_ExtLoop(rtype[type], (j > 1) ? SS2[j-1] : -1, (i<n1) ? SS1[i+1] : -1, P);
      if (E<Emin) {
        Emin=E; i_min=i; j_min=j;
      } 
    }
  }
  return Emin;
}



PUBLIC snoopT *snoop_subopt(const char *s1, const char *s2, int delta, int w, 
                            const int penalty, const int threshloop, 
                            const int threshLE, const int threshRE, const int threshDE, const int threshTE, const int threshSE, const int threshD,
                            const int distance, const int half_stem, const int max_half_stem,
                            const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2, const int fullStemEnergy) {



  /* printf("%d %d\n", min_s2, max_s2); */
  int i,j, n1, n2, E, n_subopt=0, n_max;
  char *struc;
  snoopT mfe;
  snoopT *subopt;
  int thresh;

  int Duplex_El, Duplex_Er, Loop_E;
  int Loop_D;
  Duplex_El=0; Duplex_Er=0; Loop_E=0;Loop_D=0;
  int u;
  u=0;
  n_max=16;
  subopt = (snoopT *) vrna_alloc(n_max*sizeof(snoopT));
  delay_free=1;
  mfe = snoopfold(s1, s2, penalty, threshloop, threshLE, threshRE, threshDE,threshD,
                  half_stem, max_half_stem,
                  min_s2, max_s2, min_s1, max_s1, min_d1, min_d2, fullStemEnergy);



  if(mfe.energy > 0){
          free(subopt);
        delay_free=0;
        return NULL;
  }
  thresh = MIN2((int) ((mfe.Duplex_Er + mfe.Duplex_El + mfe.Loop_E)*100+0.1 + 410) + delta, threshTE );
 /* subopt[n_subopt++]=mfe; */
  free(mfe.structure);
  
  n1 = (int)strlen(s1); n2=(int)strlen(s2);
  for (i=n1; i>0; i--) {
    for (j=1; j<=n2; j++) {
      int type, Ed;
      type = pair[S2[j]][S1[i]];
      if (!type) continue;
      E = Ed = r[i][j];
      /**
      *** if (i<n1) Ed += P->dangle3[type][SS1[i+1]]; 
      *** if (j>1)  Ed += P->dangle5[type][SS2[j-1]]; 
      *** if (type>2) Ed += P->TerminalAU;
      **/
      Ed+= E_ExtLoop(type, (j > 1) ? SS2[j-1] : -1, (i<n1) ? SS1[i+1] : -1, P);
      if (Ed>thresh) continue;
      /* too keep output small, remove hits that are dominated by a
         better one close (w) by. For simplicity we do test without
         adding dangles, which is slightly inaccurate. 
      */ 
      /* w=1; */
/*       for (ii=MAX2(i-w,1); (ii<=MIN2(i+w,n1)) && type; ii++) {  */
/*         for (jj=MAX2(j-w,1); jj<=MIN2(j+w,n2); jj++) */
/*           if (r[ii][jj]<E) {type=0; break;} */
/*       } */
      if (!type) continue;

      struc = snoop_backtrack(i,j,s2, &Duplex_El, &Duplex_Er, &Loop_E, &Loop_D, &u, penalty, threshloop,threshLE,threshRE,threshDE,threshD, 
                        half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2);
      if (Duplex_Er > threshRE || Duplex_El > threshLE || Loop_D > threshD ||
         (Duplex_Er + Duplex_El) > threshDE || 
         (Duplex_Er + Duplex_El + Loop_E) > threshTE ||
         (Duplex_Er + Duplex_El + Loop_E + Loop_D + 410) > threshSE) {
                 /* printf(" Duplex_Er %d threshRE %d Duplex_El %d threshLE %d \n" */
                /*        " Duplex_Er + Duplex_El %d  threshDE %d \n" */
                /*        " Duplex_Er + Duplex_El + Loop_E %d  threshTE %d \n" */
                /*        " Duplex_Er + Duplex_El + Loop_E + Loop_D %d  threshSE %d \n",  */
                /*          Duplex_Er , threshRE , Duplex_El ,threshLE, */
                /*          Duplex_Er + Duplex_El, threshDE, */
                /*          Duplex_Er + Duplex_El+  Loop_E , threshTE, */
                /*          Duplex_Er + Duplex_El+  Loop_E + Loop_D, threshSE);  */
                 Duplex_Er=0; 
                Duplex_El=0;
                Loop_E = 0;
                Loop_D = 0;
                u=0,
                free(struc);
                continue;
        }

      if (n_subopt+1>=n_max) {
        n_max *= 2;
        subopt = (snoopT *) vrna_realloc(subopt, n_max*sizeof(snoopT));
      }
      subopt[n_subopt].i = i-5;
      subopt[n_subopt].j = j-5;
      subopt[n_subopt].u = u-5;
      subopt[n_subopt].Duplex_Er = Duplex_Er * 0.01;
      subopt[n_subopt].Duplex_El = Duplex_El * 0.01;
      subopt[n_subopt].Loop_E = Loop_E * 0.01;
      subopt[n_subopt].Loop_D = Loop_D * 0.01;
      subopt[n_subopt].energy = (Duplex_Er +Duplex_El + Loop_E + Loop_D + 410) * 0.01 ;
      subopt[n_subopt].fullStemEnergy = (float) fullStemEnergy * 0.01;
      subopt[n_subopt++].structure = struc;

      Duplex_Er=0; Duplex_El=0; Loop_E=0; Loop_D=0;u=0;
    }
  }
  
  for (i=0; i<=n1; i++) {free(c[i]);free(r[i]);}
  free(c);free(r);
  free(S1); free(S2); free(SS1); free(SS2);
  delay_free=0;

  if (snoop_subopt_sorted) qsort(subopt, n_subopt, sizeof(snoopT), compare);
  subopt[n_subopt].i =0;
  subopt[n_subopt].j =0;
  subopt[n_subopt].structure = NULL;
  return subopt;
}

PUBLIC void snoop_subopt_XS(const char *s1, const char *s2, const int **access_s1, int delta, int w, 
                            const int penalty, const int threshloop, 
                            const int threshLE, const int threshRE, const int threshDE, const int threshTE, const int threshSE, const int threshD,
                            const int distance, const int half_stem, const int max_half_stem,
                            const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2, const int alignment_length, const char *name, const int fullStemEnergy) {



  /* printf("%d %d\n", min_s2, max_s2); */
  int i,j, E,  n_max;
  /* char *struc; */
  /* snoopT mfe; */

  int thresh;

  int Duplex_El, Duplex_Er, Loop_E;
  int Loop_D;
  Duplex_El=0; Duplex_Er=0; Loop_E=0;Loop_D=0;
  int u;
  u=0;
  n_max=16;
  delay_free=1;
  int Emin = snoopfold_XS_fill(s1, s2, access_s1,penalty, threshloop, threshLE, threshRE, threshDE,threshD,
                               half_stem, max_half_stem,
                               min_s2, max_s2, min_s1, max_s1, min_d1, min_d2);
  if(Emin > 0){
    delay_free=0;
  }
  thresh = MIN2(-100, threshTE +alignment_length*30);  
  /*   n1=(int)strlen(s1);  */
  /*   n2=(int)strlen(s2); */
  
  int n3=(int)strlen(s1);
  int n4=(int)strlen(s2);
  S1_fill = (short*)vrna_alloc(sizeof(short)*(n3+2));
  S2_fill = (short*)vrna_alloc(sizeof(short)*(n4+2));
  SS1_fill = (short*)vrna_alloc(sizeof(short)*(n3+1));
  SS2_fill = (short*)vrna_alloc(sizeof(short)*(n4+1));
  memcpy(S1_fill, S1, sizeof(short)*n3+2);
  memcpy(S2_fill, S2, sizeof(short)*n4+2);
  memcpy(SS1_fill, SS1, sizeof(short)*n3+1);
  memcpy(SS2_fill, SS2, sizeof(short)*n4+1);
  free(S1);free(S2);free(SS1);free(SS2);
  int count=0;
  for (i=n3-5; i>0; i--) {
    for (j=1; j<=n4; j++) {
      int type, Ed;
      type = pair[S2_fill[j]][S1_fill[i]];
      if (!type) continue;
      E = Ed = r_fill[i][j];
      /**
      ***if (i<n3) Ed += P->dangle3[type][SS1_fill[i+1]]; 
      ***if (j>1)  Ed += P->dangle5[type][SS2_fill[j-1]]; 
      ***if (type>2) Ed += P->TerminalAU;
      **/
      Ed+=E_ExtLoop(type, (j > 1) ? SS2[j-1] : -1, (i<n3) ? SS1[i+1] : -1, P);
      if (Ed>thresh) continue;
      
      /* to keep output small, remove hits that are dominated by a
         better one close (w) by. For simplicity we do test without
         adding dangles, which is slightly inaccurate. 
      */ 
/*       w=10;  */
/*       for (ii=MAX2(i-w,1); (ii<=MIN2(i+w,n3-5)) && type; ii++) {   */
/*         for (jj=MAX2(j-w,1); jj<=MIN2(j+w,n4-5); jj++)  */
/*           if (r_fill[ii][jj]<E) {type=0; break;}  */
/*       }  */
/*       i=ii;j=jj; */
      if (!type) continue;
      int begin=MAX2(5, i-alignment_length);
      int end  =MIN2(n3-5, i-1); 
      char *s3 = (char*) vrna_alloc(sizeof(char)*(end-begin+2)+5);
      strncpy(s3, (s1+begin), end - begin +1);
      strcat(s3,"NNNNN\0");
      int n5 = (int)strlen(s3);
      snoopT test = snoopfold_XS(s3, s2, access_s1, i, j ,penalty, 
                          threshloop, threshLE, threshRE, 
                          threshDE, threshD, half_stem, 
                          max_half_stem, min_s2, max_s2, min_s1, 
                                 max_s1, min_d1, min_d2,fullStemEnergy);
      if(test.energy==INF){
        free(s3);
        continue;
      }
      if( test.Duplex_El > threshLE * 0.01 ||test.Duplex_Er > threshRE * 0.01  || 
          test.Loop_D > threshD * 0.01 || (test.Duplex_Er + test.Duplex_El) > threshDE * 0.01 || 
          (test.Duplex_Er + test.Duplex_El + test.Loop_E) > threshTE*0.01 || (test.Duplex_Er + test.Duplex_El + test.Loop_E + test.Loop_D + 410) > threshSE*0.01) 
        { 
          free(test.structure);free(s3); 
          continue; 
        }
      char *s4; 
      s4 = (char*) vrna_alloc(sizeof(char) *(n4-9)); 
      strncpy(s4, s2+5, n4-10); 
      s4[n4-10]='\0';
      
      char *s5 = vrna_alloc(sizeof(char) * n5-test.i+2-5);
      strncpy(s5,s3+test.i-1,n5-test.i+1-5);
      s5[n5-test.i+1-5]='\0';
      float dE = ((float) (access_s1[n5-test.i+1-5][i]))*0.01;
      printf("%s %3d,%-3d;%3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f + %5.2f + %5.2f + 4.10)  (%5.2f)\n%s&%s\n" ,  
             test.structure, i  -  (n5 - test.i) ,i - 5, i - (n5 - test.u ),
             j-5, j-5 + (int)(strrchr(test.structure,'>') - strchr(test.structure,'>')), 
             test.Loop_D + test.Duplex_El + test.Duplex_Er + test.Loop_E + 4.10+dE, test.Duplex_El, 
             test.Duplex_Er, test.Loop_E, test.Loop_D,dE , test.fullStemEnergy, s5,s4);
      if(name){
        int begin_t, end_t, begin_q, end_q, and, pipe,k; 
        char *psoutput;
        begin_q=0;
        end_q=n4-10;
        begin_t=0;
        end_t=n5-test.i+ 1-5;
        and=end_t+1;
        pipe=test.u -test.i + 1;
        cut_point = end_t +1 ;
        char *catseq, *catstruct;/*  *fname;  */
        catseq = (char*) vrna_alloc(n5 + end_q -begin_q +2);
        catstruct = (char*) vrna_alloc(n5 + end_q -begin_q +2);
        strcpy(catseq, s5);
        strncpy(catstruct, test.structure, end_t);
        strcat(catseq, s4);
        strncat(catstruct, test.structure+end_t+1, end_q-begin_q+1);
        catstruct[end_t - begin_t + end_q-begin_q+2]='\0';
        catseq[end_t - begin_t + end_q-begin_q+2]='\0';
        int *relative_access;
        relative_access = vrna_alloc(sizeof(int)*strlen(s5));
        relative_access[0] = access_s1[1][i - (n5  - test.i) + 5];
        for(k=1;k<(int)strlen(s5);k++){
          relative_access[k] =  access_s1[k+1][i - (n5  - test.i) + k + 5] -  access_s1[k][i - (n5  - test.i) + k + 4];
        }

        psoutput = vrna_strdup_printf("sno_XS_%d_u_%d_%s.ps",
                                      count,
                                      i - (n5 - test.u ),
                                      name);

        PS_rna_plot_snoop_a(catseq, catstruct, psoutput,relative_access,NULL);
        free(catseq);free(catstruct);free(relative_access);
        free(psoutput);
        count++;
      }
      free(s3);free(s4);free(s5);free(test.structure);
    }
  }  
  for (i=0; i<=n3; i++) {free(c_fill[i]);free(r_fill[i]);}
  free(c_fill);free(r_fill);
  free(S1_fill); free(S2_fill); free(SS1_fill); free(SS2_fill);
  delay_free=0;
}




PRIVATE char *snoop_backtrack(int i, int j, const char* snoseq, 
                              int *Duplex_El, int *Duplex_Er, 
                              int *Loop_E, int *Loop_D, int *u, 
                              const int penalty, const int threshloop, 
                              const int threshLE, const int threshRE, const int threshDE, const int threshD,
                              const int half_stem, const int max_half_stem, 
                              const int min_s2, const int max_s2, const int min_s1, 
                              const int max_s1, const int min_d1, const int min_d2) {
  /* backtrack structure going backwards from i, and forwards from j 
     return structure in bracket notation with & as separator */
  int k, l, type, type2, E, traced, i0, j0;
  int traced_r=0; /* flag for following backtrack in c or r */
  char *st1, *st2, *struc;
  char *struc_loop;

  st1 = (char *) vrna_alloc(sizeof(char)*(n1+1));
  st2 = (char *) vrna_alloc(sizeof(char)*(n2+1));
  int *indx;
  int *mLoop;
  int *cLoop;
  folden **foldlist, **foldlist_XS;
  type=pair[S1[i]][S2[j]];
  snoexport_fold_arrays(&indx, &mLoop, &cLoop,&foldlist, &foldlist_XS ); 
  i0=i; j0=j;
  /**
  *** if (i<n1)   *Duplex_Er += P->dangle3[rtype[type]][SS1[i+1]];
  *** if (j>1)    *Duplex_Er += P->dangle5[rtype[type]][SS2[j-1]];
  *** if (type>2) *Duplex_Er += P->TerminalAU;
  **/
  *Duplex_Er += E_ExtLoop(rtype[type], (j > 1) ? SS2[j-1] : -1, (i<n1) ? SS1[i+1] : -1, P);
  while (i>0 && j<=n2-min_d2 ) {
    if(!traced_r) {
      E = r[i][j]; traced=0;
      st1[i-1] = '<';
      st2[j-1] = '>'; 
      type = pair[S1[i]][S2[j]];
      if (!type) vrna_message_error("backtrack failed in fold duplex r");
      for (k=i-1; k>0 && (i-k)<MAXLOOP_L; k--) {
        for (l=j+1; l<=n2 ; l++) {
          int LE;
          if (i-k+l-j>2*MAXLOOP_L-2) break;
          if (abs(i-k-l+j) >= ASS) continue;
          
          type2 = pair[S1[k]][S2[l]];
          if (!type2) continue;
          LE = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                         SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
          if (E == r[k][l]+LE+(i-k+l-j)*penalty) {
            traced=1; 
            i=k; j=l;
            *Duplex_Er+=LE;
            break;
          }
        }
        if (traced) break;
      }
      if(!traced){
        if(/*  pair[S1[i+1]][S2[j-1]] && */  
            j < max_s1 && j > min_s1 && 
            j > n2 - max_s2 - max_half_stem && 
            j < n2 -min_s2 -half_stem && 
            S1[i-2]==4) {
          int min_k, max_k;
          max_k = MIN2(n2-min_s2,j+max_half_stem+1);
          min_k = MAX2(j+half_stem+1, n2-max_s2);
          folden * temp;
          temp=foldlist[j+1];
          while(temp->next) {
            int k = temp->k;
            if(pair[S1[i-3]][S2[k+1]] /*&& pair[S1[i-4]][S2[k+2]]*/    ){  /* introduce structure from RNAfold */
              if(E==c[i-3][k+1]+temp->energy){
                *Loop_E=temp->energy;
                st1[i-3]= '|';
                *u=i-2;
                int a,b;
                /* int fix_ij=indx[k-1+1]+j+1; */
                for(a=0; a< MISMATCH ;a++){
                  for(b=0; b< MISMATCH ; b++){
                    int ij=indx[k-1-a+1]+j+1+b;
                    if(cLoop[ij]==temp->energy) {
                      struc_loop=snobacktrack_fold_from_pair(snoseq, j+1+b, k-a-1+1);
                    a=INF; b=INF;        
                    }
                  }
                }
                traced=1;
                traced_r=1;
                i=i-3;j=k+1;
                break;
              }
            }
            /*else*/ if (pair[S1[i-4]][S2[k+1]] /*&& pair[S1[i-5]][S2[k+2]]*/){  /* introduce structure from RNAfold */
              if(E==c[i-4][k+1]+temp->energy){
                *Loop_E=temp->energy;
                st1[i-3]= '|';
                *u=i-2;
                int a,b;
                /* int fix_ij=indx[k-1+1]+j+1; */
                for(a=0; a< MISMATCH ;a++){
                  for(b=0; b< MISMATCH ; b++){
                    int ij=indx[k-1-a+1]+j+1+b;
                    if(cLoop[ij]==temp->energy) {
                      struc_loop=snobacktrack_fold_from_pair(snoseq, j+1+b, k-a-1+1);
                      a=INF; b=INF;        
                    }
                  }
                }
                traced=1;
                traced_r=1;
                i=i-4;j=k+1;
                break;
              }
            } /* else if */
            temp=temp->next;
          } /* while temp-> next */
        } /* test on j  */
      }/* traced? */
    }/* traced_r? */
    else{
      E = c[i][j]; traced=0;
      st1[i-1] = '<';
      st2[j-1] = '>'; 
      type = pair[S1[i]][S2[j]];
      if (!type) vrna_message_error("backtrack failed in fold duplex c");
      for (k=i-1; (i-k)<MAXLOOP_L; k--) {
        for (l=j+1; l<=n2; l++) {
          int LE;
          if (i-k+l-j>2*MAXLOOP_L-2) break;
          if (abs(i-k-l+j) >= ASS) continue;
          type2 = pair[S1[k]][S2[l]];
          if (!type2) continue;
          LE = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                         SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
          if (E == c[k][l]+LE+(i-k+l-j)*penalty) {
            traced=1; 
            i=k; j=l;
            *Duplex_El+=LE;
            break;
          }
        }
        if (traced) break;
      }
    }
    if (!traced) { 
      int correction;
      correction = E_ExtLoop(type, (i>1) ? SS1[i-1] : -1, (j<n2) ? SS2[j+1] : -1, P);
      E-=correction;
      *Duplex_El+=correction;
      /**
      *** if (i>1)    {E -= P->dangle5[type][SS1[i-1]]; *Duplex_El +=P->dangle5[type][SS1[i-1]];}
      *** if (j<n2)   {E -= P->dangle3[type][SS2[j+1]]; *Duplex_El +=P->dangle3[type][SS2[j+1]];}
      *** if (type>2) {E -= P->TerminalAU;                    *Duplex_El +=P->TerminalAU;}
      **/
      if (E != P->DuplexInit) {
        vrna_message_error("backtrack failed in fold duplex end");
      } else break;
    }
  }
/*   if (i>1)  i--; */
/*   if (j<n2) j++; */  
  /* struc = (char *) vrna_alloc(i0-i+1+j-j0+1+2); */ /* declare final duplex structure */
  struc = (char *) vrna_alloc(i0-i+1+n2-1+1+2); /* declare final duplex structure */
  char * struc2;
  struc2 = (char *) vrna_alloc(n2+1);
  /* char * struct_const; */
  for (k=MAX2(i,1); k<=i0; k++) if (!st1[k-1]) st1[k-1] = '.';
  /* for (k=j0; k<=j; k++) if (!st2[k-1]) st2[k-1] = struc_loop[k-1];*/ /* '.'; normal */
  /*  char * struct_const; */
  /*  struct_const = (char *) vrna_alloc(sizeof(char)*(n2+1));   */
  for (k=1; k<=n2; k++) {
    if (!st2[k-1]) st2[k-1] = struc_loop[k-1];/* '.'; */
    struc2[k-1] = st2[k-1];/* '.'; */
    /*      if (k>=j0 && k<=j){ */
    /*              struct_const[k-1]='x'; */
    /*      } */
    /*      else{ */
    /*              if(k<j0) {struct_const[k-1]='<';} */
    /*              if(k>j) {struct_const[k-1]='>';} */
    /*      } */
  }
  char duplexseq_1[j0];
  char duplexseq_2[n2-j+2];
  if(j<n2){
    strncpy(duplexseq_1, snoseq, j0-1);
    strcpy(duplexseq_2, snoseq + j);
    duplexseq_1[j0-1]='\0';
    duplexseq_2[n2-j+1]='\0';
    duplexT temp;
    temp=duplexfold(duplexseq_1, duplexseq_2);
    *Loop_D =  MIN2(0,-410 + (int) 100 * temp.energy);
    if(*Loop_D){
      int l1,ibegin, iend, jbegin, jend;
      l1=strchr(temp.structure, '&')-temp.structure;
      ibegin=temp.i-l1;
      iend  =temp.i-1;
      jbegin=temp.j;
      jend  =temp.j+(int)strlen(temp.structure)-l1-2-1;
      for(k=ibegin+1; k<=iend+1; k++){
        struc2[k-1]=temp.structure[k-ibegin-1];
      }
      for(k=jbegin+j; k<=jend+j; k++){
        struc2[k-1]=temp.structure[l1+k-j-jbegin+1];
      } 
    }
    free(temp.structure);
  } 
  strcpy(struc, st1+MAX2(i-1,0)); strcat(struc, "&"); 
  /* strcat(struc, st2); */
  strncat(struc, struc2+5, (int)strlen(struc2)-10);
  free(struc2);
  free(struc_loop);
  free(st1); free(st2);
  /* free_arrays(); */
  return struc;
}

void Lsnoop_subopt_list_XS(const char *s1, const char *s2,  const int **access_s1, int delta, int w, 
                        const int penalty, const int threshloop, 
                        const int threshLE, const int threshRE, const int threshDE, const int threshTE,const int threshSE,const int threshD,
                        const int distance,
                        const int half_stem, const int max_half_stem,
                           const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2, const int alignment_length, const char *name, const int fullStemEnergy)
{
 
 int min_colonne=INF;
 int max_pos;
 int max;max=INF;
 /* int temp; */
 /* int nsubopt=10; */
 n1 = (int) strlen(s1);
 n2 = (int) strlen(s2);
 int *position;
 int *position_j;
 int min_j_colonne;
 int max_pos_j=INF; 
 position = (int*) vrna_alloc((n1+3)*sizeof(int));
 position_j = (int*) vrna_alloc((n1+3)*sizeof(int));

 /* int Eminj, Emin_l; */
  int i, j;/*  l1, Emin=INF, i_min=0, j_min=0; */
  /* char *struc; */
  /* snoopT mfe; */
  int *indx;
  int *mLoop;
  int *cLoop;
  folden **foldlist, **foldlist_XS;
  int Duplex_El, Duplex_Er;
  int Loop_D;
  /* int u; */
  int Loop_E;
  vrna_md_t md;

  Duplex_El=0;Duplex_Er=0;Loop_E=0, Loop_D=0;
  snoexport_fold_arrays(&indx, &mLoop, &cLoop, &foldlist, &foldlist_XS); 

  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    snoupdate_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  
  lpair = (int **) vrna_alloc(sizeof(int *) * (6));
  lc    = (int **) vrna_alloc(sizeof(int *) * (6));
  lr    = (int **) vrna_alloc(sizeof(int *) * (6));
  for (i=0; i<6; i++) {
          lc[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        lr[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        lpair[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        for(j=n2; j>-1; j--){
                lc[i][j]=INF;
                lr[i][j]=INF;
                lpair[i][j]=0;
        }
  }
  encode_seqs(s1, s2);
  int lim_maxj=n2-min_d2 ;
  int lim_minj=min_d1;
  int lim_maxi=n1-5;
  for (i=5; i<=lim_maxi; i++) {
    int idx=i%5;
    int idx_1=(i-1)%5;
    int idx_2=(i-2)%5;
    int idx_3=(i-3)%5;
    int idx_4=(i-4)%5;
    int di1, di2, di3, di4;
    di1 = access_s1[5][i]   - access_s1[4][i-1];           
    di2 =access_s1[5][i-1] - access_s1[4][i-2] + di1;
    di3 = access_s1[5][i-2] - access_s1[4][i-3] + di2;
    di4 = access_s1[5][i-3] - access_s1[4][i-4] + di3;
    di1=MIN2(di1,165);
    di2=MIN2(di2,330);
    di3=MIN2(di3,495);
    di4=MIN2(di4,660);
    for (j=lim_maxj; j>lim_minj; j--) {
      int type, type2;/*  E, k,l; */
      type = pair[S1[i]][S2[j]];
      lpair[idx][j] = type;
      lc[idx][j] = (type) ? P->DuplexInit + access_s1[1][i] : INF;
      lr[idx][j] = INF;
      if(!type) continue;
      if( /*pair[S1[i+1]][S2[j-1]] && Be sure it binds*/
          j < max_s1 && j > min_s1 &&  
          j > n2 - max_s2 - max_half_stem && 
          j < n2 -min_s2 -half_stem && S1[i-2]==4 ) { /*constraint on s2 and i*/
        int min_k, max_k;
        max_k = MIN2(n2-min_s2,j+max_half_stem+1);
        min_k = MAX2(j+half_stem+1, n2-max_s2);
        folden * temp;
        temp=foldlist[j+1];
        while(temp->next){
          int k = temp->k;
          /* if(k >= min_k-1 && k < max_k){ comment to recover normal behaviour */
          if(lpair[idx_3][k+1] && lc[idx_3][k+1] /*+di3*/ < 411 /*&& lpair[idx_4][k+2]*/){ /*  remove second condition */
            lr[idx][j]=MIN2(lr[idx][j], di3 + lc[idx_3][k+1]+temp->energy);/*--NU--*/
          }
          /*else*/ if(lpair[idx_4][k+1] && /*di4 +*/ lc[idx_4][k+1] < 411  ){/*--NUN--*/ /*  remove second condition  */
            lr[idx][j]=MIN2(lr[idx][j], di4 + lc[idx_4][k+1]+temp->energy);
          }
            /*  } */
          temp=temp->next;
        }
      }
      /* dangle 5'SIDE relative to the mRNA  */
      /**
      *** lc[idx][j] += P->dangle5[type][SS1[i-1]];
      *** lc[idx][j] += P->dangle3[type][SS2[j+1]];
      *** if (type>2) lc[idx][j] += P->TerminalAU;
      **/
      lc[idx][j]+=E_ExtLoop(type, SS1[i-1] ,  SS2[j+1] , P);
      /*       if(j<n2 && i>1){ */
      /* type2=pair[S1[i-1]][S2[j+1]]; */
        type2=lpair[idx_1][j+1];
        if(type2>0 ){
          lc[idx][j]=MIN2(lc[idx_1][j+1]+E_IntLoop(0,0,type2, rtype[type],SS1[i], SS2[j], SS1[i-1], SS2[j+1], P)+di1, lc[idx][j]);
          lr[idx][j]=MIN2(lr[idx_1][j+1]+E_IntLoop(0,0,type2, rtype[type],SS1[i], SS2[j], SS1[i-1], SS2[j+1], P)+di1, lr[idx][j]);
        }
        type2=lpair[idx_2][j+2];
        if(type2>0 ){
          lc[idx][j]=MIN2(lc[idx_2][j+2]+E_IntLoop(1,1,type2, rtype[type],SS1[i-1], SS2[j+1], SS1[i-1], SS2[j+1], P)+di2, lc[idx][j]);
          lr[idx][j]=MIN2(lr[idx_2][j+2]+E_IntLoop(1,1,type2, rtype[type],SS1[i-1], SS2[j+1], SS1[i-1], SS2[j+1], P)+di2, lr[idx][j]);
         
      }
        type2 =lpair[idx_3][j+3];
        if(type2>0 ){
          lc[idx][j]=MIN2(lc[idx_3][j+3]+E_IntLoop(2,2,type2, rtype[type],SS1[i-2], SS2[j+2], SS1[i-1], SS2[j+1], P)+di3,lc[idx][j]);
          lr[idx][j]=MIN2(lr[idx_3][j+3]+E_IntLoop(2,2,type2, rtype[type],SS1[i-2], SS2[j+2], SS1[i-1], SS2[j+1], P)+di3,lr[idx][j]);

      }
      int bla;
      int temp2;
      temp2=min_colonne;
      bla=lr[idx][j] + E_ExtLoop(rtype[type], SS2[j-1], SS1[i+1] , P);
        /**
        *** (type>2?P->TerminalAU:0)+P->dangle3[rtype[type]][SS1[i+1]]+P->dangle5[rtype[type]][SS2[j-1]];
        **/
      min_colonne=MIN2(bla, min_colonne);
      if(temp2>min_colonne){
        min_j_colonne=j;
      }
    }
    position[i]=min_colonne;
    if(max>=min_colonne){
      max=min_colonne;
      max_pos=i;
      max_pos_j=min_j_colonne;
      }
    position_j[i]=min_j_colonne;
    min_colonne=INF;
 }
  free(S1); free(S2); free(SS1); free(SS2);

  if(max<threshTE + 30*alignment_length){
    find_max_snoop_XS(s1, s2, access_s1,max,alignment_length, position, position_j,
                     delta, distance, penalty, threshloop, 
                     threshLE, threshRE, threshDE, threshTE, threshSE, threshD,
                      half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2,name,fullStemEnergy);
   }
  for (i=1; i<6; i++) {free(lc[i]);free(lr[i]);free(lpair[i]);}
  free(lc[0]);free(lr[0]);free(lpair[0]);
  free(lc);free(lr);free(lpair);
  free(position);free(position_j);
}  


PRIVATE void find_max_snoop_XS(const char *s1, const char *s2, const int **access_s1, 
                               const int max,  const int alignment_length, 
                               const int* position, const int* position_j, const int delta, 
                               const int distance, const int penalty, const int threshloop,  const int threshLE, const int threshRE, 
                               const int threshDE, const int threshTE, const int threshSE, const int threshD, 
                               const int half_stem, const int max_half_stem, const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2, const char *name, const int fullStemEnergy){
  int count=0;
  int n3=(int)strlen(s1);
  int n4=(int)strlen(s2);
  int pos=n1-4;
  int max_pos_j;
  int threshold = MIN2(threshTE + alignment_length*30, -100);
  /* printf("threshTE %d max %d\n", threshTE, max); */
  /* #pragma omp parallel for */
  /*   for(pos=n1+1;pos>distance;pos--){ */
  while(pos-->5){
    int temp_min=0;
    if(position[pos]<(threshold)){
      int search_range;
      search_range=distance+1;
      while(--search_range){
        if(position[pos-search_range]<=position[pos-temp_min]){
          temp_min=search_range;
        }
      }
      pos-=temp_min;
      max_pos_j=position_j[pos];
      int begin=MAX2(5, pos-alignment_length);
      int end  =MIN2(n3-5, pos-1); 
      char *s3 = (char*) vrna_alloc(sizeof(char)*(end-begin+2)+5);
      strncpy(s3, (s1+begin), end - begin +1);
      strcat(s3,"NNNNN\0");

      int n5 = (int)strlen(s3);
      snoopT test;
      test = snoopfold_XS(s3, s2, access_s1, pos, max_pos_j ,penalty, 
                          threshloop, threshLE, threshRE, 
                          threshDE, threshD, half_stem, 
                          max_half_stem, min_s2, max_s2, min_s1, 
                          max_s1, min_d1, min_d2, fullStemEnergy);
      if(test.energy==INF){
        free(s3);
        continue;
      }
      if( test.Duplex_El > threshLE * 0.01 ||test.Duplex_Er > threshRE * 0.01  || 
         test.Loop_D > threshD * 0.01 || (test.Duplex_Er + test.Duplex_El) > threshDE * 0.01 || 
         (test.Duplex_Er + test.Duplex_El + test.Loop_E) > threshTE*0.01 || (test.Duplex_Er + test.Duplex_El + test.Loop_E + test.Loop_D + 410) > threshSE*0.01) { 
        free(test.structure);free(s3); 
        continue; 
      }
      
      char *s4; 
      s4 = (char*) vrna_alloc(sizeof(char) *(n4-9)); 
      strncpy(s4, s2+5, n4-10); 
      s4[n4-10]='\0';

      char *s5 = vrna_alloc(sizeof(char) * n5-test.i+2-5);
      strncpy(s5,s3+test.i-1,n5-test.i+1-5);
      s5[n5-test.i+1-5]='\0';
      float dE = ((float) (access_s1[n5-test.i+1-5][pos]))*0.01;
      printf("%s %3d,%-3d;%3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f + %5.2f + %5.2f + 4.10) (%5.2f)\n%s&%s\n" ,  
             test.structure, pos  -  (n5 - test.i) ,pos - 5, pos - (n5 - test.u ),
             max_pos_j-5, max_pos_j-5 + (int)(strrchr(test.structure,'>') - strchr(test.structure,'>')), 
             test.Loop_D + test.Duplex_El + test.Duplex_Er + test.Loop_E + 4.10+dE, test.Duplex_El, 
             test.Duplex_Er, test.Loop_E, test.Loop_D,dE ,test.fullStemEnergy, s5,s4);
      if(name){
        int begin_t, end_t, begin_q, end_q, and, pipe, i; 
        char *psoutput;
        begin_q=0;
        end_q=n4-10;
        begin_t=0;
        end_t=n5-test.i+ 1-5;
        and=end_t+1;
        pipe=test.u -test.i + 1;
        cut_point = end_t +1 ;
        char *catseq, *catstruct;/*  *fname;  */
        catseq = (char*) vrna_alloc(n5 + end_q -begin_q +2);
        catstruct = (char*) vrna_alloc(n5 + end_q -begin_q +2);
        strcpy(catseq, s5);
        strncpy(catstruct, test.structure, end_t);
        strcat(catseq, s4);
        strncat(catstruct, test.structure+end_t+1, end_q-begin_q+1);
        catstruct[end_t - begin_t + end_q-begin_q+2]='\0';
        catseq[end_t - begin_t + end_q-begin_q+2]='\0';
        int *relative_access;
        relative_access = vrna_alloc(sizeof(int)*strlen(s5));

        relative_access[0] = access_s1[1][pos - (n5  - test.i) + 5];
        for(i=1;i<(int)strlen(s5);i++){
          relative_access[i] =  access_s1[i+1][pos - (n5  - test.i) + i + 5] -  access_s1[i][pos - (n5  - test.i) + i + 4];
        }

        psoutput = vrna_strdup_printf("sno_XS_%d_u_%d_%s.ps",
                                      count,
                                      pos - (n5 - test.u ),
                                      name);

        PS_rna_plot_snoop_a(catseq, catstruct, psoutput,relative_access,NULL);
        free(catseq);free(catstruct);free(relative_access);
        free(psoutput);
        count++;
      }
      free(s3);free(s4);free(s5);free(test.structure);
    }
  }
}

snoopT snoopfold_XS(const char *s1, const char *s2, const int **access_s1, const int pos_i, const int pos_j,
                 const int penalty, const int threshloop, const int threshLE, const int threshRE, const int threshDE,
                 const int threshD,
                 const int half_stem, const int max_half_stem, 
                    const int min_s2, const int max_s2, const int min_s1, const int max_s1, const int min_d1, const int min_d2, const int fullStemEnergy) {
  /*   int Eminj, Emin_l; */
  int a,b,i, j, Emin=INF, a_min=0, b_min=0;
  char *struc;
  snoopT mfe;
  int *indx;
  int *mLoop;
  int *cLoop;
  folden** foldlist, **foldlist_XS;
  int Duplex_El, Duplex_Er;
  int Loop_D;
  int u;
  int Loop_E;
  vrna_md_t md;

  Duplex_El=0;Duplex_Er=0;Loop_E=0, Loop_D=0;
  snoexport_fold_arrays(&indx, &mLoop, &cLoop,&foldlist, &foldlist_XS ); 
  n1 = (int) strlen(s1);
  n2 = (int) strlen(s2);
  
  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    snoupdate_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  
  c = (int **) vrna_alloc(sizeof(int *) * (n1+1));
  r = (int **) vrna_alloc(sizeof(int *) * (n1+1));
  for (i=0; i<=n1; i++) {
          c[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        r[i] = (int *) vrna_alloc(sizeof(int) * (n2+1));
        for(j=n2; j>-1; j--){
                c[i][j]=INF;
                r[i][j]=INF;
        }
  }
  encode_seqs(s1, s2);
  i=n1-5;
  j=pos_j;
  /* printf("tar: %s\nsno: %s\n ", s1, s2); */
  /* printf("pos_i %d pos_j %d\n", pos_i, pos_j); */
  /* printf("type %d n1 %d n2 %d S1[n1] %d S2[n2] %d", pair[S1[i]][S2[j]], n1, n2, S1[n1], S2[n2]); */
  int type, type2, E, p,q;    
  r[i][j] = P->DuplexInit; 
  /* r[i][j] += P->dangle3[rtype[type]][SS1[i+1]] + P->dangle5[rtype[type]][SS2[j-1]];  */
  
  if(pair[S1[i]][S2[j]]>2) r[i][j]+=P->TerminalAU;
  for(a=i-1; a>0;a--){ /* i-1 */
    r[a+1][0]=INF;
    for(b=j+1; b<=n2-min_d2;b++){ /* j+1 */
      r[a][b]=INF;
      type = pair[S1[a]][S2[b]]; 
       if(!type) continue; 
       if(S1[a+1]==4){ 
         folden * temp; 
         temp=foldlist_XS[b-1]; 
         while(temp->next) {     
           int k = temp->k; 
           if(pair[S1[a+3]][S2[k-1]] && k< max_s1 && k > min_s1 && k > n2 - max_s2 - max_half_stem &&  k < n2 -min_s2 -half_stem /*&& r[a+3][k-1] + access_s1[i-(a+3)+1][pos_i] < 411*/) { /* remove last condition last condition is to check if the interaction is stable enough */
             c[a][b]=MIN2(c[a][b], r[a+3][k-1]+temp->energy); 
           }
           temp=temp->next;
         }
       }
       if(S1[a+2]==4){
         folden * temp; 
         temp=foldlist_XS[b-1]; 
         while(temp->next){ 
           int k = temp->k; 
           if(pair[S1[a+4]][S2[k-1]] &&  k< max_s1 && k > min_s1 && k > n2 - max_s2 - max_half_stem &&  k < n2 -min_s2 -half_stem /*&& r[a+4][k-1] + access_s1[i-(a+4)+1][pos_i] < 411 */ ) { /* remove last condition  */
             c[a][b]=MIN2(c[a][b], r[a+4][k-1]+temp->energy); 
           }
           temp=temp->next;
         }
       }
       for(p=a+1; p<n1 && (p-a) < MAXLOOP_L; p++) { /* p < n1 */
         for (q=b-1; q > 1  ; q--) {  /* q > 1 */
           if (p-a+b-q>2*MAXLOOP_L-2) break; 
           if (abs((p-a)-(b-q)) >= ASS ) continue; 
           type2 = pair[S1[p]][S2[q]]; 
           if (!type2) continue; 
           E = E_IntLoop(p-a-1, b-q-1, type2, rtype[type],SS1[a+1], SS2[b-1], SS1[p-1], SS2[q+1],P); 
           c[a][b] = MIN2(c[a][b], c[p][q]+E); 
           r[a][b] = MIN2(r[a][b], r[p][q]+E); 
         }
       }
       E = c[a][b];  
       if (type>2) E += P->TerminalAU;  
       /*        E +=P->dangle5[rtype[type]][SS1[i+1]]; */
       /* E +=P->dangle5[rtype[type]][SS2[j-1]];  */
       E+=access_s1[i-a+1][pos_i]; 
       if (E<Emin) { 
         Emin=E; a_min=a; b_min=b; 
       }  
    }
  }
  if(Emin > 0){ 
    printf("no target found under the constraints chosen\n");
    for (i=0; i<=n1; i++) {free(r[i]);free(c[i]);}
    free(c);
    free(r); 
    free(S1); free(S2); free(SS1); free(SS2);
    mfe.energy=INF;
    return mfe;
  }  
  type2=pair[S1[a_min]][S2[b_min]];
  if(type2>2) Emin +=P->TerminalAU;
  mfe.energy = ((float) (Emin))/100;
  struc = snoop_backtrack_XS(a_min, b_min,s2, &Duplex_El, &Duplex_Er, &Loop_E, &Loop_D, 
                             &u, penalty, threshloop, threshLE, threshRE,threshDE, threshD, 
                             half_stem, max_half_stem, min_s2, max_s2, min_s1, max_s1, min_d1, min_d2); 

  mfe.i = a_min;
  mfe.j = b_min;
  mfe.u = u;
  mfe.Duplex_Er = (float) Duplex_Er/100;
  mfe.Duplex_El = (float) Duplex_El/100;
  mfe.Loop_D = (float) Loop_D/100;
  mfe.Loop_E = (float) Loop_E/100;
  mfe.energy = (float) Emin/100 ;
  mfe.fullStemEnergy = (float) fullStemEnergy/100;
  mfe.structure = struc;
  return mfe;
}

PRIVATE char *snoop_backtrack_XS(int i, int j, const char* snoseq, 
                              int *Duplex_El, int *Duplex_Er, 
                              int *Loop_E, int *Loop_D, int *u, 
                              const int penalty, const int threshloop, 
                              const int threshLE, const int threshRE, const int threshDE, const int threshD,
                              const int half_stem, const int max_half_stem, 
                              const int min_s2, const int max_s2, const int min_s1, 
                              const int max_s1, const int min_d1, const int min_d2) {
  /* backtrack structure going backwards from i, and forwards from j 
     return structure in bracket notation with & as separator */
  int k, l, type, type2, E, traced, i0, j0;
  int traced_c=0; /* flag for following backtrack in c or r */
  char *st1, *st2, *struc;
  char *struc_loop;

  st1 = (char *) vrna_alloc(sizeof(char)*(n1+1));
  st2 = (char *) vrna_alloc(sizeof(char)*(n2+1));
  int *indx;
  int *mLoop;
  int *cLoop;
  folden **foldlist, **foldlist_XS;
  type=pair[S1[i]][S2[j]];
  snoexport_fold_arrays(&indx, &mLoop, &cLoop,&foldlist, &foldlist_XS ); 
  i0=i;j0=j;
  /*   i0=MAX2(i,1); j0=MIN2(j+1,n2); */
  while (i<=n1 && j>=1 ) {
    if(!traced_c) {
      E = c[i][j]; traced=0;
      st1[i] = '<';
      st2[j-1] = '>'; 
      type = pair[S1[i]][S2[j]];
      if (!type) vrna_message_error("backtrack failed in fold duplex c");
      for (k=i+1; k>0 && (k-i)<MAXLOOP_L; k++) {
        for (l=j-1; l>=1 ; l--) {
          int LE;
          if (k-i+j-l>2*MAXLOOP_L-2) break;
          if (abs(k-i-j+l) >= ASS) continue;
          type2 = pair[S1[k]][S2[l]];
          if (!type2) continue;
          LE = E_IntLoop(k-i-1, j-l-1, type2, rtype[type],
                         SS1[i+1], SS2[j-1], SS1[k-1], SS2[l+1],P);
          if (E == c[k][l]+LE) {
            traced=1; 
            i=k; j=l;
            *Duplex_El+=LE;
            break;
          }
        }
        if (traced) break;
      }
      if(!traced){
        if(S1[i+1]==4) {
          folden * temp;
          temp=foldlist_XS[j-1];
          while(temp->next) {
            int k = temp->k;
            if(pair[S1[i+3]][S2[k-1]] && k< max_s1 && k > min_s1 && k > n2 - max_s2 - max_half_stem &&  k < n2 -min_s2 -half_stem ) {
              if(E==r[i+3][k-1]+temp->energy){
                *Loop_E=temp->energy;
                st1[i+1]= '|';
                st1[i+2]='.';
                *u=i+1;
                int a,b;
                for(a=0; a< MISMATCH ;a++){
                  for(b=0; b< MISMATCH ; b++){
                    int ij=indx[j-1-a]+k+b;
                    if(cLoop[ij]==temp->energy) {
                      struc_loop=snobacktrack_fold_from_pair(snoseq, k+b, j-1-a); 
                      a=INF; b=INF;        
                    }
                  }
                }
                traced=1;
                traced_c=1;
                i=i+3;j=k-1;
                break;
              }
            }
            temp=temp->next;
          }
        }
        if (S1[i+2]==4){  /* introduce structure from RNAfold */
          folden * temp;
          temp=foldlist_XS[j-1];
          while(temp->next) {
            int k = temp->k;
            if(pair[S1[i+4]][S2[k-1]] && k< max_s1 && k > min_s1 && k > n2 - max_s2 - max_half_stem &&  k < n2 -min_s2 -half_stem ) {
              if(E==r[i+4][k-1]+temp->energy){
                *Loop_E=temp->energy;
                st1[i+2]= '|';
                st1[i+1]=st1[i+3]='.';
                *u=i+2;
                int a,b;
                for(a=0; a< MISMATCH ;a++){
                  for(b=0; b< MISMATCH ; b++){
                    int ij=indx[j-1-a]+k+b;
                    if(cLoop[ij]==temp->energy) {
                      struc_loop=snobacktrack_fold_from_pair(snoseq, k+b, j-a-1);
                      a=INF; b=INF;        
                    }
                  }
                }
                traced=1;
                traced_c=1;
                i=i+4;j=k-1;
                break;
              }
            }
            temp=temp->next;
          }
        }
      }/* traced? */
    }/* traced_r? */
    else{
      E = r[i][j]; traced=0;
      st1[i] = '<';
      st2[j-1] = '>'; 
      type = pair[S1[i]][S2[j]];
      if (!type) vrna_message_error("backtrack failed in fold duplex r");
      for (k=i+1; k>0 && (k-i)<MAXLOOP_L; k++) {
        for (l=j-1; l>=1 ; l--) {
          int LE;
          if (k-i+j-l>2*MAXLOOP_L-2) break;
          if (abs(k-i-j+l) >= ASS) continue;
          type2 = pair[S1[k]][S2[l]];
          if (!type2) continue;
          LE = E_IntLoop(k-i-1, j-l-1, type2, rtype[type],
                         SS1[i+1], SS2[j-1], SS1[k-1], SS2[l+1],P);
          if (E == r[k][l]+LE) {
            traced=1; 
            i=k; j=l;
            *Duplex_Er+=LE;
            break;
          }
        }
        if (traced) break;
      }
    }
    if (!traced) { 
/*       if (i>1)    {E -= P->dangle5[type][SS1[i-1]]; *Duplex_El +=P->dangle5[type][SS1[i-1]];} */
/*       if (j<n2)   {E -= P->dangle3[type][SS2[j+1]]; *Duplex_El +=P->dangle3[type][SS2[j+1]];} */
      if (type>2) {E -= P->TerminalAU;        *Duplex_Er +=P->TerminalAU;}
      if (E != P->DuplexInit) {
        vrna_message_error("backtrack failed in fold duplex end");
      } else break;
    }
  }

  
  /* struc = (char *) vrna_alloc(i0-i+1+j-j0+1+2); */ /* declare final duplex structure */
  struc = (char *) vrna_alloc(i-i0+1+n2); /* declare final duplex structure */
  char * struc2;
  struc2 = (char *) vrna_alloc(n2+1);
  /* char * struct_const; */

  for (k=MIN2(i0,1); k<=i; k++) if (!st1[k-1]) st1[k-1] = '.';
  /* for (k=j0; k<=j; k++) if (!st2[k-1]) st2[k-1] = struc_loop[k-1];*/ /* '.'; normal */
  /*  char * struct_const; */
  /*  struct_const = (char *) vrna_alloc(sizeof(char)*(n2+1));   */
  for (k=1; k<=n2; k++) {
    if (!st2[k-1]) st2[k-1] = struc_loop[k-1];/* '.'; */
    struc2[k-1] = st2[k-1];/* '.'; */
    /*      if (k>=j0 && k<=j){ */
    /*              struct_const[k-1]='x'; */
    /*      } */
    /*      else{ */
    /*              if(k<j0) {struct_const[k-1]='<';} */
    /*              if(k>j) {struct_const[k-1]='>';} */
    /*      } */
  }
  char duplexseq_1[j];
  char duplexseq_2[n2-j0+2];
  if(j0<n2){
    strncpy(duplexseq_1, snoseq, j-1);
    strcpy(duplexseq_2, snoseq + j0);
    duplexseq_1[j-1]='\0';
    duplexseq_2[n2-j0+1]='\0';
    duplexT temp;
    temp=duplexfold(duplexseq_1, duplexseq_2);
    *Loop_D =  MIN2(0,-410 + (int) 100 * temp.energy);
    if(*Loop_D){
      int l1,ibegin, iend, jbegin, jend;
      l1=strchr(temp.structure, '&')-temp.structure;
      ibegin=temp.i-l1;
      iend  =temp.i-1;
      jbegin=temp.j;
      jend  =temp.j+(int)strlen(temp.structure)-l1-2-1;
      for(k=ibegin+1; k<=iend+1; k++){
        struc2[k-1]=temp.structure[k-ibegin-1];
      }
      for(k=jbegin+j0; k<=jend+j0; k++){
        struc2[k-1]=temp.structure[l1+k-j0-jbegin+1];
      } 
    }
    free(temp.structure);
  } 
  strcpy(struc, st1+MAX2(i0,1)); strcat(struc, "&"); 
  /* strcat(struc, st2); */
  strncat(struc, struc2+5, (int)strlen(struc2)-10);
  free(struc2);
  free(struc_loop);
  free(st1); free(st2);
  
    for (i=0; i<=n1; i++) {free(r[i]);free(c[i]);}
    free(c);
    free(r);
    free(S1);free(S2);free(SS1);free(SS2);
    /* free_arrays(); */
  return struc;
}

PRIVATE int covscore(const int *types, int n_seq) {
  /* calculate co-variance bonus for a pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */
#define NONE -10000 /* score for forbidden pairs */
  int k,l,s,score, pscore;
  int dm[7][7]={{0,0,0,0,0,0,0}, /* hamming distance between pairs */
                {0,0,2,2,1,2,2} /* CG */,
                {0,2,0,1,2,2,2} /* GC */,
                {0,2,1,0,2,1,2} /* GU */,
                {0,1,2,2,0,2,1} /* UG */,
                {0,2,2,1,2,0,2} /* AU */,
                {0,2,2,2,1,2,0} /* UA */};
  
  int pfreq[8]={0,0,0,0,0,0,0,0};
  for (s=0; s<n_seq; s++)
    pfreq[types[s]]++;

  if (pfreq[0]*2>n_seq) return NONE;
  for (k=1,score=0; k<=6; k++) /* ignore pairtype 7 (gap-gap) */
    for (l=k+1; l<=6; l++)
      /* scores for replacements between pairtypes    */
      /* consistent or compensatory mutations score 1 or 2  */
      score += pfreq[k]*pfreq[l]*dm[k][l];
  
  /* counter examples score -1, gap-gap scores -0.25   */
  pscore = cv_fact * 
    ((UNIT*score)/n_seq - nc_fact*UNIT*(pfreq[0] + pfreq[7]*0.25));
  return pscore;
}

/*---------------------------------------------------------------------------*/

PRIVATE short * aliencode_seq(const char *sequence) {
  unsigned int i,l;
  short *Stemp;
  l = strlen(sequence);
  Stemp = (short *) vrna_alloc(sizeof(short)*(l+2));
  Stemp[0] = (short) l;

  /* make numerical encoding of sequence */
  for (i=1; i<=l; i++)
    Stemp[i]= (short) encode_char(toupper(sequence[i-1]));

  /* for circular folding add first base at position n+1 */
  /* Stemp[l+1] = Stemp[1]; */

  return Stemp;
}

PRIVATE short * encode_seq(const char *sequence) {
  unsigned int i,l;
  short *S;
  l = strlen(sequence);
extern double nc_fact;
  S = (short *) vrna_alloc(sizeof(short)*(l+2));
  S[0] = (short) l;

  /* make numerical encoding of sequence */
  for (i=1; i<=l; i++)
    S[i]= (short) encode_char(toupper(sequence[i-1]));
  /* for circular folding add first base at position n+1 */
  S[l+1] = S[1];

  return S;
}

PRIVATE void encode_seqs(const char *s1, const char *s2) {
  unsigned int i,l;

  l = strlen(s1);
  S1 = encode_seq(s1);
  SS1= (short *) vrna_alloc(sizeof(short)*(l+1));
  /* SS1 exists only for the special X K and I bases and energy_set!=0 */
  
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    SS1[i] = alias[S1[i]];   /* for mismatches of nostandard bases */
  }

  l = strlen(s2);
  S2 = encode_seq(s2);
  SS2= (short *) vrna_alloc(sizeof(short)*(l+1));
  /* SS2 exists only for the special X K and I bases and energy_set!=0 */
  
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    SS2[i] = alias[S2[i]];   /* for mismatches of nostandard bases */
  }
}

PRIVATE int compare(const void *sub1, const void *sub2) {
  int d;
  if (((snoopT *) sub1)->energy > ((snoopT *) sub2)->energy)
    return 1;
  if (((snoopT *) sub1)->energy < ((snoopT *) sub2)->energy)
    return -1;
  d = ((snoopT *) sub1)->i - ((snoopT *) sub2)->i;
  if (d!=0) return d;
  return  ((snoopT *) sub1)->j - ((snoopT *) sub2)->j;
}



