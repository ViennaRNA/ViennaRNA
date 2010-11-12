/* Last changed Time-stamp: <2006-03-02 22:32:02 ivo> */
/*
                  minimum free energy consensus
                  RNA secondary structure prediction
                  with maximum distance base pairs

                  c Ivo Hofacker, Stephan Bernhart

                  Vienna RNA package
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "params.h"
#include "ribo.h"
#include "alifold.h"
#include "fold.h"
#include "loop_energies.h"

#ifdef _OPENMP
#include <omp.h> 
#endif


/*@unused@*/
static char rcsid[] UNUSED = "$Id: aliLfold.c,v 1.1 2007/06/23 08:49:57 ivo Exp $";


#define PAREN

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */
#define MAXSECTORS      500     /* dimension for a backtrack array */
#define LOCALITY        0.      /* locality parameter for base-pairs */
#define UNIT 100
#define MINPSCORE -2 * UNIT
#define NONE -10000 /* score for forbidden pairs */

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
PRIVATE paramT        *P = NULL;
PRIVATE int           **c;       /* energy array, given that i-j pair */
PRIVATE int           *cc;      /* linear array for calculating canonical structures */
PRIVATE int           *cc1;     /*   "     "        */
PRIVATE int           *f3;      /* energy of 5' end */
PRIVATE int           **fML;     /* multi-loop auxiliary energy array */
PRIVATE int           *Fmi;     /* holds row i of fML (avoids jumps in memory) */
PRIVATE int           *DMLi;    /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
PRIVATE int           *DMLi1;   /*             MIN(fML[i+1,k]+fML[k+1,j])  */
PRIVATE int           *DMLi2;   /*             MIN(fML[i+2,k]+fML[k+1,j])  */
PRIVATE int           **pscore; /* precomputed array of pair types */
PRIVATE unsigned int  length;
PRIVATE short         **S;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
         thus we have to initialize them before usage by a seperate function!
         OR: use copyin in the PARALLEL directive!
         e.g.:
         #pragma omp parallel for copyin(pf_params)
*/
#pragma omp threadprivate(P, c, cc, cc1, f3, fML, Fmi, DMLi, DMLi1, DMLi2, pscore, length, S)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE void  initialize_aliLfold(int length, int maxdist);
PRIVATE void  free_aliL_arrays(int maxdist);
PRIVATE void  get_arrays(unsigned int size, int maxdist);
PRIVATE void  make_pscores(const char ** AS, const char *structure,int maxdist, int start);
PRIVATE int   fill_arrays(char **strings, int maxdist, char *structure);
PRIVATE char  *backtrack(char **strings, int start, int maxdist);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE void initialize_aliLfold(int length, int maxdist){
  if (length<1) nrerror("initialize_fold: argument must be greater 0");
  get_arrays((unsigned) length, maxdist);
  make_pair_matrix();
  P = scale_parameters();
}

/*--------------------------------------------------------------------------*/

PRIVATE void get_arrays(unsigned int size, int maxdist)
{
  int i;
  c       = (int **)space(sizeof(int *)*(size+1));
  fML     = (int **)space(sizeof(int *)*(size+1));
  pscore  = (int **)space(sizeof(int *)*(size+1));
  f3      = (int *) space(sizeof(int)*(size+2));  /* has to be one longer */
  cc      = (int *) space(sizeof(int)*(maxdist+5));
  cc1     = (int *) space(sizeof(int)*(maxdist+5));
  Fmi     = (int *) space(sizeof(int)*(maxdist+5));
  DMLi    = (int *) space(sizeof(int)*(maxdist+5));
  DMLi1   = (int *) space(sizeof(int)*(maxdist+5));
  DMLi2   = (int *) space(sizeof(int)*(maxdist+5));
  for (i=size; i>(int)size-maxdist-5 && i>=0; i--) {
    c[i]      = (int *) space(sizeof(int) *(maxdist+5));
    fML[i]    = (int *) space(sizeof(int) *(maxdist+5));
    pscore[i] = (int *) space(sizeof(int )*(maxdist+5));
  }

}

/*--------------------------------------------------------------------------*/

PRIVATE void free_aliL_arrays(int maxdist) {
  int i;
  for(i=0; i<maxdist+5 && i<=length; i++){
    free(c[i]);
    free(fML[i]);
    free(pscore[i]);
  }
  free(c);
  free(fML);
  free(f3);
  free(cc);
  free(cc1);
  free(pscore);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);
}

/*--------------------------------------------------------------------------*/
PUBLIC float aliLfold(char **strings, char *structure, int maxdist) {
  int length, energy, s, n_seq, i;
  length = (int) strlen(strings[0]);
  if (maxdist>length) maxdist = length;
  initialize_aliLfold(length, maxdist);

  for (s=0; strings[s]!=NULL; s++);
  n_seq = s;
  S = (short **) space(n_seq*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if (strlen(strings[s]) != length) nrerror("uneqal seqence lengths");
    S[s] = encode_sequence(strings[s], 0);
  }

  for (i=length; i>=(int)length-(int)maxdist-4 && i>0; i--)
    make_pscores((const char **) strings,structure,maxdist,i);

  energy = fill_arrays(strings, maxdist, structure);

  free_aliL_arrays(maxdist);
  return (float) energy/100.;
}

PRIVATE int fill_arrays(char **strings, int maxdist, char *structure) {
  /* fill "c", "fML" and "f3" arrays and return  optimal energy */

  int   i, j, k, length, energy;
  int   decomp, new_fML,MLenergy ;
  int   *type, type_2, tt, s, n_seq;
  /* int   bonus=0;*/

  length = (int) strlen(strings[0]);
  for (s=0; strings[s]!=NULL; s++);
  n_seq = s;
  type = (int *) space(n_seq*sizeof(int));
  for (j=0; j<maxdist+5; j++)
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
  for (j=length; j>length-maxdist-3; j--) {
    for (i=(length-maxdist-2>0)?length-maxdist-2:1 ; i<j; i++)
      c[i][j-i] = fML[i][j-i] = INF;
  }

  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */
    for (j = i+1; j<=length && j<=i+TURN; j++) {
      c[i][j-i]=fML[i][j-i]=INF;
    }
   for (j = i+TURN+1; j <= length && j <= i+maxdist; j++) {
      int p, q, psc;
      /* bonus = 0;*/
      for (s=0; s<n_seq; s++) {
        type[s] = pair[S[s][i]][S[s][j]];
        if (type[s]==0) type[s]=7;
      }

      psc = pscore[i][j-i];

      if (psc>=cv_fact*MINPSCORE) {   /* we have a pair 2 consider */
        int new_c=0, stackEnergy=INF;
        /* hairpin ----------------------------------------------*/
        for (new_c=s=0; s<n_seq; s++)
          new_c += E_Hairpin(j-i-1,type[s],S[s][i+1],S[s][j-1],strings[s]+i-1, P);
        /*--------------------------------------------------------
          check for elementary structures involving more than one
          closing pair.
          --------------------------------------------------------*/
        for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++) {
          int minq = j-i+p-MAXLOOP-2;
          if (minq<p+1+TURN) minq = p+1+TURN;
          for (q = minq; q < j; q++) {
            if (pscore[p][q-p]<MINPSCORE) continue;


            for (energy = s=0; s<n_seq; s++) {
              type_2 = pair[S[s][q]][S[s][p]]; /* q,p not p,q! */
              if (type_2 == 0) type_2 = 7;
              energy += E_IntLoop(p-i-1, j-q-1, type[s], type_2,
                                   S[s][i+1], S[s][j-1],
                                   S[s][p-1], S[s][q+1], P);
            }
            new_c = MIN2(energy+c[p][q-p], new_c);
            if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */

          } /* end q-loop */
        } /* end p-loop */


        /* multi-loop decomposition ------------------------*/
        decomp = DMLi1[j-1-(i+1)];
        if (dangles) {
          for (s=0; s<n_seq; s++) {
            tt = rtype[type[s]];
            decomp += E_MLstem(tt, S[s][j-1], S[s][i+1], P);
          }
        }
        else{
          for(s=0; s<n_seq; s++){
            tt = rtype[type[s]];
            decomp += E_MLstem(tt, -1, -1, P);
          }
        }
        MLenergy = decomp + n_seq*P->MLclosing;
        new_c = MIN2(new_c, MLenergy);

        new_c = MIN2(new_c, cc1[j-1-(i+1)]+stackEnergy);
        cc[j-i] = new_c - psc; /* add covariance bonnus/penalty */
        if (noLonelyPairs)
          c[i][j-i] = cc1[j-1-(i+1)]+stackEnergy-psc;
        else
          c[i][j-i] = cc[j-i];

      } /* end >> if (pair) << */
      else c[i][j-i] = INF;


      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/

      new_fML = fML[i+1][j-i-1]+n_seq*P->MLbase;
      new_fML = MIN2(fML[i][j-1-i]+n_seq*P->MLbase, new_fML);
      energy = c[i][j-i]/*+P->MLintern[type]*/;
      if(dangles){
        for (s=0; s<n_seq; s++) {
          energy += E_MLstem(type[s], (i > 1) ? S[s][i-1] : -1, (j < length) ? S[s][j+1] : -1, P);
        }
      }
      else{
        for (s=0; s<n_seq; s++) {
          energy += E_MLstem(type[s], -1, -1, P);
        }
      }
      new_fML = MIN2(energy, new_fML);

      /* modular decomposition -------------------------------*/

      for (decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++)
        decomp = MIN2(decomp, Fmi[k-i]+fML[k+1][j-k-1]);

      DMLi[j-i] = decomp;               /* store for use in ML decompositon */
      new_fML = MIN2(new_fML,decomp);



      fML[i][j-i] = Fmi[j-i] = new_fML;     /* substring energy */

    } /* for (j...) */

    /* calculate energies of 5' and 3' fragments */
    {
      static int do_backtrack = 0, prev_i=0;
      static char * prev=NULL;
      char *ss;
      f3[i] = f3[i+1];
      for (j=i+TURN+1; j<length && j<=i+maxdist; j++) {
        if(c[i][j-i]<INF) {
        /*        if (c[j+1]<INF) {*/
          energy = f3[j+1]+c[i][j-i];
          if(dangles){
            for(s = 0; s < n_seq; s++){
              tt = pair[S[s][i]][S[s][j]];
              if(tt==0) tt=7;
              energy += E_ExtLoop(tt, (i>1) ? S[s][i-1] : -1, S[s][j+1], P);
            }
          }
          else{
            for(s = 0; s < n_seq; s++){
              tt = pair[S[s][i]][S[s][j]];
              if(tt==0) tt=7;
              energy += E_ExtLoop(tt, -1, -1, P);
            }
          }
          f3[i] = MIN2(f3[i], energy);
        }
      }
      if (length<=i+maxdist) {
        j=length;
        if(c[i][j-i]<INF) {
          energy = c[i][j-i];
          if(dangles){
            for (s=0; s<n_seq; s++) {
              tt = pair[S[s][i]][S[s][j]];
              if(tt==0) tt=7;
              energy += E_ExtLoop(tt, (i>1) ? S[s][i-1] : -1, -1, P);
            }
          }
          else{
            for (s=0; s<n_seq; s++) {
              tt = pair[S[s][i]][S[s][j]];
              if(tt==0) tt=7;
              energy += E_ExtLoop(tt, -1, -1, P);
            }
          }
          f3[i] = MIN2(f3[i], energy);
        }
      }
      /* backtrack partial structure */
      if (i+maxdist<length) {
        if (f3[i] != f3[i+1]) do_backtrack=1;
        else if (do_backtrack) {
          ss =  backtrack(strings, i+1 , maxdist+1);
          if ((prev_i+strlen(prev)>i+1+strlen(ss))) {
            char *outstr;
            printf("%s (%6.2f) %4d - %4d\n",
                   prev, (f3[prev_i]-f3[prev_i + strlen(prev)-1])/(100.*n_seq), prev_i,prev_i + strlen(prev)-1);
            outstr=strndup(strings[0]+prev_i-1, strlen(prev));
            printf("%s \n", outstr);
            free(outstr);
          }
          free(prev); prev=ss; prev_i = i+1;
          do_backtrack=0;
        }
      }
      if (i+maxdist==length) {
        prev =  backtrack(strings, i , maxdist);
        prev_i = i; do_backtrack=0;
      }
      if (i==1) {
        char *outstr;
        if (prev) {
          printf("%s (%6.2f) %4d - %4d\n", prev,
                         (f3[prev_i]-f3[prev_i + strlen(prev)-1])/(100.*n_seq), prev_i,prev_i + strlen(prev)-1);
          outstr=strndup(strings[0]+prev_i-1, strlen(prev));
          printf("%s \n", outstr);
        }
        if ((f3[prev_i] != f3[1]) || (!prev)) {
          ss =  backtrack(strings, i , maxdist);
          printf("%s (%6.2f) %4d\n", ss, (f3[1]-f3[1 + strlen(ss)-1])/(100.*n_seq), 1);
          outstr=strndup(strings[0], strlen(ss));
          printf("%s \n", outstr);
          free(ss);
        }
        if (prev) free(prev);
        free(outstr);
      }
    }
    {
      int ii, *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=0; j< maxdist+5; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
      if (i<=length-maxdist-4) {
        c[i-1] = c[i+maxdist+4]; c[i+maxdist+4] = NULL;
        fML[i-1] = fML[i+maxdist+4]; fML[i+maxdist+4]=NULL;
        pscore[i-1] = pscore[i+maxdist+4]; pscore[i+maxdist+4] = NULL;
        if (i>1)  make_pscores((const char**) strings, structure, maxdist, i-1);
        for (ii=0; ii<maxdist+5; ii++) {
          c[i-1][ii] = fML[i-1][ii] = INF;
        }
      }
    }
  }

  return f3[1];
}

PRIVATE char * backtrack(char **strings, int start, int maxdist) {
  /*------------------------------------------------------------------
    trace back through the "c", "f3" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This is fast, since only few structure elements are recalculated.
    ------------------------------------------------------------------*/
  sect  sector[MAXSECTORS];   /* backtracking sectors */

  int   i, j, k, energy;
  int   *type, type_2, tt, n_seq;
  /*int   bonus;*/
  int   s=0, ss;
  char *structure;
  for (s=0; strings[s]!=NULL; s++);
  n_seq = s;
  type = (int *) space(n_seq*sizeof(int));
  s=0;
  length = strlen(strings[0]);
  sector[++s].i = start;
  sector[s].j = MIN2(length, start+maxdist+1);
  sector[s].ml = (backtrack_type=='M') ? 1 : ((backtrack_type=='C')?2:0);

  structure = (char *) space((MIN2(length-start, maxdist)+3)*sizeof(char));
  for (i=0; i<=MIN2(length-start, maxdist); i++) structure[i] = '.';

  while (s>0) {
    int ml, fij, cij, traced, i1, j1, d3, d5, mm, p, q, jj=0;
    int canonical = 1;     /* (i,j) closes a canonical structure */
    i  = sector[s].i;
    j  = sector[s].j;
    ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
                              occur in the fML- (1) or in the f-array (0) */
    if (ml==2) {
      structure[i-start] = '(';
      structure[j-start] = ')';
      goto repeat1;
    }

    if (j < i+TURN+1) continue; /* no more pairs in this interval */

    fij = (ml)? fML[i][j-i] : f3[i];

    if (ml == 0) { /* backtrack in f3 */

      if (fij == f3[i+1]) {
        sector[++s].i = i+1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }
      /* i is paired. Find pairing partner */
      for (k=i+TURN+1,traced=0; k<=j; k++) {
        int cc;
        jj = k+1;
        cc = c[i][k-(i)];
        if (cc<INF) {
          if(dangles){
            for (ss=0; ss<n_seq; ss++) {
              type[ss] = pair[S[ss][i]][S[ss][k]];
              if (type[ss]==0) type[ss]=7;
              cc += E_ExtLoop(type[ss], (i>1) ? S[ss][i-1] : -1, (k<length) ? S[ss][k+1] : -1, P);
            }
          }
          else{
            for (ss=0; ss<n_seq; ss++) {
              type[ss] = pair[S[ss][i]][S[ss][k]];
              if (type[ss]==0) type[ss]=7;
              cc += E_ExtLoop(type[ss], -1, -1, P);
            }
          }
          if (fij == cc + f3[k+1]) traced=i;
        }
        if (traced) break;
      }

      if (!traced) nrerror("backtrack failed in f3");
      if (j==length) { /* backtrack only one component, unless j==length */
        sector[++s].i = jj;
        sector[s].j   = j;
        sector[s].ml  = ml;
      }
      i=traced; j=k;
      structure[i-start] = '('; structure[j-start] = ')';
      goto repeat1;
    }
    else { /* trace back in fML array */
      if (fML[i][j-1-i]+n_seq*P->MLbase == fij) {  /* 3' end is unpaired */
        sector[++s].i = i;
        sector[s].j   = j-1;
        sector[s].ml  = ml;
        continue;
      }
      if (fML[i+1][j-(i+1)]+n_seq*P->MLbase == fij) { /* 5' end is unpaired */
        sector[++s].i = i+1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }

      cij = c[i][j-i] ;
      if(dangles){
        for (ss=0; ss<n_seq; ss++) {
          tt  = pair[S[ss][i]][S[ss][j]];
          if (tt==0) tt=7;
          cij += E_ExtLoop(tt, (i>1) ? S[ss][i-1] : -1, (j<length) ? S[ss][j+1] : -1, P);
        }
      }
      else{
        for (ss=0; ss<n_seq; ss++) {
          tt  = pair[S[ss][i]][S[ss][j]];
          if (tt==0) tt=7;
          cij += E_ExtLoop(tt, -1, -1, P);
        }
      }

      if(fij==cij){
        /* found a pair */
        structure[i-start] = '('; structure[j-start] = ')';
        goto repeat1;
      }

      for (k = i+1+TURN; k <= j-2-TURN; k++)
        if (fij == (fML[i][k-i]+fML[k+1][j-(k+1)]))
          break;

      sector[++s].i = i;
      sector[s].j   = k;
      sector[s].ml  = ml;
      sector[++s].i = k+1;
      sector[s].j   = j;
      sector[s].ml  = ml;

      if (k>j-2-TURN) nrerror("backtrack failed in fML");
      continue;
    }

  repeat1:

    /*----- begin of "repeat:" -----*/
    if (canonical)  cij = c[i][j-i];

    for (ss=0; ss<n_seq; ss++) {
      type[ss] = pair[S[ss][i]][S[ss][j]];
      if (type[ss]==0) type[ss] = 7;
    }

    /*    bonus = 0;*/

    if (noLonelyPairs)
      if (cij == c[i][j-i]) {
        /* (i.j) closes canonical structures, thus
           (i+1.j-1) must be a pair                */
        for (ss=0; ss<n_seq; ss++) {
          type_2 = pair[S[ss][j-1]][S[ss][i+1]];  /* j,i not i,j */
          if (type_2==0) type_2 = 7;
          cij -= P->stack[type[ss]][type_2];
        }
        cij += pscore[i][j-i];
        structure[i+1-start] = '('; structure[j-1-start] = ')';
        i++; j--;
        canonical=0;
        goto repeat1;
      }
    canonical = 1;
    cij += pscore[i][j-i];

    {
      int cc=0;
      for (ss=0; ss<n_seq; ss++)
        cc += E_Hairpin(j-i-1, type[ss], S[ss][i+1], S[ss][j-1], strings[ss]+i-1, P);
      if (cij == cc) /* found hairpin */
        continue;
    }

    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
      int minq;
      minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+TURN) minq = p+1+TURN;
      for (q = j-1; q >= minq; q--) {
        if (c[p][q-p]>=INF) continue;
         for (ss=energy=0; ss<n_seq; ss++) {
          type_2 = pair[S[ss][q]][S[ss][p]];  /* q,p not p,q */
          if (type_2==0) type_2 = 7;
          energy += E_IntLoop(p-i-1, j-q-1, type[ss], type_2,
                               S[ss][i+1], S[ss][j-1],
                               S[ss][p-1], S[ss][q+1], P);
        }
        traced = (cij == energy+c[p][q-p]);
        if (traced) {
          structure[p-start] = '(';
          structure[q-start] = ')';
          i = p, j = q;
          goto repeat1;
        }
      }
    }

    /* end of repeat: --------------------------------------------------*/

    /* (i.j) must close a multi-loop */
    mm = n_seq*P->MLclosing;
    if(dangles){
      for (ss=0; ss<n_seq; ss++) {
        tt = rtype[type[ss]];
        mm += E_MLstem(tt, S[ss][j-1],S[ss][i+1], P);
      }
    }
    else{
      for (ss=0; ss<n_seq; ss++) {
        tt = rtype[type[ss]];
        mm += E_MLstem(tt, S[ss][j-1],S[ss][i+1], P);
      }
    }
    i1 = i+1; j1 = j-1;
    sector[s+1].ml  = sector[s+2].ml = 1;

    for (k = i+TURN+2; k < j-TURN-2; k++){
      if(cij == fML[i+1][k-(i+1)] + fML[k+1][j-1-(k+1)] + mm) break;
    }
    if (k<=j-3-TURN){ /* found the decomposition */
      sector[++s].i = i1;
      sector[s].j   = k;
      sector[++s].i = k+1;
      sector[s].j   = j1;
    } else {
        nrerror("backtracking failed in repeat");
    }

  }
  if (start+maxdist<length) {
    for (i=strlen(structure)-1; i>0 && structure[i-1] == '.'; i--)
      structure[i] = '\0';
  }
  return structure;
}

/*---------------------------------------------------------------------------*/
PRIVATE double cov_score(const char ** AS, int i, int j, float **dm) {
  int n_seq,k,l,s,score;
  int pfreq[8]={0,0,0,0,0,0,0,0};
  for (n_seq=0; AS[n_seq]!=NULL; n_seq++);
  for (s=0; s<n_seq; s++) {
    int type;
    if (S[s][i]==0 && S[s][j]==0) type = 7; /* gap-gap  */
    else {
      if ((AS[s][i] == '~')||(AS[s][j] == '~')) type = 7;
      else type = pair[S[s][i]][S[s][j]];
    }

    pfreq[type]++;
  }
  if (pfreq[0]*2>n_seq)
    return NONE;
  else
    for (k=1,score=0; k<=6; k++) /* ignore pairtype 7 (gap-gap) */
      for (l=k+1; l<=6; l++)
        /* scores for replacements between pairtypes    */
        /* consistent or compensatory mutations score 1 or 2  */
        score += pfreq[k]*pfreq[l]*dm[k][l];

  /* counter examples score -1, gap-gap scores -0.25   */
  return cv_fact * ((UNIT*score)/n_seq - nc_fact*UNIT*(pfreq[0] + pfreq[7]*0.25));
}

PRIVATE void make_pscores(const char ** AS,
                          const char *structure, int maxd, int i) {
  /* calculate co-variance bonus for each pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */
  int n,j,k,l,n_seq;
  float **dm;
  n=S[0][0];  /* length of seqs */
  for (n_seq=0; AS[n_seq]!=NULL; n_seq++);

  int olddm[7][7]={{0,0,0,0,0,0,0}, /* hamming distance between pairs */
                {0,0,2,2,1,2,2} /* CG */,
                {0,2,0,1,2,2,2} /* GC */,
                {0,2,1,0,2,1,2} /* GU */,
                {0,1,2,2,0,2,1} /* UG */,
                {0,2,2,1,2,0,2} /* AU */,
                {0,2,2,2,1,2,0} /* UA */};
  if (ribo) {
    if (RibosumFile !=NULL) dm=readribosum(RibosumFile);
    else dm=get_ribosum_slice(AS,n_seq,i-1,n);
  }
  else { /*use usual matrix*/
    dm=(float **)space(7*sizeof(float*));
    for (k=0; k<7;k++) {
      dm[k]=(float *)space(7*sizeof(float));
      for (l=0; l<7; l++)
        dm[k][l] = (float) olddm[k][l];
    }
  }

  /*first allocate space:*/
  pscore[i]=(int *)space((maxd+5)*sizeof(int));
  /*  pscore[start]-=start;*/
  /*fill pscore[start], too close*/
  for (j=i+1; (j<i+TURN+1) && (j<=n); j++) {
    pscore[i][j-i] = NONE;
  }
  for (j=i+TURN+1; ((j<=n) && (j<=i+maxd)); j++) {
    pscore[i][j-i] = cov_score(AS, i, j, dm);
  }

  if (noLonelyPairs) { /* remove unwanted lonely pairs */
    int type, otype=0, ntype=0;
    for (j=i+TURN; ((j<n)&&(j<i+maxd)); j++) {
      if ((i>1) && (j<n)) otype = cov_score(AS, i-1, j+1, dm);
      type=pscore[i][j-i];
      if (i<n) ntype=pscore[i+1][j-1-(i+1)];
      else ntype=NONE;

      if ((otype<-4*UNIT)&&(ntype<-4*UNIT))  /* worse than 2 counterex */
        pscore[i][j-i] = NONE; /* i.j can only form isolated pairs */
    }
  }

  if (fold_constrained&&(structure!=NULL)) {
    int psij, hx, *stack;
    stack = (int *) space(sizeof(int)*(n+1));
    hx=psij=0;
    /* for(hx=0, j=i+TURN; ((j<=i+maxd)&&(j<=n)); j++) {*/
    switch (structure[i-1]) {
    case 'x': /* can't pair */
      for (l=i+TURN+1; l<=i+maxd; l++) pscore[i][l-i] = NONE;
      break;
    case '(':
        hx=1;
        psij=1;
        for (l=i+1; l<=i+maxd; l++) {
          switch (structure[l-1]) {
          case '(':
            hx++;
            pscore[i][l-i] = NONE;
            break;
          case ')':
            hx--;
            if (hx!=0) pscore[i][l-i] = NONE;
            break;
          default:
            pscore[i][l-i] = NONE;
          }
          /* fallthrough */
                }
    case ')':
      for (l=i+TURN+1; l<=i+maxd; l++) pscore[i][l-i] = NONE;
      break;
    case '>':
      for (l=i+TURN+1; l<=i+maxd; l++) pscore[i][l-i] = NONE;
      break;

    }
    if (!psij) for (l=i+1; l<=i+maxd; l++) { /*no '(' constraint on i*/
      switch (structure[l-1]) {
      case '(':
        pscore[i][l-i] = NONE;
        break;
      case '<':
        pscore[i][l-i] = NONE;
        break;
      case 'x':
        pscore[i][l-i] = NONE;
        break;
      case ')':
         pscore[i][l-i] = NONE;
        break;
      }
    }
    if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in constraint string");
    }
    free(stack);
  }
}
