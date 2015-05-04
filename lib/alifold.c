/* Last changed Time-stamp: <2009-02-24 15:17:17 ivo> */
/*
                  minimum free energy folding
                  for a set of aligned sequences

                  c Ivo Hofacker

                  Vienna RNA package
*/

/**
*** \file alifold.c
**/


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "fold.h"
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "params.h"
#include "ribo.h"
#include "alifold.h"
#include "loop_energies.h"

#ifdef _OPENMP
#include <omp.h>
#endif


/*@unused@*/
static char rcsid[] UNUSED = "$Id: alifold.c,v 1.18 2009/02/27 16:25:54 ivo Exp $";

#define PAREN
#define STACK_BULGE1  1     /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1     /* new asymetry penalty */
#define MAXSECTORS    500   /* dimension for a backtrack array */
#define LOCALITY        0.  /* locality parameter for base-pairs */
#define UNIT 100
#define MINPSCORE -2 * UNIT

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
PUBLIC double cv_fact=1.; /* should be made static to not interfere with other threads */
PUBLIC double nc_fact=1.; /* should be made static to not interfere with other threads */

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE short           **S;
PRIVATE short           **S5;               /*S5[s][i] holds next base 5' of i in sequence s*/
PRIVATE short           **S3;               /*Sl[s][i] holds next base 3' of i in sequence s*/
PRIVATE char            **Ss;
PRIVATE unsigned short  **a2s;
PRIVATE paramT          *P;
PRIVATE int             *indx;              /* index for moving in the triangle matrices c[] and fMl[]*/
PRIVATE int             *c;                 /* energy array, given that i-j pair */
PRIVATE int             *cc;                /* linear array for calculating canonical structures */
PRIVATE int             *cc1;               /*   "     "        */
PRIVATE int             *f5;                /* energy of 5' end */
PRIVATE int             *fML;               /* multi-loop auxiliary energy array */
PRIVATE int             *Fmi;               /* holds row i of fML (avoids jumps in memory) */
PRIVATE int             *DMLi;              /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
PRIVATE int             *DMLi1;             /*             MIN(fML[i+1,k]+fML[k+1,j])  */
PRIVATE int             *DMLi2;             /*             MIN(fML[i+2,k]+fML[k+1,j])  */
PRIVATE int             *pscore;            /* precomputed array of pair types */
PRIVATE int             init_length = -1;
PRIVATE sect            sector[MAXSECTORS]; /* stack of partial structures for backtracking */
PRIVATE bondT           *base_pair2;
PRIVATE int             circular;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
         thus we have to initialize them before usage by a seperate function!
         OR: use copyin in the PARALLEL directive!
         e.g.:
         #pragma omp parallel for copyin(P, ...)
*/
#pragma omp threadprivate(S, S5, S3, Ss, a2s, P, indx, c, cc, cc1, f5, fML, Fmi, DMLi, DMLi1, DMLi2,\
                          pscore, init_length, sector, base_pair2)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE void  init_alifold(int length);
PRIVATE void  get_arrays(unsigned int size);
PRIVATE void  make_pscores(const short *const *S, const char **AS, int n_seq, const char *structure);
PRIVATE int   fill_arrays(const char **strings);
PRIVATE void  backtrack(const char **strings, int s);

PRIVATE void    energy_of_alistruct_pt(const char **sequences,short * ptable, int n_seq, int *energy);
PRIVATE void    stack_energy_pt(int i, const char **sequences, short *ptable,  int n_seq, int *energy);
PRIVATE int     ML_Energy_pt(int i, int n_seq, short *pt);
PRIVATE int     EL_Energy_pt(int i, int n_seq, short *pt);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

/* unsafe function that will be replaced by a threadsafe companion in the future */
PRIVATE void init_alifold(int length){

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  if (length < 1) nrerror("initialize_fold: argument must be greater 0");
  free_alifold_arrays();
  get_arrays((unsigned) length);
  init_length = length;

  indx = get_indx((unsigned)length);

  update_alifold_params();
}

PRIVATE void get_arrays(unsigned int size){
  indx    = (int *) space(sizeof(int)*(size+1));
  c       = (int *) space(sizeof(int)*((size*(size+1))/2+2));
  fML     = (int *) space(sizeof(int)*((size*(size+1))/2+2));
  pscore  = (int *) space(sizeof(int)*((size*(size+1))/2+2));
  f5      = (int *) space(sizeof(int)*(size+2));
  cc      = (int *) space(sizeof(int)*(size+2));
  cc1     = (int *) space(sizeof(int)*(size+2));
  Fmi     = (int *) space(sizeof(int)*(size+1));
  DMLi    = (int *) space(sizeof(int)*(size+1));
  DMLi1   = (int *) space(sizeof(int)*(size+1));
  DMLi2   = (int *) space(sizeof(int)*(size+1));
  if(base_pair2) free(base_pair2);
  base_pair2 = (bondT *) space(sizeof(bondT)*(1+size/2));
}

PUBLIC  void  free_alifold_arrays(void){
  if(indx)        free(indx);
  if(c)           free(c);
  if(fML)         free(fML);
  if(f5)          free(f5);
  if(cc)          free(cc);
  if(cc1)         free(cc1);
  if(pscore)      free(pscore);
  if(base_pair2)  free(base_pair2);
  if(Fmi)         free(Fmi);
  if(DMLi)        free(DMLi);
  if(DMLi1)       free(DMLi1);
  if(DMLi2)       free(DMLi2);
  if(P)           free(P);
  indx = c = fML = f5 = cc = cc1 = Fmi = DMLi = DMLi1 = DMLi2 = NULL;
  pscore      = NULL;
  base_pair2  = NULL;
  P           = NULL;
  init_length = 0;
}


PUBLIC void alloc_sequence_arrays(const char **sequences, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss, int circ){
  unsigned int s, n_seq, length;
  if(sequences[0] != NULL){
    length = strlen(sequences[0]);
    for (s=0; sequences[s] != NULL; s++);
    n_seq = s;
    *S    = (short **)          space(n_seq * sizeof(short *));
    *S5   = (short **)          space(n_seq * sizeof(short *));
    *S3   = (short **)          space(n_seq * sizeof(short *));
    *a2s  = (unsigned short **) space(n_seq * sizeof(unsigned short *));
    *Ss   = (char **)           space(n_seq * sizeof(char *));
    for (s=0; s<n_seq; s++) {
      if(strlen(sequences[s]) != length) nrerror("uneqal seqence lengths");
      (*S5)[s]  = (short *)         space((length + 2) * sizeof(short));
      (*S3)[s]  = (short *)         space((length + 2) * sizeof(short));
      (*a2s)[s] = (unsigned short *)space((length + 2) * sizeof(unsigned short));
      (*Ss)[s]  = (char *)          space((length + 2) * sizeof(char));
      (*S)[s]   = (short *)         space((length + 2) * sizeof(short));
      encode_ali_sequence(sequences[s], (*S)[s], (*S5)[s], (*S3)[s], (*Ss)[s], (*a2s)[s], circ);
    }
  }
  else nrerror("alloc_sequence_arrays: no sequences in the alignment!");
}

PUBLIC void free_sequence_arrays(unsigned int n_seq, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss){
  unsigned int s;
  for (s=0; s<n_seq; s++) {
    free((*S)[s]);
    free((*S5)[s]);
    free((*S3)[s]);
    free((*a2s)[s]);
    free((*Ss)[s]);
  }
  free(*S);
  free(*S5);
  free(*S3);
  free(*a2s);
  free(*Ss);
}

PUBLIC void update_alifold_params(void){
  P = scale_parameters();
  make_pair_matrix();
  if (init_length < 0) init_length=0;
}

PUBLIC float alifold(const char **strings, char *structure){
  int  length, energy, s, n_seq;

  circular = 0;
  length = (int) strlen(strings[0]);

#ifdef _OPENMP
  /* always init everything since all global static variables are uninitialized when entering a thread */
  init_alifold(length);
#else
  if (length>init_length) init_alifold(length);
#endif
  if (fabs(P->temperature - temperature)>1e-6)  update_alifold_params();

  for (s=0; strings[s]!=NULL; s++);
  n_seq = s;

  alloc_sequence_arrays(strings, &S, &S5, &S3, &a2s, &Ss, circular);
  make_pscores((const short **) S, strings, n_seq, structure);

  energy = fill_arrays((const char **)strings);

  backtrack((const char **)strings, 0);

#ifdef PAREN
  parenthesis_structure(structure, base_pair2, length);
#else
  letter_structure(structure, base_pair2, length);
#endif

  /*
  *  Backward compatibility:
  *  This block may be removed if deprecated functions
  *  relying on the global variable "base_pair" vanishs from within the package!
  */
  {
    if(base_pair) free(base_pair);
    base_pair = (bondT *)space(sizeof(bondT) * (1+length/2));
    memcpy(base_pair, base_pair2, sizeof(bondT) * (1+length/2));
  }

  free_sequence_arrays(n_seq, &S, &S5, &S3, &a2s, &Ss);

  if (backtrack_type=='C')
    return (float) c[indx[length]+1]/(n_seq*100.);
  else if (backtrack_type=='M')
    return (float) fML[indx[length]+1]/(n_seq*100.);
  else
    return (float) f5[length]/(n_seq*100.);
}


/**
*** the actual forward recursion to fill the energy arrays
**/
PRIVATE int fill_arrays(const char **strings) {
  int   i, j, k, p, q, length, energy, new_c;
  int   decomp, MLenergy, new_fML;
  int   s, n_seq, *type, type_2, tt;

  /* count number of sequences */
  for (n_seq=0; strings[n_seq]!=NULL; n_seq++);

  type = (int *) space(n_seq*sizeof(int));
  length = strlen(strings[0]);

  /* init energies */
  for (j=1; j<=length; j++){
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
    for (i=(j>TURN?(j-TURN):1); i<j; i++) {
      c[indx[j]+i] = fML[indx[j]+i] = INF;
    }
  }

  /* begin recursions */
  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */
    for (j = i+TURN+1; j <= length; j++) {
      int ij, psc;
      ij = indx[j]+i;

      for (s=0; s<n_seq; s++) {
        type[s] = pair[S[s][i]][S[s][j]];
        if (type[s]==0) type[s]=7;
      }

      psc = pscore[indx[j]+i];
      if (psc>=MINPSCORE) {   /* a pair to consider */
        int stackEnergy = INF;
        /* hairpin ----------------------------------------------*/


        for (new_c=s=0; s<n_seq; s++) {
          if ((a2s[s][j-1]-a2s[s][i])<3) new_c+=600;
          else  new_c += E_Hairpin(a2s[s][j-1]-a2s[s][i],type[s],S3[s][i],S5[s][j],Ss[s]+(a2s[s][i-1]), P);
       }
        /*--------------------------------------------------------
          check for elementary structures involving more than one
          closing pair.
          --------------------------------------------------------*/

        for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++) {
          int minq = j-i+p-MAXLOOP-2;
          if (minq<p+1+TURN) minq = p+1+TURN;
          for (q = minq; q < j; q++) {
            if (pscore[indx[q]+p]<MINPSCORE) continue;
            for (energy = s=0; s<n_seq; s++) {
              type_2 = pair[S[s][q]][S[s][p]]; /* q,p not p,q! */
              if (type_2 == 0) type_2 = 7;
              energy += E_IntLoop(a2s[s][p-1]-a2s[s][i], a2s[s][j-1]-a2s[s][q], type[s], type_2,
                                   S3[s][i], S5[s][j],
                                   S5[s][p], S3[s][q], P);
            }
            new_c = MIN2(new_c, energy + c[indx[q]+p]);
            if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */
          } /* end q-loop */
        } /* end p-loop */

        /* multi-loop decomposition ------------------------*/
        decomp = DMLi1[j-1];
        if(dangles){
          for(s=0; s<n_seq; s++){
            tt = rtype[type[s]];
            decomp += E_MLstem(tt, S5[s][j], S3[s][i], P);
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

        new_c = MIN2(new_c, cc1[j-1]+stackEnergy);

        cc[j] = new_c - psc; /* add covariance bonnus/penalty */
        if (noLonelyPairs)
          c[ij] = cc1[j-1]+stackEnergy-psc;
        else
          c[ij] = cc[j];

      } /* end >> if (pair) << */

      else c[ij] = INF;
      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/

      new_fML = fML[ij+1]+n_seq*P->MLbase;
      new_fML = MIN2(fML[indx[j-1]+i]+n_seq*P->MLbase, new_fML);
      energy = c[ij];
      if(dangles){
        for (s=0; s<n_seq; s++) {
          energy += E_MLstem(type[s], S5[s][i], S3[s][j], P);
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
        decomp = MIN2(decomp, Fmi[k]+fML[indx[j]+k+1]);

      DMLi[j] = decomp;               /* store for use in ML decompositon */
      new_fML = MIN2(new_fML,decomp);

      /* coaxial stacking deleted */

      fML[ij] = Fmi[j] = new_fML;     /* substring energy */
    } /* END for j */

    {
      int *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=1; j<=length; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
    }
  } /* END for i */
  /* calculate energies of 5' and 3' fragments */

  f5[TURN+1]=0;
  for (j=TURN+2; j<=length; j++) {
    f5[j] = f5[j-1];
    if (c[indx[j]+1]<INF) {
      energy = c[indx[j]+1];
      if(dangles){
        for(s = 0; s < n_seq; s++){
          tt = pair[S[s][1]][S[s][j]];
          if(tt==0) tt=7;
          energy += E_ExtLoop(tt, -1, (j<length) ? S3[s][j] : -1, P);
        }
      }
      else{
        for(s = 0; s < n_seq; s++){
          tt = pair[S[s][1]][S[s][j]];
          if(tt==0) tt=7;
          energy += E_ExtLoop(tt, -1, -1, P);
        }
      }
      f5[j] = MIN2(f5[j], energy);
    }
    for (i=j-TURN-1; i>1; i--) {
      if (c[indx[j]+i]<INF) {
        energy = f5[i-1] + c[indx[j]+i];
        if(dangles){
          for(s = 0; s < n_seq; s++){
            tt = pair[S[s][i]][S[s][j]];
            if(tt==0) tt=7;
            energy += E_ExtLoop(tt, S5[s][i], (j < length) ? S3[s][j] : -1, P);
          }
        }
        else{
          for(s = 0; s < n_seq; s++){
            tt = pair[S[s][i]][S[s][j]];
            if(tt==0) tt=7;
            energy += E_ExtLoop(tt, -1, -1, P);
          }
        }
        f5[j] = MIN2(f5[j], energy);
      }
    }
  }
  free(type);
  return(f5[length]);
}

#include "alicircfold.inc"

/**
*** backtrack in the energy matrices to obtain a structure with MFE
**/
PRIVATE void backtrack(const char **strings, int s) {
  /*------------------------------------------------------------------
    trace back through the "c", "f5" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This inverts the folding procedure, hence it's very fast.
    ------------------------------------------------------------------*/
   /* normally s=0.
     If s>0 then s items have been already pushed onto the sector stack */
  int   i, j, k, p, q, length, energy;
  int   type_2, tt, mm;
  int   b=0, cov_en = 0;
  int   n_seq;
  int *type;
  length = strlen(strings[0]);
  for (n_seq=0; strings[n_seq]!=NULL; n_seq++);
  type = (int *) space(n_seq*sizeof(int));
  if (s==0) {
    sector[++s].i = 1;
    sector[s].j = length;
    sector[s].ml = (backtrack_type=='M') ? 1 : ((backtrack_type=='C')?2:0);
  }
  while (s>0) {
    int ss, ml, fij, fi, cij, traced, i1, j1, d3, d5, jj=0;
    int canonical = 1;     /* (i,j) closes a canonical structure */
    i  = sector[s].i;
    j  = sector[s].j;
    ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
                              occur in the fML- (1) or in the f-array (0) */
    if (ml==2) {
      base_pair2[++b].i = i;
      base_pair2[b].j   = j;
      cov_en += pscore[indx[j]+i];
      goto repeat1;
    }

    if (j < i+TURN+1) continue; /* no more pairs in this interval */

    fij = (ml)? fML[indx[j]+i] : f5[j];
    fi  = (ml)?(fML[indx[j-1]+i]+n_seq*P->MLbase):f5[j-1];

    if (fij == fi) {  /* 3' end is unpaired */
      sector[++s].i = i;
      sector[s].j   = j-1;
      sector[s].ml  = ml;
      continue;
    }

    if (ml == 0) { /* backtrack in f5 */
      /* j or j-1 is paired. Find pairing partner */
      for (i=j-TURN-1,traced=0; i>=1; i--) {
        int cc, en;
        jj = i-1;
        if (c[indx[j]+i]<INF) {
          en = c[indx[j]+i] + f5[i-1];
          if(dangles){
            for(ss = 0; ss < n_seq; ss++){
              type[ss] = pair[S[ss][i]][S[ss][j]];
              if (type[ss]==0) type[ss] = 7;
              en += E_ExtLoop(type[ss], (i>1) ? S5[ss][i]: -1, (j < length) ? S3[ss][j] : -1, P);
            }
          }
          else{
            for(ss = 0; ss < n_seq; ss++){
              type[ss] = pair[S[ss][i]][S[ss][j]];
              if (type[ss]==0) type[ss] = 7;
              en += E_ExtLoop(type[ss], -1, -1, P);
            }
          }
          if (fij == en) traced=j;
        }
        if (traced) break;
      }

      if (!traced) nrerror("backtrack failed in f5");
      sector[++s].i = 1;
      sector[s].j   = jj;
      sector[s].ml  = ml;

      j=traced;
      base_pair2[++b].i = i;
      base_pair2[b].j   = j;
      cov_en += pscore[indx[j]+i];
      goto repeat1;
    }
    else { /* trace back in fML array */
      if (fML[indx[j]+i+1]+n_seq*P->MLbase == fij) { /* 5' end is unpaired */
        sector[++s].i = i+1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }

      cij = c[indx[j]+i];
      if(dangles){
        for(ss = 0; ss < n_seq; ss++){
          tt = pair[S[ss][i]][S[ss][j]];
          if(tt==0) tt=7;
          cij += E_MLstem(tt, S5[ss][i], S3[ss][j], P);
        }
      }
      else{
        for(ss = 0; ss < n_seq; ss++){
          tt = pair[S[ss][i]][S[ss][j]];
          if(tt==0) tt=7;
          cij += E_MLstem(tt, -1, -1, P);
        }
      }

      if (fij==cij){
        /* found a pair */
        base_pair2[++b].i = i;
        base_pair2[b].j   = j;
        cov_en += pscore[indx[j]+i];
        goto repeat1;
      }

      for (k = i+1+TURN; k <= j-2-TURN; k++)
        if (fij == (fML[indx[k]+i]+fML[indx[j]+k+1]))
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
    if (canonical)  cij = c[indx[j]+i];

    for (ss=0; ss<n_seq; ss++) {
      type[ss] = pair[S[ss][i]][S[ss][j]];
      if (type[ss]==0) type[ss] = 7;
    }

    if (noLonelyPairs)
      if (cij == c[indx[j]+i]) {
        /* (i.j) closes canonical structures, thus
           (i+1.j-1) must be a pair                */
        for (ss=0; ss<n_seq; ss++) {
          type_2 = pair[S[ss][j-1]][S[ss][i+1]];  /* j,i not i,j */
          if (type_2==0) type_2 = 7;
          cij -= P->stack[type[ss]][type_2];
        }
        cij += pscore[indx[j]+i];
        base_pair2[++b].i = i+1;
        base_pair2[b].j   = j-1;
        cov_en += pscore[indx[j-1]+i+1];
        i++; j--;
        canonical=0;
        goto repeat1;
      }
    canonical = 1;
    cij += pscore[indx[j]+i];

    {int cc=0;
    for (ss=0; ss<n_seq; ss++) {
        if ((a2s[ss][j-1]-a2s[ss][i])<3) cc+=600;
        else cc += E_Hairpin(a2s[ss][j-1]-a2s[ss][i], type[ss], S3[ss][i], S5[ss][j], Ss[ss]+a2s[ss][i-1], P);
      }
    if (cij == cc) /* found hairpin */
      continue;
    }
    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
      int minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+TURN) minq = p+1+TURN;
      for (q = j-1; q >= minq; q--) {

        if (c[indx[q]+p]>=INF) continue;

        for (ss=energy=0; ss<n_seq; ss++) {
          type_2 = pair[S[ss][q]][S[ss][p]];  /* q,p not p,q */
          if (type_2==0) type_2 = 7;
          energy += E_IntLoop(a2s[ss][p-1]-a2s[ss][i],a2s[ss][j-1]-a2s[ss][q],
                               type[ss], type_2,
                               S3[ss][i], S5[ss][j],
                               S5[ss][p], S3[ss][q], P);

        }
        traced = (cij == energy+c[indx[q]+p]);
        if (traced) {
          base_pair2[++b].i = p;
          base_pair2[b].j   = q;
          cov_en += pscore[indx[q]+p];
          i = p, j = q;
          goto repeat1;
        }
      }
    }

    /* end of repeat: --------------------------------------------------*/

    /* (i.j) must close a multi-loop */

    mm = n_seq*P->MLclosing;
    if(dangles){
      for(ss = 0; ss < n_seq; ss++){
        tt = rtype[type[ss]];
        mm += E_MLstem(tt, S5[ss][j], S3[ss][i], P);
      }
    }
    else{
      for(ss = 0; ss < n_seq; ss++){
        tt = rtype[type[ss]];
        mm += E_MLstem(tt, -1, -1, P);
      }
    }
    i1 = i+1;
    j1 = j-1;
    sector[s+1].ml  = sector[s+2].ml = 1;

    for (k = i1+TURN+1; k < j1-TURN-1; k++){
      if(cij == fML[indx[k]+i1] + fML[indx[j1]+k+1] + mm) break;
    }

    if (k<=j-3-TURN) { /* found the decomposition */
      sector[++s].i = i1;
      sector[s].j   = k;
      sector[++s].i = k+1;
      sector[s].j   = j1;
    } else {
        nrerror("backtracking failed in repeat");
    }

  }

  /* fprintf(stderr, "covariance energy %6.2f\n", cov_en/100.);  */

  base_pair2[0].i = b;    /* save the total number of base pairs */
  free(type);
}


PUBLIC void encode_ali_sequence(const char *sequence, short *S, short *s5, short *s3, char *ss, unsigned short *as, int circular){
  unsigned int i,l;
  unsigned short p;
  l     = strlen(sequence);
  S[0]  = (short) l;
  s5[0] = s5[1] = 0;

  /* make numerical encoding of sequence */
  for(i=1; i<=l; i++){
    short ctemp;
    ctemp=(short) encode_char(toupper(sequence[i-1]));
    S[i]= ctemp ;
  }

  if (oldAliEn){
    /* use alignment sequences in all energy evaluations */
    ss[0]=sequence[0];
    for(i=1; i<l; i++){
      s5[i] = S[i-1];
      s3[i] = S[i+1];
      ss[i] = sequence[i];
      as[i] = i;
    }
    ss[l]   = sequence[l];
    as[l]   = l;
    s5[l]   = S[l-1];
    s3[l]   = 0;
    S[l+1]  = S[1];
    s5[1]   = 0;
    if (circular) {
      s5[1]   = S[l];
      s3[l]   = S[1];
      ss[l+1] = S[1];
    }
  }
  else{
    if(circular){
      for(i=l; i>0; i--){
        char c5;
        c5 = sequence[i-1];
        if ((c5=='-')||(c5=='_')||(c5=='~')||(c5=='.')) continue;
        s5[1] = S[i];
        break;
      }
      for (i=1; i<=l; i++) {
        char c3;
        c3 = sequence[i-1];
        if ((c3=='-')||(c3=='_')||(c3=='~')||(c3=='.')) continue;
        s3[l] = S[i];
        break;
      }
    }
    else  s5[1]=s3[l]=0;

    for(i=1,p=0; i<=l; i++){
      char c5;
      c5 = sequence[i-1];
      if ((c5=='-')||(c5=='_')||(c5=='~')||(c5=='.'))
        s5[i+1]=s5[i];
      else { /* no gap */
        ss[p++]=sequence[i-1]; /*start at 0!!*/
        s5[i+1]=S[i];
      }
      as[i]=p;
    }
    for (i=l; i>=1; i--) {
      char c3;
      c3 = sequence[i-1];
      if ((c3=='-')||(c3=='_')||(c3=='~')||(c3=='.'))
        s3[i-1]=s3[i];
      else
        s3[i-1]=S[i];
    }
  }
}

PRIVATE void make_pscores(const short *const* S, const char **AS,
                          int n_seq, const char *structure) {
  /* calculate co-variance bonus for each pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */
#define NONE -10000 /* score for forbidden pairs */
  int n,i,j,k,l,s;

  int olddm[7][7]={{0,0,0,0,0,0,0}, /* hamming distance between pairs */
                  {0,0,2,2,1,2,2} /* CG */,
                  {0,2,0,1,2,2,2} /* GC */,
                  {0,2,1,0,2,1,2} /* GU */,
                  {0,1,2,2,0,2,1} /* UG */,
                  {0,2,2,1,2,0,2} /* AU */,
                  {0,2,2,2,1,2,0} /* UA */};

  float **dm;
  n=S[0][0];  /* length of seqs */
  if (ribo) {
    if (RibosumFile !=NULL) dm=readribosum(RibosumFile);
    else dm=get_ribosum(AS,n_seq,n);
  }
  else { /*use usual matrix*/
    dm=(float **)space(7*sizeof(float*));
    for (i=0; i<7;i++) {
      dm[i]=(float *)space(7*sizeof(float));
      for (j=0; j<7; j++)
        dm[i][j] = (float) olddm[i][j];
    }
  }

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
          else type = pair[S[s][i]][S[s][j]];
        }
        pfreq[type]++;
      }
      if (pfreq[0]*2+pfreq[7]>n_seq) { pscore[indx[j]+i] = NONE; continue;}
      for (k=1,score=0; k<=6; k++) /* ignore pairtype 7 (gap-gap) */
        for (l=k; l<=6; l++)
          score += pfreq[k]*pfreq[l]*dm[k][l];
      /* counter examples score -1, gap-gap scores -0.25   */
      pscore[indx[j]+i] = cv_fact *
        ((UNIT*score)/n_seq - nc_fact*UNIT*(pfreq[0] + pfreq[7]*0.25));
    }
  }

  if (noLonelyPairs) /* remove unwanted pairs */
    for (k=1; k<n-TURN-1; k++)
      for (l=1; l<=2; l++) {
        int type,ntype=0,otype=0;
        i=k; j = i+TURN+l;
        type = pscore[indx[j]+i];
        while ((i>=1)&&(j<=n)) {
          if ((i>1)&&(j<n)) ntype = pscore[indx[j+1]+i-1];
          if ((otype<cv_fact*MINPSCORE)&&(ntype<cv_fact*MINPSCORE))  /* too many counterexamples */
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

/*--------New scoring part-----------------------------------*/

PUBLIC int get_mpi(char *Alseq[], int n_seq, int length, int *mini) {
  int i, j,k;
  float ident=0;
  int pairnum=0;
  int sumident=0;
  float minimum=1.;
  for(j=0; j<n_seq-1; j++)
    for(k=j+1; k<n_seq; k++) {
      ident=0;
      for (i=1; i<=length; i++){
        if (Alseq[k][i]==Alseq[j][i]) ident++;
        pairnum++;
      }
      if ((ident/length)<minimum) minimum=ident/(float)length;
      sumident+=ident;
    }
  mini[0]=(int)(minimum*100.);
  if (pairnum>0)   return (int) (sumident*100/pairnum);
  else return 0;

}

PUBLIC float **readribosum(char *name){

  float **dm;
  char *line;
  FILE *fp;
  int i=0;
  int who=0;
  float a,b,c,d,e,f;
  int translator[7]={0,5,1,2,3,6,4};

  fp=fopen(name,"r");
  dm=(float **)space(7*sizeof(float*));
  for (i=0; i<7;i++) {
    dm[i]=(float *)space(7*sizeof(float));
  }
  while(1) { /*bisma hoit fertisch san*/
    line=get_line(fp);
    if (*line=='#') continue;
    i=0;
    i=sscanf(line,"%f %f %f %f %f %f",&a,&b,&c,&d,&e,&f);
    if (i==0) break;
    dm[translator[++who]][translator[1]]=a;
    dm[translator[who]][translator[2]]=b;
    dm[translator[who]][translator[3]]=c;
    dm[translator[who]][translator[4]]=d;
    dm[translator[who]][translator[5]]=e;
    dm[translator[who]][translator[6]]=f;
    free(line);
    if (who==6) break;
  }
  fclose(fp);
  return dm;
}


PUBLIC  float energy_of_alistruct(const char **sequences, const char *structure, int n_seq, float *energy){
  int new=0;
  unsigned int s, n;
  short *pt;

  short           **tempS;
  short           **tempS5;     /*S5[s][i] holds next base 5' of i in sequence s*/
  short           **tempS3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
  char            **tempSs;
  unsigned short  **tempa2s;

  int *tempindx;
  int *temppscore;

  int en_struct[2];

  if(sequences[0] != NULL){
    n = (unsigned int) strlen(sequences[0]);
    update_alifold_params();

    /*save old memory*/
    tempS = S; tempS3 = S3; tempS5 = S5; tempSs = Ss; tempa2s = a2s;
    tempindx = indx; temppscore = pscore;

    alloc_sequence_arrays(sequences, &S, &S5, &S3, &a2s, &Ss, 0);
    pscore  = (int *) space(sizeof(int)*((n+1)*(n+2)/2));
    indx    = get_indx(n);
    make_pscores((const short *const*)S, sequences, n_seq, NULL);
    new     = 1;

    pt = make_pair_table(structure);
    energy_of_alistruct_pt(sequences,pt, n_seq, &(en_struct[0]));

    free(pt);
    energy[0] = (float)en_struct[0]/(float)(100*n_seq);
    energy[1] = (float)en_struct[1]/(float)(100*n_seq);

    free(pscore);
    free(indx);
    free_sequence_arrays(n_seq, &S, &S5, &S3, &a2s, &Ss);

    /* restore old memory */
    S = tempS; S3 = tempS3; S5 = tempS5; Ss = tempSs; a2s = tempa2s;
    indx = tempindx; pscore = temppscore;
  }
  else nrerror("energy_of_alistruct(): no sequences in alignment!");

  return energy[0];
}

PRIVATE void energy_of_alistruct_pt(const char **sequences, short *pt, int n_seq, int *energy){
  int i, length;

  length = S[0][0];
  energy[0] =  backtrack_type=='M' ? ML_Energy_pt(0, n_seq, pt) : EL_Energy_pt(0, n_seq, pt);
  energy[1] = 0;
  for (i=1; i<=length; i++) {
    if (pt[i]==0) continue;
    stack_energy_pt(i, sequences, pt, n_seq, energy);
    i=pt[i];
  }
}

PRIVATE void stack_energy_pt(int i, const char **sequences, short *pt, int n_seq, int *energy)
{
  /* calculate energy of substructure enclosed by (i,j) */
  int ee= 0;
  int j, p, q, s;
  int *type = (int *) space(n_seq*sizeof(int));

  j = pt[i];
  for (s=0; s<n_seq; s++) {
    type[s] = pair[S[s][i]][S[s][j]];
    if (type[s]==0) {
    type[s]=7;
    }
  }
  p=i; q=j;
  while (p<q) { /* process all stacks and interior loops */
    int type_2;
    while (pt[++p]==0);
    while (pt[--q]==0);
    if ((pt[q]!=(short)p)||(p>q)) break;
    ee=0;
    for (s=0; s<n_seq; s++) {
      type_2 = pair[S[s][q]][S[s][p]];
      if (type_2==0) {
        type_2=7;
      }
      ee += E_IntLoop(a2s[s][p-1]-a2s[s][i], a2s[s][j-1]-a2s[s][q], type[s], type_2, S3[s][i], S5[s][j], S5[s][p], S3[s][q], P);
    }
    energy[0] += ee;
    energy[1] += pscore[indx[j]+i];
    i=p; j=q;
    for (s=0; s<n_seq; s++) {
      type[s] = pair[S[s][i]][S[s][j]];
      if (type[s]==0) type[s]=7;
    }
  }  /* end while */

  /* p,q don't pair must have found hairpin or multiloop */

  if (p>q) {
    ee=0;/* hair pin */
    for (s=0; s< n_seq; s++) {
      if ((a2s[s][j-1]-a2s[s][i])<3) ee+=600;
      else ee += E_Hairpin(a2s[s][j-1]-a2s[s][i], type[s], S3[s][i], S5[s][j], Ss[s]+(a2s[s][i-1]), P);
    }
    energy[0] += ee;
    energy[1] += pscore[indx[j]+i];
    free(type);
    return;
  }
  /* (i,j) is exterior pair of multiloop */
  energy[1] += pscore[indx[j]+i];
  while (p<j) {
    /* add up the contributions of the substructures of the ML */
    stack_energy_pt(p, sequences, pt, n_seq, energy);
    p = pt[p];
    /* search for next base pair in multiloop */
    while (pt[++p]==0);
  }
  energy[0] += ML_Energy_pt(i, n_seq, pt);
  free(type);
}



PRIVATE int ML_Energy_pt(int i, int n_seq, short *pt){
  /* i is the 5'-base of the closing pair */

  int   energy, tt, i1, j, p, q, u, s;
  short d5, d3;

  j = pt[i];
  i1  = i;
  p   = i+1;
  u   = 0;
  energy = 0;

  do{ /* walk around the multi-loop */
    /* hop over unpaired positions */
    while (p < j && pt[p]==0) p++;
    if(p >= j) break;
    /* get position of pairing partner */
    q  = pt[p];
    /* memorize number of unpaired positions? no, we just approximate here */
    u += p-i1-1;

    for (s=0; s< n_seq; s++) {
      /* get type of base pair P->q */
      tt = pair[S[s][p]][S[s][q]];
      if (tt==0) tt=7;
      d5 = dangles && (a2s[s][p]>1) && (tt!=0) ? S5[s][p] : -1;
      d3 = dangles && (a2s[s][q]<a2s[s][S[0][0]]) ? S3[s][q] : -1;
      energy += E_MLstem(tt, d5, d3, P);
    }
    i1  = q;
    p   = q + 1;
  }while(1);

  if(i > 0){
    energy  += P->MLclosing * n_seq;
    if(dangles){
      for (s=0; s<n_seq; s++){
        tt = pair[S[s][j]][S[s][i]];
        if (tt==0) tt=7;
        energy += E_MLstem(tt, S5[s][j], S3[s][i], P);
      }
    }
    else{
      for (s=0; s<n_seq; s++){
        tt = pair[S[s][j]][S[s][i]];
        if (tt==0) tt=7;
        energy += E_MLstem(tt, -1, -1, P);
      }
    }
  }
  u += j - i1 - 1;
  energy += u * P->MLbase * n_seq;
  return energy;
}

PRIVATE int EL_Energy_pt(int i, int n_seq, short *pt){
  int   energy, tt, i1, j, p, q, s;
  short d5, d3;

  j = pt[0];

  p   = i+1;
  energy = 0;

  do{ /* walk along the backbone */
    /* hop over unpaired positions */
    while (p < j && pt[p]==0) p++;
    if(p == j) break; /* no more stems */
    /* get position of pairing partner */
    q  = pt[p];
    for (s=0; s< n_seq; s++) {
      /* get type of base pair P->q */
      tt = pair[S[s][p]][S[s][q]];
      if (tt==0) tt=7;
      d5 = dangles && (a2s[s][p]>1) && (tt!=0) ? S5[s][p] : -1;
      d3 = dangles && (a2s[s][q]<a2s[s][S[0][0]]) ? S3[s][q] : -1;
      energy += E_ExtLoop(tt, d5, d3, P);
    }
    p   = q + 1;
  }while(p < j);

  return energy;
}

