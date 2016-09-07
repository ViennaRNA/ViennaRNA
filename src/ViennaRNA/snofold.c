/*
                  minimum free energy
                  RNA secondary structure prediction

                  c Ivo Hofacker, Chrisoph Flamm
                  original implementation by
                  Walter Fontana

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
#include "ViennaRNA/structure_utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/snofold.h"
#include "ViennaRNA/loop_energies.h"

#ifdef __GNUC__
#define INLINE inline
#else
#define INLINE
#endif

#define PAREN

#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */

/*@unused@*/
PRIVATE void  get_arrays(unsigned int size);
/* PRIVATE int   stack_energy(int i, const char *string); */
PRIVATE void  make_ptypes(const short *S, const char *structure);
PRIVATE void encode_seq(const char *sequence);
PRIVATE void  backtrack(const char *sequence, int s);
PRIVATE int   fill_arrays(const char *sequence, const int max_asymm, const int threshloop,
                          const int min_s2, const int max_s2, const int half_stem, const int max_half_stem);
/*@unused@*/


/* alifold */
PRIVATE void alisnoinitialize_fold(const int length);
PRIVATE void make_pscores(const short *const* S, const char *const* AS,int n_seq, const char *structure);
PRIVATE int   *pscore;  /* precomputed array of pair types */
PRIVATE short **Sali;
PRIVATE int alifill_arrays(const char **string, const int max_asymm, const int threshloop, 
                           const int min_s2, const int max_s2, const int half_stem, 
                           const int max_half_stem);
PRIVATE void aliget_arrays(unsigned int size);
PRIVATE short * aliencode_seq(const char *sequence);
PRIVATE int alibacktrack(const char **strings, int s);

#define UNIT 100
#define MINPSCORE -2 * UNIT
/* end alifold */

#define MAXSECTORS      500     /* dimension for a backtrack array */
#define LOCALITY        0.      /* locality parameter for base-pairs */

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define MAX2(A, B)      ((A) > (B) ? (A) : (B))
#define SAME_STRAND(I,J) (((I)>=cut_point)||((J)<cut_point))

PRIVATE vrna_param_t *P = NULL;

PRIVATE int *indx = NULL; /* index for moving in the triangle matrices c[] and fMl[]*/

PRIVATE int   *c = NULL;       /* energy array, given that i-j pair */
PRIVATE int   *cc = NULL;      /* linear array for calculating canonical structures */
PRIVATE int   *cc1 = NULL;     /*   "     "        */
PRIVATE int   *Fmi = NULL;     /* holds row i of fML (avoids jumps in memory) */
PRIVATE int   *DMLi = NULL;    /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
PRIVATE int   *DMLi1 = NULL;   /*             MIN(fML[i+1,k]+fML[k+1,j])  */
PRIVATE int   *DMLi2 = NULL;   /*             MIN(fML[i+2,k]+fML[k+1,j])  */
PRIVATE char  *ptype = NULL;   /* precomputed array of pair types */
PRIVATE short *S = NULL, *S1 = NULL;
PRIVATE int    init_length=-1;
PRIVATE int    *mLoop = NULL; /*contains the minimum of c for a xy range*/
PRIVATE folden **foldlist = NULL;
PRIVATE folden **foldlist_XS = NULL;

PRIVATE int     *BP = NULL; /* contains the structure constrainsts: BP[i]
                        -1: | = base must be paired
                        -2: < = base must be paired with j<i
                        -3: > = base must be paired with j>i
                        -4: x = base must not pair
                        positive int: base is paired with int      */


static sect sector[MAXSECTORS]; /* stack of partial structures for backtracking */

PRIVATE char  alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
/* needed by cofold/eval */
/* PRIVATE int cut_in_loop(int i); */
/* PRIVATE int min_hairpin = TURN; */

/* some definitions to take circfold into account...        */
/* PRIVATE int   *fM2 = NULL;*/        /* fM2 = multiloop region with exactly two stems, extending to 3' end        */
PUBLIC        int   Fc, FcH, FcI, FcM; /* parts of the exterior loop energies                        */
/*--------------------------------------------------------------------------*/

void snoinitialize_fold(const int length)
{
  unsigned int n;
  if (length<1) vrna_message_error("snoinitialize_fold: argument must be greater 0");
  if (init_length>0) snofree_arrays(length);
  get_arrays((unsigned) length);
  init_length=length;
  for (n = 1; n <= (unsigned) length; n++)
    indx[n] = (n*(n-1)) >> 1;        /* n(n-1)/2 */

  snoupdate_fold_params();
}

PRIVATE void alisnoinitialize_fold(const int length)
{
  unsigned int n;
  if (length<1) vrna_message_error("snoinitialize_fold: argument must be greater 0");
  if (init_length>0) snofree_arrays(length);
  aliget_arrays((unsigned) length);
  make_pair_matrix();
  init_length=length;
  for (n = 1; n <= (unsigned) length; n++)
    indx[n] = (n*(n-1)) >> 1;        /* n(n-1)/2 */
  snoupdate_fold_params();
}  


/*--------------------------------------------------------------------------*/

PRIVATE void get_arrays(unsigned int size)
{
  indx = (int *) vrna_alloc(sizeof(int)*(size+1));
  c     = (int *) vrna_alloc(sizeof(int)*((size*(size+1))/2+2));
  mLoop = (int *) vrna_alloc(sizeof(int)*((size*(size+1))/2+2));

  ptype = (char *) vrna_alloc(sizeof(char)*((size*(size+1))/2+2));
  cc    = (int *) vrna_alloc(sizeof(int)*(size+2));
  cc1   = (int *) vrna_alloc(sizeof(int)*(size+2));
  Fmi   = (int *) vrna_alloc(sizeof(int)*(size+1));
  DMLi  = (int *) vrna_alloc(sizeof(int)*(size+1));
  DMLi1  = (int *) vrna_alloc(sizeof(int)*(size+1));
  DMLi2  = (int *) vrna_alloc(sizeof(int)*(size+1));
  if (base_pair) free(base_pair);
  base_pair = (vrna_bp_stack_t *) vrna_alloc(sizeof(vrna_bp_stack_t)*(1+size/2));
  /* extra array(s) for circfold() */
}

PRIVATE void aliget_arrays(unsigned int size)
{
  indx = (int *) vrna_alloc(sizeof(int)*(size+1));
  c     = (int *) vrna_alloc(sizeof(int)*((size*(size+1))/2+2));
  mLoop = (int *) vrna_alloc(sizeof(int)*((size*(size+1))/2+2));
  pscore = (int *) vrna_alloc(sizeof(int)*((size*(size+1))/2+2));
  ptype = (char *) vrna_alloc(sizeof(char)*((size*(size+1))/2+2));
  cc    = (int *) vrna_alloc(sizeof(int)*(size+2));
  cc1   = (int *) vrna_alloc(sizeof(int)*(size+2));
  Fmi   = (int *) vrna_alloc(sizeof(int)*(size+1));
  DMLi  = (int *) vrna_alloc(sizeof(int)*(size+1));
  DMLi1  = (int *) vrna_alloc(sizeof(int)*(size+1));
  DMLi2  = (int *) vrna_alloc(sizeof(int)*(size+1));
  if (base_pair) free(base_pair);
  base_pair = (vrna_bp_stack_t *) vrna_alloc(sizeof(vrna_bp_stack_t)*(1+size/2));
  /* extra array(s) for circfold() */
}




/*--------------------------------------------------------------------------*/


void snofree_arrays(const int length)
{
  free(indx); free(c);free(cc); free(cc1);
  free(ptype);free(mLoop);
  int i;
  for(i=length;i>-1;i--){
    while(foldlist[i]!=NULL){
      folden *n = foldlist[i];
      foldlist[i] = foldlist[i]->next;
      free(n);
    }
    free(foldlist[i]);
  }
  free(foldlist);
  for(i=length;i>-1;i--){
    while(foldlist_XS[i]!=NULL){
      folden *n = foldlist_XS[i];
      foldlist_XS[i] = foldlist_XS[i]->next;
      free(n);
    }
    free(foldlist_XS[i]);
  }
  free(foldlist_XS);
  free(base_pair); base_pair=NULL; free(Fmi);
  free(DMLi); free(DMLi1);free(DMLi2);
  free(BP);
  init_length=0;
}

void alisnofree_arrays(const int length)
{
  free(indx); free(c);free(cc); free(cc1);
  free(ptype);free(mLoop);free(pscore);
  int i;
  for(i=length-1;i>-1;i--){
    while(foldlist[i]!=NULL){
      folden *n = foldlist[i];
      foldlist[i] = foldlist[i]->next;
      free(n);
    }
    free(foldlist[i]);
  }
  free(foldlist);
  free(base_pair); base_pair=NULL; free(Fmi);
  free(DMLi); free(DMLi1);free(DMLi2);
  free(BP);
  init_length=0;
}

/*--------------------------------------------------------------------------*/

void snoexport_fold_arrays(int **indx_p, int **mLoop_p, int **cLoop,  folden ***fold_p, folden ***fold_p_XS) {
  /* make the DP arrays available to routines such as subopt() */
  *indx_p = indx; *mLoop_p = mLoop;
  *cLoop = c; *fold_p = foldlist;*fold_p_XS=foldlist_XS;
}

/* void alisnoexport_fold_arrays(int **indx_p, int **mLoop_p, int **cLoop, folden ***fold_p, int **pscores) { */
/*   /\* make the DP arrays available to routines such as subopt() *\/ */
/*   *indx_p = indx; *mLoop_p = mLoop; */
/*   *cLoop = c; *fold_p = foldlist; */
/*   *pscores=pscore; */
/* } */

/*--------------------------------------------------------------------------*/










int snofold(const char *string, char *structure, const int max_assym, const int threshloop, 
              const int min_s2, const int max_s2, const int half_stem, const int max_half_stem) {
  int length, energy, bonus, bonus_cnt, s;
  
  /* Variable initialization */
  bonus = 0;
  bonus_cnt = 0;
  s     = 0;
  length = (int) strlen(string);
  
  S   = encode_sequence(string, 0);
  S1  = encode_sequence(string, 1);

  
  /* structure = (char *) vrna_alloc((unsigned) length+1); */
  
  if (length>init_length) snoinitialize_fold(length);
  else if (fabs(P->temperature - temperature)>1e-6) snoupdate_fold_params();

  

  /* encode_seq(string); */
  BP  = (int *)vrna_alloc(sizeof(int)*(length+2));
  make_ptypes(S, structure);
  energy=fill_arrays(string, max_assym, threshloop, min_s2, max_s2, half_stem, max_half_stem);
  backtrack(string, s);

  free(structure);
  free(S); free(S1); /* free(BP); */
  return energy;
}

PRIVATE void make_pscores(const short *const* S, const char *const* AS,
                          int n_seq, const char *structure) {
  /* calculate co-variance bonus for each pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */
#define NONE -10000 /* score for forbidden pairs */
  int n,i,j,k,l,s,score;
  int dm[7][7]={{0,0,0,0,0,0,0}, /* hamming distance between pairs */
                       {0,0,2,2,1,2,2} /* CG */,
                {0,2,0,1,2,2,2} /* GC */,
                {0,2,1,0,2,1,2} /* GU */,
                {0,1,2,2,0,2,1} /* UG */,
                {0,2,2,1,2,0,2} /* AU */,
                {0,2,2,2,1,2,0} /* UA */};
  n=Sali[0][0];  /* length of seqs */
  for (i=1; i<n; i++) {
    for (j=i+1; (j<i+TURN+1) && (j<=n); j++)
      pscore[indx[j]+i] = NONE;
    for (j=i+TURN+1; j<=n; j++) {
      int pfreq[8]={0,0,0,0,0,0,0,0};
      for (s=0; s<n_seq; s++) {
        int type;
        if (Sali[s][i]==0 && Sali[s][j]==0) type = 7; /* gap-gap  */
        else {
          if ((AS[s][i] == '~')||(AS[s][j] == '~')) type = 7;
          else type = pair[Sali[s][i]][Sali[s][j]];
        }

        pfreq[type]++;
      }
      if (pfreq[0]*2>n_seq) { pscore[indx[j]+i] = NONE; continue;}
      for (k=1,score=0; k<=6; k++) /* ignore pairtype 7 (gap-gap) */
        for (l=k+1; l<=6; l++)
          /* scores for replacements between pairtypes    */
          /* consistent or compensatory mutations score 1 or 2  */
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
          if ((otype<-4*UNIT)&&(ntype<-4*UNIT))  /* worse than 2 counterex */
            pscore[indx[j]+i] = NONE; /* i.j can only form isolated pairs */
          otype =  type;
          type  = ntype;
          i--; j++;
        }
      }


  if (fold_constrained&&(structure!=NULL)) {
    int psij, hx, hx2, *stack, *stack2;
    stack = (int *) vrna_alloc(sizeof(int)*(n+1));
    stack2 = (int *) vrna_alloc(sizeof(int)*(n+1));

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
          vrna_message_error("unbalanced brackets in constraints\n%s", structure);
        }
        i = stack2[--hx2];
        pscore[indx[j]+i]=NONE;
        break;
      case ')':
        if (hx<=0) {
          vrna_message_error("unbalanced brackets in constraints\n%s", structure);
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
      vrna_message_error("unbalanced brackets in constraint string\n%s", structure);
    }
    free(stack); free(stack2);
  }
}

float alisnofold(const char **strings, const int max_assym, const int threshloop, 
              const int min_s2, const int max_s2, const int half_stem, const int max_half_stem) {
  int s,n_seq, length, energy;
  char * structure;
  length = (int) strlen(strings[0]);
  /* structure = (char *) vrna_alloc((unsigned) length+1); */
  structure = NULL;
  if (length>init_length) alisnoinitialize_fold(length);
  if (fabs(P->temperature - temperature)>1e-6) snoupdate_fold_params();
  for (s=0; strings[s]!=NULL; s++);
  n_seq = s;
  Sali = (short **) vrna_alloc(n_seq*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if (strlen(strings[s]) != length) vrna_message_error("uneqal seqence lengths");
    Sali[s] = aliencode_seq(strings[s]);
  }
  make_pscores((const short **) Sali, (const char *const *) strings, n_seq, structure);
  energy=alifill_arrays(strings, max_assym, threshloop, min_s2, max_s2, half_stem, max_half_stem);
  alibacktrack((const char **)strings, 0);
  for (s=0; s<n_seq; s++) free(Sali[s]);
  free(Sali);
  /* free(structure); */
  /*  free(S)*/; free(S1); /* free(BP); */
  return (float) energy/100.;
}

PRIVATE int alifill_arrays(const char **strings, const int max_asymm, const int threshloop, 
                           const int min_s2, const int max_s2, const int half_stem, 
                           const int max_half_stem) {

  int   i, j, length, energy;
  /* int   decomp, new_fML; */
  int   *type, type_2;
  int   bonus,n_seq,s;

  
  for (n_seq=0; strings[n_seq]!=NULL; n_seq++);
  type = (int *) vrna_alloc(n_seq*sizeof(int));
  length = strlen(strings[0]);
  bonus=0;
  /*   max_separation = (int) ((1.-LOCALITY)*(double)(length-2));*/ /* not in use */
  
    /* for (i=(j>TURN?(j-TURN):1); i<j; i++) { */
    /* } */
    for (i = (length)-TURN-1; i >= 1; i--) { /* i,j in [1..length] */
      for (j = i+TURN+1; j <= length; j++) {
        int p, q, ij,psc;
        ij = indx[j]+i;
        for (s=0; s<n_seq; s++) {
          type[s] = pair[Sali[s][i]][Sali[s][j]];
          if (type[s]==0) type[s]=7;
        }
        psc = pscore[indx[j]+i];
        if (psc>=MINPSCORE) {   /* we have a pair */
        int new_c=0, stackEnergy=INF; /* seems that new_c immer den minimum von cij enthaelt */
        /* hairpin ----------------------------------------------*/
        
        for (new_c=s=0; s<n_seq; s++)
          new_c += E_Hairpin(j-i-1,type[s],Sali[s][i+1],Sali[s][j-1],strings[s]+i-1,P);
        /*--------------------------------------------------------      
          check for elementary structures involving more than one
          closing pair (interior loop).
          --------------------------------------------------------*/      
        
        for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++) {
          int minq = j-i+p-MAXLOOP-2;
          if (minq<p+1+TURN) minq = p+1+TURN;
          for (q = minq; q < j; q++) {
            if (pscore[indx[q]+p]<MINPSCORE) continue;
            if(abs((p-i) - (j-q)) > max_asymm) continue;
            for (energy = s=0; s<n_seq; s++) {
              type_2 = pair[Sali[s][q]][Sali[s][p]]; /* q,p not p,q! */
              if (type_2 == 0) type_2 = 7;
              energy += E_IntLoop(p-i-1, j-q-1, type[s], type_2,
                                  Sali[s][i+1], Sali[s][j-1],
                                  Sali[s][p-1], Sali[s][q+1],P);
            }
            new_c = MIN2(energy+c[indx[q]+p], new_c);
            if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */
            
          } /* end q-loop */
        } /* end p-loop */
        
        /* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */
        
        new_c = MIN2(new_c, cc1[j-1]+stackEnergy);
        cc[j] = new_c - psc; /* add covariance bonnus/penalty */
        c[ij]=cc[j];
        } /* end >> if (pair) << */
        else c[ij] = INF;
        /* done with c[i,j], now compute fML[i,j] */
        /* free ends ? -----------------------------------------*/
        
      }

    {
      int *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=1; j<=length; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
    }
  }
  foldlist = (folden**) vrna_alloc((length)*sizeof(folden*));

  for(i=0; i< length; i++){
    foldlist[i]=(folden*) vrna_alloc(sizeof(folden));
    foldlist[i]->next=NULL;
    foldlist[i]->k=INF+1;
    foldlist[i]->energy=INF;

  }
  folden* head; /* we save the stem loop information in a list like structure */

  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */
    int max_k, min_k;
    max_k = MIN2(length-min_s2,i+max_half_stem+1);
    min_k = MAX2(i+half_stem+1, length-max_s2);
    for (j = i+TURN+1; j <= length; j++) {
      int ij,a,b;
      ij = indx[j]+i;
      for(a=0; a< MISMATCH ;a++){
        for(b=0; b< MISMATCH ; b++){
          mLoop[ij]=MIN2(mLoop[ij],  c[indx[j-a]+i+b]);

        }
      }
      if(mLoop[ij]>=n_seq*threshloop){
        mLoop[ij]=INF;        
      }
      else{
        if(j>=min_k-1 && j < max_k){ /* comment if out to recover the known behaviour */
          head = (folden*) vrna_alloc(sizeof(folden));
          head->k=j;
          head->energy=mLoop[ij];
          head->next=foldlist[i];
          foldlist[i] = head;

        }
      }
    }
    
  }
  free(type);
  return mLoop[indx[length]+1];/* mLoop;  */
}

PRIVATE int alibacktrack(const char **strings, int s) {

  /*------------------------------------------------------------------
    trace back through the "c", "f5" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This is fast, since only few structure elements are recalculated.
    ------------------------------------------------------------------*/

  /* normally s=0.
     If s>0 then s items have been already pushed onto the sector stack */
  int   i, j, length, energy;/* , new; */
  int   type_2;
  int   bonus,n_seq,*type;  int   b=0,cov_en = 0;

  length = strlen(strings[0]);
  for (n_seq=0; strings[n_seq]!=NULL; n_seq++);
  type = (int *) vrna_alloc(n_seq*sizeof(int));
  if (s==0) {
    sector[++s].i = 1;
    sector[s].j = length;
    sector[s].ml = 2 ; 
  }
  while (s>0) {
    int ml, ss, cij, traced, i1, j1, p, q;
    int canonical = 1;     /* (i,j) closes a canonical structure */
    i  = sector[s].i;
    j  = sector[s].j;
    ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
                              occur in the fML- (1) or in the f-array (0) */
    if (ml==2) {
      base_pair[++b].i = i;
      base_pair[b].j   = j;
      goto repeat1;
    }

    if (j < i+TURN+1) continue; /* no more pairs in this interval */


  repeat1:

    /*----- begin of "repeat:" -----*/
    if (canonical)  cij = c[indx[j]+i];
    for (ss=0; ss<n_seq; ss++) {
      type[ss] = pair[Sali[ss][i]][Sali[ss][j]];
      if (type[ss]==0) type[ss] = 7;
    }
    bonus = 0;
    
    if (noLonelyPairs)
      if (cij == c[indx[j]+i]) {
        /* (i.j) closes canonical structures, thus
           (i+1.j-1) must be a pair                */
        for (ss=0; ss<n_seq; ss++) {
          type_2 = pair[Sali[ss][j-1]][Sali[ss][i+1]];  /* j,i not i,j */
          if (type_2==0) type_2 = 7;
          cij -= P->stack[type[ss]][type_2];
        }
        cij += pscore[indx[j]+i];
        base_pair[++b].i = i+1;
        base_pair[b].j   = j-1;
        cov_en += pscore[indx[j-1]+i+1];
        i++; j--;
        canonical=0;
        goto repeat1;
      }
    canonical = 1;
    cij += pscore[indx[j]+i];
    {int cc=0;
      for (ss=0; ss<n_seq; ss++)
        cc += E_Hairpin(j-i-1, type[ss], Sali[ss][i+1], Sali[ss][j-1], strings[ss]+i-1,P);
      if (cij == cc) /* found hairpin */
        continue;
    }
    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
      int minq;
      minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+TURN) minq = p+1+TURN;
      for (q = j-1; q >= minq; q--) {
        for (ss=energy=0; ss<n_seq; ss++) {
          type_2 = pair[Sali[ss][q]][Sali[ss][p]];  /* q,p not p,q */
          if (type_2==0) type_2 = 7;
          energy += E_IntLoop(p-i-1, j-q-1, type[ss], type_2,
                               Sali[ss][i+1], Sali[ss][j-1],
                              Sali[ss][p-1], Sali[ss][q+1],P);
        }
        traced = (cij == energy+c[indx[q]+p]);
        if (traced) {
          base_pair[++b].i = p;
          base_pair[b].j   = q;
          cov_en += pscore[indx[q]+p];
          i = p, j = q;
          goto repeat1;
        }
      }
    }

    /* end of repeat: --------------------------------------------------*/

    /* (i.j) must close a multi-loop */
    /* tt = rtype[type]; */
/*     mm = bonus+P->MLclosing+P->MLintern[tt]; */
/*     d5 = P->dangle5[tt][S1[j-1]]; */
/*     d3 = P->dangle3[tt][S1[i+1]]; */
    i1 = i+1; j1 = j-1;
    sector[s+1].ml  = sector[s+2].ml = 1;
    
/*      if (k<=j-3-TURN) { /\\* found the decomposition *\\/ *\/ */
/*       sector[++s].i = i1; */
/*       sector[s].j   = k; */
/*       sector[++s].i = k+1; */
/*       sector[s].j   = j1; */
/*     } /\* else { *\/ */
/*       vrna_message_error("backtracking failed in repeat"); */
/*     } */
    
  }
  base_pair[0].i = b;    /* save the total number of base pairs */
  free(type);
  return cov_en;
}

PRIVATE int fill_arrays(const char *string, const int max_asymm, const int threshloop, 
                        const int min_s2, const int max_s2, const int half_stem, const int max_half_stem) {

  int   i, j,  length, energy;
  /*   int   decomp;*/ /*, new_fML; */
  int   no_close, type, type_2;
  int   bonus;
  int min_c;
  
  min_c=INF;
  length = (int) strlen(string);
  bonus=0;
  /*   max_separation = (int) ((1.-LOCALITY)*(double)(length-2)); */ /* not in use */


  

  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */
    /* printf("i=%d\t",i);  */
    for (j = i+TURN+1; j <= length; j++) {
/*         printf("j=%d,",j); */
      int p, q, ij;
      ij = indx[j]+i;
      type = ptype[ij];
      bonus = 0;
      energy = INF;

      if ((BP[i]==j)||(BP[i]==-1)||(BP[i]==-2)) bonus -= BONUS;
      if ((BP[j]==-1)||(BP[j]==-3)) bonus -= BONUS;
      if ((BP[i]==-4)||(BP[j]==-4)) type=0;

      no_close = (((type==3)||(type==4))&&no_closingGU);

      /* if (j-i-1 > max_separation) type = 0; */ /* forces locality degree */

      if (type) {   /* we have a pair */
        int new_c=0, stackEnergy=INF; /* seems that new_c immer den minimum von cij enthaelt */
        /* hairpin ----------------------------------------------*/

        if (no_close) new_c = FORBIDDEN;
        else
          new_c = E_Hairpin(j-i-1, type, S1[i+1], S1[j-1], string+i-1,P); /* computes hair pin structure for subsequence i...j */

        /*--------------------------------------------------------      
          check for elementary structures involving more than one
          closing pair (interior loop).
          --------------------------------------------------------*/      

        for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++) {
          int minq = j-i+p-MAXLOOP-2;
          if (minq<p+1+TURN) minq = p+1+TURN;
          for (q = minq; q < j; q++) {
            
            if(abs((p-i) - (j-q)) > max_asymm) continue;
            type_2 = ptype[indx[q]+p];

            if (type_2==0) continue;
            type_2 = rtype[type_2];

            if (no_closingGU)
              if (no_close||(type_2==3)||(type_2==4))
                if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

            energy = E_IntLoop(p-i-1, j-q-1, type, type_2,
                               S1[i+1], S1[j-1], S1[p-1], S1[q+1],P);
            new_c = MIN2(energy+c[indx[q]+p], new_c);
            if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */

          } /* end q-loop */
        } /* end p-loop */




        /* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */


        new_c = MIN2(new_c, cc1[j-1]+stackEnergy);
        cc[j] = new_c;
        c[ij] = new_c;
        /*         min_c=MIN2(min_c, c[ij]); */

      } /* end >> if (pair) << */

      else c[ij] = INF;


      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/

    }

    {
      int *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=1; j<=length; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
    }
  }
  foldlist = (folden**) vrna_alloc((length+1)*sizeof(folden*));
  foldlist_XS = (folden**) vrna_alloc((length+1)*sizeof(folden*));
  /* linked list initialization*/
  for(i=0; i<=length; i++){
    foldlist[i]=(folden*) vrna_alloc(sizeof(folden));
    foldlist[i]->next=NULL;
    foldlist[i]->k=INF+1;
    foldlist[i]->energy=INF;
    foldlist_XS[i]=(folden*) vrna_alloc(sizeof(folden));
    foldlist_XS[i]->next=NULL;
    foldlist_XS[i]->k=INF+1;
    foldlist_XS[i]->energy=INF;
  }
  folden* head; /* we save the stem loop information in a list like structure */
  folden* head_XS;
  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */
    int max_k, min_k;
    max_k = MIN2(length-min_s2,i+max_half_stem+1);
    min_k = MAX2(i+half_stem+1, length-max_s2);


    for (j = i+TURN+1; j <= length; j++) {
        int ij,a,b;
              ij = indx[j]+i;
            for(a=0; a< MISMATCH ;a++){
          for(b=0; b< MISMATCH ; b++){
              mLoop[ij]=MIN2(mLoop[ij],  c[indx[j-a]+i+b]);
            /* #mLoop[ij]=MIN2(mLoop[ij], c[indx[j-2]+i]); */
            /* #mLoop[ij]=MIN2(mLoop[ij], c[indx[j]+i+1]); */
            /* #mLoop[ij]=MIN2(mLoop[ij], c[indx[j-1]+i+1]); */
            /* #mLoop[ij]=MIN2(mLoop[ij], c[indx[j-2]+i+1]); */
            /* #mLoop[ij]=MIN2(mLoop[ij], c[indx[j]+i+2]); */
            /* #mLoop[ij]=MIN2(mLoop[ij], c[indx[j-1]+i+2]); */
            /* #mLoop[ij]=MIN2(mLoop[ij], c[indx[j-2]+i+2]); */
          }
        }
        min_c = MIN2(mLoop[ij] ,min_c);
        
        if(mLoop[ij]>=threshloop){
          mLoop[ij]=INF;        
        }
        else{
          if(j>=min_k-1 && j <= max_k){ /* comment if out to recover the known behaviour */
            head = (folden*) vrna_alloc(sizeof(folden));
            head->k=j;
            head->energy=mLoop[ij];
            head->next=foldlist[i];
            foldlist[i] = head;
            head_XS = (folden*) vrna_alloc(sizeof(folden));
            head_XS->k=i;
            head_XS->energy=mLoop[ij];
            head_XS->next=foldlist_XS[j];
            foldlist_XS[j] = head_XS;            
          }
        }
    }
    
  }
/*   int count=0; */
/*    for(i=0; i< length; i++){  */
/*      folden *temp;  */
/*      temp = foldlist[i];  */
/*      while(temp->next){  */
/*        count++; */
/*        printf("count %d: i%d j%d energy %d \n", count, i, temp->k, temp->energy);  */
/*        temp=temp->next;  */
/*      }      */
/*    }  */
/*    printf("Count %d \n", count); */
/*    count=0; */
/*    for(i=length-1; i>=0; i--){  */
/*      folden *temp;  */
/*      temp = foldlist_XS[i];  */
/*      while(temp->next){  */
/*        count++; */
/*        printf("count %d: i%d j%d energy %d \n", count, temp->k,i, temp->energy);  */
/*        temp=temp->next;  */
/*      }      */
/*    }  */
/*    printf("Count %d \n", count); */
/*    return mLoop[indx[length]+1]; */ /* mLoop; */
   return min_c;
  /* printf("\nmin_array = %d\n", min_c); */
  /* return f5[length]; */
}




PRIVATE void backtrack(const char *string, int s) {

  /*------------------------------------------------------------------
    trace back through the "c", "f5" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This is fast, since only few structure elements are recalculated.
    ------------------------------------------------------------------*/

  /* normally s=0.
     If s>0 then s items have been already pushed onto the sector stack */
  int   i, j, /*k,*/ length, energy, new;
  int   no_close, type, type_2;/* , tt; */
  int   bonus;
  int   b=0;

  length = strlen(string);
  if (s==0) {
    sector[++s].i = 1;
    sector[s].j = length;
    sector[s].ml = 2 ; 
  }
  while (s>0) {
    int ml, cij, traced, i1, j1, /*d3, d5, mm,*/ p, q;
    int canonical = 1;     /* (i,j) closes a canonical structure */
    i  = sector[s].i;
    j  = sector[s].j;
    ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
                              occur in the fML- (1) or in the f-array (0) */
    if (ml==2) {
      base_pair[++b].i = i;
      base_pair[b].j   = j;
      goto repeat1;
    }

    if (j < i+TURN+1) continue; /* no more pairs in this interval */


  repeat1:

    /*----- begin of "repeat:" -----*/
    if (canonical)  cij = c[indx[j]+i];
    type = ptype[indx[j]+i];
    bonus = 0;
    if (fold_constrained) {
      if ((BP[i]==j)||(BP[i]==-1)||(BP[i]==-2)) bonus -= BONUS;
      if ((BP[j]==-1)||(BP[j]==-3)) bonus -= BONUS;
    }
    if (noLonelyPairs)
      if (cij == c[indx[j]+i]) {
        /* (i.j) closes canonical structures, thus
           (i+1.j-1) must be a pair                */
        type_2 = ptype[indx[j-1]+i+1]; type_2 = rtype[type_2];
        cij -= P->stack[type][type_2] + bonus;
        base_pair[++b].i = i+1;
        base_pair[b].j   = j-1;
        i++; j--;
        canonical=0;
        goto repeat1;
      }
    canonical = 1;
    no_close = (((type==3)||(type==4))&&no_closingGU&&(bonus==0));
    if (no_close) {
      if (cij == FORBIDDEN) continue;
    } else
      if (cij == E_Hairpin(j-i-1, type, S1[i+1], S1[j-1],string+i-1,P)+bonus)
        continue;
    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
      int minq;
      minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+TURN) minq = p+1+TURN;
      for (q = j-1; q >= minq; q--) {
        type_2 = ptype[indx[q]+p];
        if (type_2==0) continue;
        type_2 = rtype[type_2];
        if (no_closingGU)
          if (no_close||(type_2==3)||(type_2==4))
            if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */
        energy = E_IntLoop(p-i-1, j-q-1, type, type_2,
                           S1[i+1], S1[j-1], S1[p-1], S1[q+1],P);
        new = energy+c[indx[q]+p]+bonus;
        traced = (cij == new);
        if (traced) {
          base_pair[++b].i = p;
          base_pair[b].j   = q;
          i = p, j = q;
          goto repeat1;
        }
      }
    }

    /* end of repeat: --------------------------------------------------*/

    /* (i.j) must close a multi-loop */
/*     tt = rtype[type]; */
/*     mm = bonus+P->MLclosing+P->MLintern[tt]; */
/*     d5 = P->dangle5[tt][S1[j-1]]; */
/*     d3 = P->dangle3[tt][S1[i+1]]; */
    i1 = i+1; j1 = j-1;
    sector[s+1].ml  = sector[s+2].ml = 1;

/*      if (k<=j-3-TURN) { */ /* found the decomposition */
/*       sector[++s].i = i1; */
/*       sector[s].j   = k; */
/*       sector[++s].i = k+1; */
/*       sector[s].j   = j1; */
/*     } else { */
/*         vrna_message_error("backtracking failed in repeat"); */
/*     } */
/*  */
  }

  base_pair[0].i = b;    /* save the total number of base pairs */
}

char *snobacktrack_fold_from_pair(const char *sequence, int i, int j) {
  char *structure;
  sector[1].i  = i;
  sector[1].j  = j;
  sector[1].ml = 2;
  base_pair[0].i=0;
  encode_seq(sequence);
  backtrack(sequence, 1);
  structure = vrna_db_from_bp_stack(base_pair, strlen(sequence));
  free(S);free(S1);
  return structure;
}

char *alisnobacktrack_fold_from_pair(const char **strings, int i, int j, int *cov) {
  char *structure;
  int n_seq, s, length;
  length = (int) strlen(strings[0]);
  for (s=0; strings[s]!=NULL; s++);
  n_seq = s;
  sector[1].i  = i;
  sector[1].j  = j;
  sector[1].ml = 2;
  base_pair[0].i=0;
  /* encode_seq(sequence); */
  Sali = (short **) vrna_alloc(n_seq*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if (strlen(strings[s]) != length) vrna_message_error("uneqal seqence lengths");
    Sali[s] = aliencode_seq(strings[s]);
  }
  *cov=alibacktrack(strings, 1);
  structure = vrna_db_from_bp_stack(base_pair, length);
  free(S);free(S1);
  for (s=0; s<n_seq; s++) {
    free(Sali[s]);
  }
  free(Sali);
  return structure;
}



/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/


PRIVATE void encode_seq(const char *sequence) {
  unsigned int i,l;

  l = strlen(sequence);
  S = (short *) vrna_alloc(sizeof(short)*(l+2));
  S1= (short *) vrna_alloc(sizeof(short)*(l+2));
  /* S1 exists only for the special X K and I bases and energy_set!=0 */
  S[0] = (short) l;

  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    S[i]= (short) encode_char(toupper(sequence[i-1]));
    S1[i] = alias[S[i]];   /* for mismatches of nostandard bases */
  }
  /* for circular folding add first base at position n+1 and last base at
        position 0 in S1        */
  S[l+1] = S[1]; S1[l+1]=S1[1]; S1[0] = S1[l];
}

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

/*---------------------------------------------------------------------------*/

PUBLIC void snoupdate_fold_params(void)
{
  vrna_md_t md;
  if(P)
    free(P);
  set_model_details(&md);
  P = vrna_params(&md);
  make_pair_matrix();
  if (init_length < 0) init_length=0;
}

/*---------------------------------------------------------------------------*/

PRIVATE void make_ptypes(const short *S, const char *structure) {
  int n,i,j,k,l;

  n=S[0];
  for (k=1; k<n-TURN; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+TURN+l; if (j>n) continue;
      type = pair[S[i]][S[j]];
      while ((i>=1)&&(j<=n)) {
        if ((i>1)&&(j<n)) ntype = pair[S[i-1]][S[j+1]];
        if (noLonelyPairs && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */
        ptype[indx[j]+i] = (char) type;
        otype =  type;
        type  = ntype;
        i--; j++;
      }
    }

  if (fold_constrained&&(structure!=NULL)) {
    constrain_ptypes(structure, (unsigned int)n, ptype, BP, TURN, 0);
  }
}
