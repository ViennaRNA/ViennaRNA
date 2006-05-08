/* Last changed Time-stamp: <2005-10-30 21:03:40 ivo> */
/*
		  minimum free energy
		  RNA secondary structure prediction
		  with maximum distance base pairs

		  c Ivo Hofacker, Peter Stadler

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

/*@unused@*/
static char rcsid[] UNUSED = "$Id: Lfold.c,v 1.6 2006/05/08 08:29:51 ivo Exp $";


#define PAREN

#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */

PUBLIC float  Lfold(char *string, char *structure, int maxdist);
PRIVATE void   initialize_fold(int length, int maxdist);
PRIVATE void   update_fold_params(void);

/*@unused@*/
PRIVATE void  letter_structure(char *structure, int length) UNUSED;
/* PRIVATE void  parenthesis_structure(char *structure, int length); */
PRIVATE void  get_arrays(unsigned int size, int maxdist);
PRIVATE void free_arrays(int maxdist);
/* PRIVATE void  scale_parameters(void); */
/* PRIVATE int   stack_energy(int i, char *string);
   PRIVATE int   ML_Energy(int i, int is_extloop); */
PRIVATE void  make_ptypes(const short *S, int i, int maxdist, int n);
PRIVATE void  encode_seq(char *sequence);
PRIVATE char * backtrack(char *sequence, int start, int maxdist);
PRIVATE int fill_arrays(char *sequence, int maxdist);
extern  int  LoopEnergy(int n1, int n2, int type, int type_2,
			int si1, int sj1, int sp1, int sq1);
extern  int  HairpinE(int size, int type, int si1, int sj1, const char *string);

/*@unused@*/
#define MAXSECTORS      500     /* dimension for a backtrack array */
#define LOCALITY        0.      /* locality parameter for base-pairs */

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define SAME_STRAND(I,J) (((I)>=cut_point)||((J)<cut_point))

PRIVATE paramT *P = NULL;


PRIVATE int   **c;       /* energy array, given that i-j pair */
PRIVATE int   *cc;      /* linear array for calculating canonical structures */
PRIVATE int   *cc1;     /*   "     "        */
PRIVATE int   *f3;      /* energy of 5' end */
PRIVATE int   **fML;     /* multi-loop auxiliary energy array */
/* PRIVATE int   *fM1;  */   /* second ML array, only for subopt */
PRIVATE int   *Fmi;     /* holds row i of fML (avoids jumps in memory) */
PRIVATE int   *DMLi;    /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
PRIVATE int   *DMLi1;   /*             MIN(fML[i+1,k]+fML[k+1,j])  */
PRIVATE int   *DMLi2;   /*             MIN(fML[i+2,k]+fML[k+1,j])  */
PRIVATE char  **ptype;  /* precomputed array of pair types */
PRIVATE short  *S, *S1;
/* PRIVATE short  *pair_table; */
/* PRIVATE int    length; */

PRIVATE char  alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
/* needed by cofold/eval */
/* PRIVATE int cut_in_loop(int i); */
/*PRIVATE int min_hairpin = TURN; */
PRIVATE unsigned int length;
/*--------------------------------------------------------------------------*/

PRIVATE void initialize_fold(int length, int maxdist)
{
  if (length<1) nrerror("initialize_fold: argument must be greater 0");
  get_arrays((unsigned) length, maxdist);
  update_fold_params();
}

/*--------------------------------------------------------------------------*/

PRIVATE void get_arrays(unsigned int size, int maxdist)
{
  int i;
  c     = (int **) space(sizeof(int *)*(size+1));
  fML   = (int **) space(sizeof(int *)*(size+1));
  ptype = (char **) space(sizeof(char *)*(size+1));
  f3    = (int *) space(sizeof(int)*(size+2));  /* has to be one longer */
  cc    = (int *) space(sizeof(int)*(maxdist+5));
  cc1   = (int *) space(sizeof(int)*(maxdist+5));
  Fmi   = (int *) space(sizeof(int)*(maxdist+5));
  DMLi  = (int *) space(sizeof(int)*(maxdist+5));
  DMLi1  = (int *) space(sizeof(int)*(maxdist+5));
  DMLi2  = (int *) space(sizeof(int)*(maxdist+5));
  for (i=size; i>(int)size-maxdist-5 && i>=0; i--) {
    c[i] = (int *) space(sizeof(int)*(maxdist+5));
  }
  for (i=size; i>(int)size-maxdist-5 && i>=0; i--) {
    fML[i] = (int *) space(sizeof(int)*(maxdist+5));
  }
  for (i=size; i>(int)size-maxdist-5 && i>=0; i--) {
    ptype[i] = (char *) space(sizeof(char)*(maxdist+5));
  }
}

/*--------------------------------------------------------------------------*/

PRIVATE void free_arrays(int maxdist)
{
  int i;
  for (i=0; i<maxdist+5 && i<=length; i++) {
    free(c[i]); free(fML[i]); free(ptype[i]);
  }
  free(c); free(fML); free(f3); free(cc); free(cc1);
  free(ptype);

  free(base_pair); free(Fmi);
  free(DMLi); free(DMLi1);free(DMLi2);
}

/*--------------------------------------------------------------------------*/

/* PRIVATE   int   *BP; *//* contains the structure constrainsts: BP[i]
			-1: | = base must be paired
			-2: < = base must be paired with j<i
			-3: > = base must be paired with j>i
			-4: x = base must not pair
			positive int: base is paired with int      */

float Lfold(char *string, char *structure, int maxdist) {
  int i, energy, bonus=0;

  length = (int) strlen(string);
  if (maxdist>length) maxdist = length;
  initialize_fold(length, maxdist);
  if (fabs(P->temperature - temperature)>1e-6) update_fold_params();

  encode_seq(string);

  /* BP = (int *)space(sizeof(int)*(length+2)); */
  for (i=length; i>(int)length-(int)maxdist-4 && i>0; i--)
    make_ptypes(S, i, maxdist, length);

  energy = fill_arrays(string, maxdist);

  /* parenthesis_structure(structure, length); */

  free(S); free(S1);\

  energy += bonus;      /*remove bonus energies from result */
  free_arrays(maxdist);
  return (float) energy/100.;
}

PRIVATE int fill_arrays(char *string, int maxdist) {
  /* fill "c", "fML" and "f3" arrays and return  optimal energy */

  int   i, j, k, length, energy;
  int   decomp, new_fML;
  int   no_close, type, type_2, tt;
  int   bonus=0;

  length = (int) strlen(string);

  for (j=0; j<maxdist+5; j++)
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
  for (j=length; j>length-maxdist-4; j--) {
    for (i=(length-maxdist-4>0)?length-maxdist-4:1 ; i<j; i++)
      c[i][j-i] = fML[i][j-i] = INF;
  }

  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */

    for (j = i+TURN+1; j <= length && j <= i+maxdist; j++) {
      int p, q;
      bonus = 0;
      type = ptype[i][j-i];

      no_close = (((type==3)||(type==4))&&no_closingGU&&(bonus==0));

      if (type) {   /* we have a pair */
	int new_c=0, stackEnergy=INF;
	/* hairpin ----------------------------------------------*/

	if (no_close) new_c = FORBIDDEN;
	else
	  new_c = HairpinE(j-i-1, type, S1[i+1], S1[j-1], string+i-1);

	/*--------------------------------------------------------
	  check for elementary structures involving more than one
	  closing pair.
	  --------------------------------------------------------*/

	for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++) {
	  int minq = j-i+p-MAXLOOP-2;
	  if (minq<p+1+TURN) minq = p+1+TURN;
	  for (q = minq; q < j; q++) {
	    type_2 = ptype[p][q-p];

	    if (type_2==0) continue;
	    type_2 = rtype[type_2];

	    if (no_closingGU)
	      if (no_close||(type_2==3)||(type_2==4))
		if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

#if 1
	    energy = LoopEnergy(p-i-1, j-q-1, type, type_2,
				S1[i+1], S1[j-1], S1[p-1], S1[q+1]);
#else
	    /* duplicated code is faster than function call */

#endif
	    new_c = MIN2(energy+c[p][q-p], new_c);
	    if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */

	  } /* end q-loop */
	} /* end p-loop */

	/* multi-loop decomposition ------------------------*/


	if (!no_close) {
	  int MLenergy;
	  decomp = DMLi1[j-1-(i+1)];
	  if (dangles) {
	    int d3=0, d5=0;
	    tt = rtype[type];
	    d3 = P->dangle3[tt][S1[i+1]];
	    d5 = P->dangle5[tt][S1[j-1]];
	    if (dangles==2) /* double dangles */
	      decomp += d5 + d3;
	    else {          /* normal dangles */
	      decomp = MIN2(DMLi2[j-1-(i+2)]+d3+P->MLbase, decomp);
	      decomp = MIN2(DMLi1[j-2-(i+1)]+d5+P->MLbase, decomp);
	      decomp = MIN2(DMLi2[j-2-(i+2)]+d5+d3+2*P->MLbase, decomp);
	    }
	  }

	  MLenergy = P->MLclosing+P->MLintern[type]+decomp;

	  new_c = MLenergy < new_c ? MLenergy : new_c;
	}

	/* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */

	if (dangles==3) {
	  decomp = INF;
	  for (k = i+2+TURN; k < j-2-TURN; k++) {
	    type_2 = ptype[i+1][k-i-1]; type_2 = rtype[type_2];
	    if (type_2)
	      decomp = MIN2(decomp, c[i+1][k-i-1]+P->stack[type][type_2]+
			    fML[k+1][j-1-k-1]);
	    type_2 = ptype[k+1][j-1-k-1]; type_2 = rtype[type_2];
	    if (type_2)
	      decomp = MIN2(decomp, c[k+1][j-1-k-1]+P->stack[type][type_2]+
			    fML[i+1][k-i-1]);
	  }
	  /* no TermAU penalty if coax stack */
	  decomp += 2*P->MLintern[1] + P->MLclosing;
	  new_c = MIN2(new_c, decomp);
	}

	new_c = MIN2(new_c, cc1[j-1-(i+1)]+stackEnergy);
	cc[j-i] = new_c + bonus;
	if (noLonelyPairs)
	  c[i][j-i] = cc1[j-1-(i+1)]+stackEnergy+bonus;
	else
	  c[i][j-i] = cc[j-i];

      } /* end >> if (pair) << */

      else c[i][j-i] = INF;


      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/

      new_fML = fML[i+1][j-i-1]+P->MLbase;
      new_fML = MIN2(fML[i][j-1-i]+P->MLbase, new_fML);
      energy = c[i][j-i]+P->MLintern[type];
      if (dangles==2) {  /* double dangles */
	if (i>1)      energy += P->dangle5[type][S1[i-1]];
	if (j<length) energy += P->dangle3[type][S1[j+1]];
      }
      new_fML = MIN2(energy, new_fML);

      if (dangles%2==1) {  /* normal dangles */
	tt = ptype[i+1][j-i-1]; /* i+1,j */
	new_fML = MIN2(c[i+1][j-i-1]+P->dangle5[tt][S1[i]]
		       +P->MLintern[tt]+P->MLbase,new_fML);
	tt = ptype[i][j-1-i];
	new_fML = MIN2(c[i][j-1-i]+P->dangle3[tt][S1[j]]
		       +P->MLintern[tt]+P->MLbase, new_fML);
	tt = ptype[i+1][j-1-i-1];
	new_fML = MIN2(c[i+1][j-1-i-1]+P->dangle5[tt][S1[i]]+
		       P->dangle3[tt][S1[j]]+P->MLintern[tt]+2*P->MLbase, new_fML);
      }

      /* modular decomposition -------------------------------*/

      for (decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++)
	decomp = MIN2(decomp, Fmi[k-i]+fML[k+1][j-k-1]);

      DMLi[j-i] = decomp;               /* store for use in ML decompositon */
      new_fML = MIN2(new_fML,decomp);

      /* coaxial stacking */
      if (dangles==3) {
	/* additional ML decomposition as two coaxially stacked helices */
	for (decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++) {
	  type = ptype[i][k-i]; type = rtype[type];
	  type_2 = ptype[k+1][j-k-1]; type_2 = rtype[type_2];
	  if (type && type_2)
	    decomp = MIN2(decomp,
			  c[i][k-i]+c[k+1][j-k-1]+P->stack[type][type_2]);
	}

	decomp += 2*P->MLintern[1];  	/* no TermAU penalty if coax stack */
#if 0
	/* This is needed for Y shaped ML loops with coax stacking of
	   interior pairts, but backtracking will fail if activated */
	DMLi[j-i] = MIN2(DMLi[j-i], decomp);
	DMLi[j-i] = MIN2(DMLi[j-i], DMLi[j-1-i]+P->MLbase);
	DMLi[j-i] = MIN2(DMLi[j-i], DMLi1[j-(i+1)]+P->MLbase);
	new_fML = MIN2(new_fML, DMLi[j-i]);
#endif
	new_fML = MIN2(new_fML, decomp);
      }

      fML[i][j-i] = Fmi[j-i] = new_fML;     /* substring energy */

    } /* for (j...) */

    /* calculate energies of 5' and 3' fragments */
    {
      static int do_backtrack = 0, prev_i=0;
      static char * prev=NULL;
      char *ss;
      f3[i] = f3[i+1];
      for (j=i+TURN+1; j<length && j<=i+maxdist; j++) {
	type = ptype[i][j-i];
	if (type) {
	  energy = f3[j+1]+c[i][j-i];
	  if (type>2) energy += P->TerminalAU;
	  if (dangles==2) {
	    if (i>1) energy += P->dangle5[type][S1[i-1]];
	    energy += P->dangle3[type][S1[j+1]];
	  }
	  f3[i] = MIN2(f3[i], energy);
	  if (dangles%2==1) {
	    energy = c[i][j-i]+P->dangle3[type][S1[j+1]];
	    if (j+2<=length) energy += f3[j+2];
	    if (type>2) energy += P->TerminalAU;
	    f3[i] = MIN2(f3[i], energy);
	  }
	}
	type = ptype[i+1][j-i-1];
	if ((type)&&(dangles%2==1)) {
	  energy = c[i+1][j-i-1]+P->dangle5[type][S1[i]];
	  if (type>2) energy += P->TerminalAU;
	  f3[i] = MIN2(f3[i], f3[j+1]+energy);
	  energy+=P->dangle3[type][S1[j+1]];
	  if (j+1<length) energy += f3[j+2];
	  f3[i] = MIN2(f3[i], energy);
	}
      }
      if (length<=i+maxdist) {
	j=length;
	type = ptype[i][j-i];
	if (type) {
	  energy = c[i][j-i];
	  if (type>2) energy += P->TerminalAU;
	  if (dangles==2) {
	    energy += P->dangle5[type][S1[i-1]];
	  }
	  f3[i] = MIN2(f3[i], energy);
	}
	type = ptype[i+1][j-i-1];
	if ((type)&&(dangles%2==1)) {
	  energy = c[i+1][j-i-1]+P->dangle5[type][S1[i]];
	  if (type>2) energy += P->TerminalAU;
	  f3[i] = MIN2(f3[i], energy);
	}
      }
      /* backtrack partial structure */
      if (f3[i] != f3[i+1]) do_backtrack=1;
      else if (do_backtrack) {
	ss =  backtrack(string, i+1 , maxdist+1);
	if (prev) {
	  if ((i+strlen(ss)<prev_i+strlen(prev)) ||
	      strncmp(ss+prev_i-i,prev,strlen(prev))) { /* ss does not contain prev */
	    if (dangles==2)
	      printf(".%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)-1])/100., prev_i-1);
	    else
	      printf("%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)])/100., prev_i);
	  }
	  free(prev);
	}
	prev=ss; prev_i = i+1;
	do_backtrack=0;
      }
      if (i==1) {
	if (prev) {
	  if (dangles==2)
	    printf(".%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)-1])/100., prev_i-1);
	  else
	    printf("%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)])/100., prev_i);
	  free(prev); prev=NULL;
	} else do_backtrack=1;

	if (do_backtrack) {
	  ss =  backtrack(string, i , maxdist);
	  if (dangles==2)
	    printf("%s (%6.2f) %4d\n", ss, (f3[1]-f3[1+strlen(ss)-1])/100., 1);
	  else
	    printf("%s (%6.2f) %4d\n", ss, (f3[1]-f3[1+strlen(ss)])/100., 1);
	  free(ss);
	}
      }
    }
    {
      int ii, *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=0; j< maxdist+5; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
      if (i+maxdist+4<=length) {
	c[i-1] = c[i+maxdist+4]; c[i+maxdist+4] = NULL;
	fML[i-1] = fML[i+maxdist+4]; fML[i+maxdist+4]=NULL;
	ptype[i-1] = ptype[i+maxdist+4]; ptype[i+maxdist+4] = NULL;
	if (i>1) make_ptypes(S, i-1, maxdist, length);
	for (ii=0; ii<maxdist+5; ii++) {
	  c[i-1][ii] = INF;
	  fML[i-1][ii] = INF;
	}
      }
    }
  }

  return f3[1];
}

PRIVATE char * backtrack(char *string, int start, int maxdist) {
  /*------------------------------------------------------------------
    trace back through the "c", "f3" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This is fast, since only few structure elements are recalculated.
    ------------------------------------------------------------------*/
  struct sect {
    int  i;
    int  j;
    int ml;
  }
  sector[MAXSECTORS];   /* backtracking sectors */

  int   i, j, k, energy, new;
  int   no_close, type, type_2, tt;
  int   bonus;
  int   s=0;
  char *structure;

  /* length = strlen(string); */
  sector[++s].i = start;
  sector[s].j = MIN2(length, start+maxdist+1);
  sector[s].ml = (backtrack_type=='M') ? 1 : ((backtrack_type=='C')?2:0);

  structure = (char *) space((MIN2(length-start, maxdist)+3)*sizeof(char));
  for (i=0; i<=MIN2(length-start, maxdist); i++) structure[i] = '-';

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
      /* i or i+1 is paired. Find pairing partner */
      for (k=i+TURN+1,traced=0; k<=j; k++) {
	int cc, en;
	jj = k+1;
	type = ptype[i+1][k-(i+1)];
	if((type)&&(dangles%2==1)) {
	  cc = c[i+1][k-(i+1)]+P->dangle5[type][S1[i]];
	  if (type>2) cc += P->TerminalAU;
	  if (fij == cc + f3[k+1])
	    traced=i+1;
	  if (k<length)
	    if (fij == f3[k+2] + cc + P->dangle3[type][S1[k+1]]) {
	      traced=i+1; jj=k+2;
	    }
	}
	type = ptype[i][k-i];
	if (type) {
	  cc = c[i][k-i];
	  if (type>2) cc += P->TerminalAU;
	  en = cc + f3[k+1];
	  if (dangles==2) {
	    if (i>1)      en += P->dangle5[type][S1[i-1]];
	    if (k<length) en += P->dangle3[type][S1[k+1]];
	  }
	  if (fij == en) traced=i;
	  if ((dangles%2==1) && (k<length))
	    if (fij == f3[k+2]+cc+P->dangle3[type][S1[k+1]]) {
	      traced=i; jj=k+2;
	    }
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
      if ((jj==j+2) || (dangles==2)) structure[j+1-start] = '.';
      goto repeat1;
    }
    else { /* trace back in fML array */
      int cij1=INF, ci1j=INF, ci1j1=INF;

      if (fML[i][j-1-i]+P->MLbase == fij) {  /* 3' end is unpaired */
	sector[++s].i = i;
	sector[s].j   = j-1;
	sector[s].ml  = ml;
	continue;
      }
      if (fML[i+1][j-(i+1)]+P->MLbase == fij) { /* 5' end is unpaired */
	sector[++s].i = i+1;
	sector[s].j   = j;
	sector[s].ml  = ml;
	continue;
      }

      tt  = ptype[i][j-i];
      cij = c[i][j-i] + P->MLintern[tt];
      if (dangles==2) {       /* double dangles */
	if (i>1)      cij += P->dangle5[tt][S1[i-1]];
	if (j<length) cij += P->dangle3[tt][S1[j+1]];
      }
      else if (dangles%2==1) {  /* normal dangles */
	tt = ptype[i+1][j-(i+1)];
	ci1j= c[i+1][j-(i+1)] + P->dangle5[tt][S1[i]] + P->MLintern[tt]+P->MLbase;
	tt = ptype[i][j-1-i];
	cij1= c[i][j-1-i] + P->dangle3[tt][S1[j]] + P->MLintern[tt]+P->MLbase;
	tt = ptype[i+1][j-1-(i+1)];
	ci1j1=c[i+1][j-1-(i+1)] + P->dangle5[tt][S1[i]] + P->dangle3[tt][S1[j]]
	  +  P->MLintern[tt] + 2*P->MLbase;
      }

      if ((fij==cij)||(fij==ci1j)||(fij==cij1)||(fij==ci1j1)) {
	/* found a pair */
	if (fij==ci1j) i++;
	else if (fij==cij1) j--;
	else if (fij==ci1j1) {i++; j--;}
	structure[i-start] = '('; structure[j-start] = ')';
	goto repeat1;
      }

      for (k = i+1+TURN; k <= j-2-TURN; k++)
	if (fij == (fML[i][k-i]+fML[k+1][j-(k+1)]))
	  break;

      if ((dangles==3)&&(k>j-2-TURN)) { /* must be coax stack */
	ml = 2;
	for (k = i+1+TURN; k <= j-2-TURN; k++) {
	  type = ptype[i][k-i];  type= rtype[type];
	  type_2 = ptype[k+1][j-(k+1)]; type_2= rtype[type_2];
	  if (type && type_2)
	    if (fij == c[i][k-i]+c[k+1][j-(k+1)]+P->stack[type][type_2]+
		       2*P->MLintern[1])
	      break;
	}
      }

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

    type = ptype[i][j-i];

    bonus = 0;

    if (noLonelyPairs)
      if (cij == c[i][j-i]) {
	/* (i.j) closes canonical structures, thus
	   (i+1.j-1) must be a pair                */
	type_2 = ptype[i+1][j-1-(i+1)]; type_2 = rtype[type_2];
	cij -= P->stack[type][type_2] + bonus;
	structure[i+1-start] = '('; structure[j-1-start] = ')';
	i++; j--;
	canonical=0;
	goto repeat1;
      }
    canonical = 1;


    no_close = (((type==3)||(type==4))&&no_closingGU&&(bonus==0));
    if (no_close) {
      if (cij == FORBIDDEN) continue;
    } else
      if (cij == HairpinE(j-i-1, type, S1[i+1], S1[j-1],string+i-1)+bonus)
	continue;

    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
      int minq;
      minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+TURN) minq = p+1+TURN;
      for (q = j-1; q >= minq; q--) {

	type_2 = ptype[p][q-p];
	if (type_2==0) continue;
	type_2 = rtype[type_2];
	if (no_closingGU)
	  if (no_close||(type_2==3)||(type_2==4))
	    if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

	/* energy = oldLoopEnergy(i, j, p, q, type, type_2); */
	energy = LoopEnergy(p-i-1, j-q-1, type, type_2,
			    S1[i+1], S1[j-1], S1[p-1], S1[q+1]);

	new = energy+c[p][q-p]+bonus;
	traced = (cij == new);
	if (traced) {
	  structure[p-start] = '(';	structure[q-start] = ')';
	  i = p, j = q;
	  goto repeat1;
	}
      }
    }

    /* end of repeat: --------------------------------------------------*/

    /* (i.j) must close a multi-loop */
    tt = rtype[type];
    mm = bonus+P->MLclosing+P->MLintern[tt];
    d5 = P->dangle5[tt][S1[j-1]];
    d3 = P->dangle3[tt][S1[i+1]];
    i1 = i+1; j1 = j-1;
    sector[s+1].ml  = sector[s+2].ml = 1;

    for (k = i+2+TURN; k < j-2-TURN; k++) {
      int en;
      en = fML[i+1][k-(i+1)]+fML[k+1][j-1-(k+1)]+mm;
      if (dangles==2) /* double dangles */
	en += d5+d3;
      if (cij == en)
	break;
      if (dangles%2==1) { /* normal dangles */
	if (cij == (fML[i+2][k-(i+2)]+fML[k+1][j-1-(k+1)]+mm+d3+P->MLbase)) {
	  i1 = i+2;
	  break;
	}
	if (cij == (fML[i+1][k-(i+1)]+fML[k+1][j-2-(k+1)]+mm+d5+P->MLbase)) {
	  j1 = j-2;
	  break;
	}
	if (cij == (fML[i+2][k-(i+2)]+fML[k+1][j-2-(k+1)]+mm+d3+d5+P->MLbase+P->MLbase)) {
	  i1 = i+2; j1 = j-2;
	  break;
	}
      }
      /* coaxial stacking of (i.j) with (i+1.k) or (k.j-1) */
      /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
      if (dangles==3) {
	type_2 = ptype[i+1][k-(i+1)]; type_2 = rtype[type_2];
	if (type_2) {
	  en = c[i+1][k-(i+1)]+P->stack[type][type_2]+fML[k+1][j-1-(k+1)];
	  if (cij == en+2*P->MLintern[1]+P->MLclosing) {
	    ml = 2;
	    sector[s+1].ml  = 2;
	    break;
	  }
	}
	type_2 = ptype[k+1][j-1-(k+1)]; type_2 = rtype[type_2];
	if (type_2) {
	  en = c[k+1][j-1-(k+1)]+P->stack[type][type_2]+fML[i+1][k-(i+1)];
	  if (cij == en+2*P->MLintern[1]+P->MLclosing) {
	    sector[s+2].ml = 2;
	    break;
	  }
	}
      }

    }
    if (k<=j-3-TURN) { /* found the decomposition */
      sector[++s].i = i1;
      sector[s].j   = k;
      sector[++s].i = k+1;
      sector[s].j   = j1;
    } else {
#if 0
      /* Y shaped ML loops fon't work yet */
      if (dangles==3) {
	/* (i,j) must close a Y shaped ML loop with coax stacking */
	if (cij ==  fML[i+1][j-2-(i+2)] + mm + d3 + d5 + P->MLbase + P->MLbase) {
	  i1 = i+2;
	  j1 = j-2;
	} else if (cij ==  fML[i+1][j-2-(i+1)] + mm + d5 + P->MLbase)
	  j1 = j-2;
	else if (cij ==  fML[i+2][j-1-(i+2)] + mm + d3 + P->MLbase)
	  i1 = i+2;
	else /* last chance */
	  if (cij != fML[i+1][j-1-(i+1)] + mm + P->MLbase)
	    fprintf(stderr,  "backtracking failed in repeat");
	/* if we arrive here we can express cij via fML[i1,j1]+dangles */
	sector[++s].i = i1;
	sector[s].j   = j1;
      }
      else
#endif
	nrerror("backtracking failed in repeat");
    }

  }

  for (i=strlen(structure)-1; i>0 && structure[i] == '-'; i--)
    structure[i] = '\0';
  for (;i>=0; i--)
    if (structure[i]=='-') structure[i]='.';

  return structure;
}

/*---------------------------------------------------------------------------*/

PRIVATE void encode_seq(char *sequence) {
  unsigned int i,l;

  l = strlen(sequence);
  S = (short *) space(sizeof(short)*(l+1));
  S1= (short *) space(sizeof(short)*(l+1));
  /* S1 exists only for the special X K and I bases and energy_set!=0 */
  S[0] = S1[0] = (short) l;

  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    S[i]= (short) encode_char(toupper(sequence[i-1]));
    S1[i] = alias[S[i]];   /* for mismatches of nostandard bases */
  }
}

/*---------------------------------------------------------------------------*/

PRIVATE void letter_structure(char *structure, int length)
{
  int n, k, x, y;

  for (n = 0; n <= length-1; structure[n++] = ' ') ;
  structure[length] = '\0';

  for (n = 0, k = 1; k <= base_pair[0].i; k++) {
    y = base_pair[k].j;
    x = base_pair[k].i;
    if (x-1 > 0 && y+1 <= length) {
      if (structure[x-2] != ' ' && structure[y] == structure[x-2]) {
	structure[x-1] = structure[x-2];
	structure[y-1] = structure[x-1];
	continue;
      }
    }
    if (structure[x] != ' ' && structure[y-2] == structure[x]) {
      structure[x-1] = structure[x];
      structure[y-1] = structure[x-1];
      continue;
    }
    n++;
    structure[x-1] = alpha[n-1];
    structure[y-1] = alpha[n-1];
  }
}

/*---------------------------------------------------------------------------*/
#if 0
PRIVATE void parenthesis_structure(char *structure, int length)
{
  int n, k;

  for (n = 0; n <= length-1; structure[n++] = '.') ;
  structure[length] = '\0';

  for (k = 1; k <= base_pair[0].i; k++) {
    structure[base_pair[k].i-1] = '(';
    structure[base_pair[k].j-1] = ')';
  }
}
#endif
/*---------------------------------------------------------------------------*/

PRIVATE void update_fold_params(void)
{
  P = scale_parameters();
  make_pair_matrix();
}

/*---------------------------------------------------------------------------*/

PRIVATE void make_ptypes(const short *S, int i, int maxdist, int n) {
  int j,k, type;

  for (k=TURN+1; k<maxdist; k++) {
    j = i+k;
    if (j>n) continue;
    type = pair[S[i]][S[j]];
    if (noLonelyPairs && type) {
      if (!ptype[i+1][j-1-i-1])
	if (j==n || i==1 || (!pair[S[i-1]][S[j+1]])) type=0;
    }
    ptype[i][j-i]=type;
  }
}
