/* Last changed Time-stamp: <2006-02-26 00:08:50 ivo> */
/*
		  minimum free energy folding
		  for a set of aligned sequences

		  c Ivo Hofacker

		  Vienna RNA package
*/

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

/*@unused@*/
static char rcsid[] UNUSED = "$Id: alifold.c,v 1.11 2006/02/28 14:56:47 ivo Exp $";

#define PAREN

#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */

PUBLIC float  alifold(char **strings, char *structure);

PRIVATE void   init_alifold(int length);
PUBLIC  void   free_alifold_arrays(void);
PUBLIC  void   update_alifold_params(void);

PUBLIC double cv_fact=1.;
PUBLIC double nc_fact=1.;

PRIVATE void  parenthesis_structure(char *structure, int length);
PRIVATE void  get_arrays(unsigned int size);
PRIVATE void  make_pscores(const short *const *S, const char *const *AS,
			   int n_seq, const char *structure);
PRIVATE short *encode_seq(const char *sequence);
PRIVATE int fill_arrays(const char **strings);
PRIVATE void backtrack(const char **strings, int s);

/*@unused@*/
extern  int LoopEnergy(int n1, int n2, int type, int type_2,
		       int si1, int sj1, int sp1, int sq1);
extern  int HairpinE(int size, int type, int si1, int sj1, const char *string);

#define MAXSECTORS      500     /* dimension for a backtrack array */
#define LOCALITY        0.      /* locality parameter for base-pairs */

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

PRIVATE const paramT *P;

PRIVATE int *indx; /* index for moving in the triangle matrices c[] and fMl[]*/

PRIVATE int   *c;       /* energy array, given that i-j pair */
PRIVATE int   *cc;      /* linear array for calculating canonical structures */
PRIVATE int   *cc1;     /*   "     "        */
PRIVATE int   *f5;      /* energy of 5' end */
PRIVATE int   *fML;     /* multi-loop auxiliary energy array */

PRIVATE int   *Fmi;     /* holds row i of fML (avoids jumps in memory) */
PRIVATE int   *DMLi;    /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
PRIVATE int   *DMLi1;   /*             MIN(fML[i+1,k]+fML[k+1,j])  */
PRIVATE int   *DMLi2;   /*             MIN(fML[i+2,k]+fML[k+1,j])  */
PRIVATE int   *pscore;  /* precomputed array of pair types */
PRIVATE int   init_length=-1;

/*--------------------------------------------------------------------------*/

PRIVATE void init_alifold(int length)
{
  unsigned int n;
  if (length<1) nrerror("initialize_fold: argument must be greater 0");
  if (init_length>0) free_alifold_arrays();
  get_arrays((unsigned) length);
  make_pair_matrix();
  init_length=length;

  for (n = 1; n <= (unsigned) length; n++)
    indx[n] = (n*(n-1)) >> 1;        /* n(n-1)/2 */

  update_fold_params();
}

/*--------------------------------------------------------------------------*/

PRIVATE void get_arrays(unsigned int size)
{
  indx =  (int *) space(sizeof(int)*(size+1));
  c     = (int *) space(sizeof(int)*((size*(size+1))/2+2));
  fML   = (int *) space(sizeof(int)*((size*(size+1))/2+2));

  pscore = (int *) space(sizeof(int)*((size*(size+1))/2+2));
  f5    = (int *) space(sizeof(int)*(size+2));
  cc    = (int *) space(sizeof(int)*(size+2));
  cc1   = (int *) space(sizeof(int)*(size+2));
  Fmi   = (int *) space(sizeof(int)*(size+1));
  DMLi  = (int *) space(sizeof(int)*(size+1));
  DMLi1  = (int *) space(sizeof(int)*(size+1));
  DMLi2  = (int *) space(sizeof(int)*(size+1));
  base_pair = (struct bond *) space(sizeof(struct bond)*(1+size/2));
}

/*--------------------------------------------------------------------------*/

void free_alifold_arrays(void)
{
  free(indx); free(c); free(fML); free(f5); free(cc); free(cc1);
  free(pscore);
  free(base_pair); free(Fmi);
  free(DMLi); free(DMLi1);free(DMLi2);
  init_length=0;
}

static short **S;
/*--------------------------------------------------------------------------*/
#define UNIT 100
#define MINPSCORE -2 * UNIT
float alifold(char **strings, char *structure)
{
  int  length, energy, s, n_seq;

  length = (int) strlen(strings[0]);
  if (length>init_length) init_alifold(length);
  if ((P==NULL)||(fabs(P->temperature - temperature)>1e-6)) {
  	update_fold_params();  P = scale_parameters();
  }
  for (s=0; strings[s]!=NULL; s++);
  n_seq = s;
  S = (short **) space(n_seq*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if (strlen(strings[s]) != length) nrerror("uneqal seqence lengths");
    S[s] = encode_seq(strings[s]);
  }
  make_pscores((const short **) S, (const char *const *) strings, n_seq, structure);

  energy = fill_arrays((const char **)strings);

  backtrack((const char **)strings, 0);

  parenthesis_structure(structure, length);

  for (s=0; s<n_seq; s++) free(S[s]);
  free(S);

  if (backtrack_type=='C')
    return (float) c[indx[length]+1]/(n_seq*100.);
  else if (backtrack_type=='M')
    return (float) fML[indx[length]+1]/(n_seq*100.);
  else
    return (float) f5[length]/(n_seq*100.);
}

PRIVATE int fill_arrays(const char **strings) {
  int   i, j, k, p, q, length, energy, new_c;
  int   decomp, MLenergy, new_fML;
  int   s, n_seq, *type, type_2, tt;

  for (n_seq=0; strings[n_seq]!=NULL; n_seq++);
  type = (int *) space(n_seq*sizeof(int));
  length = strlen(strings[0]);
  for (j=1; j<=length; j++) {
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
  }

  for (j = 1; j<=length; j++)
    for (i=(j>TURN?(j-TURN):1); i<j; i++) {
      c[indx[j]+i] = fML[indx[j]+i] = INF;

    }

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


	for (new_c=s=0; s<n_seq; s++)
	  new_c += HairpinE(j-i-1,type[s],S[s][i+1],S[s][j-1],strings[s]+i-1);

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
	      energy += LoopEnergy(p-i-1, j-q-1, type[s], type_2,
				   S[s][i+1], S[s][j-1],
				   S[s][p-1], S[s][q+1]);
	    }
	    new_c = MIN2(energy+c[indx[q]+p], new_c);
	    if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */

	  } /* end q-loop */
	} /* end p-loop */

	/* multi-loop decomposition ------------------------*/

	decomp = DMLi1[j-1];
	if (dangles) {
	  int d3=0, d5=0;
	  for (s=0; s<n_seq; s++) {
	    tt = rtype[type[s]];
	    d3 = P->dangle3[tt][S[s][i+1]];
	    d5 = P->dangle5[tt][S[s][j-1]];
	    decomp += d5 + d3;
	  }
	}

	MLenergy = decomp + n_seq*P->MLclosing;
	for (s=0; s<n_seq; s++)
	  MLenergy += P->MLintern[type[s]];

	new_c = MLenergy < new_c ? MLenergy : new_c;

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
      for (s=0; s<n_seq; s++) {
	energy += P->MLintern[type[s]];
	if (dangles) {  /* double dangles */
	  energy += (i==1) ? /* works also for circfold */
	    P->dangle5[type[s]][S[s][length]] : P->dangle5[type[s]][S[s][i-1]];
	  /* if (j<length) */ energy += P->dangle3[type[s]][S[s][j+1]];
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

    }

    {
      int *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=1; j<=length; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
    }
  }
  /* calculate energies of 5' and 3' fragments */

  f5[TURN+1]=0;
  for (j=TURN+2; j<=length; j++) {
    f5[j] = f5[j-1];
    if (c[indx[j]+1]<INF) {
      energy = c[indx[j]+1];
      for (s=0; s<n_seq; s++) {
	int type;
	type = pair[S[s][1]][S[s][j]]; if (type==0) type=7;
	if (type>2) energy += TerminalAU;
	if ((dangles)&&(j<length))  /* double dangles */
	  energy += P->dangle3[type][S[s][j+1]];
      }
      f5[j] = MIN2(f5[j], energy);
    }
    for (i=j-TURN-1; i>1; i--) {
      if (c[indx[j]+i]<INF) {
	energy = f5[i-1]+c[indx[j]+i];
	for (s=0; s<n_seq; s++) {
	  int type;
	  type = pair[S[s][i]][S[s][j]]; if (type==0) type=7;
	  if (type>2) energy += TerminalAU;
	  if (dangles) {
	    energy += P->dangle5[type][S[s][i-1]];
	    if (j<length) energy += P->dangle3[type][S[s][j+1]];
	  }
	}
	f5[j] = MIN2(f5[j], energy);
      }
    }
  }
  free(type);
  return(f5[length]);
}

struct sect {
  int  i;
  int  j;
  int ml;
}
static sector[MAXSECTORS]; /* stack of partial structures for backtracking */

#include "alicircfold.inc"

void backtrack(const char **strings, int s) {
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
      base_pair[++b].i = i;
      base_pair[b].j   = j;
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
	  cc = c[indx[j]+i];
	  for (ss=0; ss<n_seq; ss++) {
	    type[ss] = pair[S[ss][i]][S[ss][j]];
	    if (type[ss]==0) type[ss] = 7;
	    if (type[ss]>2) cc += TerminalAU;
	  }
	  en = cc + f5[i-1];
	  if (dangles) {
	    for (ss=0; ss<n_seq; ss++) {
	      if (i>1)      en += P->dangle5[type[ss]][S[ss][i-1]];
	      if (j<length) en += P->dangle3[type[ss]][S[ss][j+1]];
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
      base_pair[++b].i = i;
      base_pair[b].j   = j;
      cov_en += pscore[indx[j]+i];
      goto repeat1;
    }
    else { /* trace back in fML array */
      int cij1=INF, ci1j=INF, ci1j1=INF;
      if (fML[indx[j]+i+1]+n_seq*P->MLbase == fij) { /* 5' end is unpaired */
	sector[++s].i = i+1;
	sector[s].j   = j;
	sector[s].ml  = ml;
	continue;
      }

      cij = c[indx[j]+i];
      for (ss=0; ss<n_seq; ss++) {
	tt  = pair[S[ss][i]][S[ss][j]];
	if (tt==0) tt=7;
	cij += P->MLintern[tt];
	if (dangles) {       /* double dangles */
	  cij += (i==1) ?
	    P->dangle5[tt][S[ss][length]] : P->dangle5[tt][S[ss][i-1]];
	  /* if (j<length) */ cij += P->dangle3[tt][S[ss][j+1]];
	}
      }

      if ((fij==cij)||(fij==ci1j)||(fij==cij1)||(fij==ci1j1)) {
	/* found a pair */
	if (fij==ci1j) i++;
	else if (fij==cij1) j--;
	else if (fij==ci1j1) {i++; j--;}
	base_pair[++b].i = i;
	base_pair[b].j   = j;
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
      cc += HairpinE(j-i-1, type[ss], S[ss][i+1], S[ss][j-1], strings[ss]+i-1);
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
	  energy += LoopEnergy(p-i-1, j-q-1, type[ss], type_2,
			       S[ss][i+1], S[ss][j-1],
			       S[ss][p-1], S[ss][q+1]);
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

    mm = n_seq*P->MLclosing;
    for (ss=d3=d5=0; ss<n_seq; ss++) {
      tt = rtype[type[ss]];
      mm += P->MLintern[tt];
      d5 += P->dangle5[tt][S[ss][j-1]];
      d3 += P->dangle3[tt][S[ss][i+1]];
    }
    i1 = i+1; j1 = j-1;
    sector[s+1].ml  = sector[s+2].ml = 1;

    for (k = i+2+TURN; k < j-2-TURN; k++) {
      int en;
      en = fML[indx[k]+i+1]+fML[indx[j-1]+k+1]+mm;
      if (dangles) /* double dangles */
	en += d5+d3;
      if (cij == en)
	break;

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

  base_pair[0].i = b;    /* save the total number of base pairs */
  free(type);
}

/*---------------------------------------------------------------------------*/

PRIVATE short * encode_seq(const char *sequence) {
  unsigned int i,l;
  short *S;
  l = strlen(sequence);
  S = (short *) space(sizeof(short)*(l+2));
  S[0] = (short) l;

  /* make numerical encoding of sequence */
  for (i=1; i<=l; i++)
    S[i]= (short) encode_char(toupper(sequence[i-1]));

  /* for circular folding add first base at position n+1 */
  S[l+1] = S[1];

  return S;
}

/*---------------------------------------------------------------------------*/

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
/*---------------------------------------------------------------------------*/

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
  n=S[0][0];  /* length of seqs */
  for (i=1; i<n; i++) {
    for (j=i+1; (j<i+TURN+1) && (j<=n); j++)
      pscore[indx[j]+i] = NONE;
    for (j=i+TURN+1; j<=n; j++) {
      int pfreq[8]={0,0,0,0,0,0,0,0};
      for (s=0; s<n_seq; s++) {
	int type;
	if (S[s][i]==0 && S[s][j]==0) type = 7; /* gap-gap  */
	else {
	  if ((AS[s][i] == '~')||(AS[s][j] == '~')) type = 7;
	  else type = pair[S[s][i]][S[s][j]];
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
    free(stack);
  }
}
