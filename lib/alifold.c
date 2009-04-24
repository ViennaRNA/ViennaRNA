/* Last changed Time-stamp: <2009-04-24 11:48:38 ivo> */
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
#include "ribo.h"
/*@unused@*/
static char rcsid[] UNUSED = "$Id: alifold.c,v 1.18 2009/02/27 16:25:54 ivo Exp $";

#define PAREN

#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */

PUBLIC float  alifold(char **strings, char *structure);

PRIVATE void   init_alifold(int length);
PUBLIC  void   free_alifold_arrays(void);

PUBLIC double cv_fact=1.;
PUBLIC double nc_fact=1.;

PRIVATE void  parenthesis_structure(char *structure, int length);
PRIVATE void  get_arrays(unsigned int size);
PRIVATE void  make_pscores(const short *const *S, const char *const *AS,
			   int n_seq, const char *structure);
PRIVATE short *encode_seq(const char *sequence, short *s5, short *s3, char *ss, unsigned short *as);
PRIVATE int fill_arrays(const char **strings);
PRIVATE void backtrack(const char **strings, int s);
PRIVATE int ML_Energy(int i,  int is_extloop,int n_seq);
PRIVATE double stack_energy(int i, char **sequences,  int n_seq, float *energy);
PRIVATE void arrays_for_energyofstruct(int n_seq, char **sequences);
PRIVATE void free_arrays_for_energyofstruct(int n_seq);

PRIVATE double energy_of_alistruct_pt(char **sequences,short * ptable, int n_seq, float *CVenergy);

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
PUBLIC float **readribosum(char *name);

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
  if (base_pair) free(base_pair);
  base_pair = (struct bond *) space(sizeof(struct bond)*(1+size/2));
}

/*--------------------------------------------------------------------------*/

void free_alifold_arrays(void)
{
  free(indx); free(c); free(fML); free(f5); free(cc); free(cc1);
  free(pscore);
  free(base_pair); base_pair=NULL; free(Fmi);
  free(DMLi); free(DMLi1);free(DMLi2);
  init_length=0;
}

static short **S;
static short **S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
static short **S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
static char **Ss;
static unsigned short **a2s;
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
  S5 = (short **) space(n_seq*sizeof(short *));
  S3 = (short **) space(n_seq*sizeof(short *));
  a2s= (unsigned short **)space(n_seq*sizeof(unsigned short *));
  Ss = (char **)space(n_seq*sizeof(char *));
  for (s=0; s<n_seq; s++) {
    if (strlen(strings[s]) != length) nrerror("uneqal seqence lengths");
    S5[s] =(short *) space ((length+2)*sizeof(short));
    S3[s] =(short *) space ((length+2)*sizeof(short));
    a2s[s]=(unsigned short *)space ((length+2)*sizeof(unsigned short));
    Ss[s]=(char *)space((length+2)*sizeof(char));
    S[s] = encode_seq(strings[s], S5[s],S3[s],Ss[s],a2s[s]);
  }
  make_pscores((const short **) S, (const char *const *) strings, n_seq, structure);

  energy = fill_arrays((const char **)strings);

  backtrack((const char **)strings, 0);

  parenthesis_structure(structure, length);

  for (s=0; s<n_seq; s++) {
    free(S[s]);
    free(S5[s]);
    free(S3[s]);
    free(a2s[s]);
    free(Ss[s]);
  }
  free(S);free(S5);free(S3);free(a2s);free(Ss);

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


	for (new_c=s=0; s<n_seq; s++) {
	  if ((a2s[s][j-1]-a2s[s][i])<3) new_c+=600;
	  else  new_c += HairpinE(a2s[s][j-1]-a2s[s][i],type[s],S3[s][i],S5[s][j],Ss[s]+(a2s[s][i-1]));
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
	      energy += LoopEnergy(a2s[s][p-1]-a2s[s][i], a2s[s][j-1]-a2s[s][q], type[s], type_2,
				   S3[s][i], S5[s][j],
				   S5[s][p], S3[s][q]);
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
	    d3 = P->dangle3[tt][S3[s][i]];
	    d5 = P->dangle5[tt][S5[s][j]];
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
	    P->dangle5[type[s]][S5[s][1]] : P->dangle5[type[s]][S5[s][i]];
	  /* if (j<length) */ energy += P->dangle3[type[s]][S3[s][j]];
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
	  energy += P->dangle3[type][S3[s][j]];
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
	    energy += P->dangle5[type][S5[s][i]];
	    if (j<length) energy += P->dangle3[type][S3[s][j]];
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
	      if (i>1)      en += P->dangle5[type[ss]][S5[ss][i]];
	      if (j<length) en += P->dangle3[type[ss]][S3[ss][j]];
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
	    P->dangle5[tt][S5[ss][1]] : P->dangle5[tt][S5[ss][i]];
	  /* if (j<length) */ cij += P->dangle3[tt][S3[ss][j]];
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
    for (ss=0; ss<n_seq; ss++) {
	if ((a2s[ss][j-1]-a2s[ss][i])<3) cc+=600;
	else cc += HairpinE(a2s[ss][j-1]-a2s[ss][i], type[ss], S3[ss][i], S5[ss][j], Ss[ss]+a2s[ss][i-1]);
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
	  energy += LoopEnergy(a2s[ss][p-1]-a2s[ss][i],a2s[ss][j-1]-a2s[ss][q],
			       type[ss], type_2,
			       S3[ss][i], S5[ss][j],
			       S5[ss][p], S3[ss][q]);

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
      d5 += P->dangle5[tt][S5[ss][j]];
      d3 += P->dangle3[tt][S3[ss][i]];
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

 PRIVATE short * encode_seq(const char *sequence, short *s5, short *s3, char *ss, unsigned short *as) {
   unsigned int i,l;
  short *S;
  unsigned short p;
  l = strlen(sequence);
  S = (short *) space(sizeof(short)*(l+2));
  S[0] = (short) l;

   s5[0]=s5[1]=0;
  /* make numerical encoding of sequence */

   for (i=1; i<=l; i++) {
     short ctemp;
     ctemp=(short) encode_char(toupper(sequence[i-1]));
     S[i]= ctemp ;
   }

   if (oldAliEn) {
     /*use alignment sequences in all energy evaluations*/
     ss[0]=sequence[0];
     for (i=1; i<l; i++) {
       s5[i]=S[i-1];
       s3[i]=S[i+1];
       ss[i]= sequence[i];
       as[i]=i;
     }
     ss[l] = sequence[l];
     as[l]=l;
     s5[l]=S[l-1];
     s3[l]=0;
     S[l+1] = S[1];
     s5[1]=0;
     if (1) {
       s5[1]=S[l];
       s3[l]=S[1];
       ss[l+1]=S[1];
     }
     return S;
   }
   else {
     if (1) {
       for (i=l; i>0; i--) {
	 char c5;
	 c5=sequence[i-1];
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
     } else  s5[1]=s3[l]=0;

     for (i=1,p=0; i<=l; i++) {
       char c5;
       c5=sequence[i-1];
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
       c3=sequence[i-1];
       if ((c3=='-')||(c3=='_')||(c3=='~')||(c3=='.'))
	 s3[i-1]=s3[i];
       else
	 s3[i-1]=S[i];
     }
   }

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

  n=S[0][0];  /* length of seqs */
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
/*Get the mean pairwise identity in steps from ?to?(ident)*/
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

/* how to chose a ribosum matrix:
ribosum matrices exist for starlike clusters of
X=45 55 60 65 70 75 80 85 90 95 100
they are further seperated by only regarding sequences with a minimal pairwise idensity of Y=25-95, step 5, not all are present.
now the question is, which matrix to use when.
the suggestion of the dr. will:
with a mpi of Z
X~Z and Y > Z ??
if we say the minimum of the pis is M,
then we may be able to use:
X~Z and Y ~ M?
I'd say we do a default matrix (e.g. 85/60) but better is try out (all 170??)
and then we use the best in average.
actually, it would be preferrable to make a very big testset and simply check it out.
(create a function to derive the best matrix)
furthermore:
default, function or user defined.

ntscd:
fijklmn
pijpklpmn

*/


PUBLIC float **readribosum(char *name) {

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







/*------------------ENERGY OF STRUCT----------------------------------*/
PRIVATE short  *pair_table;
/*write another "call function" including all the allocation, seq2num etc
  if needed*/
static int *type;

PUBLIC float energy_of_alistruct(char **sequences, const char *structure, int n_seq, float *CVenergy)
{
  double energy;
  int new=0;
  short **tempS;
  short **tempS5;     /*S5[s][i] holds next base 5' of i in sequence s*/
  short **tempS3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
  char **tempSs;
  unsigned short **tempa2s;
  int *temptype, *tempindx, *temppscore;

  if (P==NULL)  P = scale_parameters();
  /*save old memory*/
  tempS=S; tempS3=S3; tempS5=S5; tempSs=Ss; tempa2s=a2s; temptype=type;
  tempindx=indx; temppscore=pscore;

  arrays_for_energyofstruct(n_seq, sequences);
  make_pscores((const short *const*)S, (const char *const *)sequences, n_seq, NULL);
  make_pair_matrix();
  new=1;

  pair_table = make_pair_table(structure);

  energy = energy_of_alistruct_pt(sequences,pair_table, n_seq, CVenergy);

  free(pair_table);
  energy /= (100.*n_seq);
  *CVenergy /= (100.*n_seq);
  free_arrays_for_energyofstruct(n_seq);
  S=tempS;S3=tempS3; S5=tempS5; Ss=tempSs; a2s=tempa2s; type=temptype;
  indx=tempindx; pscore=temppscore;
  return energy ;
}
/*------------------------------------------------------------------*/

PRIVATE double energy_of_alistruct_pt(char **sequences,short * ptable,
				      int n_seq, float *CVenergy) {
  int   i, length;
  double energy=0;

  *CVenergy = 0;

  pair_table = ptable;
  length = S[0][0];
  energy =  backtrack_type=='M' ? (float) ML_Energy(0, 0,n_seq) :
				  (float) ML_Energy(0, 1,n_seq);

  for (i=1; i<=length; i++) {
    if (pair_table[i]==0) continue;
    energy += stack_energy(i, sequences,  n_seq, CVenergy);
    i=pair_table[i];
  }
  return energy;
}


/*---------------------------------------------------------------------------*/
PRIVATE double stack_energy(int i, char **sequences,  int n_seq, float *CVenergy)
{
  /* calculate energy of substructure enclosed by (i,j) */
  int ee= 0;
  int j, p, q, s;
  int numberofcomponents;
  double energy = 0;
  j=pair_table[i];
  for (s=0; s<n_seq; s++) {
    type[s] = pair[S[s][i]][S[s][j]];
    if (type[s]==0) {
    type[s]=7;
    }
  }
  p=i; q=j;
  while (p<q) { /* process all stacks and interior loops */
    int type_2;
    while (pair_table[++p]==0);
    while (pair_table[--q]==0);
    if ((pair_table[q]!=(short)p)||(p>q)) break;
    ee=0;
    for (s=0; s<n_seq; s++) {
      type_2 = pair[S[s][q]][S[s][p]];
      if (type_2==0) {
	type_2=7;
      }
      ee += LoopEnergy(a2s[s][p-1]-a2s[s][i], a2s[s][j-1]-a2s[s][q], type[s], type_2,
		       S3[s][i], S5[s][j], S5[s][p], S3[s][q]);
    }
    energy += ee;
    *CVenergy += pscore[indx[j]+i];
    i=p; j=q;
    for (s=0; s<n_seq; s++) {
      type[s] = pair[S[s][i]][S[s][j]];
      if (type[s]==0) {
	type[s]=7;	}
    }
  }  /* end while */

  /* p,q don't pair must have found hairpin or multiloop */

  if (p>q) {
    ee=0;/* hair pin */
    for (s=0; s< n_seq; s++) {
      if ((a2s[s][j-1]-a2s[s][i])<3) ee+=600;
      else ee += HairpinE(a2s[s][j-1]-a2s[s][i], type[s],
			  S3[s][i], S5[s][j], Ss[s]+(a2s[s][i-1]));
    }
    energy += ee;
    *CVenergy += pscore[indx[j]+i];
    return energy;
  }
  numberofcomponents=0;
  /* (i,j) is exterior pair of multiloop */
  *CVenergy += pscore[indx[j]+i];
  while (p<j) {
    /* add up the contributions of the substructures of the ML */
      energy += stack_energy(p, sequences,  n_seq, CVenergy);
    p = pair_table[p];
    numberofcomponents++;
    /* search for next base pair in multiloop */
    while (pair_table[++p]==0);
  }
  energy += ML_Energy(i,0,n_seq);

  return energy;
}

PRIVATE int ML_Energy(int i,  int is_extloop, int n_seq) {
  /* i is the 5'-base of the closing pair (or 0 for exterior loop)
     loop is scored as ML if extloop==0 else as exterior loop

     since each helix can coaxially stack with at most one of its
     neighbors we need an auxiliarry variable  cx_energy
     which contains the best energy given that the last two pairs stack.
     energy  holds the best energy given the previous two pairs do not
     stack (i.e. the two current helices may stack)
     We don't allow the last helix to stack with the first, thus we have to
     walk around the Loop twice with two starting points and take the minimum
 */

  int energy, best_energy=INF;
  int i1, j, p, q, u, x, stype, count,s;
  int mlintern[NBPAIRS+1], mlclosing, mlbase;

  if (is_extloop) {
    mlintern[0]=mlintern[1]=mlintern[2]=0;
    for (x = 3; x <= NBPAIRS; x++) /*3 od 2?*/
      mlintern[x] = P->TerminalAU; /* 0 or TerminalAU */
    mlclosing = mlbase = 0;
  } else {
    for (x = 0; x <= NBPAIRS; x++) mlintern[x] = P->MLintern[x];
    mlclosing =P->MLclosing*n_seq; mlbase = P->MLbase*n_seq;
  }
  if ( i==0 ) {
    j = pair_table[0]+1;
    stype = 0;  /* no pair */
  }
  else {
    j = pair_table[i];
    for (s=0; s<n_seq; s++) {
      type[s] = pair[S[s][j]][S[s][i]];
      if (type[s]==0) type[s]=7;
    }
  }
  i1=i; p = i+1; u=0;
  energy = 0;
  do { /* walk around the multi-loop */
    int tt = INF;

    /* hop over unpaired positions */
    while (p <= pair_table[0] && pair_table[p]==0) p++;

    /* memorize number of unpaired positions? no, we just approximate here */
    u += p-i1-1;
    /* get position of pairing partner */
    for (s=0; s< n_seq; s++) {
      if ( p == pair_table[0]+1 ) {
	q = tt = 0; /* virtual root pair */
	/*	  break; wiaso soi i da noamoi duach?*/
      }
      else {
	q  = pair_table[p];
	/* get type of base pair P->q */
	tt = pair[S[s][p]][S[s][q]]; if (tt==0) tt=7;
      }

      energy += mlintern[tt];

      if (dangles) {
	int dang5=0, dang3=0, dang;
	if ((a2s[s][p]>1) && (tt!=0)) dang5= P->dangle5[tt][S5[s][p]];
	if ((i1>0) && a2s[s][i1]<a2s[s][S[0][0]])
	  dang3 = P->dangle3[type[s]][S3[s][i1]];

	switch (p-i1-1) {
	case 0: /* adjacent helices */
	  if (dangles==2)
	    energy += dang3+dang5;
	  break;
	case 1: /* 1 unpaired base between helices */
	  dang = (dangles==2)?(dang3+dang5):MIN2(dang3, dang5);
	  energy += dang;
	  break;
	default: /* many unpaired base between helices */
	  energy += dang5 +dang3;
	}
	type[s] = tt;
      }
    }
    i1 = q; p=q+1;

  } while (q!=i);

  energy += mlclosing;
  energy += mlbase*u;

  return energy;
}

/*---------------------------------------------------------------------------*/
PRIVATE void arrays_for_energyofstruct(int n_seq, char **sequences){
  int s,n;
  n = (int) strlen(sequences[0]);
  S = (short **) space(sizeof(short *)*(n_seq+1));
  S5 = (short **) space(n_seq*sizeof(short *));
  S3 = (short **) space(n_seq*sizeof(short *));
  a2s= (unsigned short **)space(n_seq*sizeof(unsigned short *));
  Ss = (char **)space(n_seq*sizeof(char *));
  type = (int *) space(n_seq*sizeof(int));
  pscore = (int *) space(sizeof(int)*((n+1)*(n+2)/2));
   indx = (int *) space(sizeof(int)*(n+1));
   for (s = 1; s <= (unsigned) n; s++){
     indx[s] = (s*(s-1)) >> 1;
   }
   for (s=0; s<n_seq; s++) {
    if (strlen(sequences[s]) != n) nrerror("uneqal seqence lengths");
    S5[s] =(short *) space ((n+2)*sizeof(short));
    S3[s] =(short *) space ((n+2)*sizeof(short));
    a2s[s]=(unsigned short *)space ((n+2)*sizeof(unsigned short));
    Ss[s]=(char *)space((n+2)*sizeof(char));
    S[s] = encode_seq(sequences[s], S5[s], S3[s], Ss[s], a2s[s]);
  }


}
PRIVATE void free_arrays_for_energyofstruct(int n_seq)
{
  int i;
  for (i=0; i<n_seq; i++) {
    free(S[i]);
    free(S5[i]);
    free(S3[i]);
    free(Ss[i]);
    free(a2s[i]);
  }
  free(S5);
  free(S3);
  free(Ss);
  free(a2s);
  free(S);
  free(type);
  free(pscore);
  free(indx);
}
