/* Last changed Time-stamp: <2004-02-09 17:09:47 ivo> */
/*                
		  minimum free energy
		  RNA secondary structure prediction

		  c Ivo Hofacker, Chrisoph Flamm
		  original implementation by
		  Walter Fontana
		  
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
static char rcsid[] UNUSED = "$Id: cofold.c,v 1.4 2004/02/09 16:50:48 ivo Exp $";

#define PAREN

#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */
#define free_arrays free_co_arrays
#define initialize_fold initialize_cofold
#define update_fold_params update_cofold_params
PUBLIC float  cofold(const char *string, char *structure);
PRIVATE int   energy_of_struct_pt(const char *string, short *ptable,
				  short *s, short *s1);
PUBLIC void   free_arrays(void);
PUBLIC void   initialize_fold(int length);
PUBLIC void   update_fold_params(void);

extern int    logML;    /* if nonzero use logarithmic ML energy in
			     energy_of_struct */
extern int    uniq_ML;  /* do ML decomposition uniquely (for subopt) */
/*@unused@*/
PRIVATE void  letter_structure(char *structure, int length) UNUSED;
PRIVATE void  parenthesis_structure(char *structure, int length);
PRIVATE void  get_arrays(unsigned int size);
/* PRIVATE void  scale_parameters(void); */
PRIVATE int   stack_energy(int i, const char *string);
PRIVATE int   ML_Energy(int i, int is_extloop);
PRIVATE void  make_ptypes(const short *S, const char *structure);
PRIVATE void  encode_seq(const char *sequence);
PRIVATE void backtrack(const char *sequence);
PRIVATE int fill_arrays(const char *sequence);
/*@unused@*/
inline PRIVATE  int   oldLoopEnergy(int i, int j, int p, int q, int type, int type_2);
inline PRIVATE int  LoopEnergy(int n1, int n2, int type, int type_2,
			 int si1, int sj1, int sp1, int sq1);
inline PRIVATE int  HairpinE(int size, int type, int si1, int sj1, const char *string);
PRIVATE void free_end(int *array, int i, int start);

#define MAXSECTORS      500     /* dimension for a backtrack array */
#define LOCALITY        0.      /* locality parameter for base-pairs */

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define SAME_STRAND(I,J) (((I)>=cut_point)||((J)<cut_point))

PRIVATE paramT *P = NULL;

PRIVATE int *indx; /* index for moving in the triangle matrices c[] and fMl[]*/

PRIVATE int   *c;       /* energy array, given that i-j pair */
PRIVATE int   *cc;      /* linear array for calculating canonical structures */
PRIVATE int   *cc1;     /*   "     "        */
PRIVATE int   *f5;      /* energy of 5' end */
PRIVATE int   *fc;      /* energy from i to cutpoint (and vice versa if i>cut) */
PRIVATE int   *fML;     /* multi-loop auxiliary energy array */
PRIVATE int   *fM1;     /* second ML array, only for subopt */
PRIVATE int   *Fmi;     /* holds row i of fML (avoids jumps in memory) */
PRIVATE int   *DMLi;    /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
PRIVATE int   *DMLi1;   /*             MIN(fML[i+1,k]+fML[k+1,j])  */
PRIVATE int   *DMLi2;   /*             MIN(fML[i+2,k]+fML[k+1,j])  */
PRIVATE char  *ptype;   /* precomputed array of pair types */ 
PRIVATE short  *S, *S1;
PRIVATE int   init_length=-1;

PRIVATE char  alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
/* needed by cofold/eval */
PRIVATE int cut_in_loop(int i);
#undef TURN
#define TURN 0
PRIVATE int min_hairpin = TURN;
extern  int cut_point;   /* set to first pos of second seq for cofolding */
extern  int eos_debug;   /* verbose info from energy_of_struct */

/*--------------------------------------------------------------------------*/

PUBLIC void initialize_fold(int length)
{
  unsigned int n;
  if (length<1) nrerror("initialize_fold: argument must be greater 0");
  if (init_length>0) free_arrays();
  get_arrays((unsigned) length);
  init_length=length;
  
  for (n = 1; n <= (unsigned) length; n++)
    indx[n] = (n*(n-1)) >> 1;        /* n(n-1)/2 */

  update_fold_params();
}
    
/*--------------------------------------------------------------------------*/

PRIVATE void get_arrays(unsigned int size)
{
  indx = (int *) space(sizeof(int)*(size+1));
  c     = (int *) space(sizeof(int)*((size*(size+1))/2+2));
  fML   = (int *) space(sizeof(int)*((size*(size+1))/2+2));
  if (uniq_ML)
    fM1    = (int *) space(sizeof(int)*((size*(size+1))/2+2));

  ptype = (char *) space(sizeof(char)*((size*(size+1))/2+2));
  f5    = (int *) space(sizeof(int)*(size+2));  /* must have size [0..n+1] */
  fc    = (int *) space(sizeof(int)*(size+2));
  cc    = (int *) space(sizeof(int)*(size+2));
  cc1   = (int *) space(sizeof(int)*(size+2));
  Fmi   = (int *) space(sizeof(int)*(size+1));
  DMLi  = (int *) space(sizeof(int)*(size+1));
  DMLi1 = (int *) space(sizeof(int)*(size+1));
  DMLi2 = (int *) space(sizeof(int)*(size+1));
  base_pair = (struct bond *) space(sizeof(struct bond)*(1+size/2));
}

/*--------------------------------------------------------------------------*/

PUBLIC void free_arrays(void)
{
  free(indx); free(c); free(fML); free(f5); free(cc); free(cc1);
  free(fc);
  free(ptype);
  if (uniq_ML) free(fM1);

  free(base_pair); free(Fmi);
  free(DMLi); free(DMLi1);free(DMLi2);
  init_length=0;
}

/*--------------------------------------------------------------------------*/

void export_cofold_arrays(int **f5_p, int **c_p, int **fML_p, int **fM1_p, 
			int **indx_p, char **ptype_p) {
  /* make the DP arrays available to routines such as subopt() */ 
  *f5_p = f5; *c_p = c;
  *fML_p = fML; *fM1_p = fM1;
  *indx_p = indx; *ptype_p = ptype;
}

/*--------------------------------------------------------------------------*/

PRIVATE   int   *BP; /* contains the structure constrainsts: BP[i]
			-1: | = base must be paired
			-2: < = base must be paired with j<i
			-3: > = base must be paired with j>i
			-4: x = base must not pair
			positive int: base is paired with int      */

float cofold(const char *string, char *structure) {
  int i, length, energy, bonus=0, bonus_cnt=0;
  
  length = (int) strlen(string);
  if (length>init_length) initialize_fold(length);
  if (fabs(P->temperature - temperature)>1e-6) update_fold_params();
  
  encode_seq(string);
  
  BP = (int *)space(sizeof(int)*(length+2));
  make_ptypes(S, structure);
  
  energy = fill_arrays(string);
  fprintf(stderr, "energy=%d\n", energy);
  
  backtrack(string);

#ifdef PAREN
  parenthesis_structure(structure, length);
#else
  letter_structure(structure, length);
#endif

  /* check constraints */ 
  for(i=1;i<=length;i++) {
    if((BP[i]<0)&&(BP[i]>-4)) {
      bonus_cnt++;
      if((BP[i]==-3)&&(structure[i-1]==')')) bonus++;
      if((BP[i]==-2)&&(structure[i-1]=='(')) bonus++;
      if((BP[i]==-1)&&(structure[i-1]!='.')) bonus++;
    }
    
    if(BP[i]>i) {
      int l;
      bonus_cnt++;
      for(l=1; l<=base_pair[0].i; l++)
	if((i==base_pair[l].i)&&(BP[i]==base_pair[l].j)) bonus++;
    }
  }
  
  if (bonus_cnt>bonus) fprintf(stderr,"\ncould not enforce all constraints\n");
  bonus*=BONUS;

  free(S); free(S1); free(BP);

  energy += bonus;      /*remove bonus energies from result */

  if (backtrack_type=='C')
    return (float) c[indx[length]+1]/100.;
  else if (backtrack_type=='M')
    return (float) fML[indx[length]+1]/100.;
  else
    return (float) energy/100.;
}

PRIVATE int fill_arrays(const char *string) {
  /* fill "c", "fML" and "f5" arrays and return  optimal energy */

  int   i, j, k, length, energy;
  int   decomp, new_fML, max_separation;
  int   no_close, type, type_2, tt;
  int   bonus=0;

  length = (int) strlen(string);

  max_separation = (int) ((1.-LOCALITY)*(double)(length-2)); /* not in use */

  for (j=1; j<=length; j++) {
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF; fc[j]=0;
  }
   
  for (j = 1; j<=length; j++)
    for (i=1; i<=j; i++) {
      c[indx[j]+i] = fML[indx[j]+i] = INF;
      if (uniq_ML) fM1[indx[j]+i] = INF;
    }       
  
  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */
      
    for (j = i+TURN+1; j <= length; j++) {
      int p, q, ij;
      ij = indx[j]+i;
      bonus = 0;
      type = ptype[ij];

      /* enforcing structure constraints */
      if ((BP[i]==j)||(BP[i]==-1)||(BP[i]==-2)) bonus -= BONUS;
      if ((BP[j]==-1)||(BP[j]==-3)) bonus -= BONUS;
      if ((BP[i]==-4)||(BP[j]==-4)) type=0;
	 	 
      no_close = (((type==3)||(type==4))&&no_closingGU&&(bonus==0));
		 
      if (j-i-1 > max_separation) type = 0;  /* forces locality degree */

      if (type) {   /* we have a pair */
	int new_c=0, stackEnergy=INF;
	/* hairpin ----------------------------------------------*/

	if (SAME_STRAND(i,j)) {
          if (no_close) new_c = FORBIDDEN;
          else
            new_c = HairpinE(j-i-1, type, S1[i+1], S1[j-1], string+i-1);
        }
        else {
	  if (dangles) {
	    if (SAME_STRAND(i,i+1)) new_c += P->dangle3[rtype[type]][S1[i+1]];
	    if (SAME_STRAND(j-1,j)) new_c += P->dangle5[rtype[type]][S1[j-1]];
	  }
	  if (type>2) new_c += P->TerminalAU;
	}

	/*--------------------------------------------------------
	  check for elementary structures involving more than one
	  closing pair.
	  --------------------------------------------------------*/
	   
	for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++) {
	  int minq = j-i+p-MAXLOOP-2;
	  if (minq<p+1+TURN) minq = p+1+TURN;
	  for (q = minq; q < j; q++) {
	    type_2 = ptype[indx[q]+p];
	       
	    if (type_2==0) continue;
	    type_2 = rtype[type_2];

	    if (no_closingGU) 
	      if (no_close||(type_2==3)||(type_2==4))
		if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */
	    
#if 1
	    if (SAME_STRAND(i,p) && SAME_STRAND(q,j))
	      energy = LoopEnergy(p-i-1, j-q-1, type, type_2,
				  S1[i+1], S1[j-1], S1[p-1], S1[q+1]);
	    else {
	      energy = 0;
	      if (dangles) {
		int di, dj, dp, dq;
		di = (SAME_STRAND(i,i+1)) ?
		  P->dangle3[rtype[type]][S1[i+1]] : 0;
		dj = (SAME_STRAND(j-1,j)) ?
		  P->dangle5[rtype[type]][S1[j-1]] : 0;
		dp = (SAME_STRAND(p-1,p)) ?
		  P->dangle5[rtype[type_2]][S1[p-1]] : 0;
		dq = (SAME_STRAND(q,q+1)) ?
		  P->dangle3[rtype[type_2]][S1[q+1]] : 0;
		if (dangles == 2) energy = di + dj + dp + dq;
		else { /* no double dangles */
		  if (i+1 < p-1) energy = di + dp;
		  if (i+1 == p-1) energy = MIN2(di, dp);
		  if (q+1 < j-1) energy += dj + dq;
		  if (q+1 == j-1) energy += MIN2(dj, dq);
		}
	      }
	      if (type>2) energy += P->TerminalAU;
	      if (type_2>2) energy += P->TerminalAU;
	    }
#else
	    /* duplicated code is faster than function call */
	       
#endif
	    new_c = MIN2(energy+c[indx[q]+p], new_c);
	    if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */
	       
	  } /* end q-loop */
	} /* end p-loop */
	   
	/* multi-loop decomposition ------------------------*/

	
	if (!no_close) {
	  int MLenergy;
	  int d3=0, d5=0; 
	  decomp = DMLi1[j-1];
	  if (dangles) {
	    tt = rtype[type];
	    d3 = SAME_STRAND(i,i+1) ? P->dangle3[tt][S1[i+1]] : 0;
	    d5 = SAME_STRAND(j-1,j) ? P->dangle5[tt][S1[j-1]] : 0;
	    if (dangles==2) /* double dangles */
	      decomp += d5 + d3;
	    else {          /* normal dangles */
	      decomp = MIN2(DMLi2[j-1]+d3+P->MLbase, decomp);
	      decomp = MIN2(DMLi1[j-2]+d5+P->MLbase, decomp);
	      decomp = MIN2(DMLi2[j-2]+d5+d3+2*P->MLbase, decomp);
	    }
	  }
	  
	  MLenergy = P->MLclosing+P->MLintern[type]+decomp;
	  
	  new_c = MLenergy < new_c ? MLenergy : new_c;
	  
	  if (!SAME_STRAND(i,j)) { /* cut is somewhere in the multiloop*/
	    decomp = fc[i+1]+fc[j-1];
	    if (dangles==2)
	      decomp+=d5+d3;
	    else {
	      decomp = MIN2(fc[i+2]+fc[j-1]+d3, decomp);
	      decomp = MIN2(fc[i+1]+fc[j-2]+d5, decomp);
	      decomp = MIN2(fc[i+2]+fc[j-2]+d3+d5, decomp);
	    }
	    if (type>2) decomp+=P->TerminalAU;
	    new_c = MIN2(new_c, decomp);  
	  }
	} /* end >> if (!no_close) << */
	
	/* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */

	if (dangles==3) {
	  decomp = INF;
	  for (k = i+2+TURN; k < j-2-TURN; k++) {
	    type_2 = ptype[indx[k]+i+1]; type_2 = rtype[type_2]; 
	    if (type_2)
	      decomp = MIN2(decomp, c[indx[k]+i+1]+P->stack[type][type_2]+
			    fML[indx[j-1]+k+1]);
	    type_2 = ptype[indx[j-1]+k+1]; type_2 = rtype[type_2]; 
	    if (type_2)
	      decomp = MIN2(decomp, c[indx[j-1]+k+1]+P->stack[type][type_2]+
			    fML[indx[k]+i+1]);
	  }
	  /* no TermAU penalty if coax stack */
	  decomp += 2*P->MLintern[1] + P->MLclosing;
	  new_c = MIN2(new_c, decomp);
	}
	
	new_c = MIN2(new_c, cc1[j-1]+stackEnergy);
	cc[j] = new_c + bonus;
	if (noLonelyPairs)
	  c[ij] = cc1[j-1]+stackEnergy+bonus;
	else
	  c[ij] = cc[j];
	   
      } /* end >> if (pair) << */
	 
      else c[ij] = INF;


      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/

      new_fML = fML[ij+1]+P->MLbase;
      new_fML = MIN2(fML[indx[j-1]+i]+P->MLbase, new_fML);
      energy = c[ij]+P->MLintern[type];
      if (dangles==2) {  /* double dangles */
	if (i>1)      energy += P->dangle5[type][S1[i-1]];
	if (j<length) energy += P->dangle3[type][S1[j+1]];
      }
      new_fML = MIN2(energy, new_fML);
      if (uniq_ML)
	fM1[ij] = MIN2(fM1[indx[j-1]+i] + P->MLbase, energy);

      if (dangles%2==1) {  /* normal dangles */
	tt = ptype[ij+1]; /* i+1,j */
	new_fML = MIN2(c[ij+1]+P->dangle5[tt][S1[i]]
		       +P->MLintern[tt]+P->MLbase,new_fML);
	tt = ptype[indx[j-1]+i]; 
	new_fML = MIN2(c[indx[j-1]+i]+P->dangle3[tt][S1[j]]
		       +P->MLintern[tt]+P->MLbase, new_fML);
	tt = ptype[indx[j-1]+i+1];
	new_fML = MIN2(c[indx[j-1]+i+1]+P->dangle5[tt][S1[i]]+
		       P->dangle3[tt][S1[j]]+P->MLintern[tt]+2*P->MLbase, new_fML);
      }
      
      /* modular decomposition -------------------------------*/
    
      for (decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++)
	decomp = MIN2(decomp, Fmi[k]+fML[indx[j]+k+1]);
      
      DMLi[j] = decomp;               /* store for use in ML decompositon */
      new_fML = MIN2(new_fML,decomp);
      
      /* coaxial stacking */
      if (dangles==3) { 
	/* additional ML decomposition as two coaxially stacked helices */
	for (decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++) {
	  type = ptype[indx[k]+i]; type = rtype[type];
	  type_2 = ptype[indx[j]+k+1]; type_2 = rtype[type_2];
	  if (type && type_2)
	    decomp = MIN2(decomp, 
			  c[indx[k]+i]+c[indx[j]+k+1]+P->stack[type][type_2]);
	}

	decomp += 2*P->MLintern[1];  	/* no TermAU penalty if coax stack */
#if 0
	/* This is needed for Y shaped ML loops with coax stacking of
	   interior pairts, but backtracking will fail if activated */
	DMLi[j] = MIN2(DMLi[j], decomp);
	DMLi[j] = MIN2(DMLi[j], DMLi[j-1]+P->MLbase);
	DMLi[j] = MIN2(DMLi[j], DMLi1[j]+P->MLbase);
	new_fML = MIN2(new_fML, DMLi[j]);
#endif
	new_fML = MIN2(new_fML, decomp);
      }
      
      fML[ij] = Fmi[j] = new_fML;     /* substring energy */
  
    }
    
    if (i==cut_point)
      for (j=i; j<=length; j++) 
	free_end(fc, j, cut_point);
    if (i<cut_point)
      free_end(fc,i,cut_point-1);
    

    {
      int *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=1; j<=length; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
    }
  }
  
  /* calculate energies of 5' and 3' fragments */
  
  for (i=1; i<=length; i++) 
    free_end(f5, i, 1);
    
  return f5[length];
}

PRIVATE void backtrack(const char *string) {
   
  /*------------------------------------------------------------------
    trace back through the "c", "fc", "f5" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This is fast, since only few structure elements are recalculated.
    ------------------------------------------------------------------*/
  struct sect {
    int  i;
    int  j;
    int ml;
  } 
  sector[MAXSECTORS];   /* backtracking sectors */

  int   i, j, k, length, energy, new;
  int   no_close, type, type_2, tt;
  int   bonus;
  int   s=0, b=0;

  length = strlen(string);
  sector[++s].i = 1;
  sector[s].j   = length;
  sector[s].ml  = (backtrack_type=='M') ? 1 : ((backtrack_type=='C')?2:0);
   
  while (s>0) {
    int ml, fij, fi, cij, traced, i1, j1, d3, d5, mm, p, q, jj=0;
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

    
    if (ml==0) {fij = f5[j]; fi = f5[j-1];}
    else if (ml==1) {fij = fML[indx[j]+i]; fi = fML[indx[j-1]+i]+P->MLbase;}
    else /* 3 or 4 */ {
      fij = fc[j];
      fi = (ml==3) ? INF : fc[j-1];
    }
    if (fij == fi) {  /* 3' end is unpaired */
      sector[++s].i = i;
      sector[s].j   = j-1;
      sector[s].ml  = ml;
      continue;
    }
     
    if (ml==0 || ml==4) { /* backtrack in f5 or fc[i=cut,j>cut] */
      int *ff;
      ff = (ml==4) ? fc : f5;
      /* j or j-1 is paired. Find pairing partner */
      for (k=j-TURN-1,traced=0; k>=i; k--) {
	int cc, en;
	jj = k-1; 
	type = ptype[indx[j-1]+k];
	if((type)&&(dangles%2==1)&&SAME_STRAND(j-1,j)) {
	  cc = c[indx[j-1]+k]+P->dangle3[type][S1[j]];
	  if (type>2) cc += P->TerminalAU;
	  if (fij == cc + ff[k-1]) 
	    traced=j-1;
	  if (k>i) 
	    if (fij == ff[k-2] + cc + P->dangle5[type][S1[k-1]]) {
	      traced=j-1; jj=k-2;
	    }
	}
	type = ptype[indx[j]+k];
	if (type) {
	  cc = c[indx[j]+k];
	  if (type>2) cc += P->TerminalAU; 
	  en = cc + ff[k-1];
	  if (dangles==2) {
	    if (k>1 &&SAME_STRAND(k-1,k)) en += P->dangle5[type][S1[k-1]];
	    if (j<length &&SAME_STRAND(j,j+1)) en += P->dangle3[type][S1[j+1]];
	  }
	  if (fij == en) traced=j;
	  if ((dangles%2==1) && (k>1) && SAME_STRAND(k-1,k))
	    if (fij == ff[k-2]+cc+P->dangle5[type][S1[k-1]]) {
	      traced=j; jj=k-2;
	    }
	}
	if (traced) break;
      }

      if (!traced) nrerror("backtrack failed in f5 (or fc)");
      sector[++s].i = i;
      sector[s].j   = jj;
      sector[s].ml  = ml;
       
      i=k; j=traced;
      base_pair[++b].i = i;
      base_pair[b].j   = j;
      goto repeat1;
    }
    else if (ml==3) { /* backtrack in fc[i<cut,j=cut-1] */
      if (fc[i] == fc[i+1]) { /* 5' end is unpaired */
	sector[++s].i = i+1;
	sector[s].j   = j;
	sector[s].ml  = ml;
	continue;
      }
      /* i or i+1 is paired. Find pairing partner */
      for (k=i+TURN+1, traced=0; k<=j; k++) {
	jj=k+1;
	type = ptype[indx[k]+i];
	if (type) {
	  int d5, d3, en;
	  d5 = (i>1 && SAME_STRAND(i-1,i)) ? P->dangle5[type][S1[i-1]] : 0;
	  d3 = (SAME_STRAND(k,k+1)) ? P->dangle3[type][S1[k+1]] : 0;
	  en = fc[k+1]+c[indx[k]+i];
	  if (type>2) en+=P->TerminalAU;
	  
	  if (dangles==2) en += d5+d3;
	  if (fc[i]==en) traced=i;
	  else if (dangles%2==1) {
	    int tau;
	    tau = (type>2) ? P->TerminalAU : 0;
	    if (fc[i]==fc[k+2]+c[indx[k]+i]+d3+tau)
	      traced=i; jj=k+2;
	  }
	}
	if (traced) break;
	
	if (dangles%2==1) {
	  int tau;  
	  type = ptype[indx[k]+i+1];
	  tau = (type>2) ? P->TerminalAU : 0;
	  if (type) {
	    d5 = (SAME_STRAND(i, i+1)) ? P->dangle5[type][S1[i]] : 0;
	    if (fc[i] == fc[k+1]+c[indx[k]+i+1]+d5+tau)
	      traced=i+1;
	    if (k<j) {
	      d3 = (SAME_STRAND(k, k+1)) ? P->dangle3[type][S1[k+1]] : 0;
	      if (fc[i] == fc[k+2]+c[indx[k]+i+1]+d5+d3+tau)
		traced=i+1; jj=k+2;
	    }
	  }
	}
	if (traced) break;
      }
      
      if (!traced) nrerror("backtrack failed in fc[] 5' of cut");

      sector[++s].i = jj;
      sector[s].j   = j;
      sector[s].ml  = ml;

      j=k; i=traced;
      base_pair[++b].i = i;
      base_pair[b].j   = j;
      goto repeat1;
    }

    else { /* true multi-loop backtrack in fML */
      int cij1=INF, ci1j=INF, ci1j1=INF;
      if (fML[indx[j]+i+1]+P->MLbase == fij) { /* 5' end is unpaired */
	sector[++s].i = i+1;
	sector[s].j   = j;
	sector[s].ml  = ml;
	continue;
      } 

      tt  = ptype[indx[j]+i];
      cij = c[indx[j]+i] + P->MLintern[tt];
      if (dangles==2) {       /* double dangles */
	if (i>1)      cij += P->dangle5[tt][S1[i-1]];
	if (j<length) cij += P->dangle3[tt][S1[j+1]];
      }
      else if (dangles%2==1) {  /* normal dangles */
	tt = ptype[indx[j]+i+1];
	ci1j = c[indx[j]+i+1]+P->dangle5[tt][S1[i]]+P->MLintern[tt]+P->MLbase;
	tt = ptype[indx[j-1]+i];
	cij1 = c[indx[j-1]+i]+P->dangle3[tt][S1[j]]+P->MLintern[tt]+P->MLbase;
	tt = ptype[indx[j-1]+i+1];
	ci1j1 = c[indx[j-1]+i+1]+P->dangle5[tt][S1[i]]+P->dangle3[tt][S1[j]]
	  + P->MLintern[tt]+2*P->MLbase;
      }
       
      if ((fij==cij)||(fij==ci1j)||(fij==cij1)||(fij==ci1j1)) {
	/* found a pair */
	if (fij==ci1j) i++;
	else if (fij==cij1) j--;
	else if (fij==ci1j1) {i++; j--;}
	base_pair[++b].i = i;
	base_pair[b].j   = j;
	goto repeat1;
      } 

      /* find next coponent of multiloop */
      for (k = i+1+TURN; k <= j-2-TURN; k++) 
	if (fij == (fML[indx[k]+i]+fML[indx[j]+k+1])) 
	  break;
      
      if ((dangles==3)&&(k>j-2-TURN)) { /* must be coax stack */
	ml = 2;
	for (k = i+1+TURN; k <= j-2-TURN; k++) {
	  type = ptype[indx[k]+i];  type= rtype[type]; 
	  type_2 = ptype[indx[j]+k+1]; type_2= rtype[type_2];
	  if (type && type_2)
	    if (fij == c[indx[k]+i]+c[indx[j]+k+1]+P->stack[type][type_2]+
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
    if (canonical)  cij = c[indx[j]+i];
     
    type = ptype[indx[j]+i];
     
    bonus = 0;
     
    if ((BP[i]==j)||(BP[i]==-1)||(BP[i]==-2)) bonus -= BONUS;
    if ((BP[j]==-1)||(BP[j]==-3)) bonus -= BONUS;
     
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
    if (SAME_STRAND(i,j)) {
      if (no_close) {
	if (cij == FORBIDDEN)
	  continue;
      } else
	if (cij == HairpinE(j-i-1, type, S1[i+1], S1[j-1],string+i-1)+bonus)
	  continue;
    }
    else {
      int ee = 0;
      if (dangles) {
	if (SAME_STRAND(i,i+1)) ee  = P->dangle3[rtype[type]][S1[i+1]];
	if (SAME_STRAND(j-1,j)) ee += P->dangle5[rtype[type]][S1[j-1]];
      }
      if (type>2) ee += P->TerminalAU;
      if (cij == ee) continue;
    }
     
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
	 
	/* energy = oldLoopEnergy(i, j, p, q, type, type_2); */
	if (SAME_STRAND(i,p) && SAME_STRAND(q,j))
	  energy = LoopEnergy(p-i-1, j-q-1, type, type_2,
			      S1[i+1], S1[j-1], S1[p-1], S1[q+1]);
	else {
	  energy = 0;
	  if (dangles) {
	    int di, dj, dp, dq;
	    di = (SAME_STRAND(i,i+1)) ? P->dangle3[rtype[type]][S1[i+1]] : 0;
	    dj = (SAME_STRAND(j-1,j)) ? P->dangle5[rtype[type]][S1[j-1]] : 0;
	    dp = (SAME_STRAND(p-1,p)) ? P->dangle5[rtype[type_2]][S1[p-1]] : 0;
	    dq = (SAME_STRAND(q,q+1)) ? P->dangle3[rtype[type_2]][S1[q+1]] : 0;
	    if (dangles == 2) energy = di + dj + dp + dq;
	    else { /* no double dangles */
	      if (i+1 < p-1) energy = di + dp;
	      if (i+1 == p-1) energy = MIN2(di, dp);
	      if (q+1 < j-1) energy += dj + dq;
	      if (q+1 == j-1) energy += MIN2(dj, dq);
	    }
	  }
	  if (type>2) energy += P->TerminalAU;
	  if (type_2>2) energy += P->TerminalAU;
	}
 
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
     
    /* (i.j) must close a fake or true multi-loop */
    tt = rtype[type];
    mm = bonus+P->MLclosing+P->MLintern[tt];
    d5 = (SAME_STRAND(j-1,j) && dangles) ? P->dangle5[tt][S1[j-1]] : 0;
    d3 = (SAME_STRAND(i,i+1) && dangles) ? P->dangle3[tt][S1[i+1]] : 0;
    i1 = i+1; j1 = j-1;
    sector[s+1].ml  = sector[s+2].ml = 1;

    /* fake multi-loop */
    if (!SAME_STRAND(i,j)) {
      int ii=0, jj=0, decomp=0;
      decomp = fc[i+1]+fc[j-1];
      if (type>2) decomp+=P->TerminalAU;
      if (dangles==2) /* double dangles */
	decomp+=d5+d3;
      if (decomp==cij) ii=i+1, jj=j-1;
      else {
	if (dangles%2==1) { /* normal dangles */
	  int tau;
	  tau = (type>2) ? P->TerminalAU : 0;
	  if (cij == fc[i+2]+fc[j-1]+d3+tau) ii=i+2, jj=j-1;
	  if (cij == fc[i+1]+fc[j-2]+d5+tau) ii=i+1, jj=j-2;
	  if (cij == fc[i+2]+fc[j-2]+d3+d5+tau) ii=i+2, jj=j-2;
	}
      }
      if (ii) {
	sector[++s].i = ii;
	sector[s].j   = cut_point-1;
	sector[s].ml  = 3;
	sector[++s].i = cut_point;
	sector[s].j   = jj;
	sector[s].ml  = 4;
	continue;
      }
    }

    /* true multi-loop */
    for (k = i+2+TURN; k < j-2-TURN; k++) {
      int en;
      en = fML[indx[k]+i+1]+fML[indx[j-1]+k+1]+mm;
      if (dangles==2) /* double dangles */
	en += d5+d3;
      if (cij == en) 
	break;
      if (dangles%2==1) { /* normal dangles */
	if (cij == (fML[indx[k]+i+2]+fML[indx[j-1]+k+1]+mm+d3+P->MLbase)) {
	  i1 = i+2;
	  break;
	}
	if (cij == (fML[indx[k]+i+1]+fML[indx[j-2]+k+1]+mm+d5+P->MLbase)) {
	  j1 = j-2;
	  break;
	}
	if (cij == (fML[indx[k]+i+2]+fML[indx[j-2]+k+1]+mm+d3+d5+P->MLbase+P->MLbase)) {
	  i1 = i+2; j1 = j-2;
	  break;
	}
      }
      /* coaxial stacking of (i.j) with (i+1.k) or (k.j-1) */
      /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
      if (dangles==3) {
	type_2 = ptype[indx[k]+i+1]; type_2 = rtype[type_2];
	if (type_2) {
	  en = c[indx[k]+i+1]+P->stack[type][type_2]+fML[indx[j-1]+k+1];
	  if (cij == en+2*P->MLintern[1]+P->MLclosing) {
	    ml = 2;
	    sector[s+1].ml  = 2;
	    break;
	  }
	}
	type_2 = ptype[indx[j-1]+k+1]; type_2 = rtype[type_2];
	if (type_2) {
	  en = c[indx[j-1]+k+1]+P->stack[type][type_2]+fML[indx[k]+i+1];
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
      /* Y shaped ML loops don't work yet */
      if (dangles==3) {
	/* (i,j) must close a Y shaped ML loop with coax stacking */
	if (cij == fML[indx[j-2]+i+2] + mm + d3 + d5 + P->MLbase + P->MLbase) {
	  i1 = i+2;
	  j1 = j-2;
	} else if (cij ==  fML[indx[j-2]+i+1] + mm + d5 + P->MLbase) 
	  j1 = j-2;
	else if (cij ==  fML[indx[j-1]+i+2] + mm + d3 + P->MLbase) 
	  i1 = i+2;
	else /* last chance */
	  if (cij != fML[indx[j-1]+i+1] + mm + P->MLbase)
	    fprintf(stderr,  "backtracking failed in repeat");
	/* if we arrive here we can express cij via fML[i1,j1]+dangles */
	sector[++s].i = i1;
	sector[s].j   = j1;
      }
      else
#endif
	nrerror("backtracking failed in repeat");
    } 
     
  } /* end >> while (s>0) << */

  base_pair[0].i = b;    /* save the total number of base pairs */
}

#if 0
PRIVATE void free_end(int *array, int i, int start) {
  int inc, type, energy, length, j, ii, jj;
  inc = (i>start)? 1:-1;
  length = S[0];

  if (i==start) array[i]=0;
  else array[i] = array[i-inc];
  if (i>start) {
    ii = start; jj=i;
  } else {
     ii = i; jj = start;
  }
  type = ptype[indx[jj]+ii];
  if (type) {
    energy = c[indx[jj]+ii];
    if (type>2) energy += P->TerminalAU;
    if (dangles==2) { /* double dangles */
      if ((jj<length)&&SAME_STRAND(jj,jj+1))  
	energy += P->dangle3[type][S1[jj+1]];
      if ((ii>1)&&SAME_STRAND(ii-1,ii))  
	energy += P->dangle5[type][S1[ii-1]];
    }
    array[i] = MIN2(array[i], energy);
  }
  if (dangles%2==1) {
    type=ptype[indx[jj-1]+ii];
    if ((type)&&(dangles%2==1)&&SAME_STRAND(jj-1,jj)) {
      energy = c[indx[jj-1]+ii]+P->dangle3[type][S1[jj]];
      if (type>2) energy += P->TerminalAU;
      array[i] = MIN2(array[i], energy);
    }
  }
  for (j=start+inc; inc*(i-j)>TURN; j+=inc) {
    int d3, d5;
    if (i>j) { ii = j; jj = i;} /* inc>0 */
    else     { jj = i; ii = j;} /* inc<0 */
    type = ptype[indx[jj]+ii];
    d5 = (ii>1 && SAME_STRAND(ii-1,ii))? P->dangle5[type][S1[ii-1]]:0;
    d3 = (jj<length && SAME_STRAND(jj,jj+1))?P->dangle3[type][S1[jj+1]]:0;
    if (type) {
      energy = array[j-inc]+c[indx[jj]+ii];
      if (type>2) energy += P->TerminalAU;
      if (dangles==2) energy += d3 + d5;
      array[i] = MIN2(array[i], energy);
      if (dangles%2==1) { 
	if (inc>0)  energy = array[j-2] + c[indx[jj]+ii-1] + d5;
	else        energy = array[j+2] + c[indx[jj+1]+ii] + d3;
	if (type>2) energy += P->TerminalAU;
	array[i] = MIN2(array[i], energy);
      }
    }
    if (dangles%2==1) {
      int iii, jjj;
      if (inc>0) {iii = ii; jjj=jj-1;} 
      else       {iii = ii+1; jjj = jj;}
      type = ptype[indx[jjj]+iii];
      if (!type) continue;
      d5 = (iii>1 && SAME_STRAND(iii-1,iii)) ? 
	P->dangle5[type][S1[iii-1]] : 0;
      d3 = (jjj<length && SAME_STRAND(jjj,jjj+1)) ? 
	P->dangle3[type][S1[jjj+1]] : 0;
      energy = c[indx[jjj]+iii] + ((inc>0)?d3:d5);      
      if (type>2) energy += P->TerminalAU;
      array[i] = MIN2(array[i], array[j-inc]+energy);
      energy += (inc>0)?d5:d3;
      array[i] = MIN2(array[i], array[j-2*inc]+energy);
    }
  }
}
#endif



PRIVATE void free_end(int *array, int i, int start) {
  int inc, type, energy, length, j, ii, jj;
  inc = (i>start)? 1:-1;
  length = S[0];

  if (i==start) array[i]=0;
  else array[i] = array[i-inc];
  if (i>start) {
    ii = start; jj=i;
  } else {
    ii = i; jj = start;
  }
  for (j=start; inc*(i-j)>TURN; j+=inc) {
    int d3, d5;
    if (i>j) { ii = j; jj = i;} /* inc>0 */
    else     { ii = i; jj = j;} /* inc<0 */
    type = ptype[indx[jj]+ii];
    if (type) {
      d5 = (ii>1 && SAME_STRAND(ii-1,ii))? P->dangle5[type][S1[ii-1]]:0;
      d3 = (jj<length && SAME_STRAND(jj,jj+1))?P->dangle3[type][S1[jj+1]]:0;
      
      energy = c[indx[jj]+ii];
      if (type>2) energy += P->TerminalAU;
      if (dangles==2) energy += d3 + d5;
      array[i] = MIN2(array[i], array[j-inc]+energy);
      
      if (dangles%2==1) { 
	if (inc>0) {
	  if (j>1)  energy += array[j-2] + d5;
	} else
	  if (j<length)  energy += d3 + array[j+2];
	array[i] = MIN2(array[i], energy);
      }
    }
    if (dangles%2==1) { 
      int iii, jjj;
      if (inc>0) {iii = ii; jjj=jj-1;} 
      else       {iii = ii+1; jjj = jj;}
      type = ptype[indx[jjj]+iii];
      if (!type) continue;
      d5 = (iii>1 && SAME_STRAND(iii-1,iii)) ? 
	P->dangle5[type][S1[iii-1]] : 0;
      d3 = (jjj<length && SAME_STRAND(jjj,jjj+1)) ? 
	P->dangle3[type][S1[jjj+1]] : 0;
      energy = c[indx[jjj]+iii] + ((inc>0)?d3:d5);      
      if (type>2) energy += P->TerminalAU;
      array[i] = MIN2(array[i], array[j-inc]+energy);
      if (j!=start) { /* dangles on both sides */
	energy += (inc>0)?d5:d3;
	array[i] = MIN2(array[i], array[j-2*inc]+energy);
      }
    }
  }
}

/*---------------------------------------------------------------------------*/

inline int HairpinE(int size, int type, int si1, int sj1, const char *string) {
  int energy;
  energy = (size <= 30) ? P->hairpin[size] :
    P->hairpin[30]+(int)(P->lxc*log((size)/30.));
  if (tetra_loop)
    if (size == 4) { /* check for tetraloop bonus */
      char tl[7]={0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl))) 
	energy += P->TETRA_ENERGY[(ts-P->Tetraloops)/7];
    }
  if (size == 3) {
    char tl[6]={0,0,0,0,0,0}, *ts;
    strncpy(tl, string, 5);
    if ((ts=strstr(P->Triloops, tl))) 
      energy += P->Triloop_E[(ts-P->Triloops)/6];
    
    if (type>2)  /* neither CG nor GC */
      energy += P->TerminalAU; /* penalty for closing AU GU pair */
  }
  else  /* no mismatches for tri-loops */
    energy += P->mismatchH[type][si1][sj1];
  
  return energy;
}

/*---------------------------------------------------------------------------*/

inline PRIVATE int oldLoopEnergy(int i, int j, int p, int q, int type, int type_2) {
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int n1, n2, m, energy;
  n1 = p-i-1;
  n2 = j-q-1;

  if (n1>n2) { m=n1; n1=n2; n2=m; } /* so that n2>=n1 */

  if (n2 == 0)
    energy = P->stack[type][type_2];   /* stack */
  
  else if (n1==0) {                  /* bulge */
    energy = (n2<=MAXLOOP)?P->bulge[n2]:
      (P->bulge[30]+(int)(P->lxc*log(n2/30.)));

#if STACK_BULGE1
    if (n2==1) energy+=P->stack[type][type_2];
#endif
  } else {                           /* interior loop */

    if ((n1+n2==2)&&(james_rule))
      /* special case for loop size 2 */
      energy = P->int11[type][type_2][S1[i+1]][S1[j-1]];
    else {
      energy = (n1+n2<=MAXLOOP)?(P->internal_loop[n1+n2]):
	(P->internal_loop[30]+(int)(P->lxc*log((n1+n2)/30.)));
      
#if NEW_NINIO
      energy += MIN2(MAX_NINIO, (n2-n1)*P->F_ninio[2]);
#else
      m       = MIN2(4, n1);
      energy += MIN2(MAX_NINIO,((n2-n1)*P->F_ninio[m]));
#endif
      energy += P->mismatchI[type][S1[i+1]][S1[j-1]]+
	P->mismatchI[type_2][S1[q+1]][S1[p-1]];
    }
  }
  return energy;
}

/*--------------------------------------------------------------------------*/

inline int LoopEnergy(int n1, int n2, int type, int type_2,
		      int si1, int sj1, int sp1, int sq1) {
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int nl, ns, energy;
  
  if (n1>n2) { nl=n1; ns=n2;}
  else {nl=n2; ns=n1;}

  if (nl == 0)
    return P->stack[type][type_2];    /* stack */

  if (ns==0) {                       /* bulge */
    energy = (nl<=MAXLOOP)?P->bulge[nl]:
      (P->bulge[30]+(int)(P->lxc*log(nl/30.)));
    if (nl==1) energy += P->stack[type][type_2];
    else {
      if (type>2) energy += P->TerminalAU;
      if (type_2>2) energy += P->TerminalAU;
    }
    return energy;
  }
  else {                             /* interior loop */
    if (ns==1) { 
      if (nl==1)                     /* 1x1 loop */
	return P->int11[type][type_2][si1][sj1];
      if (nl==2) {                   /* 2x1 loop */
	if (n1==1)
	  energy = P->int21[type][type_2][si1][sq1][sj1];
	else
	  energy = P->int21[type_2][type][sq1][si1][sp1];
	return energy;
      }
    }
    else if (n1==2 && n2==2)         /* 2x2 loop */
      return P->int22[type][type_2][si1][sp1][sq1][sj1];
    { /* generic interior loop (no else here!)*/
      energy = (n1+n2<=MAXLOOP)?(P->internal_loop[n1+n2]):
	(P->internal_loop[30]+(int)(P->lxc*log((n1+n2)/30.)));
      
      energy += MIN2(MAX_NINIO, (nl-ns)*P->F_ninio[2]);
      
      energy += P->mismatchI[type][si1][sj1]+
	P->mismatchI[type_2][sq1][sp1];
    }
  }
  return energy;
}


/*---------------------------------------------------------------------------*/

PRIVATE void encode_seq(const char *sequence) {
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
	   
PUBLIC void update_fold_params(void)
{
  P = scale_parameters();
  make_pair_matrix();
  if (init_length < 0) init_length=0;
}

/*---------------------------------------------------------------------------*/
PRIVATE short  *pair_table;

PRIVATE float energy_of_struct(const char *string, char *structure)
{
  int   energy;
  short *ss, *ss1;

  if ((init_length<0)||(P==NULL)) update_fold_params();
  if (fabs(P->temperature - temperature)>1e-6) update_fold_params();

  if (strlen(structure)!=strlen(string))
    nrerror("energy_of_struct: string and structure have unequal length");

  /* save the S and S1 pointers in case they were already in use */
  ss = S; ss1 = S1;
  encode_seq(string);
   
  pair_table = make_pair_table(structure);

  energy = energy_of_struct_pt(string, pair_table, S, S1);
   
  free(pair_table);
  free(S); free(S1);
  S=ss; S1=ss1;
  return  (float) energy/100.;
}

PRIVATE int energy_of_struct_pt(const char *string, short * ptable, short *s, short *s1) {
  /* auxiliary function for kinfold,
     for most purposes call energy_of_struct instead */
  
  int   i, length, energy;
  
  pair_table = ptable;
  S = s;
  S1 = s1;
  
  length = S[0];
  energy =  backtrack_type=='M' ? ML_Energy(0, 0) : ML_Energy(0, 1);
  if (eos_debug>0)
    printf("External loop                           : %5d\n", energy);
  for (i=1; i<=length; i++) {
    if (pair_table[i]==0) continue;
    energy += stack_energy(i, string);          
    i=pair_table[i];
  }

  return energy;
}

/*---------------------------------------------------------------------------*/
PRIVATE int stack_energy(int i, const char *string)  
{
  /* calculate energy of substructure enclosed by (i,j) */
  int ee, energy = 0;
  int j, p, q, type;

  j=pair_table[i];
  type = pair[S[i]][S[j]];
  if (type==0) {
    type=7;
    if (eos_debug>=0)
    fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", i, j,
            string[i-1],string[j-1]);
  }
   
  p=i; q=j;
  while (p<q) { /* process all stacks and interior loops */
    int type_2;
    while (pair_table[++p]==0);
    while (pair_table[--q]==0);
    if ((pair_table[q]!=(short)p)||(p>q)) break;
    type_2 = pair[S[q]][S[p]];
    if (type_2==0) {
      type_2=7;
      if (eos_debug>=0)
      fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", p, q,
              string[p-1],string[q-1]);
    }
    /* energy += LoopEnergy(i, j, p, q, type, type_2); */
    if ( SAME_STRAND(i,p) && SAME_STRAND(q,j) )
      ee = LoopEnergy(p-i-1, j-q-1, type, type_2,
		      S1[i+1], S1[j-1], S1[p-1], S1[q+1]);
    else 
      ee = ML_Energy(cut_in_loop(i), 1);
    if (eos_debug>0)
      printf("Interior loop (%3d,%3d) %c%c; (%3d,%3d) %c%c: %5d\n",
	     i,j,string[i-1],string[j-1],p,q,string[p-1],string[q-1], ee);
    energy += ee;    
    i=p; j=q; type = rtype[type_2];
  } /* end while */
  
  /* p,q don't pair must have found hairpin or multiloop */
  
  if (p>q) {                       /* hair pin */
    if (SAME_STRAND(i,j))
      ee = HairpinE(j-i-1, type, S1[i+1], S1[j-1], string+i-1);
    else
      ee = ML_Energy(cut_in_loop(i), 1);
    energy += ee;
    if (eos_debug>0)
      printf("Hairpin  loop (%3d,%3d) %c%c              : %5d\n",
	     i, j, string[i-1],string[j-1], ee);
    
    return energy;
  }
   
  /* (i,j) is exterior pair of multiloop */
  while (p<j) {
    /* add up the contributions of the substructures of the ML */
    energy += stack_energy(p, string);
    p = pair_table[p]; 
    /* search for next base pair in multiloop */
    while (pair_table[++p]==0);
  }
  {
    int ii;
    ii = cut_in_loop(i);
    ee = (ii==0) ? ML_Energy(i,0) : ML_Energy(ii, 1);
  }
  energy += ee;
  if (eos_debug>0)
    printf("Multi    loop (%3d,%3d) %c%c              : %5d\n",
	   i,j,string[i-1],string[j-1],ee);
  
  return energy;
}

/*---------------------------------------------------------------------------*/

PRIVATE int ML_Energy(int i, int is_extloop) {
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

  int energy, cx_energy, best_energy=INF;
  int i1, j, p, q, u, x, type, count;
  int mlintern[NBPAIRS+1], mlclosing, mlbase;

  if (is_extloop) {
    for (x = 0; x <= NBPAIRS; x++) 
      mlintern[x] = P->MLintern[x]-P->MLintern[1]; /* 0 or TerminalAU */
    mlclosing = mlbase = 0; 
  } else {
    for (x = 0; x <= NBPAIRS; x++) mlintern[x] = P->MLintern[x];
    mlclosing = P->MLclosing; mlbase = P->MLbase; 
  }
 
  for (count=0; count<2; count++) { /* do it twice */
    int ld5 = 0; /* 5' dangle energy on prev pair (type) */
    if ( i==0 ) {
      j = pair_table[0]+1;
      type = 0;  /* no pair */
    }
    else {
      j = pair_table[i];
      type = pair[S[j]][S[i]]; if (type==0) type=7;
      
      if (dangles==3) { /* prime the ld5 variable */
	if (SAME_STRAND(j-1,j)) {
	  ld5 = P->dangle5[type][S1[j-1]];
	  if ((p=pair_table[j-2]) && SAME_STRAND(j-2, j-1)) 
	      if (P->dangle3[pair[S[p]][S[j-2]]][S1[j-1]]<ld5) ld5 = 0;
	}
      }
    }
    i1=i; p = i+1; u=0;
    energy = 0; cx_energy=INF;
    do { /* walk around the multi-loop */
      int tt, new_cx = INF;
      
      /* hope over unpaired positions */
      while (p <= pair_table[0] && pair_table[p]==0) p++;
            
      /* memorize number of unpaired positions */
      u += p-i1-1; 
      /* get position of pairing partner */
      if ( p == pair_table[0]+1 ) 
	q = tt = 0; /* virtual root pair */
      else {
      q  = pair_table[p];
	/* get type of base pair P->q */
      tt = pair[S[p]][S[q]]; if (tt==0) tt=7;
      }

      energy += mlintern[tt];
      cx_energy += mlintern[tt]; 
      
      if (dangles) {
        int dang5=0, dang3=0, dang;
	if ((SAME_STRAND(p-1,p))&&(p>1))
	  dang5=P->dangle5[tt][S1[p-1]];      /* 5'dangle of pq pair */
	if ((SAME_STRAND(i1,i1+1))&&(i1<S[0]))
	  dang3 = P->dangle3[type][S1[i1+1]];  /* 3' dangle of previous pair */

        switch (p-i1-1) {
        case 0: /* adjacent helices */
          if (dangles==2) 
	    energy += dang3+dang5;
          else if (dangles==3 && i1!=0) {
	    if (SAME_STRAND(i1,p)) {
	      new_cx = energy + P->stack[rtype[type]][rtype[tt]];
	      /* subtract 5'dangle and TerminalAU penalty */
	      new_cx += -ld5 - mlintern[tt]-mlintern[type]+2*mlintern[1];
	    }
            ld5=0;
            energy = MIN2(energy, cx_energy);
          }
          break;
        case 1: /* 1 unpaired base between helices */
          dang = (dangles==2)?(dang3+dang5):MIN2(dang3, dang5);
          if (dangles==3) {
	    energy = energy +dang; ld5 = dang - dang3;
	    /* may be problem here: Suppose
	       cx_energy>energy, cx_energy+dang5<energy
	       and the following helices are also stacked (i.e.
	       we'll subtract the dang5 again */
	    if (cx_energy+dang5 < energy) {
	      energy = cx_energy+dang5;
	      ld5 = dang5;
	    }
            new_cx = INF;  /* no coax stacking with mismatch for now */
          } else
            energy += dang;
          break;
        default: /* many unpaired base between helices */
          energy += dang5 +dang3;
          if (dangles==3) {
            energy = MIN2(energy, cx_energy + dang5);
            new_cx = INF;  /* no coax stacking possible */
            ld5 = dang5;
          }
        }
        type = tt;
      }
      if (dangles==3) cx_energy = new_cx;
      i1 = q; p=q+1;
    } while (q!=i);
    best_energy = MIN2(energy, best_energy); /* don't use cx_energy here */
    /* fprintf(stderr, "%6.2d\t", energy); */
    if (dangles!=3 || is_extloop) break;  /* may break cofold with co-ax */ 
    /* skip a helix and start again */
    while (pair_table[p]==0) p++;
    if (i == pair_table[p]) break;
    i = pair_table[p]; 
  }
  energy = best_energy;
  energy += mlclosing;
  /* logarithmic ML loop energy if logML */
  if ( (!is_extloop) && logML && (u>6) )
    energy += 6*mlbase+(int)(P->lxc*log((double)u/6.));
  else
    energy += mlbase*u;
  /* fprintf(stderr, "\n"); */
  return energy;
}

/*---------------------------------------------------------------------------*/

PRIVATE int cut_in_loop(int i) {
  /* walk around the loop;  return j pos of pair after cut if
     cut_point in loop else 0 */
  int  p, j;
  p = j = pair_table[i];
  do {
    i  = pair_table[p];  p = i+1;
    while ( pair_table[p]==0 ) p++;
  } while (p!=j && SAME_STRAND(i,p));
  return SAME_STRAND(i,p) ? 0 : pair_table[p];
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
    int hx, *stack;
    char type;
    stack = (int *) space(sizeof(int)*(n+1));
    
    for(hx=0, j=1; j<=n; j++) {
      switch (structure[j-1]) {
      case '|': BP[j] = -1; break;
      case 'x': /* can't pair */ 
        for (l=1; l<j-TURN; l++) ptype[indx[j]+l] = 0;
        for (l=j+TURN+1; l<=n; l++) ptype[indx[l]+j] = 0;
        break;
      case '(':
        stack[hx++]=j;
        /* fallthrough */
      case '<': /* pairs upstream */
        for (l=1; l<j-TURN; l++) ptype[indx[j]+l] = 0;
        break;
      case ')':
        if (hx<=0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced brackets in constraints");
        }
        i = stack[--hx];
        type = ptype[indx[j]+i];
	for (k=i+1; k<=n; k++) ptype[indx[k]+i] = 0;
	/* don't allow pairs i<k<j<l */
        for (l=j; l<=n; l++)
	  for (k=i+1; k<=j; k++) ptype[indx[l]+k] = 0;
	/* don't allow pairs k<i<l<j */
	for (l=i; l<=j; l++)
	  for (k=1; k<=i; k++) ptype[indx[l]+k] = 0;
	for (k=1; k<j; k++) ptype[indx[j]+k] = 0;
        ptype[indx[j]+i] = (type==0)?7:type;
        /* fallthrough */
      case '>': /* pairs downstream */
        for (l=j+TURN+1; l<=n; l++) ptype[indx[l]+j] = 0;
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
