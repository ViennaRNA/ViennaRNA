/* Last changed Time-stamp: <96/01/24 18:52:24 ivo> */
/*                
			 minimum free energy
		  RNA secondary structure prediction

			    c Ivo Hofacker
		      original implementation by
			    Walter Fontana
		  
			  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"

#define PAREN
#ifdef LETTER
#undef PAREN
#endif

#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */

PUBLIC float  fold(char *string, char *structure);
PUBLIC float  energy_of_struct(char *string, char *structure);
PUBLIC void   free_arrays(void);
PUBLIC void   initialize_fold(int length);
PUBLIC void   update_fold_params(void);
PUBLIC int    pf_dangl=0;

PRIVATE void    letter_structure(char *structure, int length);
PRIVATE void    parenthesis_structure(char *structure, int length);
PRIVATE void    get_arrays(int size);
PRIVATE void    initialize(int length);
PRIVATE void    scale_parameters(void);
PRIVATE void    make_pair_table(char *structure, short *table);
PRIVATE int     stack_energy(int i, char *string);
PRIVATE void    BP_calculate(char *structure, int *BP, int length);

#define MAXSECTORS      500     /* dimension for a backtrack array */
#define LOCALITY        0.      /* locality parameter for base-pairs */

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

PRIVATE int stack[NBPAIRS+1][NBPAIRS+1];
PRIVATE int hairpin[31];
PRIVATE int bulge[MAXLOOP+1];
PRIVATE int internal_loop[MAXLOOP+1];
PRIVATE int F_ninio[5];
PRIVATE double lxc;
PRIVATE int MLbase;
PRIVATE int MLintern[NBPAIRS+1];
PRIVATE int MLclosing[NBPAIRS+1];
PRIVATE int TETRA_ENERGY;

PRIVATE int *indx;  /* index for moving in the triangle matrices c[] and fMl[]*/

PRIVATE int   *c;       /* energy array, given that i-j pair */
PRIVATE int   *f5;      /* energy of 5' end */
PRIVATE int   *f3;      /* energy of 3' end */
PRIVATE int   *fML;     /* multi-loop auxiliary energy array */
PRIVATE int   *Fmi;     /* holds row i of fML (avoids jumps in memory) */
PRIVATE int   *DMLi;    /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
PRIVATE int   *DMLi1;   /*             MIN(fML[i+1,k]+fML[k+1,j])  */
PRIVATE int   *DMLi2;   /*             MIN(fML[i+2,k]+fML[k+1,j])  */
PRIVATE short  *S, *S1;
PRIVATE short  *pair_table;

PRIVATE char  alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

/*--------------------------------------------------------------------------*/

void initialize_fold(int length)
{
    get_arrays(length);
    initialize(length);
    scale_parameters();
    make_pair_matrix();
}
    
/*--------------------------------------------------------------------------*/

PRIVATE void get_arrays(int size)
{
   indx = (int *) space(sizeof(int)*(size+1));
   c     = (int *) space(sizeof(int)*((size*(size+1))/2+2));
   fML   = (int *) space(sizeof(int)*((size*(size+1))/2+2));
   f5    = (int *) space(sizeof(int)*(size+2));
   f3    = (int *) space(sizeof(int)*(size+2));
   Fmi   = (int *) space(sizeof(int)*(size+1));
   DMLi  = (int *) space(sizeof(int)*(size+1));
   DMLi1  = (int *) space(sizeof(int)*(size+1));
   DMLi2  = (int *) space(sizeof(int)*(size+1));
   base_pair = (struct bond *) space(sizeof(struct bond)*(1+size/2));
}

/*--------------------------------------------------------------------------*/

void free_arrays(void)
{
   free(indx); free(c); free(fML); free(f5); free(f3);
   free(base_pair); free(Fmi);
   free(DMLi); free(DMLi1);free(DMLi2);
}

/*--------------------------------------------------------------------------*/

PRIVATE void initialize(int length)
{
    register int n;

    for (n = 1; n <= length; n++)
	indx[n] = (n*(n-1)) >> 1;        /* n(n-1)/2 */
}

/*----------------------------------------------------------------------------*/

float fold(char *string, char *structure)
{
   struct sect {
      int  i;
      int  j;
      int ml;
   } 
   sector[MAXSECTORS];   /* backtracking sectors */
   
   int   i, j, k, l, p, q, length, energy, new_f, new_c, new;
   int   fij, fi, fj, cij, ci1j, cij1, ci1j1, max_separation;
   int   decomp, MLenergy, new_fML, ml;
   int   s, b, unpaired, traced, sizecorr, mm;
   int   no_close, no_close_2, type, type_2, tt;
   int   n1, n2, m, tetracorr;
   char *pos;
   int  *FF;
   int   bonus=0, bonus_cnt=0;
   int   *BP; /* contains the structure constrainsts: BP[i]
                      -1: | = base must be paired
                      -2: < = base must be paired with j<i
                      -3: > = base must be paired with j>i
                      -4: x = base must not pair
		 positive int: base is paired with int      */
   
   length = strlen(string);
   no_close_2=0;
   BP = (int *)space(sizeof(int)*(length+2));
   if (fold_constrained) BP_calculate(structure,BP,length);
    
   S = (short *) space(sizeof(short)*(length+2));
   S1= (short *) space(sizeof(short)*(length+2));
   /* S1 exists only for the special X K and I bases and energy_set!=0 */
   
   for (l=1; l<=length; l++) { /* make numerical encoding of sequence */
	 if (energy_set>0) S[l]=string[l-1]-'A'+1;
      else {
	 pos = strchr(Law_and_Order, string[l-1]);
	 if (pos==NULL) S[l]=0;
	 else S[l]= pos-Law_and_Order;
      }
      S1[l] = alias[S[l]];   /* for mismatches of nostandard bases */
   }
   S[l]=S1[l]=0;
   
   max_separation = (int) ((1.-LOCALITY)*(double)(length-2)); /* not in use */

   for (j=1; j<=length; j++) {
      Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
   }
   
   for (j = 1; j<=length; j++)
      for (i=(j>TURN?(j-TURN):1); i<j; i++)
	 c[indx[j]+i] = fML[indx[j]+i] = INF;

   for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */
      
      for (j = i+TURN+1; j <= length; j++) {

	 bonus = 0;
	 type = pair[S[i]][S[j]];

	 /* enforcing structure constraints */
	 if ((BP[i]==j) && (type==0)) type=7; /* nonstandard */
	 if ((BP[i]==j)||(BP[i]==-1)||(BP[i]==-2)) bonus -= BONUS;
	 if ((BP[j]==-1)||(BP[j]==-3)) bonus -= BONUS;
	 if ((BP[i]==-4)||(BP[j]==-4)) bonus +=BONUS;
	 	 
	 no_close = (((type==3)||(type==4))&&no_closingGU&&(bonus==0));
		 
	 /* smallest hairpin case ------------------------------------*/
      
	 if (i == j-TURN-1) {
	    c[indx[j]+i] = (type) ? hairpin[TURN] /* +
	       mismatch[S1[i]][S1[j]][S1[i+1]][S1[j-1]]*/ : INF;

	    c[indx[j]+i] += bonus;
	    
	    if (no_close) c[indx[j]+i] = FORBIDDEN;
	    
	    fML[indx[j]+i] = Fmi[j] = c[indx[j]+i]+MLintern[type];
	    continue;
	 }
	    
	 if (j-i-1 > max_separation) type = 0;  /* forces locality degree */
	    
	 if (type!=0) {

	    /* hairpin ----------------------------------------------*/
		
	    if (no_close) new_c = FORBIDDEN;
	    else {
	       new_c = (j-i-1 <= 30) ? hairpin[j-i-1] :
		  hairpin[30]+(int)(lxc*log((double)(j-i-1)/30.));
	       if (tetra_loop)
		  if (j-i-1 == 4) {
		     for (k = 0; k < N_TETRALOOPS; k++) {
			if (strncmp(string+i, Tetraloops[k], 4) == 0) 
			   new_c += TETRA_ENERGY;
		     }
		  }
	       new_c += mismatch[S1[i]][S1[j]][S1[i+1]][S1[j-1]];
	    }
	    
	    /*--------------------------------------------------------
	      check for elementary structures involving more than one
	      closing pair.
	      --------------------------------------------------------*/
	    
	    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++) {
	       
	       for (q = j-1; q >= p+1+TURN; q--) {

		  n1 = p-i-1;
		  n2 = j-q-1;

		  if (n1+n2 > MAXLOOP) break;

		  if (n1>n2) { m=n1; n1=n2; n2=m; } /* so that n2>=n1 */

		  type_2 = pair[S[p]][S[q]];
		  if ((BP[p]==q) && (type_2==0)) type_2=7; /* nonstandard */

		  if (type_2==0) continue;

		  if (no_closingGU)
		     no_close_2 = (no_close || ((type_2==3)||(type_2==4)));

		  if (n2 == 0)
		     energy = stack[type][type_2];   /* stack */

		  else if (n1==0) {                  /* bulge */

		     if (!no_close_2) 
			energy = bulge[n2];
		     else
			energy = FORBIDDEN;
#if STACK_BULGE1
		     if (n2==1) energy+=stack[type][type_2];
#endif
		  } else {                           /* interior loop */

		     if (!no_close_2) {
			energy = internal_loop[n1+n2];
			
			m       = MIN2(4, n1);
			energy += MIN2(MAX_NINIO,((n2-n1)*F_ninio[m]));
			
			energy += mismatch[S1[i]][S1[j]][S1[i+1]][S1[j-1]]+
			          mismatch[S1[p-1]][S1[q+1]][S1[p]][S1[q]];
		     }
		     else
			energy = FORBIDDEN;
		  }
			
		  new_c = MIN2(energy+c[indx[q]+p], new_c);
			
	       } /* end q-loop */
	    } /* end p-loop */

	    /* multi-loop decomposition ------------------------*/

	    if (!no_close) {
	       decomp = DMLi1[j-1];
	       if (dangles) {
		  tt = pair[S[j]][S[i]]; if ((tt==0)&&(bonus!=0)) tt=7;
		  decomp = MIN2(DMLi2[j-1]+dangle3[tt][S1[i+1]]+MLbase, decomp);
		  decomp = MIN2(DMLi1[j-2]+dangle5[tt][S1[j-1]]+MLbase, decomp);
		  decomp = MIN2(DMLi2[j-2]+dangle5[tt][S1[j-1]]+
				dangle3[tt][S1[i+1]] + 2*MLbase, decomp);
	       }

	       MLenergy = MLclosing[type]+decomp;

	       new_c = MLenergy < new_c ? MLenergy : new_c;
	    }
	    
	    c[indx[j]+i] = new_c + bonus;

	 } /* end >> if (pair) << */
	 
	 else c[indx[j]+i] = INF;

	 /* free ends ? -----------------------------------------*/

	 new_fML = fML[indx[j]+i+1]+MLbase;
	 new_fML = MIN2(fML[indx[j-1]+i]+MLbase, new_fML);
	 new_fML = MIN2(c[indx[j]+i]+MLintern[type], new_fML);
	 if (dangles) {
	    tt = pair[S[i+1]][S[j]]; if ((tt==0)&&(BP[i+1]==j)) tt=7;
	    new_fML = MIN2(c[indx[j]+i+1]+dangle5[tt][S1[i]]
			   +MLintern[tt]+MLbase,new_fML);
	    tt = pair[S[i]][S[j-1]]; if ((tt==0)&&(BP[i]==j-1)) tt=7;
	    new_fML = MIN2(c[indx[j-1]+i]+dangle3[tt][S1[j]]
			   +MLintern[tt]+MLbase, new_fML);
	 }
	 if (dangles) {
	    tt = pair[S[i+1]][S[j-1]]; if ((tt==0)&&(BP[i+1]==j-1)) tt=7;
	    new_fML = MIN2(c[indx[j-1]+i+1]+dangle5[tt][S1[i]]+
			   dangle3[tt][S1[j]]+MLintern[tt]+2*MLbase, new_fML);
	 }
		 
	 /* modular decomposition -------------------------------*/

	 for (decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++)
	    decomp = MIN2(decomp, Fmi[k]+fML[indx[j]+k+1]);

	 DMLi[j] = decomp;               /* store for use in ML decompositon */
	 new_fML = MIN2(new_fML,decomp);
	 
	 fML[indx[j]+i] = Fmi[j] = new_fML;     /* substring energy */
	 
      }
      
      FF = DMLi2;
      DMLi2 = DMLi1;
      DMLi1 = DMLi;
      DMLi = FF;
      for (j=1; j<=length; j++) { Fmi[j]=DMLi[j]=INF; }
   }

   /* calculate energies of 5' and 3' fragments */
   
   f5[TURN+1]=0;
   for (j=TURN+2; j<=length; j++) {
      f5[j] = f5[j-1];
      f5[j] = MIN2(f5[j], c[indx[j]+1]);
      type=pair[S[1]][S[j-1]]; if ((type==0)&&(BP[1]==j-1)) type=7;
      if ((type)&&(dangles))
	 f5[j] = MIN2(f5[j], c[indx[j-1]+1]+dangle3[type][S1[j]]);
      for (i=j-TURN-1; i>1; i--) {
	 type = pair[S[i]][S[j]]; if ((type==0)&&(BP[i]==j)) type=7;
	 if (type) {
	    f5[j] = MIN2(f5[j], f5[i-1]+c[indx[j]+i]);
	    if (dangles)
	       f5[j] = MIN2(f5[j], f5[i-2]+c[indx[j]+i]+dangle5[type][S1[i-1]]);
	 }
	 type = pair[S[i]][S[j-1]]; if ((type==0)&&(BP[i]==j-1)) type=7;
	 if ((type)&&(dangles)) {
	    f5[j] = MIN2(f5[j], f5[i-1]+c[indx[j-1]+i]+dangle3[type][S1[j]]);
	    f5[j] = MIN2(f5[j], f5[i-2]+c[indx[j-1]+i]+
			 dangle5[type][S1[i-1]]+dangle3[type][S1[j]]);
	 }
      }
   }

#if 0
   for (i=length-TURN-1; i>0; i--) {
      f3[i] = f3[i+1];
      for (j=i+TURN+1; j<length; j++) {
	 type = pair[S[i]][S[j]]; if ((type==0)&&(BP[i]==j)) type=7;
	 if (type) {
	    f3[i] = MIN2(f3[i], c[indx[j]+i] + f3[j+1]);
	    if(dangles)
	       f3[i] = MIN2(f3[i], c[indx[j]+i] + f3[j+2] +dangle3[type][S1[j+1]]);
	 }
	 type = pair[S[i+1]][S[j]]; if ((type==0)&&(BP[i+1]==j)) type=7;
	 if((type)&&(dangles)) {
	    f3[i] = MIN2(f3[i], c[indx[j]+i+1] + f3[j+1] +dangle5[type][S1[i]]);
	    f3[i] = MIN2(f3[i], c[indx[j]+i+1] + f3[j+2] +
			 dangle5[type][S1[i]] + dangle3[type][S1[j+1]]);
	 }
      }
      f3[i] = MIN2(f3[i], c[indx[length]+i]);
      type=pair[S[i+1]][S[length]]; if ((type==0)&&(BP[i+1]==j)) type=7;
      if ((type)&&(dangles))
	 f3[i] = MIN2(f3[i], c[indx[length]+i+1] +dangle5[type][S1[i]]);
   }
   
   if (f3[1]!=f5[length])
      fprintf(stderr, "f3[1]!=f5[n]! %d  %d\n",f3[1],f5[length]);
#endif
   
   /*------------------------------------------------------------------
     trace back through the "c", "f5" and "fML" arrays to get the
     base pairing list. No search for equivalent structures is done.
     This inverts the folding procedure, hence it's very fast.
     ------------------------------------------------------------------*/

   ci1j = cij1 = ci1j1 = INF;
   b = 0;
   s = 1;
   sector[s].i = 1;
   sector[s].j = length;
   sector[s].ml = (backtrack_type=='M') ? 1 : 0;
   
   do {
      i  = sector[s].i;
      j  = sector[s].j;
      ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to 
				occur in the fML- (1) or in the f-array (0) */
      if ((i>1)&&(!ml))
	 nrerror("Error while backtracking");
      
      if (((j-i+1)==length)&&(backtrack_type=='C')) {
	 base_pair[++b].i = i;
	 base_pair[b].j   = j;
	 goto repeat1;
      }
      
      if (j < i+TURN+1)
	 continue;

      if (!ml) { /* find outermost base pairs */
	 fij = f5[j];
	 if (fij == f5[j-1]) {
	    sector[++s].i = i;
	    sector[s].j   = j-1;
	    sector[s].ml  = 0;
	    continue;
	 } else {
	    int jj;
	    for (k=j-TURN-1,traced=0; k>1; k--) {
	       jj = k-1;
	       type = pair[S[k]][S[j-1]]; if ((type==0)&&(BP[k]==j-1)) type=7;
	       if((type)&&(dangles)) {
		  if (f5[j] == f5[k-1]+c[indx[j-1]+k]+dangle3[type][S1[j]])
		     traced=j-1;
		  if (f5[j] == f5[k-2]+c[indx[j-1]+k]+
		      dangle5[type][S1[k-1]]+dangle3[type][S1[j]]) {
		     traced=j-1; jj=k-2;
		  }
	       }
	       type = pair[S[k]][S[j]]; if ((type==0)&&(BP[k]==j)) type=7;
	       if(type) {
		  if (f5[j] == f5[k-1]+c[indx[j]+k]) traced=j;
		  if (dangles)
		     if (f5[j] == f5[k-2]+c[indx[j]+k]+dangle5[type][S1[k-1]]) {
			traced=j; jj=k-2;
		     }
	       }
	       if (traced) break;
	    }
	    if (!traced) {
	       k=1; jj=k-1; if (f5[j] == c[indx[j]+k]) traced=j;
	       type=pair[S[1]][S[j-1]];if ((type==0)&&(BP[1]==j-1)) type=7;
	       if ((type)&&(dangles))
		  if (f5[j] == c[indx[j-1]+k]+dangle3[type][S1[j]]) traced=j-1;
	    }
	    if (!traced) nrerror("backtrack failed");
	    sector[++s].i = 1;
	    sector[s].j   = jj;
	    sector[s].ml  = 0;

	    i=k; j=traced;
	    base_pair[++b].i = i;
	    base_pair[b].j   = j;
	    goto repeat1;
	 }
	 
      } else { /* trace back in fML array */
	 fij = fML[indx[j]+i];
	 fj  = fML[indx[j]+i+1]+MLbase;
	 fi  = fML[indx[j-1]+i]+MLbase;
	 tt  = pair[S[i]][S[j]]; if ((tt==0)&&(BP[i]==j)) tt=7;
	 cij = c[indx[j]+i] + MLintern[tt];
	 if (dangles) {
	    tt = pair[S[i+1]][S[j]]; if ((tt==0)&&(BP[i+1]==j)) tt=7;
	    ci1j= c[indx[j]+i+1] + dangle5[tt][S1[i]] + MLintern[tt]+MLbase;
	    tt = pair[S[i]][S[j-1]]; if ((tt==0)&&(BP[i]==j-1)) tt=7;
	    cij1= c[indx[j-1]+i] + dangle3[tt][S1[j]] + MLintern[tt]+MLbase;
	    tt = pair[S[i+1]][S[j-1]]; if ((tt==0)&&(BP[i+1]==j-1)) tt=7;
	    ci1j1=c[indx[j-1]+i+1] + dangle5[tt][S1[i]] + dangle3[tt][S1[j]]
	       +  MLintern[tt] + 2*MLbase;
	 }
	 if (fij == fj) {
	    sector[++s].i = i+1;
	    sector[s].j   = j;
	    sector[s].ml  = ml;
	    continue;
	 }
	 else if (fij == fi) {
	    sector[++s].i = i;
	    sector[s].j   = j-1;
	    sector[s].ml  = ml;
	    continue;
	 }
	 else {
	    if ((fij==cij)||(fij==ci1j)||(fij==cij1)||(fij==ci1j1)) {
	       if (fij==ci1j) i++;
	       else if (fij==cij1) j--;
	       else if (fij==ci1j1) {i++; j--;}
	       base_pair[++b].i = i;
	       base_pair[b].j   = j;
	       goto repeat1;
	    } 
	 }
      
	 for (k = i+1+TURN; k <= j-2-TURN; k++) {
	    if (fML[indx[j]+i] == (fML[indx[k]+i]+fML[indx[j]+k+1])) {
	       sector[++s].i = i;
	       sector[s].j   = k;
	       sector[s].ml  = ml;
	       sector[++s].i = k+1;
	       sector[s].j   = j;
	       sector[s].ml  = ml;
	       break;
	    }
	 }
      }
      
      if (k>j-2-TURN) fprintf(stderr,"backtrack failed\n");
      continue;
      
    repeat1:
      
      /*----- begin of "repeat:" -----*/
      
      type = pair[S[i]][S[j]];
      if ((BP[i]==j) && (type==0)) type=7; /* nonstandard */
      
      no_close = (((type==3)||(type==4))&&no_closingGU&&(bonus==0));

      unpaired = j-i-1;
      sizecorr = 0;
      tetracorr = 0;
      if (unpaired > 30) {
	 unpaired = 30;
	 sizecorr = (int)(lxc*log((double)(j-i-1)/30.));
      }
      else if (tetra_loop) if (unpaired == 4) {
	 for (k = 0; k < N_TETRALOOPS; k++) {
	    if (strncmp(string+i, Tetraloops[k], 4) == 0) 
	       tetracorr = TETRA_ENERGY;
	 }
      }
      mm = (unpaired>3) ? mismatch[S1[i]][S1[j]][S1[i+1]][S1[j-1]] : 0;
      
      bonus = 0;
      if ((BP[i]==j)||(BP[i]==-1)||(BP[i]==-2))
	 bonus -= BONUS;
      if ((BP[j]==-1)||(BP[j]==-3)) bonus -= BONUS;
      if (no_close) {
	 if (c[indx[j]+i] == FORBIDDEN) continue;
      } else
	 if (c[indx[j]+i] == hairpin[unpaired]+sizecorr+tetracorr+bonus+mm)
	    continue;
   
      for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
	 
	 for (q = j-1; q >= p+1+TURN; q--) {

	    n1 = p-i-1;
	    n2 = j-q-1;
	    if (n1+n2 > MAXLOOP) break;
	    
	    if (n1>n2) { m=n1; n1=n2; n2=m; }
	    
	    type_2 = pair[S[p]][S[q]];
	    if ((BP[p]==q) && (type_2==0)) type_2=7; /* nonstandard */
	    
	    if (type_2==0) continue;
	    
	    if (no_closingGU)
	       no_close_2 = (no_close || ((type_2==3)||(type_2==4)));

	    if (n2 == 0) energy = stack[type][type_2]; 
	    else if (n1==0) {                    /* bulge */
	       if (!no_close_2) 
		  energy = bulge[n2];
	       else
		  energy = FORBIDDEN;
#if STACK_BULGE1
	       if (n2==1) energy+=stack[type][type_2];
#endif
	    } else {                             /* interior loop */
	       if (!no_close_2) {
		  energy = internal_loop[n1+n2];
	       
		  m       = MIN2(4, n1);
		  energy += MIN2(MAX_NINIO,(abs(n1-n2)*F_ninio[m]));
		  
		  energy += mismatch[S1[i]][S1[j]][S1[i+1]][S1[j-1]]+
		            mismatch[S1[p-1]][S1[q+1]][S1[p]][S1[q]];
	       }
	       else
		  energy = FORBIDDEN;
	    }
	    
	    new = energy+c[indx[q]+p]+bonus;
	    traced = (c[indx[j]+i] == new);
	    if (traced) {
	       base_pair[++b].i = p;
	       base_pair[b].j   = q;
	       i = p, j = q;
	       goto repeat1;
	    }
	 }
      }
      
      /* end of repeat: --------------------------------------------------*/


      mm = bonus+MLclosing[type];
      tt = pair[S[j]][S[i]]; if ((tt==0)&&(BP[i]==j)) tt=7;

      for (k = i+2+TURN; k <= j-3-TURN; k++) {
	 if (c[indx[j]+i] == (fML[indx[k]+i+1]+fML[indx[j-1]+k+1]+mm)) {
	    sector[++s].i = i+1;
	    sector[s+1].j = j-1;
	    break;
	 }
	 if (dangles) {
	    if (c[indx[j]+i] == (fML[indx[k]+i+2]+fML[indx[j-1]+k+1]+mm+
				 dangle3[tt][S1[i+1]]+MLbase)) {
	       sector[++s].i = i+2;
	       sector[s+1].j = j-1;
	       break;
	    }
	    if (c[indx[j]+i] == (fML[indx[k]+i+1]+fML[indx[j-2]+k+1]+mm+
				 dangle5[tt][S1[j-1]]+MLbase)) {
	       sector[++s].i = i+1;
	       sector[s+1].j = j-2;
	       break;
	    }
	    if (c[indx[j]+i] == (fML[indx[k]+i+2]+fML[indx[j-2]+k+1]+mm+
				 dangle3[tt][S1[i+1]]+dangle5[tt][S1[j-1]]+
				 2*MLbase)) {
	       sector[++s].i = i+2;
	       sector[s+1].j = j-2;
	       break;
	    }
	 }
      }
      sector[s].ml  = 1; 
      sector[s].j   = k;
      sector[++s].i = k+1;
      sector[s].ml  = 1; 

      if (k>j-3-TURN) fprintf(stderr, "backtracking failed\n");
   }
   while (s > 0);

   base_pair[0].i = b;    /* save the total number of base pairs */

#ifdef LETTER
   letter_structure(structure, length);
#endif
#ifdef PAREN
   parenthesis_structure(structure, length);
#endif

   free(S); free(S1);
  
   bonus=0;
   bonus_cnt = 0;
   for(l=1;l<=length;l++) {
      if((BP[l]<0)&&(BP[l]>-4)) {
	 bonus_cnt++;
	 if((BP[l]==-3)&&(structure[l-1]==')')) bonus++;
	 if((BP[l]==-2)&&(structure[l-1]=='(')) bonus++;
	 if((BP[l]==-1)&&(structure[l-1]=='(')) bonus++;
	 if((BP[l]==-1)&&(structure[l-1]==')')) bonus++;
      }
      if(BP[l]>l) {
	 bonus_cnt++;
	 for(i=1;i<=b;i++)
	    if((l==base_pair[i].i)&&(BP[l]==base_pair[i].j)) bonus++;
      }
   }

   if (bonus_cnt>bonus) fprintf(stderr,"\ncould not enforce all constraints\n");
   bonus*=BONUS;

   free(BP);
   
   f5[length] += bonus;      

   if (backtrack_type=='C')
      return (float) c[indx[length]+1]/100.;
   else if (backtrack_type=='M')
      return (float) fML[indx[length]+1]/100.;
   else
      return (float) f5[length]/100.;
}

/*---------------------------------------------------------------------------*/


PRIVATE void letter_structure(char *structure, int length)
{
    int n, k, x, y;
    
    for (n = 0; n <= length-1; *(structure+n++) = ' ') ;
    *(structure+length) = '\0';

    for (n = 0, k = 1; k <= base_pair[0].i; k++) {
	y = base_pair[k].j;
	x = base_pair[k].i;
	if (x-1 > 0 && y+1 <= length) {
	    if (*(structure+x-2) != ' ' && *(structure+y) == *(structure+x-2)) {
		*(structure+x-1) = *(structure+x-2);
		*(structure+y-1) = *(structure+x-1);
		continue;
	    }
	}
	if (*(structure+x) != ' ' && *(structure+y-2) == *(structure+x)) {
	    *(structure+x-1) = *(structure+x);
	    *(structure+y-1) = *(structure+x-1);
	    continue;
	}
	n++;
	*(structure+x-1) = *(alpha+n-1);
	*(structure+y-1) = *(alpha+n-1);
    }
}

/*---------------------------------------------------------------------------*/

PRIVATE void parenthesis_structure(char *structure, int length)
{
    int n, k;
    
    for (n = 0; n <= length-1; *(structure+n++) = '.') ;
    *(structure+length) = '\0';

    for (k = 1; k <= base_pair[0].i; k++) {
	*(structure+base_pair[k].i-1) = '(';
	*(structure+base_pair[k].j-1) = ')';
    }
}
/*---------------------------------------------------------------------------*/

PRIVATE void scale_parameters(void)
{
   int i,j;
   double tempf;

   tempf = ((temperature+K0)/(37.+K0));
   for (i=0; i<31; i++) 
      hairpin[i] = (int) hairpin37[i]*(tempf);
   for (i=0; i<=MIN2(30,MAXLOOP); i++) {
      bulge[i] = (int) bulge37[i]*tempf;
      internal_loop[i]= (int) internal_loop37[i]*tempf;
   }
   lxc = lxc37*tempf;
   for (; i<=MAXLOOP; i++) {
      bulge[i] = bulge[30]+(int)(lxc*log((double)(i)/30.));
      internal_loop[i] = internal_loop[30]+(int)(lxc*log((double)(i)/30.));
   }
   for (i=0; i<5; i++)
      F_ninio[i] = (int) F_ninio37[i]*tempf;
   TETRA_ENERGY = TETRA_ENERGY37*tempf;
   
   MLbase = ML_BASE37*tempf;
   for (i=0; i<=NBPAIRS; i++) {
      MLintern[i] = ML_intern37[i]*tempf;
      MLclosing[i] = ML_closing37[i]*tempf;
   }
   if (no_closingGU) {
      MLintern[3] = MLintern[4] = FORBIDDEN;
      MLclosing[3] = MLclosing[4] = FORBIDDEN;
   }

   for (i=0; i<=NBPAIRS; i++)
      for (j=0; j<=NBPAIRS; j++)
	 stack[i][j] = enthalpies[i][j] -
			entropies[i][j]*(temperature+K0)/100.0;
}

/*---------------------------------------------------------------------------*/
	   
PUBLIC void update_fold_params(void)
{
   scale_parameters();
   make_pair_matrix();
}

/*---------------------------------------------------------------------------*/

float energy_of_struct(char *string, char *structure)
{
   int   i, l, length, energy;
   int   tt, j, lastd, ee, ld3;
   char *pos;

   energy=lastd=ld3=0;
   length = strlen(string);
   if (strlen(structure)!=length)
      nrerror("energy_of_struct: string and structure have unequal length");

   S = (short *) space(sizeof(short)*(length+1));
   S1= (short *) space(sizeof(short)*(length+1));
   pair_table = (short * ) space(sizeof(short)*(length+2)); /* [0..l+1] */

   for (l=1; l<=length; l++) {
      if (energy_set>0) S[l]=string[l-1]-'A'+1;
      else {
	 pos = strchr(Law_and_Order, string[l-1]);
	 if (pos==NULL) S[l]=0;
	 else S[l]= pos-Law_and_Order;
      }
      S1[l] = alias[S[l]];   /* for mismatches of non standard bases */
   }
   
   make_pair_table(structure, pair_table);

   for (i=1; i<=length; i++) {
      if (pair_table[i]<0) {
	 if (backtrack_type=='M') energy+=MLbase;
	 continue;
      }
      if (dangles) {      /* dangling end contributions */
	 j=pair_table[i];
	 if ((i>0)&&((pair_table[i-1]<0)||(pf_dangl))) {
	    tt = pair[S[i]][S[j]]; if (tt==0) tt=7;
	    ee = dangle5[tt][S1[i-1]];               /* 5' dangle */
	    if ((i-1==lastd)&&(!pf_dangl)) ee -= ld3;     /* subtract 3' */
	    energy += (ee<0)?ee:0;
	 }
	 if ((j<length)&&((pair_table[j+1]<0)||(pf_dangl))) {
	    tt = pair[S[i]][S[j]]; if (tt==0) tt=7;
	    ld3 = dangle3[tt][S1[j+1]]; /* 3'dangle */
	    energy += ld3;
	    lastd = j+1;                             /* store last 3'dangle */
	 }
      }
      energy += stack_energy(i, string);          
      
      if (backtrack_type=='M') energy+=MLintern[1];
      i=pair_table[i];
   }
   free(pair_table);
   free(S); free(S1);
   return  (float) energy/100.;
}

/*---------------------------------------------------------------------------*/

PRIVATE int stack_energy(int i, char *string)  
{
   /* calculate energy of substructure enclosed by (i,j) */
   
   int energy;
   int j,p,q,u,type,type_2,i1,q2,n1,n2,m,k,tetra;
   int ld3, ee, lastd, tt;

   energy = lastd = ld3 = 0;
   j=pair_table[i];
   p=i; q=j;
   while (pair_table[++p]<0);
   while (pair_table[--q]<0);
   while (pair_table[p]==q) {
      if (p>q) break;
      type = pair[S[i]][S[j]]; if (type==0) type=7;
      
      if (type==7) 
	 fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", i, j,
		 string[i-1],string[j-1]);

      type_2 = pair[S[p]][S[q]]; if (type_2==0) type_2=7;
      
      
      n1 = p-i-1;
      n2 = j-q-1;

      if (n1>n2) { m=n1; n1=n2; n2=m;}

      if (n2 == 0) {
	 if ((type)&&(type_2))
	    energy += stack[type][type_2];       /* stack */
      }
      else if (n1==0) {                          /* bulge */
	 energy += (n2<=30) ? bulge[n2] :     
	    bulge[30]+(int)(lxc*log((double)(n2)/30.));
#if STACK_BULGE1
	 if ((n2==1)&&(type)&&(type_2)) energy+=stack[type][type_2];
#endif
      } else {                                   /* interior loop */
	 energy += (n1+n2<=30) ? internal_loop[n1+n2] :
	    internal_loop[30]+(int)(lxc*log((double)(n1+n2)/30.));
	 
	 m       = MIN2(4, n1);
	 energy += MIN2(MAX_NINIO,(abs(n1-n2)*F_ninio[m]));

	 if (type)
	    energy += mismatch[S1[i]][S1[j]][S1[i+1]][S1[j-1]];
	 if (type_2)
	    energy += mismatch[S1[p-1]][S1[q+1]][S1[p]][S1[q]];
	 
      }
      i=p; j=q;
      while (pair_table[++p]<0);    /* find next pair */
      while (pair_table[--q]<0);
   } /* end while */

   /* p,q don't pair must have found hairpin or multiloop */
   
   type = pair[S[i]][S[j]];
   if (type==0) {
      type=7;
      fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", i, j,
	      string[i-1],string[j-1]);
   }
   
   if (p>q) {                       /* hair pin */
      energy += (j-i-1 <= 30) ? hairpin[j-i-1] :
	 hairpin[30]+(int)(lxc*log((double)(j-i-1)/30.));
      if (tetra_loop) 
	 if (j-i-1 == 4) {
	    for (tetra = 0; tetra < N_TETRALOOPS; tetra++) {
	       if (strncmp(string+i, Tetraloops[tetra], 4) == 0) 
		  energy += TETRA_ENERGY;
	    }
	 }
      if ((type)&&(j-i-1>3))
	 energy += mismatch[S1[i]][S1[j]][S1[i+1]][S1[j-1]];
      return energy;
   }
   
   /* (i,j) is exterior pair of multiloop */

   tt = pair[S[j]][S[i]]; if (tt==0) tt=7;
   if ((dangles) && ((pair_table[i+1]<0)||(pf_dangl))) {
      ld3 = dangle3[tt][S1[i+1]];
      energy += ld3; 
      lastd = i+1;
   }
   energy += MLclosing[1];

   k=1; u=0; i1=i;
   while (p<q) {
      u+=p-i1-1; k++;
      q2=pair_table[p];

      tt = pair[S[p]][S[q2]]; if (tt==0) tt=7;
      if (dangles) {
	 if ((pair_table[p-1]<0)||(pf_dangl)) {
	    ee = dangle5[tt][S1[p-1]];
	    if ((p-1==lastd)&&(!pf_dangl)) ee -= ld3;
	    energy += (ee<0)?ee:0;
	 }
	 if ((pair_table[q2+1]<0)||(pf_dangl)) {
	    ld3 = dangle3[tt][S1[q2+1]];
	    energy += ld3;
	    lastd = q2+1;
	 }
      }
      energy += MLintern[1];
      energy += stack_energy(p, string);
      p=q2+1; i1=q2;
      while (pair_table[p]<0) p++;
   }

   tt = pair[S[j]][S[i]]; if (tt==0) tt=7;
   if ((dangles) && ((pair_table[j-1]<0)||(pf_dangl))) {
      ee = dangle5[tt][S1[j-1]];
      if ((j-1==lastd)&&(!pf_dangl)) ee -= ld3;
      energy += (ee<0)?ee:0;
   }
   
   u+=p-i1-1;
   energy += MLbase*u; /*+ MLclosing[1] + (k-1)*MLintern[1]; */
   return energy;
}

/*---------------------------------------------------------------------------*/ 

PRIVATE void make_pair_table(char *structure, short *table)
{
   int i,j,hx;
   short *olist;
   
   hx=0;
   olist = (short *) space(strlen(structure)/2*sizeof(short)); 
             
   for (i=0; i<strlen(structure); i++) {
      switch (structure[i]) {
       case '.':
         table[i+1]= -1;
         break;
       case '(': 
         olist[hx++]=i;
         break;
       case ')':
         j = olist[--hx];
         table[i+1]=j+1;
         table[j+1]=i+1;
         break;
      }
   }
   free(olist);
}

/*---------------------------------------------------------------------------*/

PRIVATE void BP_calculate(char *structure, int *BP, int length)
{
   int i,j,ct;
   
   for(i=0;i<length;i++) {
      switch (structure[i]) {
       case '|': BP[i+1] = -1; break;
       case '<': BP[i+1] = -2; break;
       case '>': BP[i+1] = -3; break;
       case 'x': BP[i+1] = -4; break;
       case '(': 
	 ct = 1;
	 for(j=i+1;j<length;j++) {
	    if(structure[j]==')') ct--;
	    if(structure[j]=='(') ct++;
	    if(ct==0) {
	       BP[i+1]=j+1;
	       BP[j+1]=i+1;
	       break;
	    }
	 }
	 break;
      }
   }
}
