/* Last changed Time-stamp: <1998-07-07 16:43:54 ivo> */
/*
       global variables to change behaviour of folding routines
			  Vienna RNA package
*/
#include "fold_vars.h"

int  noGU = 0;           /* GU not allowed at all */
int  no_closingGU = 0;   /* GU allowed only inside stacks */
int  tetra_loop = 1;     /* Fold with specially stable 4-loops */
int  energy_set = 0;     /* 0 = BP; 1=any with GC; 2=any with AU parameters */
int  dangles = 1;	 /* use dangling end energies */
char *nonstandards = (char *)0; /* contains allowed non standard bases */
float temperature = 37.0;
int  james_rule = 1;     /* interior loops of size 2 get energy 0.8Kcal and
			    no mismatches */
struct bond  *base_pair;

FLT_OR_DBL *pr;           /* base pairing prob. matrix */
int  *iindx;              /* pr[i,j] -> pr[iindx[i]-j] */
float pf_scale=- 1;       /* scaling factor to avoid floting point overflows */
int   fold_constrained = 0; /* guess what */
int   do_backtrack=1;     /* calculate pair prob matrix in part_func() */
int    noLonelyPairs = 0; /* avoid helices of length 1 */
char backtrack_type='F';  /* 'C' require (1,N) to be bonded;
			     'M' seq is part of s multi loop */

