/* Last changed Time-stamp: <2000-10-10 11:24:26 ivo> */
/*
       global variables to change behaviour of folding routines
			  Vienna RNA package
*/
#include <string.h>
#include <stdio.h>
#include "fold_vars.h"

int  noGU = 0;           /* GU not allowed at all */
int  no_closingGU = 0;   /* GU allowed only inside stacks */
int  tetra_loop = 1;     /* Fold with specially stable 4-loops */
int  energy_set = 0;     /* 0 = BP; 1=any with GC; 2=any with AU parameters */
int  dangles = 1;	 /* use dangling end energies */
char *nonstandards = (char *)0; /* contains allowed non standard bases */
double temperature = 37.0;
int  james_rule = 1;     /* interior loops of size 2 get energy 0.8Kcal and
			    no mismatches (no longer used) */
struct bond  *base_pair;

FLT_OR_DBL *pr;          /* base pairing prob. matrix */
int  *iindx;             /* pr[i,j] -> pr[iindx[i]-j] */
double pf_scale=- 1;     /* scaling factor to avoid floating point overflows */
int   fold_constrained = 0; /* fold with constraints */
int   do_backtrack=1;     /* calculate pair prob matrix in part_func() */
int    noLonelyPairs = 0; /* avoid helices of length 1 */
char backtrack_type='F';  /* 'C' require (1,N) to be bonded;
			     'M' seq is part of s multi loop */

char * option_string(void) {
  static char options[100];
  *options = '\0';
  if (noGU) strcat(options, "-noGU ");
  if (no_closingGU) strcat(options, "-noCloseGU ");
  if (!tetra_loop) strcat(options, "-4 ");
  if (noLonelyPairs) strcat(options, "-noLP ");
  if (fold_constrained) strcat(options, "-C ");
  if (dangles!=1) sprintf(options+strlen(options), "-d%d ", dangles);
  if (temperature!=37.0) 
    sprintf(options+strlen(options), "-T %f ", temperature);
  return options;
}
