/* Last changed Time-stamp: <2008-06-27 17:21:42 ivo> */

/**
*** \file fold_vars.c
*** global variables to change behaviour of folding routines<BR>
*** Also there are some functions that make the live easier when
*** using functions of the Vienna RNA package
**/
#include <string.h>
#include <stdio.h>
#include "fold_vars.h"

int  circ = 0;
int  noGU = 0;           /* GU not allowed at all */
int  no_closingGU = 0;   /* GU allowed only inside stacks */
int  tetra_loop = 1;     /* Fold with specially stable 4-loops */
int  energy_set = 0;     /* 0 = BP; 1=any with GC; 2=any with AU parameters */
int  dangles = 2;	 /* use dangling end energies */
char *nonstandards = (char *)0; /* contains allowed non standard bases */
double temperature = 37.0;
int  james_rule = 1;     /* interior loops of size 2 get energy 0.8Kcal and
			    no mismatches (no longer used) */
int  oldAliEn   = 0;     /* use old alifold-energies (without removing gaps) */
int  ribo       = 0;     /* use ribosum instead of classic covariance term */
char *RibosumFile = NULL; /* TODO: compile ribosums into program
			     Warning: this variable will vanish */
int  csv        = 0;     /*generate comma seperated output*/
bondT  *base_pair;

FLT_OR_DBL *pr;          /* base pairing prob. matrix */
int  *iindx;             /* pr[i,j] -> pr[iindx[i]-j] */
double pf_scale=- 1;     /* scaling factor to avoid floating point overflows */
int   fold_constrained = 0; /* fold with constraints */
int   do_backtrack=1;     /* calculate pair prob matrix in part_func() */
int    noLonelyPairs = 0; /* avoid helices of length 1 */
char backtrack_type='F';  /* 'C' require (1,N) to be bonded;
			     'M' seq is part of s multi loop */

int *cut_points;
int *strand;

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
