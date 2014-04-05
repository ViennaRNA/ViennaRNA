/* Last changed Time-stamp: <2008-06-27 17:21:42 ivo> */

/**
*** \file fold_vars.c
*** global variables to change behaviour of folding routines<BR>
*** Also there are some functions that make the live easier when
*** using functions of the Vienna RNA package
**/
#include <string.h>
#include <stdio.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"

int         james_rule = 1;       /* interior loops of size 2 get energy 0.8Kcal and
                                    no mismatches (no longer used) */



char        *RibosumFile = NULL;  /* TODO: compile ribosums into program
                                    Warning: this variable will vanish */

int         csv = 0;              /*generate comma seperated output*/

bondT       *base_pair = NULL;

FLT_OR_DBL  *pr = NULL;           /* base pairing prob. matrix */

int         *iindx = NULL;        /* pr[i,j] -> pr[iindx[i]-j] */

int         fold_constrained = 0; /* fold with constraints */


int         *cut_points;

int         *strand;


PUBLIC char * option_string(void){
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

