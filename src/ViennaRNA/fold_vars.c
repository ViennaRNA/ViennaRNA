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

int         circ = 0;
int         uniq_ML   = 0;        /* do ML decomposition uniquely (for subopt) */

int         noGU = 0;             /* GU not allowed at all */

int         no_closingGU = 0;     /* GU allowed only inside stacks */

int         tetra_loop = 1;       /* Fold with specially stable 4-loops */

int         energy_set = 0;       /* 0 = BP; 1=any with GC; 2=any with AU parameters */

int         dangles = 2;          /* use dangling end energies */

char        *nonstandards = (char *)0;  /* contains allowed non standard bases */

double      temperature = 37.0;

int         james_rule = 1;       /* interior loops of size 2 get energy 0.8Kcal and
                                    no mismatches (no longer used) */

int         oldAliEn = 0;         /* use old alifold-energies (without removing gaps) */

int         ribo = 0;             /* use ribosum instead of classic covariance term */

char        *RibosumFile = NULL;  /* TODO: compile ribosums into program
                                    Warning: this variable will vanish */

int         csv = 0;              /*generate comma seperated output*/

bondT       *base_pair = NULL;

FLT_OR_DBL  *pr = NULL;           /* base pairing prob. matrix */

int         *iindx = NULL;        /* pr[i,j] -> pr[iindx[i]-j] */

double      pf_scale = -1;        /* scaling factor to avoid floating point overflows */

int         fold_constrained = 0; /* fold with constraints */

int         do_backtrack = 1;     /* calculate pair prob matrix in part_func() */

int         noLonelyPairs = 0;    /* avoid helices of length 1 */

char        backtrack_type = 'F'; /* 'C' require (1,N) to be bonded;
                                    'M' seq is part of s multi loop */

int         *cut_points;

int         canonicalBPonly = 0;  /* remove non-canonical base pairs from structure constraint */

int         *strand;

int         gquad = 0;            /* consider g-qudruplexes in the calculations */

int         max_bp_span = -1;     /* base pairs may span the entire sequence length by default */

double      cv_fact = 1.;

double      nc_fact = 1.;

PRIVATE     int rtype[8] = {0, 2, 1, 4, 3, 6, 5, 7};

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

PUBLIC void set_model_details(model_detailsT *md){
  int i = 0;

  if(md){
    md->dangles         = dangles;
    md->special_hp      = tetra_loop;
    md->noLP            = noLonelyPairs;
    md->noGU            = noGU;
    md->noGUclosure     = no_closingGU;
    md->logML           = logML;
    md->gquad           = gquad;
    md->canonicalBPonly = canonicalBPonly;
    md->circ            = circ;
    md->uniq_ML         = uniq_ML;
    md->do_backtrack    = do_backtrack;
    md->backtrack_type  = backtrack_type;
    md->energy_set      = energy_set;
    md->max_bp_span     = max_bp_span;
    md->min_loop_size   = TURN;
    md->oldAliEn        = oldAliEn;
    md->ribo            = ribo;
    md->cv_fact         = cv_fact;
    md->nc_fact         = nc_fact;
    md->temperature     = temperature;
    md->betaScale       = 1.;
    md->pf_scale        = pf_scale;

    if(nonstandards){
      memcpy(md->nonstandards, nonstandards, strlen(nonstandards)*sizeof(char));
    } else {
      md->nonstandards[0] = (char)0;
    }
    /* set default values for the pair/rtype[pair] stuff */
    memcpy(md->rtype, &(rtype[0]), 8 * sizeof(int));
    memset(md->alias, 0, (MAXALPHA + 1) * sizeof(short));
    for(i = 0;i <= MAXALPHA; i++)
      memset(md->pair[i], 0, (MAXALPHA + 1) * sizeof(int));

    fill_pair_matrices(md);

  }
}
