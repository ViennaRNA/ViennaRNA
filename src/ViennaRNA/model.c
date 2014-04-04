/*
                  Model Details structure creation/modification/destruction

                  This file contains everything which is necessary to
                  obtain, modify, and destroy the model_details datastructure
                  used in the folding recurrences throughout the ViennaRNA
                  Package

                  c Ronny Lorenx

                  Vienna RNA package
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/energy_const.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/model.h"

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/

/*  below are the evil global variables that will vanish
    as soon as we drop backward compatibility in ViennaRNA
    Package v3
*/

double      temperature = 37.0;

double      pf_scale = -1;        /* scaling factor to avoid floating point overflows */

int         dangles = 2;          /* use dangling end energies */

int         tetra_loop = 1;       /* Fold with specially stable 4-loops */

int         noLonelyPairs = 0;    /* avoid helices of length 1 */

int         noGU = 0;             /* GU not allowed at all */

int         no_closingGU = 0;     /* GU allowed only inside stacks */

int         circ = 0;

int         gquad = 0;            /* consider g-qudruplexes in the calculations */

int         canonicalBPonly = 0;  /* remove non-canonical base pairs from structure constraint */

int         uniq_ML   = 0;        /* do ML decomposition uniquely (for subopt) */

int         energy_set = 0;       /* 0 = BP; 1=any with GC; 2=any with AU parameters */

int         do_backtrack = 1;     /* calculate pair prob matrix in part_func() */

char        backtrack_type = 'F'; /* 'C' require (1,N) to be bonded;
                                    'M' seq is part of s multi loop */
char        *nonstandards = (char *)0;  /* contains allowed non standard bases */

int         max_bp_span = -1;     /* base pairs may span the entire sequence length by default */

int         oldAliEn = 0;         /* use old alifold-energies (without removing gaps) */

int         ribo = 0;             /* use ribosum instead of classic covariance term */

double      cv_fact = 1.;

double      nc_fact = 1.;

int         logML     = 0;        /* if nonzero use logarithmic ML energy in energy_of_struct */


/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE     int rtype[8] = {0, 2, 1, 4, 3, 6, 5, 7};


/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC void
set_model_details(model_detailsT *md){

  return vrna_md_set_default(md);
}

PUBLIC void
vrna_md_set_default(model_detailsT *md){

  int i = 0;

  if(md){
    md->dangles           = dangles;
    md->special_hp        = tetra_loop;
    md->noLP              = noLonelyPairs;
    md->noGU              = noGU;
    md->noGUclosure       = no_closingGU;
    md->logML             = logML;
    md->gquad             = gquad;
    md->canonicalBPonly   = canonicalBPonly;
    md->circ              = circ;
    md->uniq_ML           = uniq_ML;
    md->compute_bpp       = do_backtrack;
    md->backtrack         = 1;
    md->backtrack_type    = backtrack_type;
    md->energy_set        = energy_set;
    md->max_bp_span       = max_bp_span;
    md->min_loop_size     = TURN;
    md->oldAliEn          = oldAliEn;
    md->ribo              = ribo;
    md->cv_fact           = cv_fact;
    md->nc_fact           = nc_fact;
    md->temperature       = temperature;
    md->betaScale         = 1.;
    md->pf_scale          = pf_scale;

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
