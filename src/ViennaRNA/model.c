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

#ifdef  VRNA_BACKWARD_COMPAT

/*  below are the evil global variables that will vanish
    as soon as we drop backward compatibility in ViennaRNA
    Package v3
*/

double  temperature     = VRNA_MODEL_DEFAULT_TEMPERATURE;
double  pf_scale        = VRNA_MODEL_DEFAULT_PF_SCALE;
int     dangles         = VRNA_MODEL_DEFAULT_DANGLES;
int     tetra_loop      = VRNA_MODEL_DEFAULT_SPECIAL_HP;
int     noLonelyPairs   = VRNA_MODEL_DEFAULT_NO_LP;
int     noGU            = VRNA_MODEL_DEFAULT_NO_GU;
int     no_closingGU    = VRNA_MODEL_DEFAULT_NO_GU_CLOSURE;
int     circ            = VRNA_MODEL_DEFAULT_CIRC;
int     gquad           = VRNA_MODEL_DEFAULT_GQUAD;
int     canonicalBPonly = VRNA_MODEL_DEFAULT_CANONICAL_BP;
int     uniq_ML         = VRNA_MODEL_DEFAULT_UNIQ_ML;
int     energy_set      = VRNA_MODEL_DEFAULT_ENERGY_SET;
int     do_backtrack    = VRNA_MODEL_DEFAULT_COMPUTE_BPP;
char    backtrack_type  = VRNA_MODEL_DEFAULT_BACKTRACK_TYPE;
char    *nonstandards   = (char *)0;
int     max_bp_span     = VRNA_MODEL_DEFAULT_MAX_BP_SPAN;
int     oldAliEn        = VRNA_MODEL_DEFAULT_ALI_OLD_EN;
int     ribo            = VRNA_MODEL_DEFAULT_ALI_RIBO;
double  cv_fact         = VRNA_MODEL_DEFAULT_ALI_CV_FACT;
double  nc_fact         = VRNA_MODEL_DEFAULT_ALI_NC_FACT;
int     logML           = VRNA_MODEL_DEFAULT_LOG_ML;

#endif

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE int rtype[8] = {0, 2, 1, 4, 3, 6, 5, 7};
PRIVATE int BP_pair[NBASES][NBASES]=
/* _  A  C  G  U  X  K  I */
{{ 0, 0, 0, 0, 0, 0, 0, 0},
 { 0, 0, 0, 0, 5, 0, 0, 5},
 { 0, 0, 0, 1, 0, 0, 0, 0},
 { 0, 0, 2, 0, 3, 0, 0, 0},
 { 0, 6, 0, 4, 0, 0, 0, 6},
 { 0, 0, 0, 0, 0, 0, 2, 0},
 { 0, 0, 0, 0, 0, 1, 0, 0},
 { 0, 6, 0, 0, 5, 0, 0, 0}};


/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

/* Fill the base pair type encodings according to the model details */
PRIVATE void fill_pair_matrices(vrna_md_t *md);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC void
vrna_md_set_default(vrna_md_t *md){

  int i = 0;

  if(md){
    md->dangles           = VRNA_MODEL_DEFAULT_DANGLES;
    md->special_hp        = VRNA_MODEL_DEFAULT_SPECIAL_HP;
    md->noLP              = VRNA_MODEL_DEFAULT_NO_LP;
    md->noGU              = VRNA_MODEL_DEFAULT_NO_GU;
    md->noGUclosure       = VRNA_MODEL_DEFAULT_NO_GU_CLOSURE;
    md->logML             = VRNA_MODEL_DEFAULT_LOG_ML;
    md->gquad             = VRNA_MODEL_DEFAULT_GQUAD;
    md->canonicalBPonly   = VRNA_MODEL_DEFAULT_CANONICAL_BP;
    md->circ              = VRNA_MODEL_DEFAULT_CIRC;
    md->uniq_ML           = VRNA_MODEL_DEFAULT_UNIQ_ML;
    md->compute_bpp       = VRNA_MODEL_DEFAULT_COMPUTE_BPP;
    md->backtrack         = VRNA_MODEL_DEFAULT_BACKTRACK;
    md->backtrack_type    = VRNA_MODEL_DEFAULT_BACKTRACK_TYPE;
    md->energy_set        = VRNA_MODEL_DEFAULT_ENERGY_SET;
    md->max_bp_span       = VRNA_MODEL_DEFAULT_MAX_BP_SPAN;
    md->min_loop_size     = TURN;
    md->oldAliEn          = VRNA_MODEL_DEFAULT_ALI_OLD_EN;
    md->ribo              = VRNA_MODEL_DEFAULT_ALI_RIBO;
    md->cv_fact           = VRNA_MODEL_DEFAULT_ALI_CV_FACT;
    md->nc_fact           = VRNA_MODEL_DEFAULT_ALI_NC_FACT;
    md->temperature       = VRNA_MODEL_DEFAULT_TEMPERATURE;
    md->betaScale         = VRNA_MODEL_DEFAULT_BETA_SCALE;
    md->sfact             = 1.07;
    md->nonstandards[0]   = (char)0;

    /* set default values for the pair/rtype[pair] stuff */
    memcpy(md->rtype, &(rtype[0]), 8 * sizeof(int));
    memset(md->alias, 0, (MAXALPHA + 1) * sizeof(short));
    for(i = 0;i <= MAXALPHA; i++)
      memset(md->pair[i], 0, (MAXALPHA + 1) * sizeof(int));

    vrna_md_update(md);

  }
}

PUBLIC void
vrna_md_set_nonstandards(vrna_md_t *md, const char *ns){

  if(md)
    if(ns){
      unsigned int n = strlen(ns);
      if(n < 33){
        memcpy(md->nonstandards, ns, strlen(ns)*sizeof(char));
        md->nonstandards[n] = '\0';
      } else
        warn_user("vrna_md_set_nonstandards: list too long, dropping nonstandards!");
    }
}

PUBLIC void
vrna_md_set_dangles(vrna_md_t *md, int d){

  if(md)
    if((d >= 0) && (d <= 3))
      md->dangles = d;
}

PUBLIC int
vrna_md_get_dangles(vrna_md_t *md){

  if(md)
    return md->dangles;
  else
    return -1;
}

PUBLIC void
vrna_md_set_temperature(vrna_md_t *md, double T){

  if(md)
    if(T >= -K0)
      md->temperature = T;
}

PUBLIC double
vrna_md_get_temperature(vrna_md_t *md){

  if(md)
    return md->temperature;
  else
    return -K0 - 1.;
}

PUBLIC void
vrna_md_set_special_hp(vrna_md_t *md, int shp){

  if(md)
    md->special_hp = shp;
}

PUBLIC int
vrna_md_get_special_hp(vrna_md_t *md){

  if(md)
    return md->special_hp;
  else
    return -1;
}

PUBLIC void
vrna_md_set_gquad(vrna_md_t *md, int g){

  if(md)
    md->gquad = g;
}

PUBLIC int
vrna_md_get_gquad(vrna_md_t *md){

  if(md)
    return md->gquad;
  else
    return -1;
}

PUBLIC void
vrna_md_set_nolp(vrna_md_t *md, int nolp){

  if(md)
    md->noLP = (nolp) ? 1 : 0;
}

PUBLIC int
vrna_md_get_nolp(vrna_md_t *md){

  if(md)
    return md->noLP;
  else
    return -1;
}

PUBLIC void
vrna_md_set_betascale(vrna_md_t *md, double b){

  if(md)
    md->betaScale = b;
}

PUBLIC double
vrna_md_get_betascale(vrna_md_t *md){
  
  if(md)
    return md->betaScale;
  else
    return -1;
}

PUBLIC void
vrna_md_update(vrna_md_t *md){

  if(md)
    fill_pair_matrices(md);
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC void
vrna_md_set_globals(vrna_md_t *md){

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
    md->backtrack         = VRNA_MODEL_DEFAULT_BACKTRACK;
    md->backtrack_type    = backtrack_type;
    md->energy_set        = energy_set;
    md->max_bp_span       = max_bp_span;
    md->min_loop_size     = TURN;
    md->oldAliEn          = oldAliEn;
    md->ribo              = ribo;
    md->cv_fact           = cv_fact;
    md->nc_fact           = nc_fact;
    md->temperature       = temperature;
    md->betaScale         = VRNA_MODEL_DEFAULT_BETA_SCALE;
    md->sfact             = 1.07;

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

    vrna_md_update(md);

  }
}


PRIVATE void
fill_pair_matrices(vrna_md_t *md){

  int i,j;

  /* nullify everything */
  for(i = 0;i <= MAXALPHA; i++)
    memset(md->pair[i], 0, (MAXALPHA + 1) * sizeof(int));

  memset(md->alias, 0, (MAXALPHA + 1) * sizeof(short));

  /* start setting actual base pair type encodings */
  switch(md->energy_set){
    case  0:    for(i = 0; i < 5; i++)
                  md->alias[i] = (short) i;

                md->alias[5] = 3; /* X <-> G */
                md->alias[6] = 2; /* K <-> C */
                md->alias[7] = 0; /* I <-> default base '@' */

                for(i = 0; i < NBASES; i++)
                    for(j = 0; j < NBASES; j++)
                      md->pair[i][j] = BP_pair[i][j];

                if(md->noGU)
                  md->pair[3][4] = md->pair[4][3] = 0;

                if(md->nonstandards != NULL) {  /* allow nonstandard bp's (encoded by type=7) */
                   for(i = 0; i < (int)strlen(md->nonstandards); i += 2)
                      md->pair[vrna_nucleotide_encode(md->nonstandards[i], md)]
                        [vrna_nucleotide_encode(md->nonstandards[i+1], md)] = 7;
                }

                break;

    case 1:     for(i = 1; i < MAXALPHA;){
                  md->alias[i++] = 3;  /* A <-> G */
                  md->alias[i++] = 2;  /* B <-> C */
                }
                for(i = 1; i < MAXALPHA; i++){
                  md->pair[i][i+1] = 2;    /* AB <-> GC */
                  i++;
                  md->pair[i][i-1] = 1;    /* BA <-> CG */
                }

                break;

    case 2:     for(i = 1; i < MAXALPHA;){
                  md->alias[i++] = 1;  /* A <-> A*/
                  md->alias[i++] = 4;  /* B <-> U */
                }
                for(i = 1; i < MAXALPHA; i++){
                  md->pair[i][i+1] = 5;    /* AB <-> AU */
                  i++;
                  md->pair[i][i-1] = 6;    /* BA <-> UA */
                }

                break;

    case 3:     for(i = 1; i < MAXALPHA - 2; ){
                  md->alias[i++] = 3;  /* A <-> G */
                  md->alias[i++] = 2;  /* B <-> C */
                  md->alias[i++] = 1;  /* C <-> A */
                  md->alias[i++] = 4;  /* D <-> U */
                }
                for(i = 1; i < MAXALPHA - 2; i++){
                  md->pair[i][i+1] = 2;    /* AB <-> GC */
                  i++;
                  md->pair[i][i-1] = 1;    /* BA <-> CG */
                  i++;
                  md->pair[i][i+1] = 5;    /* CD <-> AU */
                  i++;
                  md->pair[i][i-1] = 6;    /* DC <-> UA */
                }

                break;

    default:    nrerror("Which energy_set are YOU using??");
                break;
  }

  /* set the reverse base pair types */
  for(i = 0; i <= MAXALPHA; i++){
    for(j = 0; j <= MAXALPHA; j++){
      md->rtype[md->pair[i][j]] = md->pair[j][i];
    }
  }

  /* was used for energy_set == 0
  for(i = 0; i < NBASES; i++)
      for(j = 0; j < NBASES; j++)
       md->rtype[md->pair[i][j]] = md->pair[j][i];
  */
}
