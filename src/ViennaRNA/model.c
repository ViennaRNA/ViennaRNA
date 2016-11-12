/*
                  Model Details structure creation/modification/destruction

                  This file contains everything which is necessary to
                  obtain, modify, and destroy the model_details datastructure
                  used in the folding recurrences throughout the ViennaRNA
                  Package

                  c Ronny Lorenx

                  Vienna RNA package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

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
# PRIVATE MACROS                #
#################################
*/

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
char    *nonstandards   = NULL;
int     max_bp_span     = VRNA_MODEL_DEFAULT_MAX_BP_SPAN;
int     oldAliEn        = VRNA_MODEL_DEFAULT_ALI_OLD_EN;
int     ribo            = VRNA_MODEL_DEFAULT_ALI_RIBO;
double  cv_fact         = VRNA_MODEL_DEFAULT_ALI_CV_FACT;
double  nc_fact         = VRNA_MODEL_DEFAULT_ALI_NC_FACT;
int     logML           = VRNA_MODEL_DEFAULT_LOG_ML;

/* below are some more deprecated global symbols we need to get rid off */

int         james_rule = 1;       /* interior loops of size 2 get energy 0.8Kcal and
                                    no mismatches (no longer used) */
char        *RibosumFile = NULL;  /* TODO: compile ribosums into program
                                    Warning: this variable will vanish */
int         csv = 0;              /*generate comma seperated output*/
vrna_bp_stack_t       *base_pair = NULL;
FLT_OR_DBL  *pr = NULL;           /* base pairing prob. matrix */
int         *iindx = NULL;        /* pr[i,j] -> pr[iindx[i]-j] */
int         fold_constrained = 0; /* fold with constraints */

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

PRIVATE vrna_md_t defaults = {
  VRNA_MODEL_DEFAULT_TEMPERATURE,
  1.,
  VRNA_MODEL_DEFAULT_DANGLES,
  VRNA_MODEL_DEFAULT_SPECIAL_HP,
  VRNA_MODEL_DEFAULT_NO_LP,
  VRNA_MODEL_DEFAULT_NO_GU,
  VRNA_MODEL_DEFAULT_NO_GU_CLOSURE,
  VRNA_MODEL_DEFAULT_LOG_ML,
  VRNA_MODEL_DEFAULT_CIRC,
  VRNA_MODEL_DEFAULT_GQUAD,
  VRNA_MODEL_DEFAULT_CANONICAL_BP,
  VRNA_MODEL_DEFAULT_UNIQ_ML,
  VRNA_MODEL_DEFAULT_ENERGY_SET,
  VRNA_MODEL_DEFAULT_BACKTRACK,
  VRNA_MODEL_DEFAULT_BACKTRACK_TYPE,
  VRNA_MODEL_DEFAULT_COMPUTE_BPP,
  {0},
  VRNA_MODEL_DEFAULT_MAX_BP_SPAN,
  TURN,
  VRNA_MODEL_DEFAULT_WINDOW_SIZE,
  VRNA_MODEL_DEFAULT_ALI_OLD_EN,
  VRNA_MODEL_DEFAULT_ALI_RIBO,
  VRNA_MODEL_DEFAULT_ALI_CV_FACT,
  VRNA_MODEL_DEFAULT_ALI_NC_FACT,
  1.07,
  {0, 2, 1, 4, 3, 6, 5, 7},
  {0, 1, 2, 3, 4, 3, 2, 0},
  {
    { 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 5, 0, 0, 5},
    { 0, 0, 0, 1, 0, 0, 0, 0},
    { 0, 0, 2, 0, 3, 0, 0, 0},
    { 0, 6, 0, 4, 0, 0, 0, 6},
    { 0, 0, 0, 0, 0, 0, 2, 0},
    { 0, 0, 0, 0, 0, 1, 0, 0},
    { 0, 6, 0, 0, 5, 0, 0, 0}
  }
};

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

/* Fill the base pair type encodings according to the model details */
PRIVATE void fill_pair_matrices(vrna_md_t *md);
PRIVATE void copy_nonstandards(vrna_md_t *md, const char *ns);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC vrna_md_t *
vrna_md_copy( vrna_md_t       *md_to,
              const vrna_md_t *md_from){

  int i;
  vrna_md_t *md;

  md = NULL;

  /* only process if md_from is non-NULL */
  if(md_from){
    if(!md_to){
      /* create container to be filled */
      md = (vrna_md_t *)vrna_alloc(sizeof(vrna_md_t));
    } else {
      /* or directly write to target */
      md = md_to;
    }

    /* check if not the same object */
    if(md_to != md_from){
      /* copy simple members */
      memcpy(md, md_from, sizeof(vrna_md_t));
      /* copy arrays */
      memcpy(md->rtype, &(md_from->rtype[0]), 8 * sizeof(int));
      memcpy(md->alias, &(md_from->alias[0]), (MAXALPHA + 1) * sizeof(short));
      memcpy(md->nonstandards, &(md_from->nonstandards[0]), 64 * sizeof(char));
      /* copy matrices */
      for(i = 0; i <= MAXALPHA; i++){
        memcpy(md->pair[i], (md_from->pair[i]), (MAXALPHA + 1) * sizeof(int));
      }
    }
  }

  return md;
}

PUBLIC void
vrna_md_set_default(vrna_md_t *md){

  if(md) /* copy defaults */
    vrna_md_copy(md, &defaults);
}

PUBLIC char *
vrna_md_option_string(vrna_md_t  *md){
  static char options[255];
  *options = '\0';

  if(md){
    if(md->dangles != VRNA_MODEL_DEFAULT_DANGLES)
      sprintf(options + strlen(options), "-d%d ", md->dangles);
    if(!md->special_hp)
      strcat(options, "-4 ");
    if(md->noLP)
      strcat(options, "--noLP ");
    if(md->noGU)
      strcat(options, "--noGU ");
    if(md->noGUclosure)
      strcat(options, "--noClosingGU ");
    if(md->temperature != VRNA_MODEL_DEFAULT_TEMPERATURE)
      sprintf(options + strlen(options), "-T %f ", md->temperature);
  }

  return options;
}

PRIVATE void
copy_nonstandards(vrna_md_t *md, const char *ns){

  unsigned int n = strlen(ns);
  if(n < 64){
    memcpy(md->nonstandards, ns, strlen(ns)*sizeof(char));
    md->nonstandards[n] = '\0';
  }
}

PUBLIC void
vrna_md_set_nonstandards(vrna_md_t *md, const char *ns_bases){

  const char    *c;
  unsigned int  n;
  int           i, sym;

  if(md){
    if(ns_bases){
      n = strlen(ns_bases);
      if(n < 33){ /* parse the ns_bases list */
        c = ns_bases;
        i = sym = 0;
        if(*c == '-'){
          sym=1;
          c++;
        }

        while(*c != '\0'){
          if(*c != ','){
            md->nonstandards[i++] = *c++;
            md->nonstandards[i++] = *c;
            if((sym) && (*c != *(c-1))){
              md->nonstandards[i++] = *c;
              md->nonstandards[i++] = *(c-1);
            }
          }
          c++;
        }
        md->nonstandards[i] = '\0';

#ifdef  VRNA_BACKWARD_COMPAT
        free(nonstandards);
        nonstandards = vrna_alloc(33);
        memcpy(nonstandards, &(md->nonstandards[0]), 33*sizeof(char));
#endif
      } else {
        vrna_message_warning("vrna_md_set_nonstandards: list too long, dropping nonstandards!");
      }
    } else { /* remove nonstandards */
      md->nonstandards[0] = '\0';
#ifdef  VRNA_BACKWARD_COMPAT
      free(nonstandards);
      nonstandards = NULL;
#endif
    }

    /* update pair/rtype/alias arrays accordingly */
    vrna_md_update(md);
  }
}

PUBLIC void
vrna_md_defaults_reset(vrna_md_t *md_p){

  int i = 0;

  /* first, reset to factory defaults */
  defaults.dangles           = VRNA_MODEL_DEFAULT_DANGLES;
  defaults.special_hp        = VRNA_MODEL_DEFAULT_SPECIAL_HP;
  defaults.noLP              = VRNA_MODEL_DEFAULT_NO_LP;
  defaults.noGU              = VRNA_MODEL_DEFAULT_NO_GU;
  defaults.noGUclosure       = VRNA_MODEL_DEFAULT_NO_GU_CLOSURE;
  defaults.logML             = VRNA_MODEL_DEFAULT_LOG_ML;
  defaults.gquad             = VRNA_MODEL_DEFAULT_GQUAD;
  defaults.canonicalBPonly   = VRNA_MODEL_DEFAULT_CANONICAL_BP;
  defaults.circ              = VRNA_MODEL_DEFAULT_CIRC;
  defaults.uniq_ML           = VRNA_MODEL_DEFAULT_UNIQ_ML;
  defaults.compute_bpp       = VRNA_MODEL_DEFAULT_COMPUTE_BPP;
  defaults.backtrack         = VRNA_MODEL_DEFAULT_BACKTRACK;
  defaults.backtrack_type    = VRNA_MODEL_DEFAULT_BACKTRACK_TYPE;
  defaults.energy_set        = VRNA_MODEL_DEFAULT_ENERGY_SET;
  defaults.max_bp_span       = VRNA_MODEL_DEFAULT_MAX_BP_SPAN;
  defaults.min_loop_size     = TURN;
  defaults.window_size       = VRNA_MODEL_DEFAULT_WINDOW_SIZE;
  defaults.oldAliEn          = VRNA_MODEL_DEFAULT_ALI_OLD_EN;
  defaults.ribo              = VRNA_MODEL_DEFAULT_ALI_RIBO;
  defaults.cv_fact           = VRNA_MODEL_DEFAULT_ALI_CV_FACT;
  defaults.nc_fact           = VRNA_MODEL_DEFAULT_ALI_NC_FACT;
  defaults.temperature       = VRNA_MODEL_DEFAULT_TEMPERATURE;
  defaults.betaScale         = VRNA_MODEL_DEFAULT_BETA_SCALE;
  defaults.sfact             = 1.07;
  defaults.nonstandards[0]   = '\0';

  if(md_p){ /* now try to apply user settings */
    /*
        Note that we use wrapper functions here instead of
        faster direct memory copy because we want to ensure
        that model settings always comply to the constraints
        we set in the wrappers
    */
    vrna_md_defaults_dangles(md_p->dangles);
    vrna_md_defaults_special_hp(md_p->special_hp);
    vrna_md_defaults_noLP(md_p->noLP);
    vrna_md_defaults_noGU(md_p->noGU);
    vrna_md_defaults_noGUclosure(md_p->noGUclosure);
    vrna_md_defaults_logML(md_p->logML);
    vrna_md_defaults_gquad(md_p->gquad);
    defaults.canonicalBPonly = md_p->canonicalBPonly;
    vrna_md_defaults_circ(md_p->circ);
    vrna_md_defaults_uniq_ML(md_p->uniq_ML);
    vrna_md_defaults_compute_bpp(md_p->compute_bpp);
    vrna_md_defaults_backtrack(md_p->backtrack);
    vrna_md_defaults_backtrack_type(md_p->backtrack_type);
    vrna_md_defaults_energy_set(md_p->energy_set);
    vrna_md_defaults_max_bp_span(md_p->max_bp_span);
    vrna_md_defaults_min_loop_size(md_p->min_loop_size);
    vrna_md_defaults_window_size(md_p->window_size);
    vrna_md_defaults_oldAliEn(md_p->oldAliEn);
    vrna_md_defaults_ribo(md_p->ribo);
    vrna_md_defaults_cv_fact(md_p->cv_fact);
    vrna_md_defaults_nc_fact(md_p->nc_fact);
    vrna_md_defaults_temperature(md_p->temperature);
    vrna_md_defaults_betaScale(md_p->betaScale);
    vrna_md_defaults_sfact(md_p->sfact);
    copy_nonstandards(&defaults, &(md_p->nonstandards[0]));
  }

  /* update pair/rtype/alias arrays accordingly */
  vrna_md_update(&defaults);

#ifdef  VRNA_BACKWARD_COMPAT
  temperature     = defaults.temperature;
  pf_scale        = VRNA_MODEL_DEFAULT_PF_SCALE;
  dangles         = defaults.dangles;
  tetra_loop      = defaults.special_hp;
  noLonelyPairs   = defaults.noLP;
  noGU            = defaults.noGU;
  no_closingGU    = defaults.noGUclosure;
  circ            = defaults.circ;
  gquad           = defaults.gquad;
  canonicalBPonly = defaults.canonicalBPonly;
  uniq_ML         = defaults.uniq_ML;
  energy_set      = defaults.energy_set;
  do_backtrack    = defaults.compute_bpp;
  backtrack_type  = defaults.backtrack_type;
  nonstandards    = defaults.nonstandards;
  max_bp_span     = defaults.max_bp_span;
  oldAliEn        = defaults.oldAliEn;
  ribo            = defaults.ribo;
  cv_fact         = defaults.cv_fact;
  nc_fact         = defaults.nc_fact;
  logML           = defaults.logML;
#endif
}

/* below are the setter functions for global default settings */

PUBLIC void
vrna_md_defaults_temperature(double T){

  if(T >= -K0){
    defaults.temperature = T;
#ifdef VRNA_BACKWARD_COMPAT
    temperature = T;
#endif
  } else
    vrna_message_warning("vrna_md_defaults_temperature@model.c: Temperature out of range, T must be above absolute zero. Not changing anything!");
}

PUBLIC double
vrna_md_defaults_temperature_get(void){

  return defaults.temperature;
}

PUBLIC void
vrna_md_defaults_betaScale(double b){

  defaults.betaScale = b;
}

PUBLIC double
vrna_md_defaults_betaScale_get(void){

  return defaults.betaScale;
}

PUBLIC void
vrna_md_defaults_dangles(int d){
  if((d >= 0) && (d <= 3)){
    defaults.dangles = d;
#ifdef VRNA_BACKWARD_COMPAT
    dangles = d;
#endif
  } else
    vrna_message_warning("vrna_md_defaults_dangles@model.c: Dangles out of range, must be (0 <= d <= 3). Not changing anything!");
}

PUBLIC int
vrna_md_defaults_dangles_get(void){

  return defaults.dangles;
}

PUBLIC void
vrna_md_defaults_special_hp(int flag){

  defaults.special_hp = flag ? 1 : 0;
#ifdef VRNA_BACKWARD_COMPAT
  tetra_loop = defaults.special_hp;
#endif
}

PUBLIC int
vrna_md_defaults_special_hp_get(void){

  return defaults.special_hp;
}

PUBLIC void
vrna_md_defaults_noLP(int flag){

  defaults.noLP = flag ? 1 : 0;
#ifdef VRNA_BACKWARD_COMPAT
  noLonelyPairs = defaults.noLP;
#endif
}

PUBLIC int
vrna_md_defaults_noLP_get(void){

  return defaults.noLP;
}

PUBLIC void
vrna_md_defaults_noGU(int flag){

  defaults.noGU = flag ? 1 : 0;
#ifdef VRNA_BACKWARD_COMPAT
  noGU = defaults.noGU;
#endif
  /* update pair/rtype/alias arrays accordingly */
  vrna_md_update(&defaults);
}

PUBLIC int
vrna_md_defaults_noGU_get(void){

  return defaults.noGU;
}

PUBLIC void
vrna_md_defaults_noGUclosure(int flag){

  defaults.noGUclosure = flag ? 1 : 0;
#ifdef VRNA_BACKWARD_COMPAT
  no_closingGU = defaults.noGUclosure;
#endif
}

PUBLIC int
vrna_md_defaults_noGUclosure_get(void){

  return defaults.noGUclosure;
}

PUBLIC void
vrna_md_defaults_logML(int flag){

  defaults.logML = flag ? 1 : 0;
#ifdef VRNA_BACKWARD_COMPAT
  logML = defaults.logML;
#endif
}

PUBLIC int
vrna_md_defaults_logML_get(void){

  return defaults.logML;
}

PUBLIC void
vrna_md_defaults_circ(int flag){

  defaults.circ = flag ? 1 : 0;
#ifdef VRNA_BACKWARD_COMPAT
  circ = defaults.circ;
#endif
}

PUBLIC int
vrna_md_defaults_circ_get(void){

  return defaults.circ;
}

PUBLIC void
vrna_md_defaults_gquad(int flag){

  defaults.gquad = flag ? 1 : 0;
#ifdef VRNA_BACKWARD_COMPAT
  gquad = defaults.gquad;
#endif
}

PUBLIC int
vrna_md_defaults_gquad_get(void){

  return defaults.gquad;
}

PUBLIC void
vrna_md_defaults_uniq_ML(int flag){

  defaults.uniq_ML = flag ? 1 : 0;
#ifdef VRNA_BACKWARD_COMPAT
  uniq_ML = defaults.uniq_ML;
#endif
}

PUBLIC int
vrna_md_defaults_uniq_ML_get(void){

  return defaults.uniq_ML;
}

PUBLIC void
vrna_md_defaults_energy_set(int e){

  if((e >= 0) && (e <= 3)){
    defaults.energy_set = e;
#ifdef VRNA_BACKWARD_COMPAT
    energy_set = e;
#endif
    /* update pair/rtype/alias arrays accordingly */
    vrna_md_update(&defaults);
  } else
    vrna_message_warning("vrna_md_defaults_energy_set@model.c: Energy Set out of range, must be (0 <= e <= 3). Not changing anything!");
}

PUBLIC int
vrna_md_defaults_energy_set_get(void){

  return defaults.energy_set;
}

PUBLIC void
vrna_md_defaults_backtrack(int flag){

  defaults.backtrack = flag ? 1 : 0;
}

PUBLIC int
vrna_md_defaults_backtrack_get(void){

  return defaults.backtrack;
}

PUBLIC void
vrna_md_defaults_backtrack_type(char t){

  switch(t){
    case 'M': /* fall through */
    case 'C': /* fall through */
    case 'F': defaults.backtrack_type = t;
#ifdef VRNA_BACKWARD_COMPAT
              backtrack_type = t;
#endif
              break;
    default:  vrna_message_warning("vrna_md_defaults_backtrack_type@model.c: Backtrack type must be any of 'F', 'C', or 'M'. Not changing anything!");
  }
}

PUBLIC char
vrna_md_defaults_backtrack_type_get(void){

  return defaults.backtrack_type;
}

PUBLIC void
vrna_md_defaults_compute_bpp(int flag){

  if((flag >= 0) && (flag <= 2)){
    defaults.compute_bpp = flag;
#ifdef VRNA_BACKWARD_COMPAT
    do_backtrack = flag;
#endif
  } else
    defaults.compute_bpp = 1;
}

PUBLIC int
vrna_md_defaults_compute_bpp_get(void){

  return defaults.compute_bpp;
}

PUBLIC void
vrna_md_defaults_max_bp_span(int span){

  defaults.max_bp_span = (span <= 0) ? -1 : span;
#ifdef VRNA_BACKWARD_COMPAT
  max_bp_span = defaults.max_bp_span;
#endif
}

PUBLIC int
vrna_md_defaults_max_bp_span_get(void){

  return defaults.max_bp_span;
}

PUBLIC void
vrna_md_defaults_min_loop_size(int size){

  defaults.min_loop_size = (size < 0) ? 0 : size;
}

PUBLIC int
vrna_md_defaults_min_loop_size_get(void){

  return defaults.min_loop_size;
}

PUBLIC void
vrna_md_defaults_window_size(int size){

  defaults.window_size = (size <= 0) ? -1 : size;
}

PUBLIC int
vrna_md_defaults_window_size_get(void){

  return defaults.window_size;
}

PUBLIC void
vrna_md_defaults_oldAliEn(int flag){

  defaults.oldAliEn = flag ? 1 : 0;
#ifdef VRNA_BACKWARD_COMPAT
  oldAliEn = defaults.oldAliEn;
#endif
}

PUBLIC int
vrna_md_defaults_oldAliEn_get(void){

  return defaults.oldAliEn;
}

PUBLIC void
vrna_md_defaults_ribo(int flag){

  defaults.ribo = flag ? 1 : 0;
#ifdef VRNA_BACKWARD_COMPAT
  ribo = defaults.ribo;
#endif
}

PUBLIC int
vrna_md_defaults_ribo_get(void){

  return defaults.ribo;
}

PUBLIC void
vrna_md_defaults_cv_fact(double factor){

  defaults.cv_fact = factor;
#ifdef VRNA_BACKWARD_COMPAT
  cv_fact = factor;
#endif
}

PUBLIC double
vrna_md_defaults_cv_fact_get(void){

  return defaults.cv_fact;
}

PUBLIC void
vrna_md_defaults_nc_fact(double factor){

  defaults.nc_fact = factor;
#ifdef VRNA_BACKWARD_COMPAT
  nc_fact = factor;
#endif
}

PUBLIC double
vrna_md_defaults_nc_fact_get(void){

  return defaults.nc_fact;
}

PUBLIC void
vrna_md_defaults_sfact(double factor){

  defaults.sfact = factor;
}

PUBLIC double
vrna_md_defaults_sfact_get(void){

  return defaults.sfact;
}


PUBLIC void
vrna_md_update(vrna_md_t *md){

  if(md)
    fill_pair_matrices(md);
}


/*
    This function updates the pair/alias/rtype arrays according to model settings.
    It should be called whenever there is a change in the following model settings:
    - energy_set
    - noGU
    - nonstandards
*/
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

                if(md->nonstandards[0] != '\0') {  /* allow nonstandard bp's (encoded by type=7) */
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

    default:    vrna_message_error("Which energy_set are YOU using??");
                break;
  }

  /* set the reverse base pair types */
  for(i = 0; i <= MAXALPHA; i++){
    for(j = 0; j <= MAXALPHA; j++){
      md->rtype[md->pair[i][j]] = md->pair[j][i];
    }
  }

  /* handle special cases separately */
  md->rtype[0] = 0;
  md->rtype[7] = 7;

  /* was used for energy_set == 0
  for(i = 0; i < NBASES; i++)
      for(j = 0; j < NBASES; j++)
       md->rtype[md->pair[i][j]] = md->pair[j][i];
  */
}

#ifdef  VRNA_BACKWARD_COMPAT

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC void
set_model_details(vrna_md_t *md){

  if(md){
    /* make sure there are no uninitialized data fields */
    memset(md, 0, sizeof(vrna_md_t));

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
    md->window_size       = VRNA_MODEL_DEFAULT_WINDOW_SIZE;
    md->oldAliEn          = oldAliEn;
    md->ribo              = ribo;
    md->cv_fact           = cv_fact;
    md->nc_fact           = nc_fact;
    md->temperature       = temperature;
    md->betaScale         = VRNA_MODEL_DEFAULT_BETA_SCALE;
    md->sfact             = 1.07;

    if (nonstandards)
      copy_nonstandards(md, nonstandards);

    /* set default values for the pair/rtype[pair] stuff */
    vrna_md_update(md);

  }
}


PUBLIC char *
option_string(void){

  vrna_md_t md;
  set_model_details(&md);

  return vrna_md_option_string(&md);
}


#endif


