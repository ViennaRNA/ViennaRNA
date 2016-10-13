/*
  Last changed Time-stamp: <2010-06-24 14:50:01 ivo>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: main.c,v 1.5 2008/08/28 09:40:55 ivo Exp $
*/

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>

#if HAVE_LIBRNA_API3
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/fold_vars.h> /* contains global variable cut_point */
#include <ViennaRNA/utils.h>
#include <ViennaRNA/string_utils.h>
#else
#include <fold_vars.h>
#include <fold.h>
#include <utils.h>
#endif

#include "baum.h"
#include "nachbar.h"
#include "cache_util.h"
#include "globals.h"

static char UNUSED rcsid[] ="$Id: main.c,v 1.5 2008/08/28 09:40:55 ivo Exp $";
extern void  read_parameter_file(const char fname[]);
extern void get_from_cache(cache_entry *c);

/* PRIVAT FUNCTIONS */
static void ini_energy_model(void);
static void read_data(void);
static void clean_up(void);

/**/
int main(int argc, char *argv[]) {
  int i;
  char * start, *tmp;
  
  /*
    process command-line optiones
  */
  decode_switches(argc, argv);

  /*
    initialize energy parameters
  */
  ini_energy_model();
  
  /*
    read input file
  */
  read_data();

#if HAVE_LIBRNA_API3
  /* init vrna_fold_compound_t */
  /*
    since read_data() already spliced out potential cut_point, we have to splice it
    in again, to properly init the vrna_fold_compound_t
  */
  tmp     = vrna_cut_point_insert(GAV.farbe, cut_point);
  GAV.vc  = vrna_fold_compound(tmp, &(GAV.md), VRNA_OPTION_DEFAULT);
  free(tmp);
#endif

  /*
    perform GSV.num simulations
  */
    
  start = strdup(GAV.startform); /* remember startform for next run */
  for (i = 0; i < GSV.num; i++) {

    /*
      initialize or reset ringlist to start conditions
    */
    ini_or_reset_rl();
    if (GSV.grow>0) {
      if (strlen(GAV.farbe)>GSV.glen) {
	start[GSV.glen] = '\0';
	GAV.farbe[GSV.glen] = '\0';
	strcpy(GAV.startform,start);
	strcpy(GAV.currform,start);
	GSV.len=GSV.glen;

#if HAVE_LIBRNA_API3
        GAV.vc->length = GSV.len;
#endif
      }
      clean_up_rl();
      ini_or_reset_rl();
    }

    /*
      perform simulation
    */
    for (GSV.steps = 1;; GSV.steps++) {
      cache_entry *c;

      /*
	take neighbourhood of current structure from cache if there
	else generate it from scratch
      */
      if ( (c = lookup_cache(GAV.currform)) ) get_from_cache(c);
      else move_it();
	
      /*
	select a structure from neighbourhood of current structure
	and make it to the new current structure.
	stop simulation if stop condition is met.
      */
      if ( sel_nb() > 0 ) break;

      /* if (GSV.grow>0) grow_chain(); */
    }
  }
  
  /*
    clean up memory
  */
  free(start);
  clean_up();
  return(0);
}

/**/
static void ini_energy_model(void) {

#if HAVE_LIBRNA_API3
  vrna_md_set_default(&(GAV.md));
#else
  set_model_details(&(GAV.md));
#endif

  if ( !GTV.seed ) {
    init_rand();
    GAV.subi[0] = xsubi[0];
    GAV.subi[1] = xsubi[1];
    GAV.subi[2] = xsubi[2];
  }
  else {
    xsubi[0] = GAV.subi[0];
    xsubi[1] = GAV.subi[1];
    xsubi[2] = GAV.subi[2];
  }
  GAV.md.logML        = logML = GTV.logML;
  GAV.md.dangles      = dangles = GTV.dangle;
  GAV.md.temperature  = temperature = GSV.Temp;
  GAV.md.noLP         = GTV.noLP;

  if( GTV.Par ) read_parameter_file(GAV.ParamFile);

#if HAVE_LIBRNA_API2
  GAV.params = get_scaled_parameters(temperature, GAV.md);
#endif
/*  update_fold_params_par(GAV.params); */
}

/**/
static void read_data(void) {
  char *ctmp, *c, **s;
  int i, len;

  /*
    read sequence
  */
  ctmp = get_line(stdin);
  len = strlen(ctmp);
  if (c=strchr(ctmp, '&')) {
    cut_point = (int)(c-ctmp)+1;
    for (; *c; c++) *c = *(c+1); /* splice out the & */
  }
  GAV.farbe = (char *)calloc(len+1, sizeof(char));
  assert(GAV.farbe != NULL);
  sscanf(ctmp, "%s", GAV.farbe);
  GSV.len = strlen(GAV.farbe);
  for (i = 0; i < len; i++) GAV.farbe[i] = toupper(GAV.farbe[i]);
  free (ctmp);
  /* allocate some global arrays */
  GAV.currform = (char *)calloc(GSV.len +1, sizeof(char));
  assert(GAV.currform != NULL);
  GAV.prevform = (char *)calloc(GSV.len +1, sizeof(char));
  assert(GAV.prevform != NULL);
  GAV.startform = (char *)calloc(GSV.len +1, sizeof(char));
  assert(GAV.startform != NULL);

  /*
    read start structure
  */
  if (GTV.start) {
    ctmp = get_line(stdin);
    len = strlen(ctmp);
    sscanf(ctmp, "%s", GAV.startform);
    if (c=strchr(GAV.startform, '&')) {
      for (; *c; c++) *c = *(c+1); /* splice out the & */
    }

    if (strlen(GAV.startform) != GSV.len) {
      fprintf(stderr,
	      "read_data():\n%s\n%s\n unequal length!\n",
	      GAV.farbe, GAV.startform);
      exit(1);
    }
    free (ctmp);
  }
  else { /* start structure is open chain */
    for (i = 0; i< GSV.len; i++) GAV.startform[i] = '.';
  }

  /*
    read stop structure(s)
  */
  if (GTV.stop) {
    s = GAV.stopform;
    while (( ctmp = get_line(stdin))) {
      *s = (char *)calloc(GSV.len+1, sizeof(char));
      sscanf(ctmp, "%s", *s);
      if (c=strchr(ctmp, '&')) {
	for (; *c; c++) *c = *(c+1); /* splice out the & */
      }
 
      if ( (s-GAV.stopform) >= GSV.maxS ) {
	fprintf(stderr,
		"WARNING: Can handle a maximum of %d stop structures\n",
		GSV.maxS );
        break;
      }
      
      if (strlen(*s) != GSV.len) {
	fprintf(stderr, "read_data():\n%s\n%s\n unequal length!\n",
		GAV.farbe, *s);
	exit(1);
      }

      s++;
      free (ctmp);
    }
    GSV.maxS = (s-GAV.stopform);
  }
  else {
    GSV.maxS = 1;
    GAV.stopform[0] = (char *)calloc(GSV.len+1, sizeof(char));
  }

  GAV.farbe_full = strdup(GAV.farbe);

  GAV.sE = (float *)calloc(GSV.maxS, sizeof(float)); 
}

/**/
void clean_up(void) {
  clean_up_globals();
  clean_up_rl();
  clean_up_nbList();
  kill_cache();
}
