/*
  Last changed Time-stamp: <2015-04-13 16:08:40 ivo>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: globals.c,v 1.8 2008/10/07 09:03:14 ivo Exp $
*/

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <getopt.h>

#if HAVE_LIBRNA_API3
#include <ViennaRNA/utils.h>
#include <ViennaRNA/fold_vars.h>
#else
#include <utils.h>
#include <fold_vars.h>
#endif

#include "globals.h"
#include "cmdline.h"

/* forward declarations privat functions */
static void ini_globs(void);
static void ini_gtoggles (void);
static void ini_gvars(void);
static void ini_garrays (void);
static void usage(int status);
static void display_settings(void);
static void display_fileformat(void);
static char *verbose(int optval, const char *which);
static int process_options (int argc, char *argv[]);
static void process_options_gg (int argc, char *argv[]);

static const char *costring(const char *str);

static char UNUSED rcsid[] ="$Id: globals.c,v 1.8 2008/10/07 09:03:14 ivo Exp $";
#define MAXMSG 8
static char msg[MAXMSG][60] =
{{"off"},
 {"on"},
 {"current time"},
 {"take from input file"},
 {"open chain"},
 {"mfe structure of sequence"},
 {"ViennaRNA-Package-1.4 defaults"},
 {""}
};

static struct option const long_options[] =
{ {"dangle",  required_argument, 0, 0},
  {"Temp",    required_argument, 0, 0},
  {"Par",     required_argument, 0, 0},
  {"phi",     required_argument, 0, 0},
  {"pbounds", required_argument, 0, 0},
  {"logML",   no_argument,       &GTV.logML, 0},
  {"noShift", no_argument,       &GTV.noShift, 1},
  {"noLP",    no_argument,       &GTV.noLP, 1},
  {"seed",    required_argument, 0, 0},
  {"time",    required_argument, 0, 0},
  {"num",     required_argument, 0, 0},
  {"start",   no_argument,       &GTV.start, 1},
  {"stop",    no_argument,       &GTV.stop, 1},
  {"fpt",     no_argument,       &GTV.fpt, 0},
  {"met",      no_argument,       &GTV.mc, 1},
  {"grow",    required_argument, 0,  0},
  {"glen",    required_argument, 0,  0},
  {"log",     required_argument, 0,  0},
  {"silent",  no_argument,       0,  0},
  {"lmin",    no_argument,       &GTV.lmin, 1},
  {"cut",     required_argument, 0, 0},
  {"help",    no_argument,       0, 'h'},
  {"verbose", no_argument,       0, 0},
  {NULL, 0, NULL, 0}
};

/**/
void decode_switches(int argc, char *argv[]) {
  ini_globs();
  strcpy(GAV.ProgramName, argv[0]);
  process_options_gg(argc, argv);
}

/**/
void clean_up_globals(void) {
  int i;
  free(GAV.ProgramName);
  free(GAV.BaseName);
  free(GAV.farbe);
  free(GAV.farbe_full);
  free(GAV.startform);
  free(GAV.currform);
  free(GAV.prevform);
  for (i = 0; i < GSV.maxS; i++) free(GAV.stopform[i]);
  free(GAV.stopform);
  free(GAV.sE);
#if HAVE_RNALIB_API3
  vrna_fold_compound_free(GAV.vc);
#else
# if HAVE_RNALIB2
  free(GAV.params);
# endif
#endif
}

/**/
static void usage(int status) {
  fprintf(stderr, "\n%s - Kinetic Folding Program for Nucleic Acids -\n",
	  GAV.ProgramName);
  fprintf(stderr, "Usage: %s [OPTION] < FILE\n", GAV.ProgramName);
  fprintf(stderr,
	  "Options:\n"
	  " EnergyModel\n"
	  "  --dangle <0|1|2>  set dangling end model to (non|normal|double)\n"
	  "  --Temp <float>    set simulation temperature to <float>\n"
	  "  --Par <string>    use energy-parameter-file <string>\n"
	  "  --logML           use linear multiloop-function not logarithmic\n"
	  " MoveSet\n"
	  "  --noShift         turn off shift-moves\n"
	  "  --noLP            forbit structures with isolated base-pairs\n"
	  " Simulation\n"
	  "  --seed <int=int=int>  set random seed to <int=int=int>\n"
	  "  --time <float>        set maxtime of simulation to <float>\n"
	  "  --num <int>           set number of simulations to <int>\n"
	  "  --start               set start structure\n"
	  "  --stop                set stop structure(s)\n"
	  "  --met                 use Metropolis rule not Kawasaki rule\n"
	  "  --fpt                 stop stop structure(s) is reached\n"
	  "  --grow <float>        grow chain every <float> time steps\n"
	  "  --phi <double>        set phi value to <double>\n"
	  "  --pbounds <d1=d2=d3>  set phi_min to d1\n"
	  "                            phi_inc to d2\n"
	  "                            phi_max to d2\n"
	  "                            (d? is a double value)\n"
	  " Output\n"
	  "  --log <string>  set basename of log-file to <string>\n"
	  "  --err <string>  set basename of error-log-file to <string>\n"
	  "  --silent        no output to stdout\n"
	  "  --verbose       more information to stdout\n"
	  "  --lmin          output only local minima to stdout\n"
	  "  --cut <float>   output structures with E <= <float> to stdout\n");
  display_settings();
  display_fileformat();
  exit (status);
}

/**/
static void display_fileformat(void) {
  fprintf(stderr,
	  "Input File Format:\n"
	  "1st line sequence\n"
	  "2nd line start structure (if option --start is used)\n"
	  "following lines stop structures\n\n");
}

/**/
void log_prog_params(FILE *FP) {
    fprintf( FP,
	     "#<\n#Date: %s"
	     "#EnergyModel: dangle=%d Temp=%.1f logML=%s Par=%s\n"
	     "#MoveSet: noShift=%s noLP=%s\n"
	     "#Simulation: num=%d time=%.2f seed=%s fpt=%s mc=%s\n"
	     "#Simulation: phi=%g pbounds=%s\n"
	     "#Output: log=%s silent=%s lmin=%s cut=%.2f\n",
	     time_stamp(),
	     GTV.dangle,
	     GSV.Temp,
	     verbose(GTV.logML, "logML"),
	     GAV.ParamFile,
	     verbose(GTV.noShift, NULL),
	     verbose(GTV.noLP, NULL),
	     GSV.num,
	     GSV.time,
	     verbose(GTV.seed, "seed"),
	     verbose(GTV.fpt, NULL),
	     verbose(GTV.mc, "met"),
	     GSV.phi,
	     verbose(GTV.phi, "pbounds"),
	     GAV.BaseName,
	     verbose(GTV.silent, NULL),
	     verbose(GTV.lmin, NULL),
	     GSV.cut);
    fflush(FP);
}

/**/
void log_start_stop(FILE *FP) {
  int i;
  fprintf(FP, "#%s\n#%s (%6.2f)\n", costring(GAV.farbe),
	  costring(GAV.startform), GSV.startE);
  for (i = 0; i < GSV.maxS; i++) {
    fprintf(FP, "#%s (%6.2f) X%02d\n", costring(GAV.stopform[i]), GAV.sE[i], i+1);
  }
  fprintf(FP, "(%-5hu %5hu %5hu)", GAV.subi[0], GAV.subi[1], GAV.subi[2]);
  costring(NULL);
  fflush(FP);
}

/**/
static void display_settings(void) {
  fprintf(stderr,
	  "Default Settings:\n"
	  " EnergyModel\n"
	  "  --dangle  = %d\n"
	  "  --Temp    = %.2f\n"
	  "  --Par     = %s\n"
	  "  --logML   = %s\n"
	  " MoveSet\n"
	  "  --noShift = %s\n"
	  "  --noLP    = %s\n"
	  " Simulation\n"
	  "  --seed    = %s\n"
	  "  --time    = %.1f\n"
	  "  --phi     = %g\n"
	  "  --pbounds = %s\n"
	  "  --num     = %d\n"
	  "  --start   = %s\n"
	  "  --stop    = %s\n"
	  "  --met     = %s\n"
	  "  --fpt     = %s\n"
	  " Output\n"
	  "  --log     = %s\n"
	  "  --silent  = %s\n"
	  "  --verbose = %s\n"
	  "  --lmin    = %s\n"
	  "  --cut     = %.2f\n",
	  GTV.dangle,
	  GSV.Temp,
	  GAV.ParamFile,
	  verbose(GTV.logML, "logML"),
	  verbose(GTV.noShift, NULL),
	  verbose(GTV.noLP, NULL),
	  verbose(GTV.seed, "seed"),
	  GSV.time,
	  GSV.phi,
	  verbose(GTV.phi, "pbounds"),
	  GSV.num,
	  verbose(GTV.start, "start"),
	  verbose(GTV.stop, "stop"),
	  verbose(GTV.mc, "met"),
	  verbose(GTV.fpt, NULL),
	  GAV.BaseName,
	  verbose(GTV.silent, NULL),
	  verbose(GTV.verbose, NULL),
	  verbose(GTV.lmin, NULL),
	  GSV.cut);
}

/**/
static char *verbose (int optval, const char *which) {
  if (which == NULL) return msg[optval];
  else {
    if ( strcmp(which, "seed") == 0 ) {
      if (optval == 0 ) return "clock";
      else {
	sprintf(msg[MAXMSG - 1],
		"%hu %hu %hu", GAV.subi[0], GAV.subi[1],GAV.subi[2]);
	return msg[MAXMSG - 1];
      }
    }

    if (strcmp(which, "pbounds") == 0) {
      sprintf(msg[MAXMSG - 1],
	      "%g %g %g",
	      GAV.phi_bounds[0], GAV.phi_bounds[1], GAV.phi_bounds[2]);
      return msg[MAXMSG - 1];
    }

    if ( strcmp(which, "met") == 0 ) {
      if ( optval == 0 ) return "Kawasaki";
      else return "Metropolis";
    }

    if ( strcmp(which, "logML") == 0 ) {
      if ( optval == 0 ) return "linear";
      else return "logarithmic";
    }

    if ( strcmp(which, "start") == 0 ) {
      if ( optval == 0 ) return "OpenChain";
      else return "input file";
    }

    if ( strcmp(which, "stop") == 0 ) {
      if ( optval == 0 ) return "Mfe";
      else return "input file";
    }
  }
  return "hae?";
}

/**/
static void ini_globs (void) {
  ini_gtoggles();
  ini_gvars();
  ini_garrays();
}

static void process_options_gg (int argc, char *argv[]) {
  struct gengetopt_args_info args_info;
  if (cmdline_parser (argc, argv, &args_info) != 0)
    exit(1) ;
  if (args_info.help_given) cmdline_parser_print_help();
  if (args_info.full_help_given) cmdline_parser_print_full_help();

  GTV.dangle = args_info.dangle_arg;
  GSV.Temp = args_info.Temp_arg;
  if (args_info.Par_given) 
    strncpy(GAV.ParamFile,args_info.Par_arg,255);
  GTV.Par = args_info.Par_given;
  GTV.logML = args_info.logML_flag;

  GTV.noShift = args_info.noShift_flag;
  GTV.noLP = args_info.noLP_flag;

  GTV.start= args_info.start_flag;
  GTV.stop = args_info.stop_flag;
  GTV.mc   = args_info.met_flag;
  GTV.lmin = args_info.lmin_flag;

  GTV.verbose = args_info.verbose_flag;
  GTV.silent = (GTV.verbose==0)?
    args_info.silent_flag : 0;
  if (args_info.seed_given)
    if (sscanf(args_info.seed_arg, "%hu=%hu=%hu",
	       &GAV.subi[0], &GAV.subi[1], &GAV.subi[2]) != 3)
      usage(EXIT_FAILURE);
    else GTV.seed = 1;

  /* laplace stuff */
  if (args_info.phi_given) {
    if (args_info.phi_arg>0) {
      GTV.phi = 1;
      GSV.phi = args_info.phi_arg;
    }
    else {
      fprintf(stderr, "Value of --phi must be > 0 >%lf<\n", args_info.phi_arg);
      exit(EXIT_FAILURE);
    }
  }
  if (args_info.pbounds_given) {
    if (sscanf(args_info.pbounds_arg, "%g=%g=%g",
	       &GAV.phi_bounds[0],
	       &GAV.phi_bounds[1],
	       &GAV.phi_bounds[2]) == 0)
      usage(EXIT_FAILURE);
    else
      GTV.phi = 1;
    /* check if values are proper */
    if (GAV.phi_bounds[0] > GAV.phi_bounds[2]
	|| GAV.phi_bounds[1] > GAV.phi_bounds[2]) {
      fprintf(stderr,
	      "Unmet requirements for pbounds:\n"
	      "phi_min < phi_max && phi_inc < phi_max\n"
	      "phi_min: %g phi_inc: %g phi_max: %g\n",
	      GAV.phi_bounds[0], GAV.phi_bounds[1], GAV.phi_bounds[2]);
      exit(EXIT_FAILURE);
    }
  }
  GSV.time = args_info.time_arg;
  GSV.num = args_info.num_arg;
  strncpy(GAV.BaseName, args_info.log_arg, 255);
  GSV.cut = args_info.cut_arg;
  GSV.grow = args_info.grow_arg;
  GSV.glen = args_info.glen_arg;
  GTV.lmin = args_info.lmin_flag;
  GTV.fpt  = args_info.fpt_flag;
  cmdline_parser_free(&args_info);
}
/**/
static int process_options (int argc, char *argv[]) {
  int c, itmp;
  float ftmp;
  double dtmp;
  int option_index = 0;
  while ((c = getopt_long (argc, argv, "h",
			     long_options, &option_index)) != EOF) {
      switch (c) {
      case 0:
	if (strcmp(long_options[option_index].name,"dangle")==0) {
	  itmp = -1;
	  if (sscanf(optarg, "%d", &itmp) == 0)
	    usage(EXIT_FAILURE);
	  else if (itmp == 0 || itmp == 1 || itmp == 2 || itmp == 3)
	    GTV.dangle = itmp;
	  else {
	    fprintf(stderr, "Value of --dangle must be 0|1|2 >%d<\n", itmp);
	    usage (EXIT_FAILURE);
	  }
	}

	if (strcmp(long_options[option_index].name,"Temp")==0) {
	  ftmp = -1.0;
	  if (sscanf(optarg, "%f", &ftmp) == 0)
	    usage(EXIT_FAILURE);
	  else if ( ftmp >= 0 )
	    GSV.Temp = ftmp;
	  else {
	    fprintf(stderr, "Value of --Temp must be >= 0 >%.2f<\n", ftmp);
	    usage (EXIT_FAILURE);
	  }
	}

	if (strcmp(long_options[option_index].name,"Par")==0) {
	  if (sscanf(optarg, "%s", GAV.ParamFile) == 0)
	    usage(EXIT_FAILURE);
	  else {
	    GTV.Par = 1;
	  }
	}

	if (strcmp(long_options[option_index].name,"silent")==0) {
	  GTV.silent = 1;
	  GTV.verbose = 0;
	}

	if (strcmp(long_options[option_index].name,"verbose")==0) {
	  if (GTV.silent == 0)
	    GTV.verbose = 1;
	}

	if (strcmp(long_options[option_index].name,"seed")==0) {
	  if (sscanf(optarg, "%hu=%hu=%hu",
		     &GAV.subi[0], &GAV.subi[1], &GAV.subi[2]) == 0)
	    usage(EXIT_FAILURE);
	  else GTV.seed = 1;
	}

	/* laplace stuff */
	if (strcmp(long_options[option_index].name,"pbounds")==0) {
	  if (sscanf(optarg, "%g=%g=%g",
		     &GAV.phi_bounds[0],
		     &GAV.phi_bounds[1],
		     &GAV.phi_bounds[2]) == 0)
	    usage(EXIT_FAILURE);
	  else
	    GTV.phi = 1;
	  /* check if values are proper */
	  if (GAV.phi_bounds[0] > GAV.phi_bounds[2]
	      || GAV.phi_bounds[1] > GAV.phi_bounds[2]) {
	    fprintf(stderr,
		    "Unmet requirements for pbounds:\n"
		    "phi_min < phi_max && phi_inc < phi_max\n"
		    "phi_min: %g phi_inc: %g phi_max: %g\n",
		    GAV.phi_bounds[0], GAV.phi_bounds[1], GAV.phi_bounds[2]);
	    exit(EXIT_FAILURE);
	  }
	}

	if (strcmp(long_options[option_index].name,"time")==0) {
	  dtmp = -1.0;
	  if (sscanf(optarg, "%lf", &dtmp) == 0)
	    usage(EXIT_FAILURE);
	  else if ( dtmp > 0 )
	    GSV.time = dtmp;
	  else {
	    fprintf(stderr, "Value of --time must be > 0 >%lf<\n", dtmp);
	    usage(EXIT_FAILURE);
	  }
	}

	/* laplace stuff */
	if (strcmp(long_options[option_index].name, "phi")==0) {
	  dtmp = -1.0;
	  if (sscanf(optarg, "%lf", &dtmp) == 0)
	    usage(EXIT_FAILURE);
	  else if ( dtmp > 0 ) {
	    GSV.phi = dtmp;
	    GTV.phi = 1;
	  }
	  else {
	    fprintf(stderr, "Value of --phi must be > 0 >%lf<\n", dtmp);
	    exit(EXIT_FAILURE);
	  }
	}

	if (strcmp(long_options[option_index].name,"num")==0) {
	  itmp = -1;
	  if (sscanf(optarg, "%d", &itmp) == 0)
	    usage(EXIT_FAILURE);
	  else if ( itmp > 0 )
	    GSV.num = itmp;
	  else {
	    fprintf(stderr, "Value of --num must be > 0 >%d<\n", itmp);
	    usage(EXIT_FAILURE);
	  }
	}

	if (strcmp(long_options[option_index].name,"log")==0)
	  if (sscanf(optarg, "%s", GAV.BaseName) == 0)
	    usage(EXIT_FAILURE);

	if (strcmp(long_options[option_index].name,"cut")==0)
	  if (sscanf(optarg, "%f", &GSV.cut) == 0)
	    usage(EXIT_FAILURE);

	if (strcmp(long_options[option_index].name,"grow")==0)
	  if (sscanf(optarg, "%lf", &GSV.grow) == 0)
	    usage(EXIT_FAILURE);

	if (strcmp(long_options[option_index].name,"glen")==0)
	  if (sscanf(optarg, "%d", &GSV.glen) == 0)
	    usage(EXIT_FAILURE);

	break;

      case 'h':
	usage (0);
      default:
	usage (EXIT_FAILURE);
      }
  }
  return optind;
}

/**/
static void ini_gtoggles(void) {
  GTV.Par = 0;
  GTV.seed = 0;
  GTV.dangle = 2;
  GTV.logML = 1;
  GTV.noLP = 0;
  GTV.noShift = 0;
  GTV.start = 0;
  GTV.stop = 0;
  GTV.silent = 0;
  GTV.phi = 0;
  GTV.verbose = 0;
  GTV.lmin = 0;
  GTV.fpt = 1;
  GTV.mc = 0;
}

/**/
static void ini_gvars(void) {
  GSV.len = 0;
  GSV.num = 1;
  GSV.maxS = 99;
  GSV.cut = 20;
  GSV.Temp = 37.0;
  GSV.startE = 0.0;
  GSV.stopE = 0.0;
  GSV.currE = 0.0;
  GSV.time = 500.0;
  GSV.phi = 1.0;
  GSV.simTime = 0.0;
  GSV.glen = 15;
}

/**/
static void ini_garrays(void) {
  GAV.ProgramName = (char *)calloc((size_t)256, sizeof(char));
  assert(GAV.ProgramName != NULL);
  GAV.ParamFile   = (char *)calloc((size_t)256, sizeof(char));
  assert(GAV.ParamFile != NULL);
  strcpy(GAV.ParamFile,"VRNA-1.4");
  GAV.BaseName    = (char *)calloc((size_t)256, sizeof(char));
  assert(GAV.BaseName != NULL);
  strcpy(GAV.BaseName, "kinout");
  GAV.stopform = (char **)calloc(GSV.maxS + 1, sizeof(char *));
  assert(GAV.stopform != NULL);
  GAV.farbe = NULL;
  GAV.startform = NULL;
  GAV.currform = NULL;
  GAV.prevform = NULL;
  GAV.phi_bounds[0] = 0.1;
  GAV.phi_bounds[1] = 0.1;
  GAV.phi_bounds[2] = 2.0;
}

static const char *costring(const char *str) {
  static char* buffer=NULL;
  static int size=0;
  int n;
  if ((str==NULL) && (buffer)) {
    /* make it possible to free buffer */
    free(buffer);
    return NULL;
  }
  n=strlen(str);
  if (n>size) {
    size = n+2;
    buffer = realloc(buffer, size);
  }
  if ((cut_point>0)&&(cut_point<=n)) {
    strncpy(buffer, str, cut_point-1);
    buffer[cut_point-1] = '&';
    strncpy(buffer+cut_point, str+cut_point-1, n-cut_point+1);
    buffer[n+1] = '\0';
  } else {
    strncpy(buffer, str, n+1);
  }
  return buffer;
}
/* End of file */
