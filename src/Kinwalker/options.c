/*
  Last changed Time-stamp: <2007-07-10 19:10:56 xtof>
  $Id: options.c,v 1.12 2007/11/03 16:45:58 Kinwalker Exp $
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>
#include "options.h"

static OptionS OptS;





static struct option const long_options[] =
{ 
  {"help",    no_argument,       0, 'h'},
  {"init_structure",   no_argument, 0,   0},
  {"interrupt",   no_argument, 0,   0},
  {"printfront",    no_argument,       0, 0},
  {"test", no_argument,       0, 't'},
  {"verbose", no_argument,       0, 'v'},

  {"barrier_heuristic",   required_argument, 0,   0},
  {"dangle",   required_argument, 0,   0},
  {"grouping",   required_argument, 0,   0},
  {"lookahead",    required_argument, 0,   0},

  {"maxkeep",   required_argument, 0,   0},  
  {"noLonelyPairs",   required_argument, 0,   0},
  {"transcribed",   required_argument, 0,   0},
  {"transcription_rate",   required_argument, 0,   0},
  {"windowsize",   required_argument, 0,   0},  
  {NULL, 0, NULL, 0}
};

/* forward declaration(s) of private function(s) */
static void initializeOptions (char *name);
static void processOptions(int argc, char *argv[]);
static void usage(int status);
static void warn ( char *fmt, ... );

/* ==> public function(s) <== */
/*
 */
OptionS *
decodeCML(int argc, char *argv[]) {
  initializeOptions(argv[0]);
  processOptions(argc, argv);
  return (&OptS);
}

char*
optionString(void)
{
  static char options[100];
  *options = '\0';
  /* //  sprintf(options+strlen(options), "--stack=%d ", OptS.minStackSize); */
/*   //sprintf(options+strlen(options), "--loop=%d ", OptS.minLoopSize); */
/*   //sprintf(options+strlen(options), "--model=%c ", OptS.Emodel); */

  return (options);
}


/* ==> private function(s) <== */
/*
  set options to default value(s)
*/

/*
  print usage information to stderr and exit
*/
static void
usage(int status)
{

  fprintf(stderr, "\nUsage: %s [OPTIONS] < SeqFile  > Outfile\n", OptS.ExeName);
  fprintf(stderr,
          "Options without argument:\n"
 	  " -h, --help          Print usage information for %s.\n"
          "--init_structure     Start with a structure other than the open chain.\n"
   	  "--interrupt          Allow interrupted folding trajectories when the barrier is exceeded.\n"
          "--printfront.        Creates PS plots of front progression with index i, named front_trajectory($i).ps.\n"
     	  " -t, --test          Use test sequence.\n"
	  " -v, --verbose       Verbose mode. Print debugging messages about program progress.\n\n",
	  OptS.ExeName);
  fprintf(stderr,
          "Options with argument:\n"
          "--barrier_heuristic  'M' Morgan-Higgs,'S' limits small stacks,'B' Barriers,'A' all, then take minimum. Default: >%c<\n"
          "--dangle             Dangle value of 0,1,2 as in the ViennaRNA package. Default: >%d<\n"
	  "--grouping           How barrier_heuristic 'M' treats conflict groups(\"standard\" or \"regroup\"). Default: >%s<\n"
	  "--lookahead          #BP that MorganHiggs forms its subpaths from. Default: >%d<\n",
	  OptS.barrier_heuristic,OptS.dangle, OptS.grouping, OptS.lookahead);
  fprintf(stderr,
          "--maxkeep            Breadth of breadth first seerch in barrier_heuristic='B'. Default: >%d<\n"
          "--nolonely           Value of noLonelyPairs as in ViennaRNA. Default: >%d<\n"
          "--transcribed        #bases initially transcribed, <0 means all is transcribed. Default: >%d<)\n"
          "--transcription_rate #bases transcribed per second. Default: >%f<)\n"
          "--windowsize         Max size of substructures considered for folding events during transcription, 0= all are considered. Default: >%d<)\n",
          OptS.maxkeep,OptS.noLonelyPairs,OptS.transcribed, OptS.transcription_rate,OptS.windowsize);
  exit (status);
}

/*
  call example:
  warn("calloc(): out of memory when requesting %d bytes", numbytes)
 */
static void
warn ( char *fmt, ... )
{
  va_list args;
  
  va_start( args, fmt );
  fprintf( stderr, "\nWARNING: ");
  vfprintf( stderr, fmt, args);
  fprintf( stderr,"\n");
  va_end( args );
}

/**
 * Set options to default value(s)
 */
static void
initializeOptions (char *name) {

  strcpy(OptS.ExeName, name);
  /*without argument*/
  OptS.init_structure=0;
  OptS.interrupt_trajectory=0;
  OptS.printfront=0;
  OptS.testseq      = 0;
  OptS.verbose      = 0;

  /* with argument */
  OptS.barrier_heuristic= 'M';
  OptS.dangle=0;
  strcpy(OptS.grouping,"standard");
  OptS.lookahead=1;

  OptS.maxkeep = 1;
  OptS.noLonelyPairs=2;
  OptS.transcribed=1;
  OptS.transcription_rate=200.0;
  OptS.windowsize=0;
}

/**
 * Process command-line options and update OptS structure
 */
static void
processOptions(int argc, char *argv[]) {
  int c, option_index = 0;
  while ((c = getopt_long (argc, argv, "htv", 
			   long_options, &option_index)) != EOF) {
    switch (c) {
    case 0:    
      /*options without argument*/
     if ( strcmp(long_options[option_index].name, "init_structure" ) == 0 ) {
         OptS.init_structure = 1; 
     }
     if ( strcmp(long_options[option_index].name, "interrupt" ) == 0 ) {
         OptS.interrupt_trajectory = 1; 
     }
     if ( strcmp(long_options[option_index].name, "printfront") == 0 ) {
	OptS.printfront = 1;
     }
     if ( strcmp(long_options[option_index].name, "test") == 0 ) {
	OptS.testseq = 1;
     }
     if ( strcmp(long_options[option_index].name, "verbose") == 0 ) {
	OptS.verbose = 1;
     }


      /*options with argument*/
     if ( strcmp(long_options[option_index].name, "barrier_heuristic" ) == 0 ) {
	char tmp;
	if ( sscanf(optarg, "%c", &tmp) == 0 ) {
	  warn ("Option requires argument: No barrier_heuristic specified");
	  usage (EXIT_FAILURE);
	}
	else { OptS.barrier_heuristic = tmp; }
      } 
      if ( strcmp(long_options[option_index].name, "dangle") == 0 ) {
	int tmp;
	if ( sscanf(optarg, "%d", &tmp ) == 0 ) {
	  warn ("Option requires argument: set dangle to 0, 1 or 2");
	  usage (EXIT_FAILURE);
	}
	else { OptS.dangle = tmp; }
      }
      if (strcmp(long_options[option_index].name, "grouping")==0) {
	char tmp[256];
	if (sscanf(optarg, "%s", &tmp[0]) == 0) {
	  warn ("Option requires argument: No grouping specified");
	  usage (EXIT_FAILURE);
	}
	else { strcpy(OptS.grouping,tmp); }
      }
      if ( strcmp(long_options[option_index].name, "lookahead" ) == 0 ) {
	int tmp;
	if ( sscanf(optarg, "%d", &tmp) == 0 ) {
	  warn ("Option requires argument: No lookahead values specified");
	  usage (EXIT_FAILURE);
	}
	else { OptS.lookahead = tmp; }
       }
      /* set breadth of breadth first seaerch in findpath */
      if ( strcmp(long_options[option_index].name, "maxkeep") == 0 ) {
	int tmp;
	if ( sscanf(optarg, "%d", &tmp ) == 0 ) {
	  warn ("Option requires argument: INT > 0");
	  usage (EXIT_FAILURE);
	}
	else { OptS.maxkeep = tmp; }
      }
      if ( strcmp(long_options[option_index].name, "noLonelyPairs" ) == 0 ) {
	int tmp;
	if ( sscanf(optarg, "%d", &tmp) == 0 ) {
	  warn ("Option requires argument: noLonelyPairs 0 or 1");
	  usage (EXIT_FAILURE);
	}
	else { OptS.noLonelyPairs = tmp; }
      }
      if ( strcmp(long_options[option_index].name, "transcribed") == 0 ) {
	int tmp;
	if ( sscanf(optarg, "%d", &tmp ) == 0 ) {
	  warn ("Option requires argument: No value specified for transcribed");
	  usage (EXIT_FAILURE);
	}
	else { OptS.transcribed = tmp; }
      }
      if ( strcmp(long_options[option_index].name, "transcription_rate") == 0 ) {
	float tmp;
	if ( sscanf(optarg, "%f", &tmp ) == 0 ) {
	  warn ("Option requires argument: No value specified for transcription_rate");
	  usage (EXIT_FAILURE);
	}
	else { OptS.transcription_rate = tmp; }
      }   
      if ( strcmp(long_options[option_index].name, "windowsize" ) == 0 ) {
	int tmp;
	if ( sscanf(optarg, "%d", &tmp) == 0 ) {
          fprintf(stderr,"Option requires argument: No windowsize specified");
	  usage (EXIT_FAILURE);
	}
	else { OptS.windowsize = tmp;}
      }   
      break;
    case 'h':
      usage (EXIT_SUCCESS);
      break;
    case 'v':
      OptS.verbose++;
      break;
    case 't':
      OptS.testseq = 1;
      break;
    default:
      usage (EXIT_FAILURE);
    }
  }

#if 0
  /* parse the to file names from the remaining command-line */
  if ( (optind + 2) <= argc) {
    strcpy(OptS.PhylipFile,  argv[optind++]);
    strcpy(OptS.DumpFile,    argv[optind++]);
  }
  else {
    warn ("Not enough Input Files. (Need a Phylip AND a Dump File)");
    usage(EXIT_FAILURE);
  }
#endif
}

/* End of file */
