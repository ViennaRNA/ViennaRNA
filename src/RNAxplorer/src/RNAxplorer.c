#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <regex.h>

#include <ViennaRNA/model.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/io/file_formats.h>
#include <ViennaRNA/datastructures/hash_tables.h>
#include <ViennaRNA/landscape/move.h>
#include <ViennaRNA/landscape/walk.h>
#include <ViennaRNA/landscape/neighbor.h>
#include <ViennaRNA/boltzmann_sampling.h>
#include <ViennaRNA/part_func.h>

#include "RNAwalk.h"
#include "meshpoint.h"
#include "barrier_lower_bound.h"
#include "distorted_sampling.h"
#include "distorted_samplingMD.h"
#include "repellant_sampling.h"
#include "PathFinder.h"
#include "gradient_walker.h"

#include "RNAxplorer_cmdl.h"

enum strategies_avail_e {
  MOVE_GRADIENT_WALK,
  PATHFINDER_SADDLE_GRADIENT_WALK,  /* perform gradient walks from saddle point to find meshpoint(s) */
  PATHFINDER_SADDLE_MONTE_CARLO,    /* perform rejection-less monte carlo (Gillespie) simulation away from saddle point */
  PATHFINDER_SADDLE_MONTE_CARLO_SA, /* perform rejection-less monte carlo (Gillespie) simulation away from saddle point while cooling down the system (simulated annealing) */
  PATHFINDER_TWO_D_REPRESENTATIVES, /* use 2D representatives as meshpoints */
  FINDPATH,                         /* use findpath heuristic to obtain optimal refolding path */
  TWO_D_LOWER_BOUND,                /* compute optimal folding path based on 2D representatives using Dijkstra algorithm on andscape projection */
  REPELLENT_SAMPLING,
  ATTRACTION_SAMPLING,
  TEMPERATURE_SCALING_SAMPLING,
  REPELLENT_SAMPLING_HEURISTIC,
  RETRIEVE_LOCAL_MINIMA             /* perform gradient walks for and arbitrary set of structures */
};

struct options_s {
  enum strategies_avail_e strategy;
  vrna_md_t               md;

  /* findpath option(s) */
  int                     max_keep;

  /* common options */
  int                     iterations;
  int                     samples;

  /* PathFinder options */
  int                     max_storage;

  /* 2D fold options */
  int                     max_d1;
  int                     max_d2;

  /* simulated annealing options */
  int                     simulated_annealing;
  float                   t_start;
  float                   t_end;
  float                   cooling_rate;

  rnax_path_finder_opt_t  pathfinder;

  /* retrieve local minima options */
  double temperature_celsius;
  int shift_moves;
  char *parameter_file;

  /* repulsive sampling options*/
  char *sequence;
  int penalize_structures;
  char *struc1;
  char *struc2;
  int granularity;
  int num_samples;
  float exploration_factor;
  float min_exploration_percent;
  int cluster; //flag
  char *lmin_file;
  char *TwoD_file;
  int non_red; //flag
  char *non_red_file;
  int explore_two_neighborhood; //flag
  int post_filter_two; //flag
  int ediff_penalty; //flag
  float mu;
  int verbose;
};

typedef int (xplorer_func)(const char       *rec_id,
                           const char       *orig_sequence,
                           char             **structures,
                           struct options_s *opt);

typedef struct {
  enum strategies_avail_e strategy;
  xplorer_func            *f;
  const char              *name;
} strategies;


struct options_s *
process_arguments(int   argc,
                  char  *argv[]);


char **
extract_structures(unsigned int n,
                   int          maybe_multiline,
                   char         **rec_rest);


int
moves_gradient_descent(const char       *rec_id,
                       const char       *orig_sequence,
                       char             **structures,
                       struct options_s *opt);


int
paths_findpath(const char       *rec_id,
               const char       *orig_sequence,
               char             **structures,
               struct options_s *opt);


int
paths_pathfinder_gd(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt);


int
paths_pathfinder_mc(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt);


int
paths_pathfinder_mcsa(const char        *rec_id,
                      const char        *orig_sequence,
                      char              **structures,
                      struct options_s  *opt);


int
paths_pathfinder_db(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt);


int
paths_pathfinder_dbba(const char        *rec_id,
                      const char        *orig_sequence,
                      char              **structures,
                      struct options_s  *opt);


int
sampling_repulsion(const char       *rec_id,
                   const char       *orig_sequence,
                   char             **structures,
                   struct options_s *opt);


int
sampling_attraction(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt);


int
sampling_temperature(const char       *rec_id,
                     const char       *orig_sequence,
                     char             **structures,
                     struct options_s *opt);

int
sampling_repellent_heuristic(const char       *rec_id,
                   const char       *orig_sequence,
                   char             **structures,
                   struct options_s *opt);

int
retrieve_local_minima(const char       *rec_id,
                      const char       *orig_sequence,
                      char             **structures,
                      struct options_s *opt);


/**
*** \file RNAxplorer.c
**/


#define NUM_STRATEGIES    12

static strategies known_strategies[NUM_STRATEGIES] = {
  /* code, function, name */
  {
    MOVE_GRADIENT_WALK,
    &moves_gradient_descent,
    "Gradient Descent Moves"
  },
  {
    PATHFINDER_SADDLE_GRADIENT_WALK,
    &paths_pathfinder_gd,
    "PathFinder - Gradient Descent from Saddle"
  },
  {
    PATHFINDER_SADDLE_MONTE_CARLO,
    &paths_pathfinder_mc,
    "PathFinder - MCMC from Saddle"
  },
  {
    PATHFINDER_SADDLE_MONTE_CARLO_SA,
    &paths_pathfinder_mcsa,
    "PathFinder - MCMC from Saddle (Simulated Annealing)"
  },
  {
    PATHFINDER_TWO_D_REPRESENTATIVES,
    &paths_pathfinder_db,
    "PathFinder - 2D Representatives"
  },
  {
    FINDPATH,
    &paths_findpath,
    "Findpath"
  },
  {
    TWO_D_LOWER_BOUND,
    &paths_pathfinder_dbba,
    "Energy Barrier Lower Bound from 2D Representation"
  },
  {
    REPELLENT_SAMPLING,
    &sampling_repulsion,
    "Repellent Sampling Scheme"
  },
  {
    ATTRACTION_SAMPLING,
    &sampling_attraction,
    "Directed Sampling Scheme"
  },
  {
    TEMPERATURE_SCALING_SAMPLING,
    &sampling_temperature,
    "Temperature Scaling Sampling Scheme"
  },
  {
    REPELLENT_SAMPLING_HEURISTIC,
    &sampling_repellent_heuristic,
    "Repellent Sampling Heuristic"
  },
  {
      RETRIEVE_LOCAL_MINIMA,
      &retrieve_local_minima,
      "Retrieve local minima for an arbitrary set of secondary structures"
  }
};

static char       *extended_options = NULL;

//percentage of distortion per reference.
size_t            length_indicesAndPercentages  = 0;
double            *indicesAndPercentages        = NULL;


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
int
main(int  argc,
     char *argv[])
{
  FILE              *input_stream     = stdin;
  xplorer_func      *processing_func  = NULL;

  struct options_s  *options = process_arguments(argc, argv);

  for (int i = 0; i < NUM_STRATEGIES; i++)
    if (options->strategy == known_strategies[i].strategy) {
      processing_func = known_strategies[i].f;
      break;
    }

  // read records from stdin only if not sufficient parameteres are given for the repellant sampling heuristic!
  if(processing_func == sampling_repellent_heuristic && options->sequence != NULL && strcmp(options->sequence,"") != 0){
    processing_func("", options->sequence, NULL, options);

    /* free options */
    free(options->TwoD_file);
    free(options->non_red_file);
    free(options->lmin_file);
    free(options->parameter_file);
    free(options->sequence);
    free(options->struc1);
    free(options->struc2);
    free(options);
    exit(EXIT_SUCCESS);
  }


  int           istty_in  = isatty(fileno(input_stream));

  unsigned int  read_opt = 0;

  if (istty_in)
    vrna_message_input_seq("Input sequence (upper or lower case) followed by structures");

  /* set options we wanna pass to vrna_file_fasta_read_record() */
  if (istty_in)
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;

  /* initialize random number generator from libRNA */
  vrna_init_rand();


  /* main loop that processes each record obtained from input stream */
  int read_fasta_records = 0;
  do {
    char          *rec_sequence, *rec_id, **rec_rest;
    unsigned int  rec_type;
    int           maybe_multiline;

    rec_id          = NULL;
    rec_rest        = NULL;
    maybe_multiline = 0;

    rec_type = vrna_file_fasta_read_record(&rec_id,
                                           &rec_sequence,
                                           &rec_rest,
                                           input_stream,
                                           read_opt);

    if (rec_type & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT)){
        if(read_fasta_records == 0){
            vrna_message_input_seq("Input at least one valid fasta record! (a header line '> header' and a sequence line (upper or lower case) followed by structures)");
        }
        free(rec_sequence);
        free(rec_id);
        int i;
        for(i=0; rec_rest[i]; i++)
          free(rec_rest[i]);
        free(rec_rest);
        break;
    }

    /*
     ########################################################
     # init everything according to the data we've read
     ########################################################
     */
    if (rec_id) {
      maybe_multiline = 1;
      /* remove '>' from FASTA header */
      rec_id = memmove(rec_id, rec_id + 1, strlen(rec_id));
    }

    unsigned int  n = strlen(rec_sequence);

    char          **structures = extract_structures(n, maybe_multiline, rec_rest);

    processing_func(rec_id, rec_sequence, structures, options);
    read_fasta_records++;

    free(rec_sequence);
    free(rec_id);

    /* free the rest of current dataset */
    if (structures) {
      for (int i = 0; structures[i]; i++)
        free(structures[i]);
      free(structures);
    }
    else{
      if(rec_rest){
        int i;
        for(i=0; rec_rest[i]; i++)
          free(rec_rest[i]);
        free(rec_rest);
      }
    }
    if (istty_in)
      vrna_message_input_seq("Input sequence (upper or lower case) followed by structures");
  } while (1);

  /* free options */
  free(options->TwoD_file);
  free(options->non_red_file);
  free(options->lmin_file);
  free(options->parameter_file);
  free(options->sequence);
  free(options->struc1);
  free(options->struc2);
  free(options);

  return EXIT_SUCCESS;
}


struct options_s *
default_options(void)
{
  struct options_s *options = (struct options_s *)vrna_alloc(sizeof(struct options_s));

  /* default strategy */
  options->strategy = REPELLENT_SAMPLING_HEURISTIC;

  /* default energy model settings */
  vrna_md_set_default(&(options->md));
  options->md.uniq_ML = 1; /* we certainly require unique multibranch loop decomposition in any case */

  /* findpath option(s) */
  options->max_keep = 10;

  /* common options to many methods */
  options->samples    = 1000;
  options->iterations = 1;

  /* PathFinder options */
  options->max_storage = 10;

  /* 2D fold options */
  options->max_d1 = 5;
  options->max_d2 = 5;

  /* simulated annealing options */
  options->simulated_annealing  = 0;
  options->t_start              = 37.0 + K0;
  options->t_end                = 0. + K0;
  options->cooling_rate         = 0.9998;

  return options;
}


struct options_s *
process_arguments(int   argc,
                  char  *argv[])
{
  struct RNAxplorer_args_info args_info;

  struct options_s            *options = default_options();

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNAxplorer_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* temperature */
  options->temperature_celsius = 37;
  if (args_info.temp_given){
    temperature = args_info.temp_arg;
    options->temperature_celsius = args_info.temp_arg;
  }
  options->shift_moves = 0;
  if (args_info.shift_moves_flag){
      options->shift_moves = 1;
  }
  options->parameter_file = NULL;
  if (args_info.parameter_file_given){
      options->parameter_file = args_info.parameter_file_arg;
  }

  /* method */
  if (args_info.method_given) {
    char *m = args_info.method_arg;
    if (!strcmp(m, "MC")) {
      options->strategy = PATHFINDER_SADDLE_MONTE_CARLO;
    } else if (!strcmp(m, "MC-SA")) {
      options->strategy   = PATHFINDER_SADDLE_MONTE_CARLO_SA;
      simulatedAnnealing  = 1;
    } else if (!strcmp(m, "GW")) {
      options->strategy = PATHFINDER_SADDLE_GRADIENT_WALK;
    } else if (!strcmp(m, "DB-MFE")) {
      options->strategy = PATHFINDER_TWO_D_REPRESENTATIVES;
    } else if (!strcmp(m, "BLUBB")) {
      options->strategy = TWO_D_LOWER_BOUND;
    } else if (!strcmp(m, "SM")) {
      options->strategy = ATTRACTION_SAMPLING;
    } else if (!strcmp(m, "RS")) {
      /* Repellant Sampling */
      options->strategy = REPELLENT_SAMPLING;
    } else if (!strcmp(m, "RSH")) {
      /* Repellant Sampling Heuristic*/
      options->strategy = REPELLENT_SAMPLING_HEURISTIC;
    } else if (!strcmp(m, "RL")) {
      options->strategy = RETRIEVE_LOCAL_MINIMA;
    }
  }

  /* maximum number of simulations / iterations */
  if (args_info.extended_opt_given)
    extended_options = strdup(args_info.extended_opt_arg);

  /* maximum number of simulations / iterations */
  if (args_info.iterations_given)
    options->iterations = args_info.iterations_arg;

  /* maxkeep for Flamm et al. direct path heuristics */
  if (args_info.maxKeep_given)
    options->max_keep = args_info.maxKeep_arg;

  /* Amount of best solutions to store per iteration */
  if (args_info.maxStore_given)
    options->max_storage = args_info.maxStore_arg;

  if (args_info.circ_given)
    options->md.circ = 1;

  if (args_info.cooling_rate_given)
    treduction = args_info.cooling_rate_arg;

  if (args_info.tstart_given)
    tstart = args_info.tstart_arg + K0;

  if (args_info.tstop_given)
    tstop = args_info.tstop_arg + K0;

  if (args_info.penalizeBackWalks_given)
    backWalkPenalty = 1;

  if (args_info.basinStructure_given)
    options->strategy = MOVE_GRADIENT_WALK;

  if (args_info.maxD_given)
    options->max_d1 = options->max_d2 = args_info.maxD_arg;

  if (args_info.maxD1_given)
    options->max_d1 = args_info.maxD1_arg;

  if (args_info.maxD2_given)
    options->max_d2 = args_info.maxD2_arg;

  if (args_info.betaScale_given)
    options->md.betaScale = args_info.betaScale_arg;

  if (args_info.p0_given) {
    int     i, j = 1, lmintmp;
    double  poptmp = 0.;

    length_indicesAndPercentages  = args_info.p0_given;
    indicesAndPercentages         = (double *)calloc(2 * args_info.p0_given + 1, sizeof(double));
    *indicesAndPercentages        = 1;
    for (i = 0; i < args_info.p0_given; i++) {
      if (sscanf(args_info.p0_arg[i], "%d=%lg", &lmintmp, &poptmp) == 0)
        exit(EXIT_FAILURE);

      if (lmintmp < 1) {
        fprintf(stderr, "States in --p0 must be >=1\n");
        exit(EXIT_FAILURE);
      } else {
        *(indicesAndPercentages + j)      = (double)lmintmp;
        *(indicesAndPercentages + j + 1)  = poptmp;
        *indicesAndPercentages            += 2;
        j                                 += 2;
      }
    }
  }

  /* repulsive sampling args */
  options->penalize_structures = 0;
  options->granularity = 1;
  options->num_samples = 1;
  options->exploration_factor = 1;
  options->cluster = 0;
  options->TwoD_file = NULL;
  options->ediff_penalty = 0;
  options->verbose = 0;
  options->lmin_file = NULL;
  options->non_red_file = NULL;
  options->non_red = 0;
  options->mu = 0.1;
  options->post_filter_two = 0;

  if(args_info.sequence_given){
      options->sequence = vrna_alloc(strlen(args_info.sequence_arg)+1);
      strcpy(options->sequence, args_info.sequence_arg);
  }
  options->penalize_structures = args_info.penalize_structures_flag;
  if(args_info.struc1_given){
      options->struc1 = vrna_alloc(strlen(args_info.struc1_arg)+1);
      strcpy(options->struc1, args_info.struc1_arg);
  }
  if(args_info.struc2_given){
      options->struc2 = vrna_alloc(strlen(args_info.struc2_arg)+1);
      strcpy(options->struc2, args_info.struc2_arg);
  }

  options->granularity = args_info.granularity_arg;
  options->num_samples = args_info.num_samples_arg;

  if(args_info.exploration_factor_given){
      options->exploration_factor = args_info.exploration_factor_arg;
  }
  if(args_info.min_exploration_percent_given){
      options->min_exploration_percent = args_info.min_exploration_percent_arg;
  }
  options->cluster = args_info.cluster_flag;
  if(args_info.lmin_file_given){
      options->lmin_file = vrna_alloc(strlen(args_info.lmin_file_arg)+1);
      strcpy(options->lmin_file, args_info.lmin_file_arg);
  }
  if(args_info.TwoD_file_given){
      options->TwoD_file = vrna_alloc(strlen(args_info.TwoD_file_arg)+1);
      strcpy(options->TwoD_file, args_info.TwoD_file_arg);
  }
  options->non_red = args_info.nonred_flag;
  if(args_info.nonred_file_given){
      options->non_red_file = vrna_alloc(strlen(args_info.nonred_file_arg)+1);
      strcpy(options->non_red_file, args_info.nonred_file_arg);
  }
  options->explore_two_neighborhood = args_info.explore_two_neighborhood_flag;
  options->post_filter_two = args_info.post_filter_two_flag;
  options->ediff_penalty = args_info.ediff_penalty_flag;
  options->mu = args_info.mu_arg;
  options->verbose = args_info.verbose_flag;
  /* free allocated memory of command line data structure */
  RNAxplorer_cmdline_parser_free(&args_info);

  return options;
}


char **
extract_structures(unsigned int n,
                   int          maybe_multiline,
                   char         **rec_rest)
{
  char    **structures = NULL;
  size_t  num_structures = 0;
  size_t  size_structures = 10;
  size_t  l, l_prev;
  int     i, read_on;

  if ((rec_rest) && (rec_rest[0])) {
    structures = (char **)vrna_alloc(sizeof(char *) * size_structures);

    read_on = 0;
    l_prev  = 0;

    for (i = 0; rec_rest[i]; i++) {
      switch (rec_rest[i][0]) {
        case '\0':  /* fall-through */
        case '#':   /* fall-through */
        case '%':   /* fall-through */
        case ';':   /* fall-through */
        case '/':   /* fall-through */
        case '*':   /* fall-through */
        case ' ':   /* fall-through */
        case '\t':
          break;

        case '(':
        case ')':
        case '.':
        case '+':
          l = strlen(rec_rest[i]);
          if (l + l_prev < n) {
            if (maybe_multiline) {
              read_on = 1;
            } else {
              vrna_message_error(
                "sequence and structure (at line %d) have unequal lengths (%u vs. %u)",
                i,
                n,
                l + l_prev);
            }
          } else if (l + l_prev > n) {
            vrna_message_error(
              "sequence and structure (at line %d) have unequal lengths (%u vs. %u)",
              i,
              n,
              l + l_prev);
          }

          if (l_prev == 0)
            structures[num_structures] = (char *)vrna_alloc(sizeof(char) * (n + 1));

          memcpy(structures[num_structures] + l_prev, &(rec_rest[i][0]), sizeof(char) * l);
          structures[num_structures][l_prev + l] = '\0';

          if (l + l_prev == n) {
            num_structures++;
            l_prev  = 0;
            read_on = 0;
            if (num_structures == size_structures - 1) {
              size_structures *= 1.4;
              structures      = (char **)vrna_realloc(structures, sizeof(char *) * size_structures);
            }
          } else if (read_on) {
            l_prev = l;
          }

          break;
      }

      free(rec_rest[i]);
    }

    free(rec_rest);

    if (num_structures > 0) {
      structures =
        (char **)vrna_realloc(structures, sizeof(char *) * (num_structures + 1));
      structures[num_structures] = NULL;
    } else {
      free(structures);
      structures = NULL;
    }
  }

  return structures;
}


/*
 *  #############################################################
 *  # Below are the wrappers for the different operation modes  #
 *  #############################################################
 */

/*
 *  *******************
 *  * 1. Move Modes   *
 *  *******************
 */
int
moves_gradient_descent(const char       *rec_id,
                       const char       *orig_sequence,
                       char             **structures,
                       struct options_s *opt)
{
  char *rec_sequence = strdup(orig_sequence);

  vrna_seq_toRNA(rec_sequence);
  vrna_seq_toupper(rec_sequence);

  vrna_fold_compound_t *fc = vrna_fold_compound(rec_sequence,
                                                &(opt->md),
                                                VRNA_OPTION_EVAL_ONLY);

  for (int i = 0; structures[i]; i++) {
    short *pt = vrna_ptable(structures[i]);

    (void)vrna_path(fc, pt, 0, VRNA_PATH_DEFAULT | VRNA_PATH_NO_TRANSITION_OUTPUT);

    char *loc_min = vrna_db_from_ptable(pt);
    fprintf(stdout, "%4d\t%s\n", i, loc_min);

    free(loc_min);
    free(pt);
  }

  vrna_fold_compound_free(fc);
  free(rec_sequence);

  return 1; /* success */
}


/*
 *  ****************************
 *  * 2. Path / Barrier Modes  *
 *  ****************************
 */
int
paths_findpath(const char       *rec_id,
               const char       *orig_sequence,
               char             **structures,
               struct options_s *opt)
{
  char        *rec_sequence = strdup(orig_sequence);
  vrna_path_t *foldingPath, *Saddle, *r;

  vrna_seq_toRNA(rec_sequence);
  vrna_seq_toupper(rec_sequence);

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute direct folding path");
    free(rec_sequence);
    return 0; /* failure */
  }

  vrna_fold_compound_t *fc = vrna_fold_compound(rec_sequence,
                                                &(opt->md),
                                                VRNA_OPTION_EVAL_ONLY);
  
  foldingPath = vrna_path_findpath(fc,
                                   structures[0],
                                   structures[1],
                                   opt->max_keep);

  Saddle = getSaddlePoint(foldingPath);

  fprintf(stdout, "# direct Path:\n# barrier: %6.2f\n\n", Saddle->en - foldingPath->en);
  for (r = foldingPath; r->s; r++)
    fprintf(stdout, "%s %6.2f\n", r->s, r->en);

  free_path(foldingPath);
  vrna_fold_compound_free(fc);
  free(rec_sequence);

  return 1; /* success */
}


int
paths_pathfinder_gd(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt)
{
  char *sequence;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute folding path");
    return 0; /* failure */
  }

  sequence = strdup(orig_sequence);
  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

#if 1
  vrna_path_t *foldingPath, *r;
  rnax_path_finder_opt_t *options = rnax_path_finder_options();
  options->md          = opt->md;
  options->iterations  = opt->iterations;
  options->max_keep    = opt->max_keep;

  foldingPath = rnax_path_finder(sequence,
                                 structures[0],
                                 structures[1],
                                 options);

  fprintf(stdout,
          "# barrier: %6.2f\n\n",
          getSaddlePoint(foldingPath)->en);

  for (r = foldingPath; r->s; r++)
    fprintf(stdout, "%s %6.2f\n", r->s, r->en);

  free_path(foldingPath);
  free(options);
  fflush(stdout);
#else
  levelSaddlePoint(sequence,
                   structures[0],
                   structures[1],
                   opt->iterations,
                   opt->max_keep,
                   GRADIENT_WALK,
                   opt->max_storage);
#endif

  free(sequence);

  return 1; /* success */
}


int
paths_pathfinder_mc(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt)
{
  char *sequence;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute folding path");
    return 0; /* failure */
  }

  sequence = strdup(orig_sequence);
  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  levelSaddlePoint(sequence,
                   structures[0],
                   structures[1],
                   opt->iterations,
                   opt->max_keep,
                   MC_METROPOLIS,
                   opt->max_storage);

  free(sequence);

  return 1; /* success */
}


int
paths_pathfinder_mcsa(const char        *rec_id,
                      const char        *orig_sequence,
                      char              **structures,
                      struct options_s  *opt)
{
  char *sequence;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute folding path");
    return 0; /* failure */
  }

  sequence = strdup(orig_sequence);
  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  levelSaddlePoint(sequence,
                   structures[0],
                   structures[1],
                   opt->iterations,
                   opt->max_keep,
                   MC_METROPOLIS,
                   opt->max_storage);

  free(sequence);

  return 1; /* success */
}


int
paths_pathfinder_db(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt)
{
  char        *sequence;
  vrna_path_t *foldingPath, *r;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute folding path");
    return 0; /* failure */
  }

  sequence = strdup(orig_sequence);
  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  foldingPath = levelSaddlePoint2(sequence,
                                  structures[0],
                                  structures[1],
                                  0,
                                  opt->iterations,
                                  opt->max_keep,
                                  opt->max_storage,
                                  opt->max_d1,
                                  opt->max_d2);

  fprintf(stdout,
          "\n# done\n\n# Path with detours:\n# barrier: %6.2f\n\n",
          getSaddlePoint(foldingPath)->en);

  for (r = foldingPath; r->s; r++)
    fprintf(stdout, "%s %6.2f\n", r->s, r->en);

  free_path(foldingPath);
  free(sequence);
  fflush(stdout);

  return 1; /* success */
}


int
paths_pathfinder_dbba(const char        *rec_id,
                      const char        *orig_sequence,
                      char              **structures,
                      struct options_s  *opt)
{
  char *sequence;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute energy barrier approximation");
    return 0; /* failure */
  }

  sequence = strdup(orig_sequence);
  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  barrier_estimate_2D(sequence,
                      &(opt->md),
                      structures[0],
                      structures[1],
                      opt->max_d1,
                      opt->max_d2);

  free(sequence);

  return 1; /* success */
}


/*
 *  *********************
 *  * 3. Sampling Modes *
 *  *********************
 */
int
sampling_repulsion(const char       *rec_id,
                   const char       *orig_sequence,
                   char             **structures,
                   struct options_s *opt)
{
  char                  *sequence = strdup(orig_sequence);

  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  vrna_fold_compound_t  *fc = vrna_fold_compound(sequence,
                                                 &(opt->md),
                                                 VRNA_OPTION_DEFAULT);

  repellant_sampling(fc);

  vrna_fold_compound_free(fc);
  free(sequence);

  return 1; /* success */
}


int
sampling_attraction(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt)
{
  char                  *sequence;
  int                   num_ref = 0;
  vrna_fold_compound_t  *fc;
  gridLandscapeT        *grid;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to perform directed sampling");
    return 0; /* failure */
  }

  /* count number of reference structures */
  for (num_ref = 0; structures[num_ref]; num_ref++);

  sequence = strdup(orig_sequence);

  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  fc = vrna_fold_compound(sequence,
                          &(opt->md),
                          VRNA_OPTION_DEFAULT);

#if 0
  grid = estimate_landscape(fc,
                            structures,
                            num_ref,
                            opt->iterations,
                            extended_options);
#endif
  grid = estimate_landscapeMD(fc,
                              (const char **)structures,
                              num_ref,
                              opt->iterations,
                              extended_options,
                              indicesAndPercentages,
                              length_indicesAndPercentages);

  printLandscape(grid, fc);

  free_gridLandscape(grid);

  vrna_fold_compound_free(fc);
  free(sequence);

  return 1; /* success */
}


int
sampling_temperature(const char       *rec_id,
                     const char       *orig_sequence,
                     char             **structures,
                     struct options_s *opt)
{
  vrna_message_warning("Not implemented yet!");
  return 1; /* success */
}


struct sc_data{
    vrna_hash_table_t base_pairs;
    vrna_hash_table_t weights;
};

typedef struct key_value_ {
  vrna_move_t *key;
  int value;
} key_value;

static unsigned
hash_function_base_pairs(void           *hash_entry,
                   unsigned long  hashtable_size)
{
  key_value *kv = (key_value *)hash_entry;
  unsigned int  hash_value = ((unsigned int)kv->key->pos_5) << 16;
  hash_value = hash_value | kv->key->pos_3;
  return hash_value % hashtable_size;
}


static int
hash_comparison_base_pairs(void *x,
                     void *y)
{
  key_value *hem_x  = ((key_value *)x);
  key_value *hem_y  = ((key_value *)y);

  if ((x == NULL) ^ (y == NULL))
    return 1;

  return !(hem_x->key->pos_5 == hem_y->key->pos_5 && hem_x->key->pos_3 == hem_y->key->pos_3);
}


static int
free_base_pairs(void *x)
{
    //do nothing (free whole array)
  return 0;
}

static vrna_hash_table_t
create_hashtable(int hashbits)
{
  vrna_callback_ht_free_entry       *my_free          = free_base_pairs;
  vrna_callback_ht_compare_entries  *my_comparison    = hash_comparison_base_pairs;
  vrna_callback_ht_hash_function    *my_hash_function = hash_function_base_pairs;
  vrna_hash_table_t                 ht                = vrna_ht_init(hashbits,
                                                                     my_comparison,
                                                                     my_hash_function,
                                                                     my_free);

  return ht;
}

typedef struct hashtable_list_ {
  unsigned long     length;
  unsigned long     allocated_size;
  float      *list_weights;
  int        *list_counts;
  key_value         **list_key_value_pairs;
  vrna_hash_table_t ht_pairs; // lookup table;
} hashtable_list;

static hashtable_list
create_hashtable_list(int hashbits)
{
  hashtable_list ht_list;

  ht_list.allocated_size          = 10;
  ht_list.length                  = 0;
  ht_list.list_weights = vrna_alloc(sizeof(float) * ht_list.allocated_size);
  ht_list.list_counts    = vrna_alloc(sizeof(int) * ht_list.allocated_size);
  ht_list.list_key_value_pairs    = vrna_alloc(sizeof(key_value *) * ht_list.allocated_size);
  ht_list.ht_pairs         = create_hashtable(hashbits);
  return ht_list;
}


static
void
free_hashtable_list(hashtable_list *ht_list)
{
  vrna_ht_free(ht_list->ht_pairs);
  free(ht_list->list_weights);
  free(ht_list->list_counts);

  int i = 0;
  for (; i < ht_list->length; i++){
    free(ht_list->list_key_value_pairs[i]->key);
    free(ht_list->list_key_value_pairs[i]);
  }
  free(ht_list->list_key_value_pairs);
}

static
void
hashtable_list_add_weight_and_count(hashtable_list *htl, vrna_move_t *key, float weight)
{
  if (htl->ht_pairs != NULL) {
    key_value to_check;
    to_check.key = key;
    //to_check->value = 0; //not checked anyways --> not set

    //to_check.key = energy;
    //to_check.value = count;
    key_value *lookup_result = NULL;
    lookup_result = vrna_ht_get(htl->ht_pairs, (void *)&to_check);
    if (lookup_result == NULL) {
      //value is not in list.
      if (htl->length >= htl->allocated_size) {
        htl->allocated_size           += 10;
        htl->list_weights  =
          vrna_realloc(htl->list_weights, sizeof(float) * htl->allocated_size);
        htl->list_counts  =
          vrna_realloc(htl->list_counts, sizeof(int) * htl->allocated_size);
        htl->list_key_value_pairs = vrna_realloc(htl->list_key_value_pairs,
                                                 sizeof(key_value *) * htl->allocated_size);
      }

      int           list_index = (int)htl->length;
      htl->list_weights[list_index]  = weight;
      htl->list_counts[list_index]  = 1;
      key_value *to_insert = vrna_alloc(sizeof(key_value));
      to_insert->key = vrna_alloc(sizeof(vrna_move_t));
      to_insert->key->pos_5 = key->pos_5;
      to_insert->key->pos_3 = key->pos_3;
      to_insert->value = list_index;
      htl->list_key_value_pairs[list_index] = to_insert;
      htl->length++;
      int           res         = vrna_ht_insert(htl->ht_pairs, (void *)to_insert);
      if (res != 0)
        fprintf(stderr, "dos.c: hash table insert failed!");
    } else {
      // the energy-index pair is already in the list.
      int list_index = lookup_result->value;
      htl->list_counts[list_index] += 1;
      htl->list_weights[list_index] += weight;
    }
  }
}

void store_basepair_sc(vrna_fold_compound_t *fc, hashtable_list *data, char *structure, float weight, int distance_based /*= False */){
    if(distance_based == 1){
        return;
    }
    else{
        short *pt = vrna_ptable(structure);
        // count number of pairs in structure to repell
        int cnt = 0;
        for(int i = 1; i < pt[0]; i++){
            if(pt[i] > i){
                cnt = cnt + 1;
            }
        }
        if(cnt > 0){
            weight = weight / (float)cnt;
        }
        // add repulsion
        int i;
        for(i = 1; i < pt[0]; i++){
            if(pt[i] > i){
                vrna_move_t key;
                key.pos_5 = i;
                key.pos_3 = pt[i];
                hashtable_list_add_weight_and_count(data, &key, weight);
                /*
                if key not in data['base_pairs']:
                    data['base_pairs'][key] = 1
                    data['weights'][key] = weight
                else:
                    data['base_pairs'][key] = data['base_pairs'][key] + 1
                    data['weights'][key] = data['weights'][key] + weight
                 */
            }
        }
        // remove previous soft constraints
        //fc->sc_remove();
        vrna_sc_init(fc);

        // add latest penalties for base pairs
        int j;
        //int list_index;
        for(int k = 0; k < (int)data->length; k++){ // in data['weights'].keys():
            i = data->list_key_value_pairs[k]->key->pos_5;
            j = data->list_key_value_pairs[k]->key->pos_3;
            //list_index = data->list_key_value_pairs[k]->value;
            float pair_weight = data->list_weights[k];
            //fc->sc_add_bp(i, j, pair_weight);
            vrna_sc_add_bp(fc, i, j, pair_weight, VRNA_OPTION_DEFAULT);
        }
        free(pt);
    }
}


/* ----------------------------------------------------------------- */

/*
 * --------------------------------------------------------------------
 * mix -- mix 3 32-bit values reversibly.
 * For every delta with one or two bits set, and the deltas of all three
 * high bits or all three low bits, whether the original value of a,b,c
 * is almost all zero or is uniformly distributed,
 * If mix() is run forward or backward, at least 32 bits in a,b,c
 * have at least 1/4 probability of changing.
 * If mix() is run forward, every bit of c will change between 1/3 and
 * 2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
 * mix() takes 36 machine instructions, but only 18 cycles on a superscalar
 * machine (like a Pentium or a Sparc).  No faster mixer seems to work,
 * that's the result of my brute-force search.  There were about 2^^68
 * hashes to choose from.  I only tested about a billion of those.
 * --------------------------------------------------------------------
 */
static void mix(unsigned int a, unsigned int b, unsigned int c)
  {
    a -= b; a -= c; a ^= (c >> 13);
    b -= c; b -= a; b ^= (a << 8);
    c -= a; c -= b; c ^= (b >> 13);
    a -= b; a -= c; a ^= (c >> 12);
    b -= c; b -= a; b ^= (a << 16);
    c -= a; c -= b; c ^= (b >> 5);
    a -= b; a -= c; a ^= (c >> 3);
    b -= c; b -= a; b ^= (a << 10);
    c -= a; c -= b; c ^= (b >> 15);
  }

typedef struct key_value_structure_ {
  char *key;
  int value;
} key_value_structure;

unsigned int
hash_function_string(void           *x,
                     unsigned long  hashtable_size)
{
  register char  *k;           /* the key */
  register unsigned int   length;       /* the length of the key */
  register unsigned int   initval = 0;  /* the previous hash, or an arbitrary value */
  register unsigned int   a, b, c, len;

  /* Set up the internal state */
  k   = ((key_value_structure *)x)->key;
  len = length = (unsigned int)strlen(k);
  a   = b = 0x9e3779b9; /* the golden ratio; an arbitrary value */
  c   = initval;        /* the previous hash value */

  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a +=
      ((unsigned int)k[0] + ((unsigned int)k[1] << 8) + ((unsigned int)k[2] << 16) + ((unsigned int)k[3] << 24));
    b +=
      ((unsigned int)k[4] + ((unsigned int)k[5] << 8) + ((unsigned int)k[6] << 16) + ((unsigned int)k[7] << 24));
    c +=
      ((unsigned int)k[8] + ((unsigned int)k[9] << 8) + ((unsigned int)k[10] << 16) +
       ((unsigned int)k[11] << 24));
    mix(a, b, c);
    k   += 12;
    len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch (len) {
    /* all the case statements fall through */
    case 11:
      c += ((unsigned int)k[10] << 24);
    case 10:
      c += ((unsigned int)k[9] << 16);
    case 9:
      c += ((unsigned int)k[8] << 8);
    /* the first byte of c is reserved for the length */
    case 8:
      b += ((unsigned int)k[7] << 24);
    case 7:
      b += ((unsigned int)k[6] << 16);
    case 6:
      b += ((unsigned int)k[5] << 8);
    case 5:
      b += k[4];
    case 4:
      a += ((unsigned int)k[3] << 24);
    case 3:
      a += ((unsigned int)k[2] << 16);
    case 2:
      a += ((unsigned int)k[1] << 8);
    case 1:
      a += (unsigned int)k[0];
      /* case 0: nothing left to add */
  }
  mix(a, b, c);
  /*-------------------------------------------- report the result */
  return c % hashtable_size;
}


/* ----------------------------------------------------------------- */
int
hash_compare_string(void  *x,
                void  *y)
{
  return strcmp(((key_value_structure *)x)->key,
                ((key_value_structure *)y)->key);
}


int
hash_free_string(void *hash_entry)
{
  //free(((key_value_structure *)hash_entry)->key);
  return 0;
}

static vrna_hash_table_t
create_hashtable_string(int hashbits)
{
  vrna_callback_ht_free_entry       *my_free          = hash_free_string;
  vrna_callback_ht_compare_entries  *my_comparison    = hash_compare_string;
  vrna_callback_ht_hash_function    *my_hash_function = hash_function_string;
  vrna_hash_table_t                 ht                = vrna_ht_init(hashbits,
                                                                     my_comparison,
                                                                     my_hash_function,
                                                                     my_free);

  return ht;
}

typedef struct hashtable_list_index_weight_ {
  unsigned long     length;
  unsigned long     allocated_size;
  double      *list_weights;
  int        *list_index;
  key_value_structure         **list_key_value_pairs;
  vrna_hash_table_t ht_pairs; // lookup table;
} hashtable_list_index_weight;

static hashtable_list_index_weight
create_hashtable_list_index_weight(int hashbits)
{
  hashtable_list_index_weight ht_list;

  ht_list.allocated_size          = 10;
  ht_list.length                  = 0;
  ht_list.list_weights = vrna_alloc(sizeof(double) * ht_list.allocated_size);
  ht_list.list_index    = vrna_alloc(sizeof(int) * ht_list.allocated_size);
  ht_list.list_key_value_pairs    = vrna_alloc(sizeof(key_value_structure *) * ht_list.allocated_size);
  ht_list.ht_pairs         = create_hashtable_string(hashbits);
  return ht_list;
}


static
void
free_hashtable_list_index_weight(hashtable_list_index_weight *ht_list)
{
  vrna_ht_free(ht_list->ht_pairs);
  free(ht_list->list_weights);
  free(ht_list->list_index);

  int i = 0;
  for (; i < ht_list->length; i++){
    free(ht_list->list_key_value_pairs[i]->key);
    free(ht_list->list_key_value_pairs[i]);
  }
  free(ht_list->list_key_value_pairs);
}


static
key_value_structure* hashtable_list_index_weight_lookup(hashtable_list_index_weight *htl, char *structure_key){
  if (htl->ht_pairs != NULL) {
    key_value_structure to_check;
    to_check.key = structure_key;
    key_value_structure *lookup_result = NULL;
    lookup_result = vrna_ht_get(htl->ht_pairs, (void *)&to_check);
    return lookup_result;
  }
  else{
    return NULL;
  }
}

static
void
hashtable_list_add_weight_and_index(hashtable_list_index_weight *htl, char *structure_key, int index, double weight)
{
  if (htl->ht_pairs != NULL) {
    key_value_structure *lookup_result = hashtable_list_index_weight_lookup(htl, structure_key);
    if (lookup_result == NULL) {
      //value is not in list.
      if (htl->length >= htl->allocated_size) {
        htl->allocated_size           += 10;
        htl->list_weights  =
          vrna_realloc(htl->list_weights, sizeof(double) * htl->allocated_size);
        htl->list_index  =
          vrna_realloc(htl->list_index, sizeof(int) * htl->allocated_size);
        htl->list_key_value_pairs = vrna_realloc(htl->list_key_value_pairs,
                                                 sizeof(key_value_structure *) * htl->allocated_size);
      }

      int           list_index = (int)htl->length;
      htl->list_weights[list_index]  = weight;
      htl->list_index[list_index]  = index;
      key_value_structure *to_insert = vrna_alloc(sizeof(key_value_structure));
      to_insert->key = vrna_alloc(sizeof(char)*(strlen(structure_key)+1));
      to_insert->key = strcpy(to_insert->key, structure_key);
      to_insert->value = list_index;
      htl->list_key_value_pairs[list_index] = to_insert;
      htl->length++;
      int           res         = vrna_ht_insert(htl->ht_pairs, (void *)to_insert);
      if (res != 0)
        fprintf(stderr, "dos.c: hash table insert failed!");
    } else {
      // the energy-index pair is already in the list.
      int list_index = lookup_result->value;
      htl->list_index[list_index] = index;
      htl->list_weights[list_index] = weight;
    }
  }
}

short * detect_local_minimum(vrna_fold_compound_t *fc, short *structure_pt){
    short * result_pt = vrna_ptable_copy(structure_pt);
    vrna_path_gradient(fc, result_pt, VRNA_PATH_DEFAULT | VRNA_PATH_NO_TRANSITION_OUTPUT);
    return result_pt;
}

/**
 * Take a local minimum and detect nearby local minimum via extended gradient walks.
 * If a lower energy structure within radius 2 is detected (lower than the local minimum
 * of a normal gradient walk), then another gradient walk is applied for the lower structure.
 * @brief detect_local_minimum_two
 * @param fc
 * @param structure_pt
 * @return
 */
char * detect_local_minimum_two(vrna_fold_compound_t *fc, char *structure){
    short *structure_pt = vrna_ptable(structure);
    short *deepest_neighbor = vrna_ptable_copy(structure_pt);
    int deepest_neighbor_ddG  = 1;
    //pt                    = RNA.IntVector(RNA.ptable(structure))
    vrna_move_t *neigh                 = vrna_neighbors(fc, structure_pt, VRNA_MOVESET_DELETION | VRNA_MOVESET_INSERTION);

    vrna_move_t *nb;
    for(nb = neigh; nb->pos_5 != 0 && nb->pos_3 != 0; nb++){
        int dG_nb   = vrna_eval_move_pt(fc, structure_pt, nb->pos_5, nb->pos_3);
        short *pt_nb = vrna_ptable_copy(structure_pt);
        vrna_move_apply(pt_nb, nb);

        vrna_move_t *neigh_2                = vrna_neighbors(fc, pt_nb, VRNA_MOVESET_DELETION | VRNA_MOVESET_INSERTION);

         vrna_move_t *nb_2;
        for(nb_2 = neigh_2; nb_2->pos_5 != 0 && nb_2->pos_3 != 0; nb_2++){
            int dG_nb2   = vrna_eval_move_pt(fc, pt_nb, nb_2->pos_5, nb_2->pos_3);

            int ddG     = dG_nb + dG_nb2;
            if(ddG < 0 && ddG < deepest_neighbor_ddG){
                deepest_neighbor_ddG  = ddG;
                memcpy(deepest_neighbor, pt_nb, sizeof(short)*(pt_nb[0]+1));
                vrna_move_apply(deepest_neighbor, nb_2);

            }
        }
        free(neigh_2);
        free(pt_nb);
    }
    free(neigh);

    if(memcmp(deepest_neighbor, structure_pt, sizeof(short)*(structure_pt[0]+1)) != 0){
        short *tmp_deepest = detect_local_minimum(fc, deepest_neighbor);
        free(deepest_neighbor);
        char * deepest_neighbor_string = vrna_db_from_ptable(tmp_deepest);
        free(tmp_deepest);
        char * result = detect_local_minimum_two(fc, deepest_neighbor_string);
        free(deepest_neighbor_string);
        free(structure_pt);
        return result;
    }
    else{
        free(deepest_neighbor);
        char *result = vrna_db_from_ptable(structure_pt);
        free(structure_pt);
        return result;
    }
}

typedef struct structure_and_index_{
    char *structure;
    int index;
} structure_and_index;


/*
 * --------------------------------------------------------------------
 * hash() -- hash a variable-length key into a 32-bit value
 * k       : the key (the unaligned variable-length array of bytes)
 * len     : the length of the key, counting by bytes
 * initval : can be any 4-byte value
 * Returns a 32-bit value.  Every bit of the key affects every bit of
 * the return value.  Every 1-bit and 2-bit delta achieves avalanche.
 * About 6*len+35 instructions.
 *
 * The best hash table sizes are powers of 2.  There is no need to do
 * mod a prime (mod is sooo slow!).  If you need less than 32 bits,
 * use a bitmask.  For example, if you need only 10 bits, do
 * h = (h & hashmask(10));
 * In which case, the hash table should have hashsize(10) elements.
 *
 * If you are hashing n strings (char **)k, do it like this:
 * for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);
 *
 * By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
 * code any way you wish, private, educational, or commercial.  It's free.
 *
 * See http://burtleburtle.net/bob/hash/evahash.html
 * Use for hash table lookup, or anything where one collision in 2^^32 is
 * acceptable.  Do NOT use for cryptographic purposes.
 * --------------------------------------------------------------------
 */
unsigned int
ht_db_hash_func_strings(void           *x,
                     unsigned long  hashtable_size)
{
  register char  *k;           /* the key */
  register unsigned int   length;       /* the length of the key */
  register unsigned int   initval = 0;  /* the previous hash, or an arbitrary value */
  register unsigned int   a, b, c, len;

  /* Set up the internal state */
  k   = ((structure_and_index *)x)->structure;
  len = length = (unsigned int)strlen(k);
  a   = b = 0x9e3779b9; /* the golden ratio; an arbitrary value */
  c   = initval;        /* the previous hash value */

  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a +=
      ((unsigned int)k[0] + ((unsigned int)k[1] << 8) + ((unsigned int)k[2] << 16) + ((unsigned int)k[3] << 24));
    b +=
      ((unsigned int)k[4] + ((unsigned int)k[5] << 8) + ((unsigned int)k[6] << 16) + ((unsigned int)k[7] << 24));
    c +=
      ((unsigned int)k[8] + ((unsigned int)k[9] << 8) + ((unsigned int)k[10] << 16) +
       ((unsigned int)k[11] << 24));
    mix(a, b, c);
    k   += 12;
    len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch (len) {
    /* all the case statements fall through */
    case 11:
      c += ((unsigned int)k[10] << 24);
    case 10:
      c += ((unsigned int)k[9] << 16);
    case 9:
      c += ((unsigned int)k[8] << 8);
    /* the first byte of c is reserved for the length */
    case 8:
      b += ((unsigned int)k[7] << 24);
    case 7:
      b += ((unsigned int)k[6] << 16);
    case 6:
      b += ((unsigned int)k[5] << 8);
    case 5:
      b += k[4];
    case 4:
      a += ((unsigned int)k[3] << 24);
    case 3:
      a += ((unsigned int)k[2] << 16);
    case 2:
      a += ((unsigned int)k[1] << 8);
    case 1:
      a += k[0];
      /* case 0: nothing left to add */
  }
  mix(a, b, c);
  /*-------------------------------------------- report the result */
  return c % hashtable_size;
}


/* ----------------------------------------------------------------- */
int
ht_db_comp_strings(void  *x,
                void  *y)
{
  return strcmp(((structure_and_index *)x)->structure,
                ((structure_and_index *)y)->structure);
}


int
ht_db_free_entry_strings(void *hash_entry)
{
  return 0;
}


static vrna_hash_table_t
create_string_hashtable(int hashbits)
{
  vrna_callback_ht_free_entry       *my_free          = ht_db_free_entry_strings;
  vrna_callback_ht_compare_entries  *my_comparison    = ht_db_comp_strings;
  vrna_callback_ht_hash_function    *my_hash_function = ht_db_hash_func_strings;
  vrna_hash_table_t                 ht                = vrna_ht_init(hashbits,
                                                                     my_comparison,
                                                                     my_hash_function,
                                                                     my_free);

  return ht;
}

typedef struct hashtable_list_strings_ {
  unsigned long     length;
  unsigned long     allocated_size;
  int        *list_counts;
  float        *list_energies;
  structure_and_index         **list_key_value_pairs; // structure and index
  vrna_hash_table_t ht_pairs; // lookup table;
} hashtable_list_strings;

static hashtable_list_strings
create_hashtable_list_strings(int hashbits)
{
  hashtable_list_strings ht_list;

  ht_list.allocated_size          = 10;
  ht_list.length                  = 0;
  ht_list.list_counts    = vrna_alloc(sizeof(int) * ht_list.allocated_size);
  ht_list.list_energies    = vrna_alloc(sizeof(float) * ht_list.allocated_size);
  ht_list.list_key_value_pairs    = vrna_alloc(sizeof(structure_and_index *) * ht_list.allocated_size);
  ht_list.ht_pairs         = create_string_hashtable(hashbits);
  return ht_list;
}


static
void
free_hashtable_list_strings(hashtable_list_strings *ht_list)
{
    if(ht_list){
      vrna_ht_free(ht_list->ht_pairs);
      ht_list->ht_pairs = NULL;
      free(ht_list->list_counts);
      free(ht_list->list_energies);

      int i = 0;
      for (; i < (int)ht_list->length; i++){
        if(ht_list->list_key_value_pairs[i] != NULL){
          free(ht_list->list_key_value_pairs[i]->structure);
          free(ht_list->list_key_value_pairs[i]);
        }
      }
      free(ht_list->list_key_value_pairs);
    }
}

static
void
hashtable_list_strings_add_structure_and_count(hashtable_list_strings *htl, short *pt_structure, float energy, int count)
{
  if (htl->ht_pairs != NULL) {
    structure_and_index to_check;
    to_check.structure = vrna_db_from_ptable(pt_structure);
    //to_check->value = 0; //not checked anyways --> not set

    //to_check.key = energy;
    //to_check.value = count;
    structure_and_index *lookup_result = NULL;
    lookup_result = vrna_ht_get(htl->ht_pairs, (void *)&to_check);
    if (lookup_result == NULL) {
      //value is not in list.
      if (htl->length >= htl->allocated_size) {
        htl->allocated_size           += 10;
        htl->list_counts  =
          vrna_realloc(htl->list_counts, sizeof(int) * htl->allocated_size);
        htl->list_energies  =
          vrna_realloc(htl->list_energies, sizeof(float) * htl->allocated_size);
        htl->list_key_value_pairs = vrna_realloc(htl->list_key_value_pairs,
                                                 sizeof(structure_and_index *) * htl->allocated_size);
      }

      int           list_index = (int)htl->length;
      htl->list_counts[list_index]  = count;
      htl->list_energies[list_index]  = energy;
      structure_and_index *to_insert = vrna_alloc(sizeof(structure_and_index));
      to_insert->structure = vrna_db_from_ptable(pt_structure);
      to_insert->index = list_index;
      htl->list_key_value_pairs[list_index] = to_insert;
      htl->length++;
      int           res         = vrna_ht_insert(htl->ht_pairs, (void *)to_insert);
      if (res != 0)
        fprintf(stderr, "Error: hash table insert failed!");
    } else {
      // the energy-index pair is already in the list.
      int list_index = lookup_result->index;
      htl->list_counts[list_index] += count;
    }
    free(to_check.structure);
  }
}

void reduce_lm_two_neighborhood(vrna_fold_compound_t *fc, hashtable_list_strings *lm, int verbose /*= False*/){
    int cnt       = 1;
    int cnt_max   = (int)lm->length; //len(lm)
    int lm_remove_allocated = 10;
    int lm_remove_length = 0;
    int* lm_remove = vrna_alloc(sizeof(int)*(lm_remove_allocated)); // = list()
    hashtable_list_strings lm_novel = create_hashtable_list_strings(13); //  = dict()

    if(verbose == 1){
        fprintf(stderr, "Applying 2-Neighborhood Filter...\n");
    }
    int i;
    for(i=0; i < (int)lm->length; i++){ //s in lm:
        char *s = lm->list_key_value_pairs[i]->structure;
        int s_count = lm->list_counts[i];
        if(verbose){
            fprintf(stderr, "\rApplying 2-Neighborhood Filter...%6d / %6d\n", cnt, cnt_max);
        }

        char *ss = detect_local_minimum_two(fc, s);

        if(strcmp(s, ss) != 0){ //(s != ss){
            // store structure for removal
            //lm_remove.append(s)
            if(lm_remove_length >= lm_remove_allocated){
                lm_remove_allocated += 100;
                lm_remove = vrna_realloc(lm_remove, sizeof(int)*lm_remove_allocated);
            }
            lm_remove[lm_remove_length++] = i;
            //if ss not in lm:
            structure_and_index to_check;
            to_check.structure = ss;
            structure_and_index *lookup_result = NULL;
            lookup_result = vrna_ht_get(lm->ht_pairs, (void *)&to_check);
            if (lookup_result == NULL) {
              //value is not in list.
                // store structure for novel insertion
                //if ss not in lm_novel:
                lookup_result = vrna_ht_get(lm_novel.ht_pairs, (void *)&to_check);
                if (lookup_result == NULL) {
                    //lm_novel[ss] = { 'count' : lm[s]['count'], 'energy' : fc.eval_structure(ss) }
                    float energy_kcal = vrna_eval_structure(fc,ss);
                    int count = 1;
                    short *ss_pt = vrna_ptable(ss);
                    hashtable_list_strings_add_structure_and_count(&lm_novel, ss_pt, energy_kcal, count);
                    free(ss_pt);
                }
                else{
                    int lm_ss_count = lm_novel.list_counts[lookup_result->index];
                    //lm_novel[ss]['count'] = lm_novel[ss]['count'] + lm[s]['count']
                    lm_novel.list_counts[lookup_result->index] = lm_ss_count + s_count;
                }
            }
            else{
                // update current minima list
                //lm[ss]['count'] = lm[ss]['count'] + lm[s]['count']
                lm->list_counts[lookup_result->index] += s_count;
            }
        }
        free(ss);
        cnt = cnt + 1;
    }
    // remove obsolete local minima
    //for a in lm_remove:
    //    del lm[a]
    for(i = 0; i < lm_remove_length; i++){
        int lm_index_to_remove = lm_remove[i];
        //remove at index and fill the empty locations in the subsequent merge process.
        lm->list_counts[lm_index_to_remove] = -1;
        lm->list_energies[lm_index_to_remove] = -1;
        structure_and_index* key_to_remove = lm->list_key_value_pairs[lm_index_to_remove];
        vrna_ht_remove(lm->ht_pairs, (void*)key_to_remove);
        free(key_to_remove->structure);
        free(key_to_remove);
        lm->list_key_value_pairs[lm_index_to_remove] = NULL;
    }
    free(lm_remove);

    // add newly detected local minima
    //lm.update(lm_novel)
    //merge lm and lm_novel

    //get empty places
    int length_empty_places = 0;
    int allocated_empty_places = 10;
    int *lm_empty_places = vrna_alloc(sizeof(int)*allocated_empty_places);
    for(i=0; i < lm->length; i++){
        if(lm->list_key_value_pairs[i] == NULL){
            if(length_empty_places >= allocated_empty_places){
                allocated_empty_places += 100;
                lm_empty_places = vrna_realloc(lm_empty_places, sizeof(int)*allocated_empty_places);
            }
            lm_empty_places[length_empty_places++] = i;
        }
    }

    int insert_index = 0;
    for(i=0; i < lm_novel.length; i++){
        structure_and_index *to_check = lm_novel.list_key_value_pairs[i];
        structure_and_index *lookup_result = vrna_ht_get(lm->ht_pairs, (void *)to_check);
        if (lookup_result == NULL) {
            if(insert_index < length_empty_places){
                int empty_place = lm_empty_places[insert_index];
                lm->list_counts[empty_place] = lm_novel.list_counts[i];
                lm->list_energies[empty_place] = lm_novel.list_energies[i];
                structure_and_index *to_insert = vrna_alloc(sizeof(structure_and_index));
                to_insert->structure = vrna_alloc(sizeof(char)*(strlen(lm_novel.list_key_value_pairs[i]->structure)+1));
                strcpy(to_insert->structure, lm_novel.list_key_value_pairs[i]->structure);
                to_insert->index = empty_place;
                lm->list_key_value_pairs[empty_place] = to_insert;
            }
            else{
                //insert at the end
                short *pt = vrna_ptable(lm_novel.list_key_value_pairs[i]->structure);
                float energy = lm_novel.list_energies[i];
                int count = lm_novel.list_counts[i];
                hashtable_list_strings_add_structure_and_count(lm,pt, energy, count);
                free(pt);
            }
            insert_index++;
        }
    }
    free(lm_empty_places);
    free_hashtable_list_strings(&lm_novel);

    if(verbose){
       fprintf(stderr, "\rApplying 2-Neighborhood Filter...done             \n");
    }
}

char ** generate_samples(vrna_fold_compound_t *fc, int number, int non_redundant /*=False */){
    char **samples;

    if(non_redundant){
        samples = vrna_pbacktrack_num(fc, number, VRNA_PBACKTRACK_NON_REDUNDANT);
    }
    else{
        samples = vrna_alloc(sizeof(char*)*(number+1));
        int i;
        for(i =0; i < number; i++){
            char *structure = vrna_pbacktrack(fc);
            samples[i] = structure;
        }
        samples[number] = NULL;
    }
    return samples;
}


int find_max_count(hashtable_list_strings *htl){
    int max_count = (int)htl->length > 0 ? htl->list_counts[0] : 0;
    int res_index = -1;
    int i;
    for(i=0; i < (int)htl->length; i++){
        if(htl->list_key_value_pairs[i] != NULL && htl->list_counts[i] >= max_count){
            max_count = htl->list_counts[i];
            res_index = i;
        }
    }
    return res_index;
}

int find_min_count(hashtable_list_strings *htl){
    int min_count = (int)htl->length > 0 ? htl->list_counts[0] : 0;
    int res_index = -1;
    int i;
    for(i=0; i < (int)htl->length; i++){
        if(htl->list_key_value_pairs[i] != NULL && htl->list_counts[i] <= min_count){
            min_count = htl->list_counts[i];
            res_index = i;
        }
    }
    return res_index;
}

int find_max_energy(hashtable_list_strings *htl){
    float max_energy = (int)htl->length > 0 ? htl->list_energies[0] : 0.f;
    int res_index = -1;
    int i;
    for(i=0; i < (int)htl->length; i++){
        if(htl->list_key_value_pairs[i] != NULL && htl->list_energies[i] >= max_energy){
            max_energy = htl->list_energies[i];
            res_index = i;
        }
    }
    return res_index;
}

int find_min_energy(hashtable_list_strings *htl){
    float min_energy = (int)htl->length > 0 ? htl->list_energies[0] : 0.f;
    int res_index = -1;
    int i;
    for(i=0; i < (int)htl->length; i++){
        if(htl->list_key_value_pairs[i] != NULL && htl->list_energies[i] <= min_energy){
            min_energy = htl->list_energies[i];
            res_index = i;
        }
    }
    return res_index;
}

typedef struct index_energy_ {
  int index;
  float energy;
} index_energy;

int compare_energy_index(const void * elem1, const void * elem2)
{
  index_energy* f = (index_energy*)elem1;
  index_energy* s = (index_energy*)elem2;
    if (f->energy > s->energy) return  1;
    if (f->energy < s->energy) return -1;
    return 0;
}

index_energy *sort_key_values_by_energy(hashtable_list_strings htl){

  index_energy *index_energy_list = vrna_alloc(sizeof(index_energy)*(htl.length+1));
  int i;
  for(i=0;i<htl.length;i++){
    index_energy_list[i].index = i;
    index_energy_list[i].energy = htl.list_energies[i];
  }
  qsort(index_energy_list, htl.length, sizeof(index_energy), compare_energy_index);
  index_energy_list[htl.length].index = -1;
  return index_energy_list;
}

void RNAlocmin_output(char *sequence, hashtable_list_strings htl, char *filename /*= None*/){
    FILE *f;
    if(filename)
        f = fopen(filename, "w");
    else
        f = stdout;

    fprintf(f, "     %s\n", sequence);
    //for i,s in enumerate(sorted(m.keys(), key=lambda x: m[x]['energy'])):
    index_energy * sorted_indices = sort_key_values_by_energy(htl);

    int i;
    int si;
    for(si=0; sorted_indices[si].index > -1; si++){
        i = sorted_indices[si].index;
        char *s = htl.list_key_value_pairs[i]->structure;
        int count = htl.list_counts[i];
        float energy = htl.list_energies[i];
        fprintf(f, "%4d %s %6.2f %6d\n", si, s, energy, count);
    }
    if(filename)
        fclose(f);
    free(sorted_indices);
}

void RNA2Dfold_output(char *sequence, char *s_mfe, float mfe, char *s1, char *s2, hashtable_list_strings m, char *filename /*= None*/){
    int n         =strlen(sequence); //= len(sequence)
    //int *distances = vrna_alloc(sizeof(int)*n*n); //= [ [ None for j in range(0, n) ] for i in range(0, n) ];
    char **structures = vrna_alloc(sizeof(char*)*n*n);
    float *energies = vrna_alloc(sizeof(float)*n*n);
    float inf = 1000000;
    int i;
    for(i=0; i < n*n; i++)
        energies[i] = inf;

    //for s in m.keys():
    for(i=0; i < m.length; i++){
        if(m.list_key_value_pairs[i] == NULL)
            continue;
        int d1 = vrna_bp_distance(s1, m.list_key_value_pairs[i]->structure);
        int d2 = vrna_bp_distance(s2, m.list_key_value_pairs[i]->structure);
        float energy = m.list_energies[i];
        int index = d1*n + d2;
        if(energies[index] != inf || energy < energies[index]){
            energies[index] = energy;
            structures[index] = m.list_key_value_pairs[i]->structure;
        }

        //if not distances[d1][d2] or m[s]['energy'] < distances[d1][d2]['e']:
        //    distances[d1][d2] = { 's': s, 'e': m[s]['energy'] }
    }
    FILE *f;
    if(filename)
        f = fopen(filename, "w");
    else
        f = stdout;

    fprintf(f, "%s\n%s (%6.2f)\n%s\n%s\n\n\n", sequence, s_mfe, mfe, s1, s2);

    int j;
    for(i=0; i<n; i++){ // i in range(0, n):
        for(j=0; j < n; j++){ // j in range(0, n):
            int index = i*n + j;
            if(structures[index] != NULL)  //distances[i][j] != None:
                fprintf(f, "%d\t%d\ta\tb\tc\t%6.2f\td\t%s\n", i, j, energies[index], structures[index]);

        }
    }
    if(filename)
        fclose(f);
    free(structures);
    free(energies);
}


int
sampling_repellent_heuristic(const char       *rec_id,
                             const char       *orig_sequence,
                             char             **structures,
                             struct options_s *opt){
    if(opt->sequence && strlen(opt->sequence) == 0){
        fprintf(stderr, "Error: the input sequence is not given!");
        exit(1);
    }
    //if(strlen(opt->sequence) != strlen(opt->struc1) || strlen(opt->sequence) != strlen(opt->struc2)){
    //    fprintf(stderr, "Error: the input sequence has to have the same length as structure 1 and structure 2!");
    //    exit(1);
    //}
    if(rec_id)
        printf("%s\n",rec_id);


    char *sequence = NULL;
    char *structure1 = NULL;
    char *structure2 = NULL;

    if(structures){
        int i;
        for(i = 0; structures[i]; i++){
            // count only
        }
        if(i>=2 && orig_sequence && strlen(orig_sequence) == strlen(structures[0]) && strlen(orig_sequence) == strlen(structures[1])){
            sequence = vrna_alloc(sizeof(char)*(strlen(orig_sequence)+1));
            strcpy(sequence, orig_sequence);
            structure1 = structures[0];
            structure2 = structures[1];
        }
    }
    if(sequence == NULL && orig_sequence != NULL){
      sequence = vrna_alloc(sizeof(char)*(strlen(orig_sequence)+1));
      strcpy(sequence, orig_sequence);
    }
    // if sequence was not in stdinput, read it from parameters
    if(sequence == NULL && opt->sequence != NULL){
      sequence = vrna_alloc(sizeof(char)*(strlen(opt->sequence)+1));
      strcpy(sequence, opt->sequence);
    }

    if(sequence == NULL){
        fprintf(stderr, "Error: please enter an RNA sequence as input!\n");
        exit(EXIT_FAILURE);
    }

    if(opt->TwoD_file)
    {
      if(structure1 == NULL || structure2 == NULL){
        if(sequence != NULL && opt->struc1 != NULL && opt->struc2 != NULL){
          if(strlen(sequence) == strlen(opt->struc1) && strlen(opt->sequence) == strlen(opt->struc2))
          {
            sequence = vrna_alloc(sizeof(char)*(strlen(opt->sequence)+1));
            strcpy(sequence, opt->sequence);
            structure1 = opt->struc1;
            structure2 = opt->struc2;
          }
          else{
            //error
            fprintf(stderr, "Error: sequence and structures have different lengths!\n");
            exit(EXIT_FAILURE);
          }
        }
        else{
          //error
          fprintf(stderr, "Error: provide at least a fasta file with structures as standard input or use the command line parameters for sequence, structure 1 and structure 2!\n");
          exit(EXIT_FAILURE);
        }
      }
    }
    printf("%s\n%s\n%s\n", sequence, structure1, structure2);

    // table for penalizing base pairs
    hashtable_list sc_data;
    if(!opt->penalize_structures)
      sc_data = create_hashtable_list(13);


    /*
    Do main stuff
    */
    // init random number generator in RNAlib
    vrna_init_rand();

    // prepare RNAlib fold_compound
    vrna_md_t md;
    vrna_md_set_default(&md);
    md.uniq_ML = 1;
    md.compute_bpp = 0;

    vrna_exp_param_t *exp_param = vrna_exp_params(&md);
    double kT =  exp_param->kT;

    vrna_fold_compound_t *fc = vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT | VRNA_OPTION_MFE | VRNA_OPTION_PF);
    vrna_fold_compound_t *fc_base = vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT | VRNA_OPTION_MFE | VRNA_OPTION_PF);

    // compute MFE structure
    char *ss_mfe = vrna_alloc(sizeof(char)*(strlen(sequence)+1));
    float mfe = vrna_mfe(fc,ss_mfe);
    vrna_pf(fc,ss_mfe);

    hashtable_list_index_weight penalized_structures;
    if(opt->penalize_structures)
      penalized_structures = create_hashtable_list_index_weight(13);

    hashtable_list_strings minima = create_hashtable_list_strings(13);// = dict()

    int sample_list_length=0;
    int sample_list_allocated = 100;
    char **sample_list = vrna_alloc(sizeof(char *) * (sample_list_allocated+1)); // = []

    hashtable_list_strings pending_lm = create_hashtable_list_strings(13);// = dict()
    //int num_sc = 1;

    int num_iter      = (int)(ceil(opt->num_samples / (float)opt->granularity));
    int samples_left  = opt->num_samples;

    int current_num_samples = 0;
    int it;
    for(it = 0; it < num_iter; it++){ // in range(0, num_iter):
        // determine number of samples for this round
        // usually, this is 'granularity'
        if(samples_left < opt->granularity)
            current_num_samples = samples_left;
        else
            current_num_samples = opt->granularity;

        samples_left = samples_left - current_num_samples;

        // generate samples through stocastic backtracing
        char **sample_set = generate_samples(fc, current_num_samples, opt->non_red);

        fprintf(stderr, "\rsamples so far: %6d / %6d", opt->num_samples - samples_left, opt->num_samples);

        // store samples of this round to global list of samples
        //sample_list = sample_list + sample_set;
        int i;
        for(i = 0; sample_set[i]; i++){}
        if(sample_list_length + i >= sample_list_allocated){
                        sample_list_allocated += i;
                        sample_list = vrna_realloc(sample_list, sizeof(char*)*(sample_list_allocated+1));
                    }
        memcpy(sample_list + sample_list_length, sample_set, sizeof(char*)*(i));
        sample_list_length += i;
        sample_list[sample_list_length] = NULL;

        hashtable_list_strings current_lm = create_hashtable_list_strings(13); // = dict()

        // go through list of sampled structures and determine corresponding local minima
        for(i=0; sample_set[i]; i++){ // s in sample_set:
            char *s = sample_set[i];
            short *s_pt = vrna_ptable(s);
            short *ss = detect_local_minimum(fc_base, s_pt);
            free(s_pt);

            char *ss_string = vrna_db_from_ptable(ss);
            structure_and_index to_check;
            to_check.structure = ss_string;
            structure_and_index *lookup_result = vrna_ht_get(current_lm.ht_pairs, (void *)&to_check);
            if (lookup_result == NULL) {
                float energy_kcal = vrna_eval_structure(fc_base, ss_string);
                int count = 1;
                hashtable_list_strings_add_structure_and_count(&current_lm, ss, energy_kcal, count);
            }
            else{
                int index = lookup_result->index;
                current_lm.list_counts[index] += 1;
            }
            free(ss_string);
            free(ss);
        }

        free(sample_set);

        // explore the 2-neighborhood of current local minima to reduce total number of local minima
        if(opt->explore_two_neighborhood)
            reduce_lm_two_neighborhood(fc_base, &current_lm, opt->verbose);


        // transfer local minima obtained in this iteration to list of pending local minima
        for(i=0; i < (int)current_lm.length; i++){ // ss in current_lm:
          if(current_lm.list_key_value_pairs[i] == NULL)
            continue; // maybe it was removed in 2-neighborhood filter.

          char *ss_string = current_lm.list_key_value_pairs[i]->structure;
          structure_and_index to_check;
          to_check.structure = ss_string;
          structure_and_index *lookup_result = vrna_ht_get(pending_lm.ht_pairs, (void *)&to_check);
          if (lookup_result == NULL) {
              short *ss = vrna_ptable(ss_string);
              float energy_kcal = current_lm.list_energies[i];
              int count = current_lm.list_counts[i];
              hashtable_list_strings_add_structure_and_count(&pending_lm, ss, energy_kcal, count);
              free(ss);
          }
          else{
              pending_lm.list_counts[lookup_result->index] += current_lm.list_counts[i];
          }
        }
        free_hashtable_list_strings(&current_lm);

        if(it < num_iter - 1){
            // find out which local minima we've seen the most in this sampling round
            int struct_cnt_max_index = find_max_count(&pending_lm); // max(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['count']))
            //int struct_cnt_min_index = find_min_count(&pending_lm); //min(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['count']))
            //int struct_en_max_index = find_max_energy(&pending_lm); // max(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['energy']))
            int struct_en_min_index = find_min_energy(&pending_lm); // min(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['energy']))

            int cnt_once  = 0;
            int cnt_other = 0;

            //for key, value in sorted(pending_lm.iteritems(), key=lambda (k, v): v['energy']):
            index_energy * sorted_indices = sort_key_values_by_energy(pending_lm);

            int i;
            int si;
            for(si=0; sorted_indices[si].index > -1; si++){
                i = sorted_indices[si].index;
                if(pending_lm.list_key_value_pairs[i] == NULL)
                    continue;
                int count = pending_lm.list_counts[i];
                if(count == 1)
                    cnt_once = cnt_once + 1;
                else
                    cnt_other = cnt_other + 1;

                // check whether we've seen other local minim way more often than those we've seen just once
                if(cnt_other > 0 && cnt_once > 0){
                    if(cnt_once <= (opt->mu * cnt_other)){
                        float repell_en = 0;
                        if(opt->ediff_penalty){
                            //repell_en = kt_fact * (value['energy'] - pending_lm[struct_en_min]['energy'])
                            float energy = pending_lm.list_energies[i];
                            float en_min = pending_lm.list_energies[struct_en_min_index];
                            repell_en = opt->exploration_factor * (energy - en_min);
                        }
                        else{
                            //repell_en = kt_fact * kT / 1000.
                            repell_en = (float)(opt->exploration_factor * kT / 1000.0);
                        }
                        char *struct_cnt_max = pending_lm.list_key_value_pairs[struct_cnt_max_index]->structure;
                        if(opt->penalize_structures){
                          int is_penalized = 0;
                          key_value_structure *kv = hashtable_list_index_weight_lookup(&penalized_structures, struct_cnt_max);
                          int last_reference_id;
                          double last_reference_weight;
                          if(kv != NULL){
                            is_penalized = 1;
                            last_reference_id = penalized_structures.list_index[kv->value];
                            last_reference_weight = penalized_structures.list_weights[kv->value];
                          }

                          if(is_penalized == 0){
                            last_reference_id = rnax_add_repulsion(fc, struct_cnt_max, repell_en);
                            last_reference_weight = repell_en;
                            //TODO: hashmap add: {strucutre: id, weight}
                            hashtable_list_add_weight_and_index(&penalized_structures, struct_cnt_max, last_reference_id, last_reference_weight);
                          }
                          else{
                            last_reference_weight = last_reference_weight * opt->exploration_factor;
                            rnax_change_repulsion(fc, last_reference_id, last_reference_weight);
                          }
                        }
                        else{
                          store_basepair_sc(fc, &sc_data, struct_cnt_max, repell_en, 0);
                        }

                        vrna_pf(fc, NULL);

                        //for cmk in pending_lm.keys():
                        int j;
                        for(j = 0; j < (int)pending_lm.length; j++){
                            structure_and_index to_check;
                            to_check.structure = pending_lm.list_key_value_pairs[j]->structure;
                            structure_and_index *lookup_result = vrna_ht_get(minima.ht_pairs, (void *)&to_check);
                            if (lookup_result == NULL) {
                                int counts = pending_lm.list_counts[j];
                                float energy = pending_lm.list_energies[j];
                                short *s_pt = vrna_ptable(to_check.structure);
                                hashtable_list_strings_add_structure_and_count(&minima, s_pt,energy, counts);
                                free(s_pt);
                            }
                            else{
                                minima.list_counts[lookup_result->index] += pending_lm.list_counts[j];
                            }
                        }
                        if(pending_lm.ht_pairs){
                            free_hashtable_list_strings(&pending_lm);
                            pending_lm = create_hashtable_list_strings(13);
                        }
                        break;
                    }
                }
            }
            free(sorted_indices);
        }
        else{
            //for cmk in pending_lm.keys():
            for(i = 0; i < (int)pending_lm.length; i++){
                structure_and_index to_check;
                to_check.structure = pending_lm.list_key_value_pairs[i]->structure;
                structure_and_index *lookup_result = vrna_ht_get(minima.ht_pairs, (void *)&to_check);
                if (lookup_result == NULL) {
                    int counts = pending_lm.list_counts[i];
                    float energy = pending_lm.list_energies[i];
                    short *s_pt = vrna_ptable(to_check.structure);
                    hashtable_list_strings_add_structure_and_count(&minima, s_pt,energy, counts);
                    free(s_pt);
                }
                else{
                    minima.list_counts[lookup_result->index] += pending_lm.list_counts[i];
                }
            }
            if(pending_lm.ht_pairs)
                free_hashtable_list_strings(&pending_lm);
        }
    }
    fprintf(stderr," ... done\n");

    if(opt->post_filter_two)
        reduce_lm_two_neighborhood(fc_base, &minima, opt->verbose);

    /* write local minima file and samples file */
    if(opt->lmin_file){
        RNAlocmin_output(sequence, minima, opt->lmin_file);


        char *sample_file= vrna_alloc(sizeof(char)*1);
        sample_file[0] = '\0';
        int rind = 0;
        char *pch = strrchr(opt->lmin_file, '.');
        rind = pch - opt->lmin_file;

        //lmin_file.rfind(".")
        if(pch != NULL){
            sample_file = vrna_realloc(sample_file, sizeof(char)*(strlen(sample_file)+1+rind+8));
            sample_file = strncat(sample_file, opt->lmin_file, rind);
            sample_file = strcat(sample_file, ".samples");
        }
        else{
            sample_file = vrna_realloc(sample_file, sizeof(char)*(strlen(opt->lmin_file)+1+8));
            sample_file = strcat(sample_file, opt->lmin_file);
            sample_file = strcat(sample_file, ".samples");
        }
        FILE *f = fopen(sample_file, "w");
        fprintf(f, "     %s\n", sequence);
        int i;
        if ((sample_list) && (sample_list[0]))
          for(i=0; sample_list[i]; i++){ // s in sample_list:
              fprintf(f,"%s\n",sample_list[i]);
          }
        fclose(f);
        free(sample_file);
    }
    else{
      // print samples and local minima to stdout
      fprintf(stdout, "Samples: \n");
      if ((sample_list) && (sample_list[0])){
        int i;
        for(i=0; sample_list[i]; i++){ // s in sample_list:
            fprintf(stdout,"%s\n",sample_list[i]);
        }
      }
      fprintf(stdout, "Local minima: \n");
      RNAlocmin_output(sequence, minima, NULL);
    }


    if(opt->TwoD_file){
        //(ss, mfe) = fc_base.mfe()
        char *ss = vrna_alloc(sizeof(char)*(strlen(sequence)+1));
        float mfe = vrna_mfe(fc_base, ss);

        RNA2Dfold_output(sequence, ss, mfe, structure1, structure2, minima, opt->TwoD_file);
        free(ss);
    }

    // read a list of sample structures and produce list of non redundant local minima for it
    if(opt->non_red_file){
        char *lmin_nonred_file = "local_minima_nonred.txt";
        unsigned int nr_samples_count = 0;
        unsigned int nr_samples_allocated = 100;
        char **nonredundant_samples = vrna_alloc(sizeof(char*)*nr_samples_allocated); // = []
        FILE *f = fopen(opt->non_red_file, "r");
        if (f == NULL){
            fprintf(stderr, "Error: cannot open non-redundant file!\n");
            exit(EXIT_FAILURE);
        }
        regex_t re;
        regmatch_t rm[1];
        int length = strlen(sequence);
        int buffer_size = ceil(log10(length));
        char * buffer = vrna_alloc(sizeof(char)*(buffer_size+1));
        sprintf(buffer,"%d",length);
        buffer[buffer_size] = '\0';
        /* create structure pattern with exact sequene length */
        char *pattern_start = "([\\.\\(\\)]{";
        char *tofind = vrna_alloc(sizeof(char)*(sizeof(pattern_start) + (buffer_size+1)*2 + 1 + 2 + 1));
        strcpy(tofind,pattern_start);
        strcat(tofind,buffer);
        strcat(tofind,",");
        strcat(tofind,buffer);
        strcat(tofind,"})");
        //fprintf(stderr, "pattern %s\n", tofind);
        if (regcomp(&re, tofind, REG_EXTENDED) != 0)
        {
            fprintf(stderr, "Failed to compile regex '%s'\n", tofind);
            exit(EXIT_FAILURE);
        }
        free(buffer);
        free(tofind);

        char * line = NULL;
        size_t len = 0;
        ssize_t read;
        unsigned int n_lines = 0;
        while ((read = getline(&line, &len, f)) != -1) {
            if (n_lines > 0){
              //printf("Retrieved line of length %zu:\n", read);
              //printf("%s", line);
              int retval = 0;
              if ((retval = regexec(&re, line, 1, rm, REG_EXTENDED)) == 0)
              {
                  //int length = rm[0].rm_eo - rm[0].rm_so;
                  //if(length == strlen(sequence)){
                /* copy the structure match (which is as long as the sequence) into a new string */
                    char *structure = vrna_alloc(sizeof(char)*(length+1));
                    strncpy(structure, line + rm[0].rm_so, sizeof(char) * (length));
                    structure[length] = '\0';

                    if(nr_samples_count >= nr_samples_allocated-1){
                      nr_samples_allocated += 100;
                      nonredundant_samples = vrna_realloc(nonredundant_samples, sizeof(char*)*nr_samples_allocated);
                    }
                    nonredundant_samples[nr_samples_count++] = structure;
                  //}
              }
            }
            n_lines++;
        }
        regfree(&re);
        nonredundant_samples[nr_samples_count] = NULL;
        fclose(f);
        if (line)
            free(line);
        //with open(nonredundant_sample_file) as f:
        //    nonredundant_samples = f.readlines()
        //nonredundant_samples = [x.strip() for x in nonredundant_samples]
        //nonredundant_samples.pop(0)

        hashtable_list_strings nonredundant_minima = create_hashtable_list_strings(13); //= dict()

        //for s in nonredundant_samples:
        int i;
        for(i=0; nonredundant_samples[i]; i++){
            char *s = nonredundant_samples[i];
            short *pt =vrna_ptable(s); //= RNA.IntVector(RNA.ptable(s))
            //fc_base.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
            vrna_path_gradient(fc_base, pt, VRNA_PATH_DEFAULT | VRNA_PATH_NO_TRANSITION_OUTPUT);

            char *ss = vrna_db_from_ptable(pt); //RNA.db_from_ptable(list(pt))

            structure_and_index to_check;
            to_check.structure = ss;
            structure_and_index *lookup_result = vrna_ht_get(nonredundant_minima.ht_pairs, (void *)&to_check);
            int count = 1;
            if (lookup_result == NULL) {
                float energy = vrna_eval_structure(fc_base,ss);
                hashtable_list_strings_add_structure_and_count(&nonredundant_minima, pt, energy, count);
                //new_minima = new_minima + 1;
            }
            else{
                nonredundant_minima.list_counts[lookup_result->index] += 1;
            }
            free(pt);
            free(s);
            free(ss);
            /*
            if ss not in nonredundant_minima:
                 nonredundant_minima[ss] = { 'count' : 1, 'energy' : fc_base.eval_structure(ss) }
                 new_minima = new_minima + 1
            else:
                 nonredundant_minima[ss]['count'] = nonredundant_minima[ss]['count'] + 1
            */
        }
        free(nonredundant_samples);

        //f = open(lmin_nonred_file, 'w')
        f =fopen(lmin_nonred_file, "w");
        //f.write("     %s\n" % sequence)
        fprintf(f, "     %s\n", sequence);

        //for i,s in enumerate(sorted(nonredundant_minima.keys(), key=lambda x: nonredundant_minima[x]['energy'])):
        //    f.write("%4d %s %6.2f %6d\n" % (i, s, nonredundant_minima[s]['energy'], nonredundant_minima[s]['count']))
        //f.close()
        index_energy * sorted_indices = sort_key_values_by_energy(nonredundant_minima);

        int si;
        for(si=0; sorted_indices[si].index > -1; si++){
            i = sorted_indices[si].index;
            char *s = nonredundant_minima.list_key_value_pairs[i]->structure;
            int count = nonredundant_minima.list_counts[i];
            float energy = nonredundant_minima.list_energies[i];
            fprintf(f, "%4d %s %6.2f %6d\n", si, s, energy, count);
        }
        fclose(f);
        free(sorted_indices);

        if(opt->TwoD_file){
            //distances = [ [ None for j in range(0, 200) ] for i in range(0, 200) ];
            int n = strlen(sequence);
            float *energies = vrna_alloc(sizeof(float)*n*n);
            int inf = 1000000;
            int i;
            for(i=0; i < n*n; i++)
                energies[i] = inf;

            for(i=0; i < nonredundant_minima.length; i++){
                char *s = nonredundant_minima.list_key_value_pairs[i]->structure;
                int d1 = vrna_bp_distance(structure1, s);
                int d2 = vrna_bp_distance(structure2, s);
                int index = d1 * n + d2;
                int energy = nonredundant_minima.list_energies[i];
                if(energies[index] != inf || energy < energies[index])
                    energies[index] = energy;

                //if not distances[d1][d2] or nonredundant_minima[s]['energy'] < distances[d1][d2]:
                //    distances[d1][d2] = nonredundant_minima[s]['energy']
            }
            f = fopen("sv11_fake_nonred.2D.out", "w");
            fprintf(f, "%s\n%s (%6.2f)\n%s\n%s\n\n\n", sequence, structure1, mfe, structure1, structure2);
            int j;
            for(i = 0; i < n; i++){ // in range(0, 200):
                for(j = 0; j < n; j++){ // in range(0, 200):
                    //if distances[i][j] != None:
                    //    f.write("%d\t%d\ta\tb\tc\t%6.2f\n" % (i, j, distances[i][j]))
                    if(energies[i*n + j] != inf)
                        fprintf(f,"%d\t%d\ta\tb\tc\t%6.2f\n", i, j, energies[i*n + j]);
                }
            }
            fclose(f);
            free(energies);
        }
        free_hashtable_list_strings(&nonredundant_minima);
    }

    /* free */
    free(sequence);
    free(ss_mfe);
    if(!opt->penalize_structures){
      if(sc_data.ht_pairs)
        free_hashtable_list(&sc_data);
    }
    else{
      if(penalized_structures.ht_pairs)
        free_hashtable_list_index_weight(&penalized_structures);
    }
    if(minima.ht_pairs)
      free_hashtable_list_strings(&minima);
    if(pending_lm.ht_pairs)
      free_hashtable_list_strings(&pending_lm);
    int i;
    if ((sample_list) && (sample_list[0]))
      for(i=0; sample_list[i]; i++)
        free(sample_list[i]);
    free(sample_list);
    vrna_fold_compound_free(fc);
    vrna_fold_compound_free(fc_base);
    free(exp_param);

    return EXIT_SUCCESS;
}

int
retrieve_local_minima(const char       *rec_id,
                      const char       *orig_sequence,
                      char             **structures,
                      struct options_s *opt){
  double temperature_celsius = opt->temperature_celsius;
  int shift_moves = opt->shift_moves;
  char *parameter_file = opt->parameter_file;
  return gradient_walker(temperature_celsius, shift_moves, parameter_file, orig_sequence, structures);
}


