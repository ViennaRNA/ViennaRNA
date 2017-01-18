/*
 *              Distances of Secondary Structure Ensembles
 *        Peter F Stadler, Ivo L Hofacker, Sebastian Bonhoeffer
 *                      Vienna RNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/dist_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/profiledist.h"
#include "ViennaRNA/ProfileAln.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/PS_dot.h"
#include "RNApaln_cmdl.h"

#define MAXLENGTH  10000
#define MAXSEQ      1000


static double gapo    = 1.5, gape = 0.666, seqw = 0.5;
static int    endgaps = 0;

PRIVATE void command_line(int   argc,
                          char  *argv[]);


PRIVATE void print_aligned_lines(FILE *somewhere);


PRIVATE char  task;
PRIVATE char  outfile[FILENAME_MAX_LENGTH];
PRIVATE char  ruler[] = "....,....1....,....2....,....3....,....4"
                        "....,....5....,....6....,....7....,....8";
static int    noconv = 0;

int
main(int  argc,
     char *argv[])

{
  float *T[MAXSEQ];
  char  *seq[MAXSEQ];
  int   i, j, istty, n = 0;
  int   type, length, taxa_list = 0;
  float dist;
  FILE  *somewhere = NULL;
  char  *structure;
  char  *line = NULL, fname[FILENAME_MAX_LENGTH], *list_title = NULL;
  plist *pr_pl, *mfe_pl;


  command_line(argc, argv);

  if ((outfile[0] == '\0') && (task == 'm') && (edit_backtrack))
    strcpy(outfile, "backtrack.file");

  if (outfile[0] != '\0')
    somewhere = fopen(outfile, "w");

  if (somewhere == NULL)
    somewhere = stdout;

  istty = (isatty(fileno(stdout)) && isatty(fileno(stdin)));

  while (1) {
    if ((istty) && (n == 0)) {
      printf("\nInput sequence;  @ to quit\n");
      printf("%s\n", ruler);
    }

    type = 0;
    do {
      /* get sequence to fold */
      if (line != NULL)
        free(line);

      *fname = '\0';
      if ((line = vrna_read_line(stdin)) == NULL) {
        type = 999;
        break;
      }

      if (line[0] == '@')
        type = 999;

      if (line[0] == '*') {
        if (taxa_list == 0) {
          if (task == 'm')
            taxa_list = 1;

          printf("%s\n", line);
          type = 0;
        } else {
          list_title  = strdup(line);
          type        = 888;
        }
      }

      if (line[0] == '>') {
        if (sscanf(line, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname) != 0)
          strcat(fname, "_dp.ps");

        if (taxa_list)
          printf("%d : %s\n", n + 1, line + 1);
        else
          printf("%s\n", line);

        type = 0;
      }

      if (isalpha(line[0])) {
        char *cp;
        cp = strchr(line, ' ');
        if (cp)
          *cp = '\0';

        type = 1;
      }
    } while (type == 0);

    if ((task == 'm') && (type > 800)) {
      if (taxa_list)
        printf("* END of taxa list\n");

      printf("> p %d (pdist)\n", n);
      for (i = 1; i < n; i++) {
        for (j = 0; j < i; j++) {
          printf("%g ", profile_aln(T[i], seq[i], T[j], seq[j]));
          if (edit_backtrack)
            fprintf(somewhere, "> %d %d\n", i + 1, j + 1);

          print_aligned_lines(somewhere);
        }
        printf("\n");
      }
      if (type == 888) {
        /* do another distance matrix */
        n = 0;
        printf("%s\n", list_title);
        free(list_title);
      }
    }

    if (type > 800) {
      for (i = 0; i < n; i++)
        free_profile(T[i]);
      if (type == 888)
        continue;

      if (outfile[0] != '\0')
        (void)fclose(somewhere);

      if (line != NULL)
        free(line);

      return 0; /* finito */
    }

    length = (int)strlen(line);
    for (i = 0; i < length; i++) {
      line[i] = toupper(line[i]);
      if (!noconv && line[i] == 'T')
        line[i] = 'U';
    }

    pr_pl = mfe_pl = NULL;
    {
      double mfe, kT;
      kT        = (temperature + 273.15) * 1.98717 / 1000.; /* in Kcal */
      structure = (char *)vrna_alloc((length + 1) * sizeof(char));
      mfe       = fold(line, structure);
      /* get pairlist from dot-bracket string */
      mfe_pl    = vrna_plist(structure, 0.95 * 0.95);
      pf_scale  = exp(-(1.07 * mfe) / kT / length);
      /* init_pf_fold(length); <- obsolete */
      (void)pf_fold(line, structure);
      /* get pairlist of probability matrix */
      assign_plist_from_pr(&pr_pl, pr, length, 1e-5);
    }

    if (*fname == '\0')
      sprintf(fname, "%d_dp.ps", n + 1);

    /* PS_dot_plot(line, fname); <- NOT THREADSAFE and obsolete function! */

    /* call threadsafe dot plot printing function */
    PS_dot_plot_list(line, fname, pr_pl, mfe_pl, "");

    T[n]    = Make_bp_profile_bppm(pr, length);
    seq[n]  = strdup(line);
    if ((istty) && (task == 'm'))
      printf("%s\n", structure);

    free(structure);
    free(mfe_pl);
    free(pr_pl);
    free_arrays();
    free_pf_arrays();

    n++;
    switch (task) {
      case 'p':
        if (n == 2) {
          dist = profile_aln(T[0], seq[0], T[1], seq[1]);
          printf("%g\n", dist);
          print_aligned_lines(somewhere);
          free_profile(T[0]);
          free_profile(T[1]);
          free(seq[0]);
          free(seq[1]);
          n = 0;
        }

        break;
      case 'f':
        if (n > 1) {
          dist = profile_aln(T[1], seq[1], T[0], seq[0]);
          printf("%g\n", dist);
          print_aligned_lines(somewhere);
          free_profile(T[1]);
          free(seq[1]);
          n = 1;
        }

        break;
      case 'c':
        if (n > 1) {
          dist = profile_aln(T[1], seq[1], T[0], seq[0]);
          printf("%g\n", dist);
          print_aligned_lines(somewhere);
          free_profile(T[0]);
          free(seq[0]);
          T[0]    = T[1];
          seq[0]  = seq[1];
          n       = 1;
        }

        break;

      case 'm':
        break;

      default:
        vrna_message_error("This can't happen.");
    }   /* END switch task */
    (void)fflush(stdout);
  }     /* END while */
  if (line != NULL)
    free(line);

  return 0;
}


/* ----------------------------------------------------------------- */

PRIVATE void
command_line(int  argc,
             char *argv[])
{
  int                       i, sym;
  char                      *ns_bases   = NULL, *c;
  char                      *ParamFile  = NULL;
  struct  RNApaln_args_info args_info;

  task = 'p';

  if (RNApaln_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* temperature */
  if (args_info.temp_given)
    temperature = args_info.temp_arg;

  /* do not take special tetra loop energies into account */
  if (args_info.noTetra_given)
    tetra_loop = 0;

  /* set dangle model */
  if (args_info.dangles_given)
    dangles = args_info.dangles_arg;

  /* do not allow weak pairs */
  if (args_info.noLP_given)
    noLonelyPairs = 1;

  /* do not allow wobble pairs (GU) */
  if (args_info.noGU_given)
    noGU = 1;

  /* do not allow weak closing pairs (AU,GU) */
  if (args_info.noClosingGU_given)
    no_closingGU = 1;

  /* set energy model */
  if (args_info.energyModel_given)
    energy_set = args_info.energyModel_arg;

  /* take another energy parameter set */
  if (args_info.paramFile_given)
    ParamFile = strdup(args_info.paramFile_arg);

  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if (args_info.nsp_given)
    ns_bases = strdup(args_info.nsp_arg);

  /* set the mode */
  if (args_info.mode_given) {
    if (strlen(args_info.mode_arg) != 1) {
      RNApaln_cmdline_parser_print_help();
      exit(EXIT_FAILURE);
    } else {
      task = *args_info.mode_arg;
    }

    switch (task) {
      case 'p':
      case 'm':
      case 'f':
      case 'c':
        break;
      default:
        RNApaln_cmdline_parser_print_help();
        exit(EXIT_FAILURE);
    }
  }

  if (args_info.printAlignment_given) {
    if (strcmp(args_info.printAlignment_arg, "stdout")) {
      strncpy(outfile, args_info.printAlignment_arg, FILENAME_MAX_LENGTH - 1);
      outfile[FILENAME_MAX_LENGTH - 1] = '\0';
    } else {
      outfile[0] = '\0';
    }

    edit_backtrack = 1;
  }

  /* gap opening penalty */
  if (args_info.gapo_given)
    gapo = args_info.gapo_arg;

  /* gap extension penalty */
  if (args_info.gape_given)
    gape = args_info.gape_arg;

  /* sequence weight */
  if (args_info.seqw_given)
    seqw = args_info.seqw_arg;

  /* endgaps */
  if (args_info.endgaps_given)
    endgaps = 1;

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (args_info.noconv_given)
    noconv = 1;

  /* free allocated memory of command line data structure */
  RNApaln_cmdline_parser_free(&args_info);

  /* fprintf(stderr, "%f %f %f %d\n", gapo, gape, seqw, -endgaps); */
  set_paln_params(gapo, gape, seqw, 1 - endgaps);

  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (ns_bases != NULL) {
    nonstandards  = vrna_alloc(33);
    c             = ns_bases;
    i             = sym = 0;
    if (*c == '-') {
      sym = 1;
      c++;
    }

    while (*c != '\0') {
      if (*c != ',') {
        nonstandards[i++] = *c++;
        nonstandards[i++] = *c;
        if ((sym) && (*c != *(c - 1))) {
          nonstandards[i++] = *c;
          nonstandards[i++] = *(c - 1);
        }
      }

      c++;
    }
  }
}


/*--------------------------------------------------------------------------*/

PRIVATE void
print_aligned_lines(FILE *somewhere)
{
  if (edit_backtrack)
    fprintf(somewhere, "%s\n%s\n%s\n%s\n",
            aligned_line[2], aligned_line[0],
            aligned_line[3], aligned_line[1]);
}


/*--------------------------------------------------------------------------*/
