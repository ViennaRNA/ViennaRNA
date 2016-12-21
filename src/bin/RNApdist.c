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
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/dist_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/profiledist.h"
#include "RNApdist_cmdl.h"


#define MAXLENGTH  10000
#define MAXSEQ      1000

PRIVATE void command_line(int       argc,
                          char      *argv[],
                          vrna_md_t *md);


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
  float     *T[MAXSEQ];
  int       i, j, istty, n = 0;
  int       type, length, taxa_list = 0;
  float     dist;
  FILE      *somewhere = NULL;
  char      *structure;
  char      *line = NULL, fname[FILENAME_MAX_LENGTH], *list_title = NULL;
  plist     *pr_pl, *mfe_pl;
  vrna_md_t md;

  /* assign globally stored model details */
  set_model_details(&md);

  pr_pl = mfe_pl = NULL;

  command_line(argc, argv, &md);

  if ((outfile[0] == '\0') && (task == 'm') && edit_backtrack)
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

    type = 0; /* what is this type anyway? */
    do {
      /* get sequence to fold */
      *fname = '\0';

      if (line != NULL)
        free(line);

      /* end do...while() loop if no more input */
      if ((line = vrna_read_line(stdin)) == NULL) {
        type = 999;
        break;
      }

      /* end do...while() loop if user requested abort */
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

      /* we currently read a fasta header */
      if (line[0] == '>') {
        if (sscanf(line, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname) != 0)
          strcat(fname, "_dp.ps");

        if (taxa_list)
          printf("%d : %s\n", n + 1, line + 1);
        else
          printf("%s\n", line);

        type = 0;
      }

      /* this seems to be a crude way to identify an RNA sequence */
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
          printf("%g ", profile_edit_distance(T[i], T[j]));
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

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if (!noconv)
      vrna_seq_toRNA(line);

    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(line);

    if (*fname == '\0')
      sprintf(fname, "%d_dp.ps", n + 1);

    vrna_fold_compound_t *vc = vrna_fold_compound(line, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);

    structure = (char *)vrna_alloc((vc->length + 1) * sizeof(char));

    (void)vrna_pf(vc, structure);

    pr_pl = vrna_plist_from_probs(vc, 1e-5);
    /* fake plist for lower part, since it stays empty */
    mfe_pl      = (plist *)vrna_alloc(sizeof(plist));
    mfe_pl[0].i = mfe_pl[0].j = 0;

    /* call threadsafe dot plot printing function */
    PS_dot_plot_list(line, fname, pr_pl, mfe_pl, "");

    T[n] = Make_bp_profile_bppm(vc->exp_matrices->probs, vc->length);

    if ((istty) && (task == 'm'))
      printf("%s\n", structure);

    free(structure);
    free(mfe_pl);
    free(pr_pl);
    vrna_fold_compound_free(vc);

    n++;
    switch (task) {
      case 'p':
        if (n == 2) {
          dist = profile_edit_distance(T[0], T[1]);
          printf("%g\n", dist);
          print_aligned_lines(somewhere);
          free_profile(T[0]);
          free_profile(T[1]);
          n = 0;
        }

        break;
      case 'f':
        if (n > 1) {
          dist = profile_edit_distance(T[1], T[0]);
          printf("%g\n", dist);
          print_aligned_lines(somewhere);
          free_profile(T[1]);
          n = 1;
        }

        break;
      case 'c':
        if (n > 1) {
          dist = profile_edit_distance(T[1], T[0]);
          printf("%g\n", dist);
          print_aligned_lines(somewhere);
          free_profile(T[0]);
          T[0]  = T[1];
          n     = 1;
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
command_line(int        argc,
             char       *argv[],
             vrna_md_t  *md)
{
  int                         i, sym;
  char                        *ns_bases   = NULL, *c;
  char                        *ParamFile  = NULL;
  struct  RNApdist_args_info  args_info;

  task = 'p';

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNApdist_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* temperature */
  if (args_info.temp_given)
    md->temperature = temperature = args_info.temp_arg;

  /* do not take special tetra loop energies into account */
  if (args_info.noTetra_given)
    md->special_hp = tetra_loop = 0;

  /* set dangle model */
  if (args_info.dangles_given) {
    md->dangles = dangles = args_info.dangles_arg;
    if (dangles)
      md->dangles = dangles = 2;
  }

  /* set energy model */
  if (args_info.energyModel_given)
    md->energy_set = energy_set = args_info.energyModel_arg;

  /* do not allow weak pairs */
  if (args_info.noLP_given)
    md->noLP = noLonelyPairs = 1;

  /* do not allow wobble pairs (GU) */
  if (args_info.noGU_given)
    md->noGU = noGU = 1;

  /* do not allow weak closing pairs (AU,GU) */
  if (args_info.noClosingGU_given)
    md->noGUclosure = no_closingGU = 1;

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (args_info.noconv_given)
    noconv = 1;

  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if (args_info.nsp_given)
    ns_bases = strdup(args_info.nsp_arg);

  /* take another energy parameter set */
  if (args_info.paramFile_given)
    ParamFile = strdup(args_info.paramFile_arg);

  if (args_info.compare_given) {
    switch (args_info.compare_arg[0]) {
      case 'p':
      case 'm':
      case 'f':
      case 'c':
        task = args_info.compare_arg[0];
        break;
      default:
        RNApdist_cmdline_parser_print_help();
        exit(EXIT_FAILURE);
    }
  }

  if (args_info.backtrack_given) {
    if (strcmp(args_info.backtrack_arg, "none"))
      strncpy(outfile, args_info.backtrack_arg, FILENAME_MAX_LENGTH - 1);

    edit_backtrack = 1;
  }

  /* free allocated memory of command line data structure */
  RNApdist_cmdline_parser_free(&args_info);

  /* do some preparations */
  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (ns_bases != NULL)
    vrna_md_set_nonstandards(md, ns_bases);
}


/* ------------------------------------------------------------------------- */

PRIVATE void
print_aligned_lines(FILE *somewhere)
{
  if (edit_backtrack)
    fprintf(somewhere, "%s\n%s\n", aligned_line[0], aligned_line[1]);
}


/*--------------------------------------------------------------------------*/
