/*
 * Ineractive Access to cofolding routines
 * c Ivo L Hofacker
 * Vienna RNA package
 */

/**
*** \file RNAup.c
**/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include "ViennaRNA/fold.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/part_func_up.h"
#include "ViennaRNA/duplex.h"
#include "ViennaRNA/energy_const.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/constraints.h"
#include "RNAup_cmdl.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define EQUAL(A, B) (fabs((A)-(B)) < 1000 * DBL_EPSILON)

PRIVATE void    tokenize(char *line,
                         char **seq1,
                         char **seq2);


PRIVATE void    seperate_bp(char  **inter,
                            int   len1,
                            char  **intra_l,
                            char  **intra_s);


PRIVATE void    print_interaction(interact    *Int,
                                  char        *s1,
                                  char        *s2,
                                  pu_contrib  *p_c,
                                  pu_contrib  *p_c2,
                                  int         w,
                                  int         incr3,
                                  int         incr5);


PRIVATE void    print_unstru(pu_contrib *p_c,
                             int        w);


PRIVATE int     compare_unpaired_values(const void  *p1,
                                        const void  *p2);


PRIVATE int     move_useless_unpaired_values(const void *p1,
                                             const void *p2);


PRIVATE void    adjustUnpairedValues(int ***unpaired_values); /* this function sorts and cleans up the unpaired values given at command line */


/* defaults for -u and -w */
PRIVATE int     default_u; /* -u options for plotting: plot pr_unpaired for 4 nucleotides */
PRIVATE double  RT;

/*--------------------------------------------------------------------------*/
int
main(int  argc,
     char *argv[])
{
  struct RNAup_args_info  args_info;
  unsigned int            input_type, up_mode;
  char                    my_contrib[10], *up_out, *name, fname1[FILENAME_MAX_LENGTH],
                          fname2[FILENAME_MAX_LENGTH], fname_target[FILENAME_MAX_LENGTH], *ParamFile,
                          *ns_bases, *c, *structure, *head, *input_string, *s1, *s2, *s3, *s_target, *cstruc1,
                          *cstruc2, *cstruc_target, *cstruc_combined, *cmdl_parameters, *orig_s1, *orig_s2,
                          *orig_target;
  int                     i, j, length1, length2, length_target, sym, istty,
                          rotated, noconv, max_u, **unpaired_values, ulength_num;
  double                  energy, min_en, sfact;

  /* variables for output */
  pu_contrib              *unstr_out, *unstr_short, *unstr_target, *contrib1, *contrib2;
  interact                *inter_out;
  /* pu_out *longer; */

  /* commandline parameters */
  int                     w       = 25; /* length of region of interaction */
  int                     incr3   = 0;  /* add x unpaired bases after 3'end of short RNA*/
  int                     incr5   = 0;  /* add x unpaired bases after 5'end of short RNA*/
  int                     header  = 1;  /* print header in output file */
  int                     output  = 1;  /* create output  file */

  /* more default settings for RNAup */
  up_mode       = RNA_UP_MODE_1; /* default RNAup mode, single sequence unpaired probabilities */
  my_contrib[0] = 'S';
  my_contrib[1] = '\0';

  default_u = 4;

  /* early initializing */
  noconv          = 0;
  max_u           = 0;
  ulength_num     = 0;      /* number of ulength values given on commandline */
  sfact           = 1.07;
  dangles         = 2;
  do_backtrack    = 1;
  rotated         = 0;
  input_string    = s1 = s2 = s3 = s_target = cstruc1 = cstruc2 = cstruc_target = cstruc_combined = NULL;
  length1         = length2 = length_target = 0;
  inter_out       = NULL;
  unstr_out       = unstr_short = unstr_target = contrib1 = contrib2 = NULL;
  structure       = ParamFile = ns_bases = head = orig_s1 = orig_s2 = orig_target = NULL;
  up_out          = NULL;
  fname_target[0] = '\0';
  /* allocate init length for commandline parameter string */

  cmdl_parameters = NULL;
  vrna_strcat_printf(&cmdl_parameters, "RNAup ");

  /*
   #############################################
   # check the command line prameters
   #############################################
   */
  if (RNAup_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* do not create header */
  if (args_info.no_header_given)
    header = 0;

  /* temperature */
  if (args_info.temp_given) {
    temperature = args_info.temp_arg;
    /* collect parameter if header is needed */
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "-T %f ", temperature);
  }

  /* structure constraint */
  if (args_info.constraint_given) {
    fold_constrained = 1;
    /* collect parameter if header is needed */
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "-C ");
  }

  /* do not take special tetra loop energies into account */
  if (args_info.noTetra_given) {
    tetra_loop = 0;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "-4 ");
  }

  /* set dangle model */
  if (args_info.dangles_given) {
    if ((args_info.dangles_arg != 0) && (args_info.dangles_arg != 2))
      vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    else
      dangles = args_info.dangles_arg;

    if (header)
      vrna_strcat_printf(&cmdl_parameters, "-d %d ", dangles);
  }

  /* do not allow weak pairs */
  if (args_info.noLP_given) {
    noLonelyPairs = 1;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "--noLP ");
  }

  /* do not allow wobble pairs (GU) */
  if (args_info.noGU_given) {
    noGU = 1;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "--noGU ");
  }

  /* do not allow weak closing pairs (AU,GU) */
  if (args_info.noClosingGU_given) {
    no_closingGU = 1;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "--noClosingGU ");
  }

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (args_info.noconv_given) {
    noconv = 1;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "--noconv ");
  }

  /* set energy model */
  if (args_info.energyModel_given) {
    energy_set = args_info.energyModel_arg;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "-e %d ", energy_set);
  }

  /* take another energy parameter set */
  if (args_info.paramFile_given) {
    ParamFile = strdup(args_info.paramFile_arg);
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "-P %s ", ParamFile);
  }

  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if (args_info.nsp_given) {
    ns_bases = strdup(args_info.nsp_arg);
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "--nsp %s ", ns_bases);
  }

  /* set pf scaling factor */
  if (args_info.pfScale_given) {
    sfact = args_info.pfScale_arg;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "-S %f ", sfact);
  }

  /* set the maximal length of interaction region */
  if (args_info.window_given) {
    w = args_info.window_arg;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "-w %d ", w);
  }

  /* do not make an output file */
  if (args_info.no_output_file_given)
    output = 0;

  /* set mode to unpaired regions in both RNAs */
  if (args_info.include_both_given) {
    up_mode = RNA_UP_MODE_3;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "--include_both ");
  }

  /* set interaction mode 1 (pairwise interaction) */
  if (args_info.interaction_pairwise_given) {
    up_mode = RNA_UP_MODE_2;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "--interaction_pairwise ");
  }

  /* set interaction mode 2 (first sequence interacts with all others) */
  if (args_info.interaction_first_given) {
    up_mode = RNA_UP_MODE_3;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "--interaction_first ");
  }

  /* extend unpaired region 5' */
  if (args_info.extend5_given) {
    incr5 = args_info.extend5_arg;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "-5 %d ", incr5);
  }

  /* extend unpaired region 3' */
  if (args_info.extend3_given) {
    incr3 = args_info.extend3_arg;
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "-3 %d ", incr3);
  }

  /* set contribution output */
  if (args_info.contributions_given) {
    strncpy(my_contrib, args_info.contributions_arg, 10);
    if (header)
      vrna_strcat_printf(&cmdl_parameters, "-c %s ", my_contrib);
  }

  /* set length(s) of unpaired (unstructured) region(s) */
  int min, max, tmp;
  i = (args_info.ulength_given == 0) ? 1 : args_info.ulength_given;

  /* here's the very new way of treating multiple ulength values/ranges */
  unpaired_values = (int **)vrna_alloc(sizeof(int *) * (i + 1));

  if (header && args_info.ulength_given)
    vrna_strcat_printf(&cmdl_parameters, "-u ");

  for (i = 0; i < args_info.ulength_given; i++) {
    unpaired_values[++ulength_num] = (int *)vrna_alloc(2 * sizeof(int));
    /* we got a ulength range... */
    if (sscanf(args_info.ulength_arg[i], "%d-%d", &min, &max) == 2) {
      if (min > max) {
        tmp = min;
        min = max;
        max = tmp;
      }

      if (min == max) {
        unpaired_values[ulength_num][0] = min;
        unpaired_values[ulength_num][1] = -1;
      } else {
        unpaired_values[ulength_num][0] = min;
        unpaired_values[ulength_num][1] = max;
      }

      max_u = MAX2(max_u, max);
    } else if (sscanf(args_info.ulength_arg[i], "%d", &max) == 1) {
      unpaired_values[ulength_num][0] = max;
      unpaired_values[ulength_num][1] = -1;
      max_u                           = MAX2(max_u, max);
    }

    if (header) {
      if (i < args_info.ulength_given - 1)
        vrna_strcat_printf(&cmdl_parameters, "%s,", args_info.ulength_arg[i]);
      else
        vrna_strcat_printf(&cmdl_parameters, "%s ", args_info.ulength_arg[i]);
    }
  }
  if (i == 0) {
    /* use default settings */
    unpaired_values[++ulength_num]  = (int *)vrna_alloc(2 * sizeof(int));
    unpaired_values[ulength_num][0] = 4;
    unpaired_values[ulength_num][1] = -1;
    max_u                           = 4;
  }

  /* store number of entries at position [0][0] */
  unpaired_values[0]    = (int *)vrna_alloc(2 * sizeof(int));
  unpaired_values[0][0] = ulength_num;
  unpaired_values[0][1] = -1;

  /* adjust ranges and values such that everything is non-redundant
   * WARNING: after this step ulength_num may not reflect actual number of ulength values anymore */
  adjustUnpairedValues(&unpaired_values);

  /* free allocated memory of command line data structure */
  RNAup_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */

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

  istty = isatty(fileno(stdout)) && isatty(fileno(stdin));
  if ((fold_constrained) && (istty)) {
    vrna_message_constraint_options(VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_X | VRNA_CONSTRAINT_DB_RND_BRACK);
    printf("constraints for intramolecular folding only:\n");
    vrna_message_constraint_options(VRNA_CONSTRAINT_NO_HEADER | VRNA_CONSTRAINT_DB_ANG_BRACK);
    printf("constraints for cofolding (intermolecular folding) only:\n");
    vrna_message_constraint_options(VRNA_CONSTRAINT_NO_HEADER | VRNA_CONSTRAINT_DB_PIPE);
  }

  RT = ((temperature + K0) * GASCONST / 1000.0);
  /*
   #############################################
   # main loop: continue until end of file
   #############################################
   */
  do {
    rotated   = 0;
    cut_point = -1;
    fname1[0] = '\0';
    fname2[0] = '\0';
    /*
     ########################################################
     # handle user input from 'stdin'
     ########################################################
     */
    if (istty) {
      switch (up_mode) {
        case RNA_UP_MODE_1:   /* just calculate the probability of beeing unpaired for the given interval(s) */
          vrna_message_input_seq_simple();
          break;
        case RNA_UP_MODE_2:   /* pairwise interaction mode, former -Xp mode */
          vrna_message_input_seq("Use either '&' to connect the 2 sequences or give each sequence on an extra line.");
          break;
        case RNA_UP_MODE_3:   /* consecutive multi interaction mode ;) first sequence pairs with all following, former -Xf mode */
                              /* either we wait for the first two sequences */
          if (s_target == NULL)
            vrna_message_input_seq("Give each sequence on an extra line. "
                                   "The first seq. is stored, every other seq. is compared to the first one.");
          /* or we already have them and wait for the next sequence */
          else
            vrna_message_input_seq("Enter another sequence.");

          break;
        default:
          vrna_message_error("This should never happen (again)");
          break;
      }
    }

    /* extract filename from fasta header if available */
    while ((input_type = get_input_line(&input_string, 0)) & VRNA_INPUT_FASTA_HEADER) {
      (void)sscanf(input_string, "%" XSTR(FILENAME_ID_LENGTH) "s", fname1);
      printf(">%s\n", input_string); /* print fasta header if available */
      free(input_string);
    }

    /* break on any error, EOF or quit request */
    if (input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR))
      break;                            /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else
      tokenize(input_string, &s1, &s2); /* this also frees the input_string */

    length1 = (int)strlen(s1);
    length2 = (s2) ? (int)strlen(s2) : 0;

    /* now we got the first and maybe a second sequence */

    /* check if we have to change the mode we are operating in */
    if ((cut_point != -1) && (up_mode & RNA_UP_MODE_1)) {
      vrna_message_warning("Two concatenated sequences given, switching to pairwise interaction mode!");
      up_mode = RNA_UP_MODE_2;
    }

    int read_again = 0;

    switch (up_mode) {
      case RNA_UP_MODE_2:
        if (cut_point == -1)
          read_again = 1;

        break;
      case RNA_UP_MODE_3:
        if (cut_point == -1) {
          if (s_target == NULL)
            read_again = 1;
        } else if (s_target != NULL) {
          vrna_message_error(
            "After the first sequence (pair): Input a single sequence (no &)!\n"
            "Each input seq. is compared to the very first seq. given.\n"
            );
        }

        break;
      default:
        break;
    }

    if (read_again) {
      /* we are in this block only if we just have 1 sequence yet but need a second, too */

      /* extract filename from fasta header if available */
      while ((input_type = get_input_line(&input_string, 0)) & VRNA_INPUT_FASTA_HEADER) {
        (void)sscanf(input_string, "%" XSTR(FILENAME_ID_LENGTH) "s", fname2);
        printf(">%s\n", input_string); /* print fasta header if available */
        free(input_string);
      }
      /* break on any error, EOF or quit request */
      if (input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR))
        break;                            /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
      else
        tokenize(input_string, &s2, &s3); /* this also frees the input_string */

      if (cut_point != -1)
        vrna_message_error("Don't confuse me by mixing concatenated (&) with single sequences! Go, have some sleep and then check your input again...");

      length2 = (int)strlen(s2);
    }

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if (!noconv) {
      vrna_seq_toRNA(s1);
      vrna_seq_toRNA(s2);
    }

    /* store case-unmodified sequence */
    orig_s1 = strdup(s1);
    if (s2)
      orig_s2 = strdup(s2);

    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(s1);
    vrna_seq_toupper(s2);

    /** read structure constraint(s) if necessary */
    if (fold_constrained) {
      char  *cstruc_tmp = NULL;
      int   old_cut     = -1;

      input_type = get_input_line(&input_string, VRNA_INPUT_NOSKIP_COMMENTS);
      if (input_type & VRNA_INPUT_QUIT) {
        break;
      } else if ((input_type & VRNA_INPUT_MISC) && (strlen(input_string) > 0)) {
        old_cut   = cut_point;
        cut_point = -1;
        tokenize(input_string, &cstruc1, &cstruc2);
      } else {
        vrna_message_error("constraints missing");
      }

      /* now that we've got the constraining structure(s) check if the input was valid */
      if (old_cut != cut_point)
        vrna_message_error("RNAup -C: mixed single/dual sequence or constraint strings or different cut points");

      read_again = 0;

      if (cut_point == -1) {
        if (up_mode & RNA_UP_MODE_2)
          read_again = 1;
        else if ((up_mode & RNA_UP_MODE_3) && (s_target == NULL))
          read_again = 1;
      }

      if (read_again) {
        input_type = get_input_line(&input_string, VRNA_INPUT_NOSKIP_COMMENTS);
        if (input_type & VRNA_INPUT_QUIT)
          break;
        else if ((input_type & VRNA_INPUT_MISC) && (strlen(input_string) > 0))
          tokenize(input_string, &cstruc2, &cstruc_tmp);
        else
          vrna_message_error("constraints missing");

        if (cut_point != -1)
          vrna_message_error("Don't confuse me by mixing concatenated (&) with single sequences! Go, have some sleep and then check your input again...");
      }

      /* check length(s) of input sequence(s) and constraint(s) */
      if (strlen(cstruc1) != length1) {
        fprintf(stderr, "%s\n%s\n", s1, cstruc1);
        vrna_message_error("RNAup -C: constraint string and structure have unequal length");
      }

      if (s2 != NULL) {
        if (strlen(cstruc2) != length2) {
          fprintf(stderr, "%s\n%s\n", s2, cstruc2);
          vrna_message_error("RNAup -C: constraint string and structure have unequal length");
        }
      }
    } /* thats all for constraint folding */

    /* rotate input sequences if upmode>=2 to ensure first sequence is the longer one */
    if (up_mode & (RNA_UP_MODE_2 | RNA_UP_MODE_3)) {
      if (length1 < length2) {
        rotated = 1;
        /* rotate the sequences such that the longer is the first */
        int   l = length2;
        length2 = length1;
        length1 = l;
        char  *s = s2;
        s2      = s1;
        s1      = s;
        s       = orig_s2;
        orig_s2 = orig_s1;
        orig_s1 = s;
        /* also rotate the file names */
        char f[FILENAME_MAX_LENGTH];
        strncpy(f, fname2, FILENAME_ID_LENGTH);
        strncpy(fname2, fname1, FILENAME_ID_LENGTH);
        strncpy(fname1, f, FILENAME_ID_LENGTH);
        /* rotate constraint strings as well */
        if (fold_constrained) {
          s       = cstruc2;
          cstruc2 = cstruc1;
          cstruc1 = s;
        }
      }
    }

    /* check ulength values against sequences given */
    if (max_u > length1)
      vrna_message_error("maximum unpaired region exceeds sequence length");

    if (up_mode & RNA_UP_MODE_3) {
      /* if we haven't seen the target yet, store it now */
      if (s_target == NULL) {
        if (rotated) {
          s_target      = s2;
          orig_target   = orig_s2;
          s2            = NULL;
          orig_s2       = NULL;
          length_target = length2;
          strcpy(fname_target, fname2);
          if (fold_constrained) {
            cstruc_target = cstruc2;
            cstruc2       = NULL;
          }
        } else {
          s_target      = s1;
          orig_target   = orig_s1;
          length_target = length1;
          s1            = s2;
          orig_s1       = orig_s2;
          s2            = NULL;
          orig_s2       = NULL;
          length1       = length2;
          strcpy(fname_target, fname1);
          strcpy(fname1, fname2);
          if (fold_constrained) {
            cstruc_target = cstruc1;
            cstruc1       = cstruc2;
            cstruc2       = NULL;
          }
        }

        fname2[0] = '\0';
      }
    }

    /*
     ########################################################
     # done with 'stdin' handling
     ########################################################
     */

    /* compose file names */

    /* first file name */
    if (fname1[0] != '\0') {
      up_out = vrna_strdup_printf("%s", fname1);
      if (up_mode & (RNA_UP_MODE_2 | RNA_UP_MODE_3)) {
        if (fname2[0] != '\0')
          vrna_strcat_printf(&up_out, "_%s", fname2);
        else if (fname_target[0] != '\0')
          vrna_strcat_printf(&up_out, "_%s", fname_target);
      }
    } else {
      up_out = vrna_strdup_printf("RNA");
    }

    if (!(up_mode & RNA_UP_MODE_1))
      vrna_strcat_printf(&up_out, "_w%d", w);

    structure = (char *)vrna_alloc(sizeof(char) * (MAX2(length_target, MAX2(length1, length2)) + 1));


    /* begin actual computations */
    update_fold_params();

    /* calc mfe of first sequence */
    if (cstruc1 != NULL)
      strncpy(structure, cstruc1, length1 + 1);

    min_en = fold(s1, structure);

    (void)fflush(stdout);

    /* calc probability to be unstructured for 1st sequence (in upmode=3 this is not the target!) */

    int wplus = w;
    if (!(up_mode & RNA_UP_MODE_3)) {
      wplus += incr3 + incr5;
      /* reset window size if maximum unstructured region is exceeds it */
      if (max_u > wplus)
        wplus = max_u;
    }

    /* reset window size if sequence length is shorter */
    if (length1 < wplus)
      wplus = length1;

    /* calc mfe for first sequence (2nd if upmode = 3) */
    if (cstruc1 != NULL)
      strncpy(structure, cstruc1, length1 + 1);

    min_en    = fold(s1, structure);
    pf_scale  = exp(-(sfact * min_en) / RT / length1);
    if (length1 > 2000)
      vrna_message_info(stderr, "scaling factor %f", pf_scale);

    if (cstruc1 != NULL)
      strncpy(structure, cstruc1, length1 + 1);

    energy    = pf_fold(s1, structure);
    unstr_out = pf_unstru(s1, wplus);
    free_pf_arrays();

    if (fold_constrained && !(up_mode & RNA_UP_MODE_1)) {
      cstruc_combined = (char *)vrna_alloc(sizeof(char) * (length1 + length2 + 1));
      strncpy(cstruc_combined, cstruc1, length1 + 1);
      strcat(cstruc_combined, cstruc2);
    }

    contrib1  = contrib2 = NULL;
    inter_out = NULL;

    switch (up_mode) {
      case RNA_UP_MODE_1:
        for (i = 1; i <= unpaired_values[0][0]; i++) {
          j = unpaired_values[i][0];
          do
            print_unstru(unstr_out, j);
          while (++j <= unpaired_values[i][1]);
        }
        if (output && header)
          head = vrna_strdup_printf("# %s\n# %d %s\n# %s", cmdl_parameters, length1, fname1, orig_s1);

        contrib1 = unstr_out;
        break;
      case RNA_UP_MODE_2:
        inter_out = pf_interact(s1, s2, unstr_out, NULL, w, cstruc_combined, incr3, incr5);
        print_interaction(inter_out, orig_s1, orig_s2, unstr_out, NULL, w, incr3, incr5);
        if (output && header)
          head = vrna_strdup_printf("# %s\n# %d %s\n# %s\n# %d %s\n# %s", cmdl_parameters, length1, fname1, orig_s1, length2, fname2, orig_s2);

        contrib1 = unstr_out;
        break;
      case RNA_UP_MODE_3:   /* calculate prob. unstruct. for target seq */
        if (unstr_target == NULL) {
          wplus = w + incr3 + incr5;
          if (max_u > wplus)
            wplus = max_u;

          if (length_target < wplus)
            wplus = length_target;

          if (cstruc_target != NULL)
            strncpy(structure, cstruc_target, length_target + 1);

          min_en    = fold(s_target, structure);
          pf_scale  = exp(-(sfact * min_en) / RT / length_target);
          if (length_target > 2000)
            vrna_message_info(stderr, "scaling factor %f", pf_scale);

          if (cstruc_target != NULL)
            strncpy(structure, cstruc_target, length_target + 1);

          energy        = pf_fold(s_target, structure);
          unstr_target  = pf_unstru(s_target, wplus);
          free_pf_arrays();                     /* for arrays for pf_fold(...) */
        }

        /* check if target sequence is actually longer than query, if not rotate both sequences */
        if (length_target < length1) {
          inter_out = pf_interact(s1, s_target, unstr_out, unstr_target, w, cstruc_combined, incr3, incr5);
          print_interaction(inter_out, orig_s1, orig_target, unstr_out, unstr_target, w, incr3, incr5);
          contrib1  = unstr_out;
          contrib2  = unstr_target;
        } else {
          inter_out = pf_interact(s_target, s1, unstr_target, unstr_out, w, cstruc_combined, incr3, incr5);
          print_interaction(inter_out, orig_target, orig_s1, unstr_target, unstr_out, w, incr3, incr5);
          contrib1  = unstr_target;
          contrib2  = unstr_out;
        }

        if (output && header)
          head = vrna_strdup_printf("# %s\n# %d %s\n# %s\n# %d %s\n# %s", cmdl_parameters, length_target, fname_target, orig_target, length1, fname1, orig_s1);

        break;
    }

    /* create additional output */
    if (output) {
      /* how to best compose a reasonable filename */
      printf("RNAup output in file: ");
      /* since we do not limit the amount of ulength values anymore we just put
       * the maximum length into the filename, the actual printed lengths
       * should be somewhere in the output itself */
      name = vrna_strdup_printf("%s_u%d.out", up_out, unpaired_values[0][0]);
      printf("%s\n", name);

      Up_plot(contrib1, contrib2, inter_out, name, unpaired_values, my_contrib, head, up_mode);
    }

    /*
     ########################################################
     # clean up
     ########################################################
     */

    /* we can save the pu contribution structure of the target sequence for the next run */
    if (unstr_target != NULL) {
      if (length_target < length1)
        free_pu_contrib_struct(contrib1);
      else
        free_pu_contrib_struct(contrib2);
    } else {
      if (contrib1 != NULL)
        free_pu_contrib_struct(contrib1);

      if (contrib2 != NULL)
        free_pu_contrib_struct(contrib2);
    }

    if (inter_out != NULL)
      free_interact(inter_out);

    /* free all unnecessary character arrays */
    free(structure);
    free(s1);
    free(s2);
    free(orig_s1);
    free(orig_s2);
    free(cstruc1);
    free(cstruc2);
    free(head);
    free(cstruc_combined);
    free(up_out);
    free(name);

    structure = s1 = s2 = orig_s1 = orig_s2 = cstruc1 = cstruc2 = head = cstruc_combined = NULL;
    up_out    = name = NULL;

    free_arrays(); /* for arrays for fold(...) */
  } while (1);
  free(cmdl_parameters);

  return EXIT_SUCCESS;
}


PRIVATE int
compare_unpaired_values(const void  *p1,
                        const void  *p2)
{
  if ((*(int **)p1)[0] > (*(int **)p2)[0])
    return 1;

  if ((*(int **)p1)[0] < (*(int **)p2)[0])
    return -1;

  return 0;
}


PRIVATE int
move_useless_unpaired_values(const void *p1,
                             const void *p2)
{
  if ((*(int **)p1)[1] < (*(int **)p2)[1])
    return 1;

  if ((*(int **)p1)[0] > (*(int **)p2)[0])
    return -1;

  return 0;
}


PRIVATE void
adjustUnpairedValues(int ***unpaired_values)
{
  int i, last_max, real_count;

  if (*unpaired_values == NULL)
    return;

  if ((*unpaired_values)[0][0] <= 0)
    return;

  /* sort the ranges array */
  qsort(&((*unpaired_values)[1]), (*unpaired_values)[0][0], sizeof(int **), compare_unpaired_values);

  last_max    = (*unpaired_values)[1][1] != -1 ? (*unpaired_values)[1][1] : (*unpaired_values)[1][0];
  real_count  = 1;
  for (i = 2; i <= (*unpaired_values)[0][0]; i++) {
    if ((*unpaired_values)[i][1] == -1) {
      /* we just have a single value */
      if ((*unpaired_values)[i][0] <= last_max) {
        (*unpaired_values)[i][1] = -2; /* mark this entry to be removed */
      } else {
        last_max = (*unpaired_values)[i][0];
        real_count++;
      }
    } else {
      /* we have a range of values */
      if (((*unpaired_values)[i][0] <= last_max) && ((*unpaired_values)[i][1] <= last_max)) {
        (*unpaired_values)[i][1] = -2; /* mark this entry to be removed as the whole range is already covered */
      } else if (((*unpaired_values)[i][0] <= last_max) && ((*unpaired_values)[i][1] > last_max)) {
        (*unpaired_values)[i][0]  = last_max + 1;
        last_max                  = last_max + 1;
        if ((*unpaired_values)[i][1] == last_max)
          (*unpaired_values)[i][1] = -1;        /* range reduced to single value */
        else
          last_max = (*unpaired_values)[i][1];  /* maximum of range */

        real_count++;
      } else {
        last_max = (*unpaired_values)[i][1];
        real_count++;
      }
    }
  }

  /* sort entries again to get rid of useless ones */
  qsort(&((*unpaired_values)[1]), (*unpaired_values)[0][0], sizeof(int **), move_useless_unpaired_values);

  /* free memory we dont need anymore */
  for (i = real_count + 1; i <= (*unpaired_values)[0][0]; i++)
    free((*unpaired_values)[i]);
  (*unpaired_values) = (int **)realloc((*unpaired_values), (real_count + 1) * sizeof(int **));

  (*unpaired_values)[0][0] = real_count;
  /* sort the array again */
  qsort(&((*unpaired_values)[1]), (*unpaired_values)[0][0], sizeof(int **), compare_unpaired_values);
}


/* call:  tokenize(line,&seq1,&seq2); the sequence string is split at the "&"
 * and the first seq is written in seq1, the second into seq2  */
/* using sscanf instead of strcpy get's rid of trainling junk on the input line */
void
tokenize(char *line,
         char **seq1,
         char **seq2)
{
  char  *pos;
  int   cut = -1;

  pos = strchr(line, '&');
  if (pos) {
    cut     = (int)(pos - line) + 1;
    (*seq1) = (char *)vrna_alloc((cut + 1) * sizeof(char));
    (*seq2) = (char *)vrna_alloc(((strlen(line) - cut) + 2) * sizeof(char));

    if (strchr(pos + 1, '&'))
      vrna_message_error("more than one cut-point in input");

    *pos = '\0';
    (void)sscanf(line, "%s", *seq1);
    (void)sscanf(pos + 1, "%s", *seq2);
  } else {
    (*seq1) = (char *)vrna_alloc((strlen(line) + 1) * sizeof(char));
    (*seq2) = NULL;
    sscanf(line, "%s", *seq1);
  }

  if (cut > -1) {
    if (cut_point == -1) {
      cut_point = cut;
    } else if (cut_point != cut) {
      fprintf(stderr, "cut_point = %d cut = %d\n", cut_point, cut);
      vrna_message_error("Sequence and Structure have different cut points.");
    }
  }

  free(line);
  return;
}


/* divide the constraints string in intermolecular constrains (inter)
 * and intramolecular constrains within both sequences */
/* len1 is the length of the LONGER input seq ! */
void
seperate_bp(char  **inter,
            int   len1,
            char  **intra_l,
            char  **intra_s)
{
  int   i, j, len;
  short *pt = NULL;
  char  *temp_inter, *pt_inter;

  len = strlen((*inter));
  /* printf("inter\n%s\n",(*inter)); */
  i           = len + 1;
  temp_inter  = (char *)vrna_alloc(sizeof(char) * i);
  /* to make a pair_table convert <|> to (|) */
  pt_inter = (char *)vrna_alloc(sizeof(char) * i);
  /* if shorter seq is first seq in constrained string, write the
   * longer one as the first one */
  temp_inter[strlen((*inter))]  = '\0';
  pt_inter[strlen((*inter))]    = '\0';
  if (cut_point < len1) {
    /* write the constrain for the longer seq first */
    for (j = 0, i = cut_point - 1; i < len; i++, j++) {
      switch ((*inter)[i]) {
        case '(':
          temp_inter[j] = ')';
          pt_inter[j]   = ')';
          break;
        case ')':
          temp_inter[j] = '(';
          pt_inter[j]   = '(';
          break;
        default:
          temp_inter[j] = (*inter)[i];
          pt_inter[j]   = '.';
      }
    }
    /* then add the constrain for the shorter seq */
    for (i = 0; i < cut_point - 1; i++, j++) {
      switch ((*inter)[i]) {
        case '(':
          temp_inter[j] = ')';
          pt_inter[j]   = ')';
          break;
        case ')':
          temp_inter[j] = '(';
          pt_inter[j]   = '(';
          break;
        default:
          temp_inter[j] = (*inter)[i];
          pt_inter[j]   = '.';
      }
    }
    cut_point = len1 + 1;
    strcpy((*inter), temp_inter);
  } else {
    for (i = 0; i < strlen((*inter)); i++) {
      switch ((*inter)[i]) {
        case '(':
          pt_inter[i] = '(';
          break;
        case ')':
          pt_inter[i] = ')';
          break;
        default:
          pt_inter[i] = '.';
      }
    }
  }

  pt = vrna_ptable(pt_inter);

  /* intramolecular structure in longer (_l) and shorter (_s) seq */
  (*intra_l)                              = (char *)vrna_alloc(sizeof(char) * (len1 + 1));
  (*intra_s)                              = (char *)vrna_alloc(sizeof(char) * (strlen((*inter)) - len1 + 2));
  (*intra_l)[len1]                        = '\0';
  (*intra_s)[strlen((*inter)) - len1 + 1] = '\0';
  /* now seperate intermolecular from intramolecular bp */
  for (i = 1; i <= pt[0]; i++) {
    if (pt[i] == 0) {
      temp_inter[i - 1] = (*inter)[i - 1];
      if (i < cut_point) {
        (*intra_l)[i - 1] = (*inter)[i - 1];
        if ((*inter)[i - 1] == '|')
          (*intra_l)[i - 1] = '.';
      } else {
        (*intra_s)[i - cut_point] = (*inter)[i - 1];
        if ((*inter)[i - 1] == '|')
          (*intra_s)[i - cut_point] = '.';
      }
    } else {
      if (i < cut_point) {
        /* intermolekular bp */
        if (pt[i] >= cut_point) {
          temp_inter[i - 1]             = (*inter)[i - 1];
          (*intra_l)[i - 1]             = '.';
          (*intra_s)[pt[i] - cut_point] = '.';
        } else {
          /* intramolekular bp */
          (*intra_l)[i - 1] = (*inter)[i - 1];
          temp_inter[i - 1] = '.';
        }
      } else {
        /* i>=cut_point */
        /* intermolekular bp */
        if (pt[i] < cut_point) {
          temp_inter[i - 1] = (*inter)[i - 1];
          /* (*intra_s)[i-1] = '.'; */
        } else {
          /* intramolekular bp */
          (*intra_s)[i - cut_point] = (*inter)[i - 1];
          temp_inter[i - 1]         = '.';
        }
      }
    }
  }

  /* printf("%s -1\n%s -2\n%s -3\n%s -4\n",(*inter),temp_inter,(*intra_l),(*intra_s)); */
  strcpy((*inter), temp_inter);
  free(temp_inter);
  free(pt_inter);
  free(pt);
}


PRIVATE void
print_interaction(interact    *Int,
                  char        *s1,
                  char        *s2,
                  pu_contrib  *p_c,
                  pu_contrib  *p_c2,
                  int         w,
                  int         incr3,
                  int         incr5)
{
  char    *i_long, *i_short;
  int     i, l_l, l_s, len1, end5, end3, i_min, j_min, l1, add_a, add_b, nix_up;
  double  p_c_S;
  double  G_min, Gi_min, Gul, G_sum, Gus, diff;
  duplexT mfe;
  char    *struc;

  G_min   = Int->Gikjl;
  Gi_min  = Int->Gikjl_wo;
  len1    = Int->length;

  /* use duplexfold() to fold the interaction site */
  l_l     = (Int->i - Int->k + 1);
  i_long  = (char *)vrna_alloc(sizeof(char) * (l_l + 1));
  l_s     = (Int->l - Int->j + 1);
  i_short = (char *)vrna_alloc(sizeof(char) * (l_s + 1));

  strncpy(i_long, &s1[Int->k - 1], l_l);
  i_long[l_l] = '\0';
  strncpy(i_short, &s2[Int->j - 1], l_s);
  i_short[l_s] = '\0';

  mfe = duplexfold(i_long, i_short);

  i_min = mfe.i;
  j_min = mfe.j;
  l1    = strchr(mfe.structure, '&') - mfe.structure;

  /* printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", mfe.structure, i_min+1-l1,
   * i_min, j_min, j_min+strlen(mfe.structure)-l1-2, mfe.energy ); */

  /* structure by duplexfold is shorter than structure by RNAup:*/

  add_a   = add_b = 0; /* length difference in longer / shorter sequence*/
  nix_up  = 0;
  if (((i_min + 1 - l1) - i_min) != (Int->k - Int->i))
    add_a = Int->i - Int->k + 2;

  if (((j_min + strlen(mfe.structure) - l1 - 2) - j_min) != (Int->l - Int->j))
    add_b = Int->l - Int->j + 2;

  /* printf("add_a %d   add_b %d\n",add_a,add_b); */
  if (add_a || add_b) {
    nix_up = 1;
    if (add_a && add_b == 0)
      add_b = Int->l - Int->j + 2;

    if (add_a == 0 && add_b)
      add_a = Int->i - Int->k + 2;

    struc = (char *)vrna_alloc(sizeof(char) * (add_a + add_b + 3));
    for (i = 0; i < (add_a + add_b - 1); i++) {
      if (i != l_l)
        struc[i] = '.';

      if (i == l_l)
        struc[i] = '&';
    }
    struc[i] = '\0';
  } else {
    l1    = strlen(mfe.structure);
    struc = (char *)vrna_alloc(sizeof(char) * (l1 + 1));
    strcpy(struc, mfe.structure);
  }

  end5  = MAX(1, Int->k - incr5);
  end3  = MIN(MIN(l_l - 1 + incr3, w + incr3 + incr5), len1);
  p_c_S = p_c->H[end5][end3] + p_c->I[end5][end3] + p_c->M[end5][end3] + p_c->E[end5][end3];
  Gul   = -RT *log(p_c_S);


  if (p_c2 == NULL) {
    G_sum = Gi_min + Gul;

    /* printf("dG = dGint + dGu_l\n"); */
    printf("%s %3d,%-3d : %3d,%-3d (%.2f = %.2f + %.2f)\n",
           struc, Int->k, Int->i, Int->j, Int->l, G_min, Gi_min, Gul);
    printf("%s&%s\n", i_long, i_short);
  } else {
    p_c_S = p_c2->H[Int->j][(Int->l) - (Int->j)] +
            p_c2->I[Int->j][(Int->l) - (Int->j)] +
            p_c2->M[Int->j][(Int->l) - (Int->j)] +
            p_c2->E[Int->j][(Int->l) - (Int->j)];
    Gus = -RT *log(p_c_S);


    G_sum = Gi_min + Gul + Gus;
    /* printf("dG = dGint + dGu_l + dGu_s\n"); */
    printf("%s %3d,%-3d : %3d,%-3d (%.2f = %.2f + %.2f + %.2f)\n",
           struc, Int->k, Int->i, Int->j, Int->l, G_min, Gi_min, Gul, Gus);
    printf("%s&%s\n", i_long, i_short);
  }

  if (!EQUAL(G_min, G_sum)) {
    printf("ERROR\n");
    diff = fabs((G_min) - (G_sum));
    printf("diff %.18f\n", diff);
  }

  if (nix_up)
    vrna_message_warning("RNAduplex structure doesn't match any structure of RNAup structure ensemble");

  free(i_long);
  free(i_short);
  free(mfe.structure);
  free(struc);
}


/* print coordinates and free energy for the region of highest accessibility */
PRIVATE void
print_unstru(pu_contrib *p_c,
             int        w)
{
  int     i, j, len, min_i, min_j;
  double  dG_u, min_gu;

  if (p_c != NULL) {
    min_gu  = 1000.0;
    len     = p_c->length;

    for (i = 1; i <= len; i++) {
      for (j = i; j < MIN((i + w), len + 1); j++) {
        double blubb;
        if ((j - i + 1) == w && i + w - 1 <= len) {
          blubb = p_c->H[i][j - i] + p_c->I[i][j - i] + p_c->M[i][j - i] + p_c->E[i][j - i];
          dG_u  = -RT *log(blubb);


          if (dG_u < min_gu) {
            min_gu  = dG_u;
            min_i   = i;
            min_j   = i + w - 1;
          }
        }
      }
    }
    printf("%4d,%4d \t (%.3f) \t for u=%3d\n", min_i, min_j, min_gu, w);
  } else {
    vrna_message_error("error with prob unpaired");
  }
}
