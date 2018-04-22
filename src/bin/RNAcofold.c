/*
 *                c Ivo L Hofacker, Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/commands.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/constraints_SHAPE.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func_co.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/structure_utils.h"
#include "ViennaRNA/read_epars.h"
#include "RNAcofold_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"

#include "ViennaRNA/color_output.inc"

struct options {
  int             filename_full;
  char            *filename_delim;
  int             pf;
  int             doT;
  int             doC;
  int             noPS;
  int             noconv;
  double          bppmThreshold;
  int             verbose;
  vrna_md_t       md;
  vrna_cmd_t                        *commands;

  dataset_id                        id_control;

  char            *concentration_file;

  char            *constraint_file;
  int             constraint_batch;
  int             constraint_enforce;
  int             constraint_canonical;

  int             shape;
  char            *shape_file;
  char            *shape_method;
  char            *shape_conversion;

  int             csv_output;
  int             csv_header;
  char            csv_output_delim;
};


PRIVATE vrna_dimer_pf_t do_partfunc(char      *string,
                                    int       length,
                                    int       Switch,
                                    plist     **tpr,
                                    plist     **mf,
                                    double    kT,
                                    struct options *opt);


PRIVATE double *read_concentrations(FILE *fp);


PRIVATE void
print_concentrations(vrna_dimer_conc_t *result,
                     double            *startconc);

static int
process_input(FILE            *input_stream,
              const char      *input_filename,
              struct options  *opt);


static void
write_csv_header(struct options *opt);

void
postscript_layout(vrna_fold_compound_t *fc,
                  const char  *orig_sequence,
                  const char  *structure,
                  const char *SEQ_ID,
                  struct options *opt);


static char *
get_filename(const char     *id,
             const char     *suffix,
             const char     *filename_default,
             struct options *opt);


/*--------------------------------------------------------------------------*/

void
init_default_options(struct options *opt)
{
  opt->filename_full = 0;
  opt->filename_delim = NULL;
  opt->pf = 0;
  opt->doT = 0; /* compute dimer free energies etc. */
  opt->noPS = 0;
  opt->noconv = 0;
  opt->bppmThreshold = 1e-5;
  opt->verbose = 0;
  opt->commands = NULL;
  opt->id_control = NULL;
  set_model_details(&(opt->md));

  opt->doC = 0; /* toggle to compute concentrations */
  opt->concentration_file     = NULL;

  opt->constraint_file = NULL;
  opt->constraint_batch = 0;
  opt->constraint_enforce = 0;
  opt->constraint_canonical = 0;

  opt->shape = 0;
  opt->shape_file = NULL;
  opt->shape_method = NULL;
  opt->shape_conversion = NULL;

  opt->csv_output = 0;  /* flag indicating whether we produce one-line outputs, a.k.a. CSV */
  opt->csv_header = 1;  /* print header for one-line output */
  opt->csv_output_delim = ',';  /* delimiting character for one-line output */
}


static char **
collect_unnamed_options(struct RNAcofold_args_info  *ggostruct,
                        int                         *num_files)
{
  char  **input_files = NULL;
  int   i;

  *num_files = 0;

  /* collect all unnamed options */
  if (ggostruct->inputs_num > 0) {
    input_files = (char **)vrna_realloc(input_files, sizeof(char *) * ggostruct->inputs_num);
    for (i = 0; i < ggostruct->inputs_num; i++)
      input_files[(*num_files)++] = strdup(ggostruct->inputs[i]);
  }

  return input_files;
}


int
main(int  argc,
     char *argv[])
{
  struct        RNAcofold_args_info args_info;
  char                      **input_files;
  int                       num_input;
  char                              *command_file;
  struct options                    opt;

  /*
   #############################################
   # init variables and parameter options
   #############################################
   */
  init_default_options(&opt);

  num_input = 0;

  command_file      = NULL;

  /*
   #############################################
   # check the command line prameters
   #############################################
   */
  if (RNAcofold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* get basic set of model details */
  ggo_get_md_eval(args_info, opt.md);
  ggo_get_md_fold(args_info, opt.md);
  ggo_get_md_part(args_info, opt.md);

  /* check dangle model */
  if ((opt.md.dangles < 0) || (opt.md.dangles > 3)) {
    vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    opt.md.dangles = dangles = 2;
  }

  /* SHAPE reactivity data */
  ggo_get_SHAPE(args_info, opt.shape, opt.shape_file, opt.shape_method, opt.shape_conversion);

  ggo_get_id_control(args_info, opt.id_control, "Sequence", "sequence", "_", 4, 1);

  ggo_get_constraints_settings(args_info,
                               fold_constrained,
                               opt.constraint_file,
                               opt.constraint_enforce,
                               opt.constraint_batch);

  /* enforce canonical base pairs in any case? */
  if (args_info.canonicalBPonly_given)
    opt.constraint_canonical = 1;

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (args_info.noconv_given)
    opt.noconv = 1;

  /*  */
  if (args_info.noPS_given)
    opt.noPS = 1;

  /* concentrations from stdin */
  if (args_info.concentrations_given)
    opt.doC = opt.doT = opt.pf = 1;

  /* set the bppm threshold for the dotplot */
  if (args_info.bppmThreshold_given)
    opt.bppmThreshold = MIN2(1., MAX2(0., args_info.bppmThreshold_arg));

  /* concentrations in file */
  if (args_info.concfile_given) {
    opt.concentration_file  = strdup(args_info.concfile_arg);
    opt.doC       = opt.doT = opt.pf = 1;
  }

  /* partition function settings */
  if (args_info.partfunc_given) {
    opt.pf = 1;
    if (args_info.partfunc_arg != -1)
      opt.md.compute_bpp = args_info.partfunc_arg;
  }

  if (args_info.all_pf_given) {
    opt.doT = opt.pf = 1;
    if (args_info.all_pf_arg != 1)
      opt.md.compute_bpp = args_info.all_pf_arg;
    else
      opt.md.compute_bpp = 1;
  }

  if (args_info.verbose_given)
    opt.verbose = 1;

  if (args_info.commands_given)
    command_file = strdup(args_info.commands_arg);

  /* filename sanitize delimiter */
  if (args_info.filename_delim_given)
    opt.filename_delim = strdup(args_info.filename_delim_arg);
  else if (get_id_delim(opt.id_control))
    opt.filename_delim = strdup(get_id_delim(opt.id_control));

  if ((opt.filename_delim) && isspace(*(opt.filename_delim))) {
    free(opt.filename_delim);
    opt.filename_delim = NULL;
  }

  /* full filename from FASTA header support */
  if (args_info.filename_full_given)
    opt.filename_full = 1;

  /* output format changes */
  if (args_info.output_format_given) {
    switch (*(args_info.output_format_arg)) {
      case 'D':
      /* fall-through */
      case 'd':
        opt.csv_output = 1;
        break;
      case 'V':
      /* fall-through */
      case 'v':
        opt.csv_output = 0;
        break;
      default:
        vrna_message_warning("unknown output format \"%c\", using defaults!",
                             *(args_info.output_format_arg));
        break;
    }
  }

  /* one-line output delimiter */
  if (args_info.csv_delim_given) {
    opt.csv_output_delim = *(args_info.csv_delim_arg);
    if (!opt.csv_output_delim) {
      vrna_message_warning("Delimiting character for One-Line output is missing, using defaults!");
      opt.csv_output_delim = ',';
    }
  }

  /* one-line output header */
  if (args_info.csv_noheader_given)
    opt.csv_header = 0;

  input_files = collect_unnamed_options(&args_info, &num_input);

  /* free allocated memory of command line data structure */
  RNAcofold_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */
  if (opt.pf && opt.md.gquad)
    vrna_message_error(
      "G-Quadruplex support is currently not available for partition function computations");

  if (command_file != NULL)
    opt.commands = vrna_file_commands_read(command_file, VRNA_CMD_PARSE_HC | VRNA_CMD_PARSE_SC);

  /*
   ################################################
   # process input files or handle input from stdin
   ################################################
   */
  if (num_input > 0) {
    int i, skip;
    for (skip = i = 0; i < num_input; i++) {
      if (!skip) {
        FILE *input_stream = fopen((const char *)input_files[i], "r");

        if (!input_stream)
          vrna_message_error("Unable to open %d. input file \"%s\" for reading", i + 1,
                             input_files[i]);

        if (opt.verbose) {
          vrna_message_info(stderr,
                            "Processing %d. input file \"%s\"",
                            i + 1,
                            input_files[i]);
        }

        if (process_input(input_stream, (const char *)input_files[i], &opt) == 0)
          skip = 1;

        fclose(input_stream);
      }

      free(input_files[i]);
    }
  } else {
    (void)process_input(stdin, NULL, &opt);
  }

  free(input_files);

  free(command_file);
  vrna_commands_free(opt.commands);

  free(opt.filename_delim);
  free(opt.concentration_file);

  free_id_data(opt.id_control);

  return EXIT_SUCCESS;
}


static int
process_input(FILE            *input_stream,
              const char      *input_filename,
              struct options  *opt)
{
  char                              *structure, *rec_sequence, *orig_sequence, *rec_id, **rec_rest,
                                    *tmp_string, *csv_output_string;
  unsigned int                      rec_type, read_opt;
  int i, length, cl, istty, istty_in, istty_out, ret;
  double                            min_en, kT, *concentrations;
  plist                             *prAB, *prAA, *prBB, *prA, *prB, *mfAB, *mfAA, *mfBB, *mfA,
                                    *mfB;

  ret               = 1;
  rec_type          = read_opt = 0;
  structure         = NULL;
  rec_id            = rec_sequence = orig_sequence = NULL;
  rec_rest          = NULL;

  csv_output_string = NULL; /* string holding output in case we produce one-line outputs */


  istty_in  = isatty(fileno(stdin));
  istty_out = isatty(fileno(stdout));
  istty     = isatty(fileno(stdout)) && isatty(fileno(stdin));

  /* print user help if we get input from tty */
  if (istty) {
    if (fold_constrained) {
      vrna_message_constraint_options(
        VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_X | VRNA_CONSTRAINT_DB_ANG_BRACK |
        VRNA_CONSTRAINT_DB_RND_BRACK);
      vrna_message_input_seq("Input sequence (upper or lower case) followed by structure constraint\n"
                             "Use '&' to connect 2 sequences that shall form a complex.");
    } else {
      vrna_message_input_seq("Use '&' to connect 2 sequences that shall form a complex.");
    }
  }

  /* set options we wanna pass to vrna_file_fasta_read_record() */
  if (istty)
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;

  if (!fold_constrained)
    read_opt |= VRNA_INPUT_NO_REST;

  if ((opt->csv_output) && (opt->csv_header))
    write_csv_header(opt);

  /*
   #############################################
   # main loop: continue until end of file
   #############################################
   */
  while (
    !((rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, NULL, read_opt))
      & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))) {
    /*
     ########################################################
     # init everything according to the data we've read
     ########################################################
     */
    char  *SEQ_ID         = NULL;
    int   maybe_multiline = 0;

    if (rec_id) {
      maybe_multiline = 1;
      /* remove '>' from FASTA header */
      rec_id = memmove(rec_id, rec_id + 1, strlen(rec_id));
    }

    /* construct the sequence ID */
    set_next_id(&rec_id, opt->id_control);
    SEQ_ID = fileprefix_from_id(rec_id, opt->id_control, opt->filename_full);

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if (!opt->noconv)
      vrna_seq_toRNA(rec_sequence);

    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(rec_sequence);

    vrna_fold_compound_t *vc = vrna_fold_compound(rec_sequence,
                                                  &(opt->md),
                                                  VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID);
    length    = vc->length;
    structure = (char *)vrna_alloc(sizeof(char) * (length + 1));

    /* parse the rest of the current dataset to obtain a structure constraint */
    if (fold_constrained) {
      if (opt->constraint_file) {
        vrna_constraints_add(vc, opt->constraint_file, VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID);
      } else {
        int           cp        = -1;
        char          *cstruc   = NULL;
        unsigned int  coptions  = (maybe_multiline) ? VRNA_OPTION_MULTILINE : 0;

        cstruc  = vrna_extract_record_rest_structure((const char **)rec_rest, 0, coptions);
        cstruc  = vrna_cut_point_remove(cstruc, &cp);
        if (vc->cutpoint != cp) {
          vrna_message_error("Sequence and Structure have different cut points.\n"
                             "sequence: %d, structure: %d",
                             vc->cutpoint, cp);
        }

        cl = (cstruc) ? (int)strlen(cstruc) : 0;

        if (cl == 0)
          vrna_message_warning("Structure constraint is missing");
        else if (cl < length)
          vrna_message_warning("Structure constraint is shorter than sequence");
        else if (cl > length)
          vrna_message_error("Structure constraint is too long");

        if (cstruc) {
          unsigned int constraint_options = VRNA_CONSTRAINT_DB_DEFAULT;

          if (opt->constraint_enforce)
            constraint_options |= VRNA_CONSTRAINT_DB_ENFORCE_BP;

          if (opt->constraint_canonical)
            constraint_options |= VRNA_CONSTRAINT_DB_CANONICAL_BP;

          vrna_constraints_add(vc, (const char *)cstruc, constraint_options);
        }

        free(cstruc);
      }
    }

    if (opt->shape) {
      vrna_constraints_add_SHAPE(vc,
                                 opt->shape_file,
                                 opt->shape_method,
                                 opt->shape_conversion,
                                 opt->verbose,
                                 VRNA_OPTION_DEFAULT| VRNA_OPTION_HYBRID);
    }

    if (opt->commands)
      vrna_commands_apply(vc, opt->commands, VRNA_CMD_PARSE_HC | VRNA_CMD_PARSE_SC);

    if (istty) {
      if (cut_point == -1) {
        vrna_message_info(stdout, "length = %d", length);
      } else {
        vrna_message_info(stdout,
                          "length1 = %d\nlength2 = %d",
                          cut_point - 1,
                          length - cut_point + 1);
      }
    }

    if (opt->doC) {
      FILE *fp;
      if (opt->concentration_file) {
        /* read from file */
        fp = fopen(opt->concentration_file, "r");
        if (fp == NULL)
          vrna_message_error("could not open concentration file %s", opt->concentration_file);

        concentrations = read_concentrations(fp);
        fclose(fp);
      } else {
        printf("Please enter concentrations [mol/l]\n format: ConcA ConcB\n return to end\n");
        concentrations = read_concentrations(stdin);
      }
    }

    /*
     ########################################################
     # begin actual computations
     ########################################################
     */

    /* compute mfe of AB dimer */
    min_en  = vrna_mfe_dimer(vc, structure);
    mfAB    = vrna_plist(structure, 0.95);

    /* check whether the constraint allows for any solution */
    if ((fold_constrained && opt->constraint_file) || (opt->commands)) {
      if (min_en == (double)(INF / 100.)) {
        vrna_message_error(
          "Supplied structure constraints create empty solution set for sequence:\n%s",
          orig_sequence);
        exit(EXIT_FAILURE);
      }
    }

    {
      char *pstruct, *msg = NULL;
      pstruct = vrna_cut_point_insert(structure, vc->cutpoint);

      if (opt->csv_output) {
        csv_output_string = vrna_strdup_printf(
          "%ld%c"                         /* sequence number */
          "%s%c"                          /* sequence id */
          "\"%s\"%c"                      /* sequence */
          "\"%s\"%c"                      /* MFE structure */
          "%6.2f",                        /* MFE */
          get_current_id(opt->id_control), opt->csv_output_delim,
          (rec_id) ? rec_id : "", opt->csv_output_delim,
          orig_sequence, opt->csv_output_delim,
          (pstruct) ? pstruct : "", opt->csv_output_delim,
          min_en);
      } else {
        print_fasta_header(stdout, rec_id);
        fprintf(stdout, "%s\n", orig_sequence);

        if (istty)
          msg = vrna_strdup_printf("\n minimum free energy = %6.2f kcal/mol", min_en);
        else
          msg = vrna_strdup_printf(" (%6.2f)", min_en);

        print_structure(stdout, pstruct, msg);
      }

      (void)fflush(stdout);

      if (!opt->noPS)
        postscript_layout(vc, orig_sequence, pstruct, SEQ_ID, opt);

      free(pstruct);
      free(msg);
    }

    if (length > 2000)
      vrna_mx_mfe_free(vc);

    /* compute partition function */
    if (opt->pf) {
      int     Blength, Alength;
      char    *Astring, *Bstring, *orig_Astring, *orig_Bstring;
      vrna_dimer_pf_t AB, AA, BB;
      vrna_dimer_conc_t *conc_result;

      prAB  = NULL;
      prAA  = NULL;
      prBB  = NULL;
      prA   = NULL;
      prB   = NULL;

      if (opt->md.dangles == 1) {
        vc->params->model_details.dangles = dangles = 2;   /* recompute with dangles as in pf_fold() */
        min_en                            = vrna_eval_structure(vc, structure);
        vc->params->model_details.dangles = dangles = 1;
      }

      vrna_exp_params_rescale(vc, &min_en);
      kT = vc->exp_params->kT / 1000.;

      if (length > 2000)
        vrna_message_info(stderr, "scaling factor %f", vc->exp_params->pf_scale);

      /* compute partition function */
      AB = vrna_pf_dimer(vc, structure);

      if (opt->md.compute_bpp)
        prAB = vrna_plist_from_probs(vc, opt->bppmThreshold);

      if (opt->doT) {
        if (vc->cutpoint <= 0) {
          vrna_message_warning(
            "Sorry, i cannot do that with only one molecule, please give me two or leave it");
          free(mfAB);
          free(prAB);
          continue;
        }

        if (opt->md.dangles == 1)
          opt->md.dangles = 2;


        Alength = vc->cutpoint - 1;                                 /* length of first molecule */
        Blength = length - vc->cutpoint + 1;                        /* length of 2nd molecule   */

        Astring = (char *)vrna_alloc(sizeof(char) * (Alength + 1)); /*Sequence of first molecule*/
        Bstring = (char *)vrna_alloc(sizeof(char) * (Blength + 1)); /*Sequence of second molecule*/
        strncat(Astring, rec_sequence, Alength);
        strncat(Bstring, rec_sequence + Alength + 1, Blength);

        orig_Astring  = (char *)vrna_alloc(sizeof(char) * (Alength + 1)); /*Sequence of first molecule*/
        orig_Bstring  = (char *)vrna_alloc(sizeof(char) * (Blength + 1)); /*Sequence of second molecule*/
        strncat(orig_Astring, orig_sequence, Alength);
        strncat(orig_Bstring, orig_sequence + Alength + 1, Blength);

        /* compute AA dimer */
        AA = do_partfunc(Astring, Alength, 2, &prAA, &mfAA, kT, opt);
        /* compute BB dimer */
        BB = do_partfunc(Bstring, Blength, 2, &prBB, &mfBB, kT, opt);
        /*free_co_pf_arrays();*/

        /* compute A monomer */
        do_partfunc(Astring, Alength, 1, &prA, &mfA, kT, opt);

        /* compute B monomer */
        do_partfunc(Bstring, Blength, 1, &prB, &mfB, kT, opt);

        if (opt->md.compute_bpp) {
          vrna_pf_dimer_probs(AB.F0AB, AB.FA, AB.FB, prAB, prA, prB, Alength, vc->exp_params);
          vrna_pf_dimer_probs(AA.F0AB, AA.FA, AA.FA, prAA, prA, prA, Alength, vc->exp_params);
          vrna_pf_dimer_probs(BB.F0AB, BB.FA, BB.FA, prBB, prA, prB, Blength, vc->exp_params);
        }
      }

      if (opt->doC)
        conc_result = vrna_pf_dimer_concentrations(AB.FcAB, AA.FcAB, BB.FcAB, AB.FA, AB.FB, concentrations, vc->exp_params);

      /* done with computations, let's produce output */
      if (opt->md.compute_bpp) {
        char *costruc, *msg = NULL;
        costruc = vrna_cut_point_insert(structure, vc->cutpoint);
        if (opt->csv_output) {
          msg = vrna_strdup_printf(
            "%s%c"            /* mfe output we've collected earlier */
            "\"%s\"%c"        /* pairing propensity */
            "%6.2f",          /* free energy of ensemble */
            csv_output_string, opt->csv_output_delim,
            costruc, opt->csv_output_delim,
            AB.FAB);
          free(csv_output_string);
          csv_output_string = msg;
        } else {
          if (istty_in)
            msg = vrna_strdup_printf("\n free energy of ensemble = %6.2f kcal/mol", AB.FAB);
          else
            msg = vrna_strdup_printf(" [%6.2f]", AB.FAB);

          print_structure(stdout, costruc, msg);
          free(msg);
        }

        free(costruc);
      } else if (opt->csv_output) {
        char *msg = vrna_strdup_printf(
          "%s%c"                  /* mfe output we've collected earlier */
          "%6.2f",                /* free energy of ensemble */
          csv_output_string, opt->csv_output_delim,
          AB.FAB);
        free(csv_output_string);
        csv_output_string = msg;
      } else {
        char *msg = vrna_strdup_printf(" free energy of ensemble = %6.2f kcal/mol", AB.FAB);
        print_structure(stdout, NULL, msg);
        free(msg);
      }

      if (opt->csv_output) {
        char *msg = vrna_strdup_printf(
          "%s%c"                /* output we've collected earlier */
          "%g%c"                /* probability of MFE structure */
          "%6.2f",              /* delta G binding */
          csv_output_string, opt->csv_output_delim,
          exp((AB.FAB - min_en) / kT), opt->csv_output_delim,
          AB.FcAB - AB.FA - AB.FB);
        free(csv_output_string);
        csv_output_string = msg;
      } else {
        char *msg = vrna_strdup_printf(" frequency of mfe structure in ensemble %g"
                                       "; delta G binding=%6.2f",
                                       exp((AB.FAB - min_en) / kT),
                                       AB.FcAB - AB.FA - AB.FB);
        print_structure(stdout, NULL, msg);
        free(msg);
      }

      if (opt->doT) {
       if (opt->csv_output) {
          char *msg = vrna_strdup_printf(
            "%s%c"                /* output we've collected earlier */
            "%.6f%c"              /* AB */
            "%.6f%c"              /* AA */
            "%.6f%c"              /* BB */
            "%.6f%c"              /* A */
            "%.6f",               /* B */
            csv_output_string, opt->csv_output_delim,
            AB.FcAB, opt->csv_output_delim,
            AA.FcAB, opt->csv_output_delim,
            BB.FcAB, opt->csv_output_delim,
            AB.FA, opt->csv_output_delim,
            AB.FB);
          free(csv_output_string);
          csv_output_string = msg;
        } else {
          print_comment(stdout, "Free Energies:");
          char *thead = NULL, *tline = NULL;
          thead = strdup("AB\t\tAA\t\tBB\t\tA\t\tB");
          tline = vrna_strdup_printf("%.6f\t%6f\t%6f\t%6f\t%6f",
                                     AB.FcAB, AA.FcAB, BB.FcAB, AB.FA, AB.FB);
          print_table(stdout, thead, tline);
          free(thead);
          free(tline);
        }
      }

      /* produce EPS dot plot(s) */
      if (opt->md.compute_bpp) {
        if (opt->doT) {
          char  *seq    = NULL;
          char  *comment      = NULL;
          char  *fname_dot    = NULL;
          char  *filename_dot = NULL;

          filename_dot  = get_filename(SEQ_ID, "dp5.ps", "dot5.ps", opt);

          /*AB dot_plot*/
          fname_dot = vrna_strdup_printf("AB%s", filename_dot);
          seq       = vrna_strdup_printf("%s%s", orig_Astring, orig_Bstring);
          comment   = vrna_strdup_printf("\n%%Heterodimer AB FreeEnergy= %.9f\n", AB.FcAB);
          (void)vrna_plot_dp_PS_list(seq, Alength + 1, fname_dot, prAB, mfAB, comment);
          free(comment);
          free(seq);
          free(fname_dot);

          /*AA dot_plot*/
          fname_dot = vrna_strdup_printf("AA%s", filename_dot);
          seq       = vrna_strdup_printf("%s%s", orig_Astring, orig_Astring);
          comment   = vrna_strdup_printf("\n%%Homodimer AA FreeEnergy= %.9f\n", AA.FcAB);
          (void)vrna_plot_dp_PS_list(seq, Alength + 1, fname_dot, prAA, mfAA, comment);
          free(comment);
          free(seq);
          free(fname_dot);

          /*BB dot_plot*/
          fname_dot = vrna_strdup_printf("BB%s", filename_dot);
          seq       = vrna_strdup_printf("%s%s", orig_Bstring, orig_Bstring);
          comment   = vrna_strdup_printf("\n%%Homodimer BB FreeEnergy= %.9f\n", BB.FcAB);
          (void)vrna_plot_dp_PS_list(seq, Blength + 1, fname_dot, prBB, mfBB, comment);
          free(comment);
          free(seq);
          free(fname_dot);

          /*A dot plot*/
          fname_dot = vrna_strdup_printf("A%s", filename_dot);
          comment   = vrna_strdup_printf("\n%%Monomer A FreeEnergy= %.9f\n", AB.FA);
          (void)vrna_plot_dp_PS_list(orig_Astring, -1, fname_dot, prA, mfA, comment);
          free(fname_dot);
          free(comment);

          /*B monomer dot plot*/
          fname_dot = vrna_strdup_printf("B%s", filename_dot);
          comment   = vrna_strdup_printf("\n%%Monomer B FreeEnergy= %.9f\n", AB.FB);
          (void)vrna_plot_dp_PS_list(orig_Bstring, -1, fname_dot, prB, mfB, comment);
          free(fname_dot);
          free(comment);

          free(filename_dot);
        } else {
          char  *filename_dot = get_filename(SEQ_ID, "dp.ps", "dot.ps", opt);

          if (filename_dot)
            (void)vrna_plot_dp_PS_list(orig_sequence, vc->cutpoint, filename_dot, prAB, mfAB, "doof");

          free(filename_dot);
        }
      }

      /* print concentrations table */
      if (opt->doC) {
        print_concentrations(conc_result, concentrations);
        free(conc_result);
        free(concentrations);
      }

      free(Astring);
      free(Bstring);
      free(orig_Astring);
      free(orig_Bstring);
      free(prAB);
      free(prAA);
      free(prBB);
      free(prA);
      free(prB);
      free(mfAB);
      free(mfAA);
      free(mfBB);
      free(mfA);
      free(mfB);

    }   /*end if(pf)*/

    if (opt->csv_output)
      fprintf(stdout, "%s\n", csv_output_string);

    free(csv_output_string);

    if (!opt->doT)
      vrna_mx_pf_free(vc);

    (void)fflush(stdout);

    /* clean up */
    free(rec_id);
    free(rec_sequence);
    free(orig_sequence);
    free(structure);
    /* free the rest of current dataset */
    if (rec_rest) {
      for (i = 0; rec_rest[i]; i++)
        free(rec_rest[i]);
      free(rec_rest);
    }

    rec_id    = rec_sequence = orig_sequence = structure = NULL;
    rec_rest  = NULL;
    vrna_fold_compound_free(vc);

    free(SEQ_ID);

    if (opt->shape || (opt->constraint_file && (!opt->constraint_batch))) {
      ret = 0;
      break;
    }

    /* print user help for the next round if we get input from tty */
    if (istty) {
      printf("Use '&' to connect 2 sequences that shall form a complex.\n");
      if (fold_constrained) {
        vrna_message_constraint_options(
          VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_X | VRNA_CONSTRAINT_DB_ANG_BRACK |
          VRNA_CONSTRAINT_DB_RND_BRACK);
        vrna_message_input_seq(
          "Input sequence (upper or lower case) followed by structure constraint\n");
      } else {
        vrna_message_input_seq_simple();
      }
    }
  }

  return ret;
}


static void
write_csv_header(struct options *opt)
{
  /* compose header line for CSV output */
  char *header = vrna_strdup_printf(
    "seq_num%c"
    "seq_id%c"
    "seq%c"
    "mfe_struct%c"
    "mfe",
    opt->csv_output_delim,
    opt->csv_output_delim,
    opt->csv_output_delim,
    opt->csv_output_delim);

  if (opt->pf) {
    if (opt->md.compute_bpp) {
      char *tmp = vrna_strdup_printf(
        "%s%c"
        "bpp_string%c"
        "ensemble_energy",
        header,
        opt->csv_output_delim,
        opt->csv_output_delim);
      free(header);
      header = tmp;
    } else {
      char *tmp = vrna_strdup_printf(
        "%s%c"
        "ensemble_energy",
        header,
        opt->csv_output_delim);
      free(header);
      header = tmp;
    }

    if (opt->doT) {
      char *tmp = vrna_strdup_printf(
        "%s%c"
        "AB%c"
        "AA%c"
        "BB%c"
        "A%c"
        "B",
        header,
        opt->csv_output_delim,
        opt->csv_output_delim,
        opt->csv_output_delim,
        opt->csv_output_delim,
        opt->csv_output_delim);
      free(header);
      header = tmp;
    }
  }

  fprintf(stdout, "%s\n", header);
  free(header);
}


static char *
get_filename(const char     *id,
             const char     *suffix,
             const char     *filename_default,
             struct options *opt)
{
  char *tmp_string, *filename = NULL;

  if (id) {
    filename = vrna_strdup_printf("%s%s%s", id, opt->filename_delim, suffix);
    /* sanitize file names */
    tmp_string = vrna_filename_sanitize(filename, opt->filename_delim);
    free(filename);
    filename = tmp_string;
  } else {
    filename = vrna_strdup_printf("%s", filename_default);
  }

  return filename;
}


PRIVATE vrna_dimer_pf_t
do_partfunc(char      *string,
            int       length,
            int       Switch,
            plist     **tpr,
            plist     **mfpl,
            double    kT,
            struct options *opt)
{
  /*compute mfe and partition function of dimer or monomer*/
  char                  *Newstring;
  char                  *tempstruc;
  double                min_en;
  vrna_md_t *md;
  vrna_dimer_pf_t       X;
  vrna_fold_compound_t  *vc;

  md  = &(opt->md);

  switch (Switch) {
    case 1:   /* monomer */
      tempstruc = (char *)vrna_alloc((unsigned)length + 1);
      vc        = vrna_fold_compound(string,
                                     md,
                                     VRNA_OPTION_MFE | VRNA_OPTION_PF);
      min_en  = vrna_mfe(vc, tempstruc);
      *mfpl   = vrna_plist(tempstruc, 0.95);
      vrna_mx_mfe_free(vc);

      vrna_exp_params_rescale(vc, &min_en);

      X = vrna_pf_dimer(vc, NULL);
      if (md->compute_bpp)
        *tpr = vrna_plist_from_probs(vc, opt->bppmThreshold);

      vrna_fold_compound_free(vc);
      free(tempstruc);
      break;

    case 2:   /* dimer */
      tempstruc = (char *)vrna_alloc((unsigned)length * 2 + 2);
      Newstring = (char *)vrna_alloc(sizeof(char) * (length * 2 + 2));
      strcat(Newstring, string);
      strcat(Newstring, "&");
      strcat(Newstring, string);

      vc = vrna_fold_compound(Newstring,
                              md,
                              VRNA_OPTION_MFE | VRNA_OPTION_PF | VRNA_OPTION_HYBRID);

      min_en  = vrna_mfe_dimer(vc, tempstruc);
      *mfpl   = vrna_plist(tempstruc, 0.95);
      vrna_mx_mfe_free(vc);

      vrna_exp_params_rescale(vc, &min_en);

      X = vrna_pf_dimer(vc, NULL);
      if (md->compute_bpp)
        *tpr = vrna_plist_from_probs(vc, opt->bppmThreshold);

      vrna_fold_compound_free(vc);

      free(Newstring);
      free(tempstruc);
      break;

    default:
      printf("Error in get_partfunc\n, computing neither mono- nor dimere!\n");
      exit(42);
  }

  return X;
}


void
postscript_layout(vrna_fold_compound_t *fc,
                  const char  *orig_sequence,
                  const char  *structure,
                  const char *SEQ_ID,
                  struct options *opt)
{
  char *filename_plot, *annot, *tmp_string;

  filename_plot = NULL;
  annot = NULL;

  if (SEQ_ID) {
    filename_plot = vrna_strdup_printf("%s%sss.ps", SEQ_ID, opt->filename_delim);
    tmp_string    = vrna_filename_sanitize(filename_plot, opt->filename_delim);
    free(filename_plot);
    filename_plot = tmp_string;
  } else {
    filename_plot = strdup("rna.ps");
  }

  if (fc->cutpoint >= 0) {
    annot = vrna_strdup_printf("1 %d 9  0 0.9 0.2 omark\n"
                               "%d %d 9  1 0.1 0.2 omark\n",
                               fc->cutpoint - 1,
                               fc->cutpoint + 1,
                               fc->length + 1);
  }

  if (filename_plot)
    (void)vrna_file_PS_rnaplot_a(orig_sequence, structure, filename_plot, annot, NULL, &(opt->md));

  free(filename_plot);
  free(annot);
}


PRIVATE void
print_concentrations(vrna_dimer_conc_t *result,
                     double            *startconc)
{
  /* compute and print concentrations out of free energies, calls get_concentrations */
  ;
  int               i, n;

  print_table(stdout, "Initial concentrations\t\trelative Equilibrium concentrations", NULL);
  print_table(stdout, "A\t\tB\t\tAB\t\tAA\t\tBB\t\tA\t\tB", NULL);
  for (n = 0; (startconc[2 * n] > 0) || (startconc[2 * n + 1] > 0); n++); /* count */
  for (i = 0; i < n; i++) {
    double  tot     = result[i].Ac_start + result[i].Bc_start;
    char    *tline  = vrna_strdup_printf("%-10g\t%-10g\t%.5f \t%.5f \t%.5f \t%.5f \t%.5f",
                                         result[i].Ac_start,
                                         result[i].Bc_start,
                                         result[i].ABc / tot,
                                         result[i].AAc / tot,
                                         result[i].BBc / tot,
                                         result[i].Ac / tot,
                                         result[i].Bc / tot);
    print_table(stdout, NULL, tline);
    free(tline);
  }
}


PRIVATE double *
read_concentrations(FILE *fp)
{
  /* reads concentrations, returns list of double, -1. marks end */
  char    *line;
  double  *startc;
  int     i = 0, n = 2;

  startc = (double *)vrna_alloc((2 * n + 1) * sizeof(double));

  while ((line = vrna_read_line(fp)) != NULL) {
    int c;
    if (i + 4 >= 2 * n) {
      n       *= 2;
      startc  = (double *)vrna_realloc(startc, (2 * n + 1) * sizeof(double));
    }

    c = sscanf(line, "%lf %lf", &startc[i], &startc[i + 1]);
    free(line);
    if (c < 2)
      break;

    i += 2;
  }
  startc[i] = startc[i + 1] = 0;
  return startc;
}
