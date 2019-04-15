/*
 *                Local version of RNAalifold
 *
 *                c Ivo L Hofacker, Stephan Bernhart
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/plotting/alignments.h"
#include "ViennaRNA/plotting/structures.h"
#include "ViennaRNA/plotting/utils.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/alifold.h"
#include "ViennaRNA/Lfold.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/params/io.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/io/file_formats_msa.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/constraints/basic.h"
#include "ViennaRNA/constraints/SHAPE.h"
#include "RNALalifold_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"

#include "ViennaRNA/color_output.inc"

#define DEFAULT_SPAN  70;

typedef struct {
  char      **names;
  char      **strings;
  char      **strings_orig;
  char      *prefix;
  int       columns;
  vrna_md_t *md;
  int       ss_eps;
  int       msa_eps;
  int       msa_stk;
  int       csv;
  int       mis;
  float     threshold;
  int       n_seq;
  int       dangle_model;
  int       split_energies;
  int       with_shapes;
} hit_data;


PRIVATE void
print_hit_cb(int        start,
             int        end,
             const char *structure,
             float      en,
             void       *data);


int
main(int  argc,
     char *argv[])
{
  FILE                          *clust_file;
  struct RNALalifold_args_info  args_info;
  char                          *string, *structure, *prefix, *tmp_string,
                                **AS, **names, *filename_in, *filename_delim, **input_files,
                                *tmp_id, *tmp_structure, **shape_files, *shape_method;
  unsigned int                  input_format_options, longest_string, aln_options;
  int                           n_seq, i, maxdist, unchangednc, unchangedcv, quiet, mis, istty,
                                alnPS, aln_columns, aln_out, ssPS, input_file_num, with_shapes,
                                *shape_file_association, verbose, s, tmp_number,
                                split_contributions;
  long int                      first_alignment_number;
  float                         e_max;
  vrna_md_t                     md;
  vrna_fold_compound_t          *fc;
  dataset_id                    id_control;

  clust_file              = stdin;
  string                  = structure = prefix = NULL;
  mis                     = 0;
  maxdist                 = DEFAULT_SPAN;
  do_backtrack            = unchangednc = unchangedcv = 1;
  ribo                    = 0;
  alnPS                   = 0;
  ssPS                    = 0;
  aln_out                 = 0;
  aln_columns             = 60;
  input_format_options    = VRNA_FILE_FORMAT_MSA_CLUSTAL; /* default to ClustalW format */
  aln_options             = VRNA_ALN_UPPERCASE;           /* we always require uppercase sequence letters internally */
  filename_in             = NULL;
  input_files             = NULL;
  input_file_num          = 0;
  shape_files             = NULL;
  shape_file_association  = 0;
  shape_method            = NULL;
  with_shapes             = 0;
  verbose                 = 0;
  quiet                   = 0;
  e_max                   = -0.1; /* threshold in kcal/mol per nucleotide in a hit */
  split_contributions     = 0;

  vrna_md_set_default(&md);

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNALalifold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* get basic set of model details */
  ggo_get_md_eval(args_info, md);
  ggo_get_md_fold(args_info, md);

  /* temperature */
  ggo_get_temperature(args_info, md.temperature);

  /* check dangle model */
  if ((md.dangles < 0) || (md.dangles > 3)) {
    vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    md.dangles = 2;
  }

  ggo_get_id_control(args_info, id_control, "Alignment", "alignment", "_", 4, 1);

  /* set cfactor */
  if (args_info.cfactor_given) {
    md.cv_fact  = args_info.cfactor_arg;
    unchangedcv = 0;
  }

  /* set nfactor */
  if (args_info.nfactor_given) {
    md.nc_fact  = args_info.nfactor_arg;
    unchangednc = 0;
  }

  /* calculate most informative sequence */
  if (args_info.mis_given)
    mis = 1;

  if (args_info.csv_given)
    csv = 1;

  if (args_info.aln_given) {
    alnPS = aln_out = ssPS = 1;
    if (args_info.aln_arg)
      prefix = strdup(args_info.aln_arg);
  }

  if (args_info.aln_EPS_ss_given)
    ssPS = 1;

  if (args_info.aln_EPS_given)
    alnPS = 1;

  if (args_info.aln_EPS_cols_given)
    aln_columns = args_info.aln_EPS_cols_arg;

  if (args_info.aln_stk_given) {
    aln_out = 1;
    if (args_info.aln_stk_arg) {
      if (prefix) {
        vrna_message_info(stdout,
                          "multiple output prefixes detected, using \"%s\"",
                          args_info.aln_stk_arg);
        free(prefix);
      }

      prefix = strdup(args_info.aln_stk_arg);
    }
  }

  if (args_info.ribosum_file_given) {
    RibosumFile = strdup(args_info.ribosum_file_arg);
    md.ribo     = ribo = 1;
  }

  if (args_info.ribosum_scoring_given) {
    RibosumFile = NULL;
    md.ribo     = ribo = 1;
  }

  if (args_info.verbose_given)
    verbose = 1;

  if (args_info.quiet_given) {
    if (verbose)
      vrna_message_warning(
        "Can not be verbose and quiet at the same time! I keep on being chatty...");
    else
      quiet = 1;
  }

  /* SHAPE reactivity data */
  if (args_info.shape_given) {
    if (verbose)
      vrna_message_info(stderr, "SHAPE reactivity data correction activated");

    with_shapes             = 1;
    shape_files             = (char **)vrna_alloc(sizeof(char *) * (args_info.shape_given + 1));
    shape_file_association  = (int *)vrna_alloc(sizeof(int *) * (args_info.shape_given + 1));

    /* find longest string in argument list */
    longest_string = 0;
    for (s = 0; s < args_info.shape_given; s++)
      if (strlen(args_info.shape_arg[s]) > longest_string)
        longest_string = strlen(args_info.shape_arg[s]);

    tmp_string  = (char *)vrna_alloc(sizeof(char) * (longest_string + 1));
    tmp_number  = 0;

    for (s = 0; s < args_info.shape_given; s++) {
      /* check whether we have int=string style that specifies a SHAPE file for a certain sequence number in the alignment */
      if (sscanf(args_info.shape_arg[s], "%d=%s", &tmp_number, tmp_string) == 2) {
        shape_files[s]            = strdup(tmp_string);
        shape_file_association[s] = tmp_number - 1;
      } else {
        shape_files[s]            = strdup(args_info.shape_arg[s]);
        shape_file_association[s] = s;
      }

      if (verbose) {
        vrna_message_info(stderr,
                          "Using SHAPE reactivity data provided in file %s for sequence %d",
                          shape_files[s],
                          shape_file_association[s] + 1);
      }
    }

    shape_file_association[s] = -1;

    free(tmp_string);
  }

  if (with_shapes)
    shape_method = strdup(args_info.shapeMethod_arg);

  /* alignment file name given as unnamed option? */
  if (args_info.inputs_num == 1) {
    filename_in = strdup(args_info.inputs[0]);
    clust_file  = fopen(filename_in, "r");
    if (clust_file == NULL) {
      vrna_message_warning("unable to open %s", filename_in);
      vrna_message_error("Input file can't be read!");
    }

    /*
     *  Use default alignment file formats.
     *  This may be overridden when we parse the
     *  --input-format parameter below
     */
    input_format_options = VRNA_FILE_FORMAT_MSA_DEFAULT;
  }

  /* get all input file name(s) */
  if (args_info.inputs_num > 0) {
    input_files = (char **)vrna_realloc(input_files, sizeof(char *) * args_info.inputs_num);
    for (i = 0; i < args_info.inputs_num; i++)
      input_files[input_file_num++] = strdup(args_info.inputs[i]);
  }

  if (args_info.input_format_given) {
    switch (args_info.input_format_arg[0]) {
      case 'C': /* ClustalW format */
        input_format_options = VRNA_FILE_FORMAT_MSA_CLUSTAL;
        break;

      case 'S': /* Stockholm 1.0 format */
        input_format_options = VRNA_FILE_FORMAT_MSA_STOCKHOLM;
        break;

      case 'F': /* FASTA format */
        input_format_options = VRNA_FILE_FORMAT_MSA_FASTA;
        break;

      case 'M': /* MAF format */
        input_format_options = VRNA_FILE_FORMAT_MSA_MAF;
        break;

      default:
        vrna_message_warning("Unknown input format specified");
        break;
    }
  }

  /* filename sanitize delimiter */
  if (args_info.filename_delim_given)
    filename_delim = strdup(args_info.filename_delim_arg);
  else if (get_id_delim(id_control))
    filename_delim = strdup(get_id_delim(id_control));
  else
    filename_delim = NULL;

  if ((filename_delim) && isspace(*filename_delim)) {
    free(filename_delim);
    filename_delim = NULL;
  }

  if (args_info.threshold_given)
    e_max = (float)args_info.threshold_arg;

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (!(args_info.noconv_given))
    aln_options |= VRNA_ALN_RNA;

  if (args_info.split_contributions_given)
    split_contributions = 1;

  /* free allocated memory of command line data structure */
  RNALalifold_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */
  if ((ribo == 1) && (unchangednc))
    md.nc_fact = 0.5;

  if ((ribo == 1) && (unchangedcv))
    md.cv_fact = 0.6;

  istty = isatty(fileno(stdout)) && isatty(fileno(stdin));

  /* check whether the user wants a non-default base pair span */
  if (md.max_bp_span != -1)
    maxdist = md.window_size = md.max_bp_span;
  else
    md.max_bp_span = md.window_size = maxdist;

  /*
   #############################################
   # begin calculations
   #############################################
   */
  if (filename_in) {
    unsigned int format_guess = vrna_file_msa_detect_format(filename_in, input_format_options);
    if (format_guess == VRNA_FILE_FORMAT_MSA_UNKNOWN) {
      char *format = NULL;
      switch (input_format_options) {
        case VRNA_FILE_FORMAT_MSA_CLUSTAL:
          format = strdup("Clustal");
          break;
        case VRNA_FILE_FORMAT_MSA_STOCKHOLM:
          format = strdup("Stockholm");
          break;
        case VRNA_FILE_FORMAT_MSA_FASTA:
          format = strdup("FASTA");
          break;
        case VRNA_FILE_FORMAT_MSA_MAF:
          format = strdup("MAF");
          break;
        default:
          format = strdup("Unknown");
          break;
      }
      vrna_message_error(
        "Your input file is missing sequences! Either your file is empty, or not in %s format!",
        format);
      free(format);
    }

    input_format_options = format_guess;
  }

  first_alignment_number = get_current_id(id_control);

  while (!feof(clust_file)) {
    char  **MSA_orig  = NULL;
    char  *MSA_ID     = NULL;
    fflush(stdout);
    if (istty && (clust_file == stdin)) {
      switch (input_format_options & (~VRNA_FILE_FORMAT_MSA_QUIET)) {
        case VRNA_FILE_FORMAT_MSA_CLUSTAL:
          vrna_message_input_msa("Input aligned sequences in ClustalW format\n"
                                 "(press Ctrl+d when finished to indicate the end of your input)");
          break;

        case VRNA_FILE_FORMAT_MSA_STOCKHOLM:
          vrna_message_input_msa("Input aligned sequences in Stockholm format (Insert one alignment at a time!)\n"
                                 "(Note, Stockholm entries always end with a line that only contains '//'");
          break;

        case VRNA_FILE_FORMAT_MSA_FASTA:
          vrna_message_input_msa("Input aligned sequences in FASTA format\n"
                                 "(press Ctrl+d when finished to indicate the end of your input)");
          break;

        case VRNA_FILE_FORMAT_MSA_MAF:
          vrna_message_input_msa("Input aligned sequences in MAF format (Insert one alignment at a time!)\n"
                                 "(Note, a MAF alignment always ends with an empty line)");
          break;

        default:
          vrna_message_error("Which input format are you using?");
          break;
      }
    }

    if (quiet)
      input_format_options |= VRNA_FILE_FORMAT_MSA_QUIET;

    /* read the first record from input file */
    n_seq = vrna_file_msa_read_record(clust_file,
                                      &names,
                                      &AS,
                                      &tmp_id,
                                      &tmp_structure,
                                      input_format_options);
    fflush(stdout);
    fflush(stderr);

    if (n_seq <= 0) {
      /* skip empty alignments */
      free(names);
      free(AS);
      free(tmp_id);
      free(tmp_structure);
      names         = NULL;
      AS            = NULL;
      tmp_id        = NULL;
      tmp_structure = NULL;
      continue;
    }

    /* construct alignment ID */
    MSA_ID = fileprefix_from_id_alifold(tmp_id, id_control, 0);

    print_fasta_header(stdout, MSA_ID);

    /*
     *  store original alignment and create a new one for internal computations
     *  which require all-uppercase letters, and RNA/DNA alphabet
     */
    MSA_orig  = AS;
    AS        = vrna_aln_copy((const char **)MSA_orig, aln_options);
    fc        = vrna_fold_compound_comparative((const char **)AS,
                                               &md,
                                               VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

    hit_data data;

    data.names          = names;
    data.strings        = AS;
    data.strings_orig   = MSA_orig;
    data.prefix         = MSA_ID;
    data.columns        = aln_columns;
    data.md             = &md;
    data.ss_eps         = ssPS;
    data.msa_eps        = alnPS;
    data.msa_stk        = aln_out;
    data.csv            = csv;
    data.mis            = mis;
    data.threshold      = e_max;
    data.n_seq          = n_seq;
    data.dangle_model   = md.dangles;
    data.split_energies = split_contributions;
    data.with_shapes    = with_shapes;

    if (with_shapes) {
      for (s = 0; shape_file_association[s] != -1; s++);

      if ((s != n_seq) && (!quiet))
        vrna_message_warning(
          "Number of sequences in alignment does not match number of provided SHAPE reactivity data files! ");

      shape_files             = (char **)vrna_realloc(shape_files, (n_seq + 1) * sizeof(char *));
      shape_file_association  = (int *)vrna_realloc(shape_file_association,
                                                    (n_seq + 1) * sizeof(int));
    }

    if (with_shapes) {
      vrna_constraints_add_SHAPE_ali(fc, \
                                     shape_method, \
                                     (const char **)shape_files, \
                                     shape_file_association, \
                                     verbose, \
                                     VRNA_OPTION_MFE);
    }

    (void)vrna_mfe_window_cb(fc, &print_hit_cb, (void *)&data);

    string =
      (mis) ? vrna_aln_consensus_mis((const char **)AS,
                                     &md) : vrna_aln_consensus_sequence((const char **)AS, &md);
    printf("%s\n", string);

    (void)fflush(stdout);
    free(string);
    free(structure);

    vrna_aln_free(AS);
    vrna_aln_free(MSA_orig);
    vrna_aln_free(names);

    vrna_fold_compound_free(fc);
    free(MSA_ID);
    free(tmp_id);
    free(tmp_structure);

    /* break after first record if constraint folding and not explicitly instructed otherwise */
    if (with_shapes)
      break;
  } /* end of input */

  if (first_alignment_number == get_current_id(id_control)) {
    char *format = NULL;
    switch (input_format_options) {
      case VRNA_FILE_FORMAT_MSA_CLUSTAL:
        format = strdup("Clustal");
        break;
      case VRNA_FILE_FORMAT_MSA_STOCKHOLM:
        format = strdup("Stockholm");
        break;
      case VRNA_FILE_FORMAT_MSA_FASTA:
        format = strdup("FASTA");
        break;
      case VRNA_FILE_FORMAT_MSA_MAF:
        format = strdup("MAF");
        break;
      default:
        format = strdup("Unknown");
        break;
    }
    vrna_message_error(
      "Your input file is missing sequences! Either your file is empty, or not in %s format!",
      format);
    free(format);
  }

  if (clust_file != stdin)
    fclose(clust_file);

  free(shape_files);

  free(filename_in);

  for (i = 0; i < input_file_num; i++)
    free(input_files[i]);
  free(input_files);

  free(filename_delim);

  free_id_data(id_control);

  return EXIT_SUCCESS;
}


PRIVATE void
print_hit_cb(int        start,
             int        end,
             const char *structure,
             float      en,
             void       *data)
{
  char      **sub;
  char      *fname, *tmp_string;
  char      *id;
  char      **A, *cons;
  char      **names;
  char      **strings;
  char      **strings_orig;
  char      *prefix;
  char      *msg, *ss;
  int       columns, dangle_model;
  vrna_md_t *md;
  int       with_ss, with_msa, with_stk, with_csv, with_mis, with_shapes, split_energies;
  float     threshold;

  names           = ((hit_data *)data)->names;
  strings         = ((hit_data *)data)->strings;
  strings_orig    = ((hit_data *)data)->strings_orig;
  prefix          = ((hit_data *)data)->prefix;
  columns         = ((hit_data *)data)->columns;
  md              = ((hit_data *)data)->md;
  with_ss         = ((hit_data *)data)->ss_eps;
  with_msa        = ((hit_data *)data)->msa_eps;
  with_stk        = ((hit_data *)data)->msa_stk;
  with_csv        = ((hit_data *)data)->csv;
  with_mis        = ((hit_data *)data)->mis;
  threshold       = ((hit_data *)data)->threshold;
  dangle_model    = ((hit_data *)data)->dangle_model;
  split_energies  = ((hit_data *)data)->split_energies;

  if ((dangle_model == 2) && (start > 1)) {
    ss = vrna_strdup_printf(".%s", structure);
    start--;
  } else {
    ss = vrna_strdup_printf("%s", structure);
  }

  if ((en / (float)(end - start + 1)) <= threshold) {
    sub   = vrna_aln_slice((const char **)strings, (unsigned int)start, (unsigned int)end);
    cons  =
      (with_mis) ? vrna_aln_consensus_mis((const char **)sub,
                                          md) : vrna_aln_consensus_sequence((const char **)sub, md);
    A = vrna_annotate_covar_db((const char **)sub, ss, md);


    if (split_energies) {
      vrna_fold_compound_t  *sub_fc = vrna_fold_compound_comparative((const char **)sub,
                                                                     md,
                                                                     VRNA_OPTION_DEFAULT);
      float                 real_en = vrna_eval_structure(sub_fc, ss);
      with_shapes = ((hit_data *)data)->with_shapes;
      if (with_shapes) {
        float cov_en = vrna_eval_covar_structure(sub_fc, ss);
        if (with_csv == 1) {
          msg = vrna_strdup_printf(",%6.2f,%6.2f,%6.2f,%d,%d",
                                   real_en,
                                   -cov_en,
                                   en - real_en + cov_en,
                                   start,
                                   end);
        } else {
          msg = vrna_strdup_printf(" (%6.2f = %6.2f + %6.2f + %6.2f) %4d - %4d",
                                   en,
                                   real_en,
                                   -cov_en,
                                   en - real_en + cov_en,
                                   start,
                                   end);
        }
      } else {
        if (with_csv == 1) {
          msg = vrna_strdup_printf(",%6.2f,%6.2f,%d,%d", real_en, en - real_en, start, end);
        } else {
          msg = vrna_strdup_printf(" (%6.2f = %6.2f + %6.2f) %4d - %4d",
                                   en,
                                   real_en,
                                   en - real_en,
                                   start,
                                   end);
        }
      }

      vrna_fold_compound_free(sub_fc);
    } else {
      if (with_csv == 1)
        msg = vrna_strdup_printf(",%6.2f,%d,%d", en, start, end);
      else
        msg = vrna_strdup_printf(" (%6.2f) %4d - %4d", en, start, end);
    }

    print_structure(stdout, ss, msg);
    free(msg);

    if (with_ss) {
      if (prefix)
        fname = vrna_strdup_printf("%s_ss_%d_%d.eps", prefix, start, end);
      else
        fname = vrna_strdup_printf("ss_%d_%d.eps", start, end);

      tmp_string = vrna_filename_sanitize(fname, "_");
      free(fname);
      fname = tmp_string;

      (void)vrna_file_PS_rnaplot_a(cons, ss, fname, A[0], A[1], md);
      free(fname);
    }

    free(A[0]);
    free(A[1]);
    free(A);
    free(cons);

    if (with_msa) {
      if (prefix)
        fname = vrna_strdup_printf("%s_aln_%d_%d.eps", prefix, start, end);
      else
        fname = vrna_strdup_printf("aln_%d_%d.eps", start, end);

      tmp_string = vrna_filename_sanitize(fname, "_");

      vrna_file_PS_aln_slice(tmp_string,
                             (const char **)strings_orig,
                             (const char **)names,
                             ss - start + 1,
                             start,
                             end,
                             0,
                             columns);
      free(fname);
      free(tmp_string);
    }

    if (with_stk) {
      char          **sub_orig = vrna_aln_slice((const char **)strings_orig,
                                                (unsigned int)start,
                                                (unsigned int)end);
      unsigned int  options = VRNA_FILE_FORMAT_MSA_STOCKHOLM |
                              VRNA_FILE_FORMAT_MSA_APPEND;
      if (with_mis)
        options |= VRNA_FILE_FORMAT_MSA_MIS;

      if (prefix) {
        id    = vrna_strdup_printf("%s_aln_%d_%d", prefix, start, end);
        fname = vrna_strdup_printf("%s.stk", prefix);

        tmp_string = vrna_filename_sanitize(fname, "_");
        free(fname);
        fname = tmp_string;

        vrna_file_msa_write(fname,
                            (const char **)names,
                            (const char **)sub_orig,
                            id,
                            ss,
                            "RNALalifold prediction",
                            options);
        free(fname);
      } else {
        id = vrna_strdup_printf("aln_%d_%d", start, end);
        vrna_file_msa_write("RNALalifold_results.stk",
                            (const char **)names,
                            (const char **)sub_orig,
                            id,
                            ss,
                            "RNALalifold prediction",
                            options);
      }

      free(id);
      vrna_aln_free(sub_orig);
    }

    vrna_aln_free(sub);
  }

  free(ss);
}
