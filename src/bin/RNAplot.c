/*
 * Plot RNA structures using different layout algorithms
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/plotting/probabilities.h"
#include "ViennaRNA/plotting/structures.h"
#include "ViennaRNA/plotting/alignments.h"
#include "ViennaRNA/plotting/utils.h"
#include "ViennaRNA/plotting/layouts.h"
#include "ViennaRNA/constraints/basic.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/io/file_formats_msa.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/utils/alignments.h"

#include "RNAplot_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"
#include "parallel_helpers.h"

#define PRIVATE static

struct options {
  unsigned int  msa_format;
  int           filename_full;
  char          *filename_delim;
  int           msa;
  int           mis;
  int           annotate_covar;
  int           aln_PS;
  int           aln_PS_cols;
  vrna_md_t     md;
  dataset_id    id_control;

  char          *pre;
  char          *post;
  char          format[5];
  int           plot_type;

  int           jobs;
  unsigned int  next_record_number;
  int           drawArcs;
  int           checkAncestorIntersections;
  int           checkSiblingIntersections;
  int           checkExteriorIntersections;
  int           allowFlipping;
  int           optimize;
};


struct record_data {
  unsigned int    number;
  char            *id;
  char            *sequence;
  char            *SEQ_ID;
  char            **rest;
  char            *input_filename;
  int             multiline_input;
  struct options  *options;
  int             tty;
};


struct record_data_msa {
  unsigned int    number;
  char            *MSA_ID;
  char            **alignment;
  char            **names;
  char            *structure;
  unsigned int    n_seq;
  struct options  *options;
  int             tty;
};


static int
process_input(FILE            *input_stream,
              const char      *filename,
              struct options  *opt);


static int
process_alignment_input(FILE            *input_stream,
                        const char      *filename,
                        struct options  *opt);


static void
process_record(struct record_data *record);


static void
process_alignment_record(struct record_data_msa *record);


void
init_default_options(struct options *opt)
{
  opt->msa_format     = VRNA_FILE_FORMAT_MSA_STOCKHOLM;
  opt->filename_full  = 0;
  opt->filename_delim = NULL;
  opt->msa            = 0;
  opt->mis            = 0;
  opt->annotate_covar = 0;
  opt->aln_PS         = 0;
  opt->aln_PS_cols    = 60;
  vrna_md_set_default(&(opt->md));

  opt->pre        = NULL;
  opt->post       = NULL;
  opt->plot_type  = 1;
  strncpy(&(opt->format[0]), "ps", 4);
  opt->format[4] = '\0';

  opt->jobs                       = 1;
  opt->next_record_number         = 0;
  opt->drawArcs                   = 1;
  opt->checkAncestorIntersections = 1;
  opt->checkSiblingIntersections  = 1;
  opt->checkExteriorIntersections = 1;
  opt->allowFlipping              = 0;
  opt->optimize                   = 1;
}


static char **
collect_unnamed_options(struct RNAplot_args_info  *ggostruct,
                        int                       *num_files)
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


static char **
append_input_files(struct RNAplot_args_info *ggostruct,
                   char                     **files,
                   int                      *numfiles)
{
  int i;

  if (ggostruct->infile_given) {
    files = (char **)vrna_realloc(files, sizeof(char *) * (*numfiles + ggostruct->infile_given));
    for (i = 0; i < ggostruct->infile_given; i++)
      files[(*numfiles)++] = strdup(ggostruct->infile_arg[i]);
  }

  return files;
}


int
main(int  argc,
     char *argv[])
{
  struct RNAplot_args_info  args_info;
  char                      **input_files;
  int                       num_input;
  struct  options           opt;

  num_input = 0;

  init_default_options(&opt);

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNAplot_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  if (args_info.msa_given) {
    opt.msa = 1;

    if (args_info.mis_given)
      opt.mis = 1;

    if (args_info.covar_given)
      opt.annotate_covar = 1;

    if (args_info.aln_given)
      opt.aln_PS = 1;

    if (args_info.aln_EPS_cols_given)
      opt.aln_PS_cols = args_info.aln_EPS_cols_arg;
  }

  if (opt.msa)
    ggo_get_id_control(args_info, opt.id_control, "Alignment", "alignment", "_", 4, 1);
  else
    ggo_get_id_control(args_info, opt.id_control, "Sequence", "sequence", "_", 4, 1);

  if (args_info.layout_type_given)
    rna_plot_type = args_info.layout_type_arg;

  if (args_info.pre_given)
    opt.pre = strdup(args_info.pre_arg);

  if (args_info.post_given)
    opt.post = strdup(args_info.post_arg);

  if (args_info.output_format_given) {
    strncpy(&(opt.format[0]), args_info.output_format_arg, 4);
    opt.format[4] = '\0';
  }

  /* filename sanitize delimiter */
  if (args_info.filename_delim_given)
    opt.filename_delim = strdup(args_info.filename_delim_arg);
  else if (get_id_delim(opt.id_control))
    opt.filename_delim = strdup(get_id_delim(opt.id_control));
  else
    opt.filename_delim = NULL;

  if ((opt.filename_delim) && isspace(*(opt.filename_delim))) {
    free(opt.filename_delim);
    opt.filename_delim = NULL;
  }

  /* full filename from FASTA header support */
  if (args_info.filename_full_given)
    opt.filename_full = 1;

  if (args_info.jobs_given) {
#if VRNA_WITH_PTHREADS
    int thread_max = max_user_threads();
    if (args_info.jobs_arg == 0) {
      /* use maximum of concurrent threads */
      int proc_cores, proc_cores_conf;
      if (num_proc_cores(&proc_cores, &proc_cores_conf)) {
        opt.jobs = MIN2(thread_max, proc_cores_conf);
      } else {
        vrna_message_warning("Could not determine number of available processor cores!\n"
                             "Defaulting to serial computation");
        opt.jobs = 1;
      }
    } else {
      opt.jobs = MIN2(thread_max, args_info.jobs_arg);
    }

    opt.jobs = MAX2(1, opt.jobs);
#else
    vrna_message_warning(
      "This version of RNAplot has been built without parallel input processing capabilities");
#endif
  }

  input_files = collect_unnamed_options(&args_info, &num_input);
  input_files = append_input_files(&args_info, input_files, &num_input);

  /* free allocated memory of command line data structure */
  RNAplot_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */

  int (*processing_func)(FILE           *stream,
                         const char     *filename,
                         struct options *opt);

  if (opt.msa)
    processing_func = &process_alignment_input;
  else
    processing_func = &process_input;

  /*
   #############################################
   # main loop: continue until end of file
   #############################################
   */
  INIT_PARALLELIZATION(opt.jobs);

  if (num_input > 0) {
    int i, skip;
    for (skip = i = 0; i < num_input; i++) {
      if (!skip) {
        FILE *input_stream = fopen((const char *)input_files[i], "r");

        if (!input_stream)
          vrna_message_error("Unable to open %d. input file \"%s\" for reading", i + 1,
                             input_files[i]);

        if (processing_func(input_stream, (const char *)input_files[i], &opt) == 0)
          skip = 1;

        fclose(input_stream);
      }

      free(input_files[i]);
    }
  } else {
    (void)processing_func(stdin, NULL, &opt);
  }

  UNINIT_PARALLELIZATION
  /*
   ################################################
   # post processing
   ################################################
   */
  free(opt.filename_delim);


  free(opt.pre);
  free(opt.post);

  free_id_data(opt.id_control);

  return EXIT_SUCCESS;
}


static int
process_input(FILE            *input_stream,
              const char      *input_filename,
              struct options  *opt)
{
  int           ret       = 1;
  int           istty_in  = isatty(fileno(input_stream));
  int           istty_out = isatty(fileno(stdout));

  unsigned int  read_opt = 0;

  /* set options we wanna pass to vrna_file_fasta_read_record() */
  if (istty_in && istty_out) {
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
    vrna_message_input_seq("Input sequence (upper or lower case) followed by structure");
  }

  /*
   #############################################
   # main loop: continue until end of file
   #############################################
   */
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

    if (rec_type & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))
      break;

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

    /* construct the sequence ID */
    set_next_id(&rec_id, opt->id_control);

    struct record_data *record = (struct record_data *)vrna_alloc(sizeof(struct record_data));

    record->number          = opt->next_record_number;
    record->sequence        = rec_sequence;
    record->SEQ_ID          = fileprefix_from_id(rec_id, opt->id_control, opt->filename_full);
    record->id              = rec_id;
    record->rest            = rec_rest;
    record->multiline_input = maybe_multiline;
    record->options         = opt;
    record->tty             = istty_in && istty_out;
    record->input_filename  = (input_filename) ? strdup(input_filename) : NULL;

    RUN_IN_PARALLEL(process_record, record);

    /* print user help for the next round if we get input from tty */
    if (istty_in && istty_out)
      vrna_message_input_seq("Input sequence (upper or lower case) followed by structure");
  } while (1);

  return ret;
}


static int
process_alignment_input(FILE            *input_stream,
                        const char      *input_filename,
                        struct options  *opt)
{
  int           ret           = 1;
  unsigned int  input_format  = opt->msa_format;
  int           istty_in      = isatty(fileno(input_stream));

  /* detect input file format if reading from file */
  if (input_filename) {
    unsigned int format_guess = vrna_file_msa_detect_format(input_filename, opt->msa_format);

    if (format_guess == VRNA_FILE_FORMAT_MSA_UNKNOWN) {
      char *msg = "Your input file is missing sequences! "
                  "Either your file is empty, or not in %s format!";

      switch (opt->msa_format) {
        case VRNA_FILE_FORMAT_MSA_CLUSTAL:
          vrna_message_error(msg, "Clustal");
          break;

        case VRNA_FILE_FORMAT_MSA_STOCKHOLM:
          vrna_message_error(msg, "Stockholm");
          break;

        case VRNA_FILE_FORMAT_MSA_FASTA:
          vrna_message_error(msg, "FASTA");
          break;

        case VRNA_FILE_FORMAT_MSA_MAF:
          vrna_message_error(msg, "MAF");
          break;

        default:
          vrna_message_error(msg, "Unknown");
          break;
      }
    }

    input_format = format_guess;
  }

  /* process input stream */
  while (!feof(input_stream)) {
    char  **alignment, **names, *tmp_id, *tmp_structure;
    int   n_seq;

    names         = NULL;
    alignment     = NULL;
    tmp_id        = NULL;
    tmp_structure = NULL;

    if (istty_in) {
      switch (input_format & (~VRNA_FILE_FORMAT_MSA_QUIET)) {
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

    input_format |= VRNA_FILE_FORMAT_MSA_QUIET;

    /* read record from input file */
    n_seq = vrna_file_msa_read_record(input_stream,
                                      &names,
                                      &alignment,
                                      &tmp_id,
                                      &tmp_structure,
                                      input_format);

    if (n_seq <= 0) {
      /* skip empty alignments */
      free(names);
      free(alignment);
      free(tmp_id);
      free(tmp_structure);
      names         = NULL;
      alignment     = NULL;
      tmp_id        = NULL;
      tmp_structure = NULL;
      continue;
    }

    /* prepare record data structure */
    struct record_data_msa *record =
      (struct record_data_msa *)vrna_alloc(sizeof(struct record_data_msa));

    record->number  = opt->next_record_number;
    record->MSA_ID  = fileprefix_from_id_alifold(tmp_id,
                                                 opt->id_control,
                                                 1);
    record->alignment = alignment;
    record->names     = names;
    record->structure = tmp_structure;
    record->n_seq     = (unsigned int)n_seq;

    record->tty = istty_in;

    record->options = opt;

    /* process the record we've just read */
    RUN_IN_PARALLEL(process_alignment_record, record);

    free(tmp_id);
  }

  return ret;
}


static void
process_record(struct record_data *record)
{
  struct options  *opt;
  char            *rec_sequence, *structure, *tmp_string, *ffname;

  opt           = record->options;
  rec_sequence  = strdup(record->sequence);
  ffname        = NULL;

  structure = vrna_extract_record_rest_structure(
    (const char **)record->rest,
    0,
    (record->multiline_input) ? VRNA_OPTION_MULTILINE : 0);

  if (!structure)
    vrna_message_error("structure missing for record %d\n", record->number);

  if (strlen(rec_sequence) != strlen(structure))
    vrna_message_error("sequence and structure have unequal length");

  if (record->SEQ_ID)
    ffname = vrna_strdup_printf("%s%sss", record->SEQ_ID, opt->filename_delim);
  else
    ffname = vrna_strdup_printf("rna");

  switch (opt->format[0]) {
    case 'p':
      tmp_string = vrna_strdup_printf("%s.ps", ffname);
      free(ffname);
      ffname = vrna_filename_sanitize(tmp_string, opt->filename_delim);
      free(tmp_string);

      if (rna_plot_type == VRNA_PLOT_TYPE_PUZZLER) {
        /* RNA puzzler behavior specification */
        vrna_plot_options_puzzler_t *puzzler;
        vrna_plot_layout_t          *layout;

        puzzler           = vrna_plot_options_puzzler();
        puzzler->filename = ffname;
        puzzler->drawArcs = 1;

        puzzler->checkAncestorIntersections = opt->checkAncestorIntersections;
        puzzler->checkSiblingIntersections  = opt->checkSiblingIntersections;
        puzzler->checkExteriorIntersections = opt->checkExteriorIntersections;
        puzzler->allowFlipping              = opt->allowFlipping;
        puzzler->optimize                   = opt->optimize;

        layout = vrna_plot_layout_puzzler(structure,
                                          puzzler);

        THREADSAFE_FILE_OUTPUT(
          vrna_file_PS_rnaplot_layout(rec_sequence,
                                      structure,
                                      ffname,
                                      opt->pre,
                                      opt->post,
                                      &(opt->md),
                                      layout));

        vrna_plot_layout_free(layout);
        vrna_plot_options_puzzler_free(puzzler);
      } else {
        THREADSAFE_FILE_OUTPUT(
          vrna_file_PS_rnaplot_a(rec_sequence,
                                 structure,
                                 ffname,
                                 opt->pre,
                                 opt->post,
                                 &(opt->md)));
      }

      break;

    case 'g':
      tmp_string = vrna_strdup_printf("%s.gml", ffname);
      free(ffname);
      ffname = vrna_filename_sanitize(tmp_string, opt->filename_delim);
      free(tmp_string);

      THREADSAFE_FILE_OUTPUT(
        gmlRNA(rec_sequence, structure, ffname, 'x')
        );
      break;

    case 'x':
      tmp_string = vrna_strdup_printf("%s.ss", ffname);
      free(ffname);
      ffname = vrna_filename_sanitize(tmp_string, opt->filename_delim);
      free(tmp_string);

      THREADSAFE_FILE_OUTPUT(
        xrna_plot(rec_sequence, structure, ffname)
        );
      break;

    case 's':
      tmp_string = vrna_strdup_printf("%s.svg", ffname);
      free(ffname);
      ffname = vrna_filename_sanitize(tmp_string, opt->filename_delim);
      free(tmp_string);

      THREADSAFE_FILE_OUTPUT(
        svg_rna_plot(rec_sequence, structure, ffname)
        );
      break;

    default:
      RNAplot_cmdline_parser_print_help();
      exit(EXIT_FAILURE);
  }

  /* clean up */
  free(record->id);
  free(record->SEQ_ID);
  free(record->sequence);
  free(rec_sequence);
  free(structure);
  free(ffname);

  /* free the rest of current dataset */
  if (record->rest) {
    for (int i = 0; record->rest[i]; i++)
      free(record->rest[i]);
    free(record->rest);
  }

  free(record->input_filename);

  free(record);
}


static void
process_alignment_record(struct record_data_msa *record)
{
  char            *consensus_sequence, *structure, *ffname, *tmp_string, **A,
                  *pre, *post;
  struct options  *opt;

  if (!record->structure)
    vrna_message_error("structure missing for record %d\n", record->number);

  opt       = record->options;
  structure = vrna_db_from_WUSS(record->structure);
  pre       = opt->pre;
  post      = opt->post;

  /*
   ########################################################
   # begin actual calculations
   ########################################################
   */
  /* generate consensus sequence */
  consensus_sequence = (opt->mis) ?
                       vrna_aln_consensus_mis((const char **)record->alignment, &(opt->md)) :
                       vrna_aln_consensus_sequence((const char **)record->alignment, &(opt->md));

  if (record->MSA_ID)
    ffname = vrna_strdup_printf("%s%sss", record->MSA_ID, opt->filename_delim);
  else
    ffname = vrna_strdup_printf("alirna");

  if (opt->annotate_covar) {
    A = vrna_annotate_covar_db((const char **)record->alignment, structure, &(opt->md));

    if (pre)
      pre = vrna_strdup_printf("%s\n%s", pre, A[0]);
    else
      pre = vrna_strdup_printf("%s", A[0]);

    if (post)
      post = vrna_strdup_printf("%s\n%s", post, A[1]);
    else
      post = vrna_strdup_printf("%s", A[1]);

    free(A[0]);
    free(A[1]);
    free(A);
  }

  switch (opt->format[0]) {
    case 'p':
      tmp_string = vrna_strdup_printf("%s.ps", ffname);
      free(ffname);
      ffname = vrna_filename_sanitize(tmp_string, opt->filename_delim);
      free(tmp_string);

      THREADSAFE_FILE_OUTPUT(
        vrna_file_PS_rnaplot_a(consensus_sequence,
                               structure,
                               ffname,
                               pre,
                               post,
                               &(opt->md)));

      break;

    case 'g':
      tmp_string = vrna_strdup_printf("%s.gml", ffname);
      free(ffname);
      ffname = vrna_filename_sanitize(tmp_string, opt->filename_delim);
      free(tmp_string);

      THREADSAFE_FILE_OUTPUT(
        gmlRNA(consensus_sequence, structure, ffname, 'x')
        );
      break;

    case 'x':
      tmp_string = vrna_strdup_printf("%s.ss", ffname);
      free(ffname);
      ffname = vrna_filename_sanitize(tmp_string, opt->filename_delim);
      free(tmp_string);

      THREADSAFE_FILE_OUTPUT(
        xrna_plot(consensus_sequence, structure, ffname)
        );
      break;

    case 's':
      tmp_string = vrna_strdup_printf("%s.svg", ffname);
      free(ffname);
      ffname = vrna_filename_sanitize(tmp_string, opt->filename_delim);
      free(tmp_string);

      THREADSAFE_FILE_OUTPUT(
        svg_rna_plot(consensus_sequence, structure, ffname)
        );
      break;

    default:
      RNAplot_cmdline_parser_print_help();
      exit(EXIT_FAILURE);
  }

  if (opt->aln_PS) {
    free(ffname);

    if (record->MSA_ID)
      ffname = vrna_strdup_printf("%s%saln.ps", record->MSA_ID, opt->filename_delim);
    else
      ffname = vrna_strdup_printf("aln.ps");

    THREADSAFE_FILE_OUTPUT(
      vrna_file_PS_aln(ffname,
                       (const char **)record->alignment,
                       (const char **)record->names,
                       structure,
                       opt->aln_PS_cols)
      );
  }

  /* clean up */
  if (pre != opt->pre)
    free(pre);

  if (post != opt->post)
    free(post);

  vrna_aln_free(record->alignment);
  vrna_aln_free(record->names);
  free(consensus_sequence);
  free(record->MSA_ID);
  free(record->structure);
  free(structure);
  free(ffname);

  free(record);
}
