/*
 *
 *        Calculate Energy of given Sequences and Structures
 *                         c Ivo L Hofacker
 *                        Vienna RNA Pckage
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/params/io.h"
#include "ViennaRNA/constraints/basic.h"
#include "ViennaRNA/constraints/SHAPE.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/io/file_formats_msa.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/datastructures/char_stream.h"
#include "ViennaRNA/datastructures/stream_output.h"
#include "ViennaRNA/color_output.inc"

#include "RNAeval_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"
#include "parallel_helpers.h"

#define DBL_ROUND(a, digits) (round((a) * pow(10., (double)(digits))) / pow(10., (double)(digits)))

struct options {
  unsigned int    msa_format;
  int             filename_full;
  int             noconv;
  int             verbose;
  int             aln;
  int             mis;
  vrna_md_t       md;
  dataset_id      id_control;

  int             shape;
  char            *shape_file;
  char            *shape_method;
  char            *shape_conversion;

  int             jobs;
  int             keep_order;
  unsigned int    next_record_number;
  vrna_ostream_t  output_queue;
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


struct output_stream {
  vrna_cstr_t data;
  vrna_cstr_t err;
};


static int
process_input(FILE            *input_stream,
              const char      *input_filename,
              struct options  *opt);


static int
process_alignment_input(FILE            *input_stream,
                        const char      *input_filename,
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
  opt->noconv         = 0;
  opt->verbose        = 0;
  opt->aln            = 0;
  opt->mis            = 0;
  vrna_md_set_default(&(opt->md));

  opt->shape            = 0;
  opt->shape_file       = NULL;
  opt->shape_method     = NULL;
  opt->shape_conversion = NULL;

  opt->jobs               = 1;
  opt->keep_order         = 1;
  opt->next_record_number = 0;
  opt->output_queue       = NULL;
}


static char **
collect_unnamed_options(struct RNAeval_args_info  *ggostruct,
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
append_input_files(struct RNAeval_args_info *ggostruct,
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


void
flush_cstr_callback(void          *auxdata,
                    unsigned int  i,
                    void          *data)
{
  struct output_stream *s = (struct output_stream *)data;

  /* flush/free errors first */
  vrna_cstr_free(s->err);

  /* flush/free data[k] */
  vrna_cstr_free(s->data);

  free(s);
}


int
main(int  argc,
     char *argv[])
{
  struct RNAeval_args_info  args_info;
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
  if (RNAeval_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* get basic set of model details */
  ggo_get_md_eval(args_info, opt.md);
  ggo_get_circ(args_info, opt.md.circ);

  /* temperature */
  ggo_get_temperature(args_info, opt.md.temperature);

  /* check dangle model */
  if ((opt.md.dangles < 0) || (opt.md.dangles > 3)) {
    vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    opt.md.dangles = dangles = 2;
  }

  /* SHAPE reactivity data */
  ggo_get_SHAPE(args_info, opt.shape, opt.shape_file, opt.shape_method, opt.shape_conversion);

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (args_info.noconv_given)
    opt.noconv = 1;

  /* logarithmic multiloop energies */
  if (args_info.logML_given)
    opt.md.logML = logML = 1;

  /* be verbose */
  if (args_info.verbose_given)
    opt.verbose = 1;

  if (args_info.msa_given) {
    opt.aln = 1;

    /* set cfactor */
    if (args_info.cfactor_given)
      opt.md.cv_fact = cv_fact = args_info.cfactor_arg;

    /* set nfactor */
    if (args_info.nfactor_given)
      opt.md.nc_fact = nc_fact = args_info.nfactor_arg;

    if (args_info.mis_given)
      opt.mis = 1;

    if (args_info.old_given)
      opt.md.oldAliEn = 1;

    if (args_info.ribosum_file_given) {
      RibosumFile = strdup(args_info.ribosum_file_arg);
      opt.md.ribo = 1;
    }

    if (args_info.ribosum_scoring_given) {
      RibosumFile = NULL;
      opt.md.ribo = 1;
    }
  }

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
      "This version of RNAeval has been built without parallel input processing capabilities");
#endif

    if (args_info.unordered_given)
      opt.keep_order = 0;
  }

  input_files = collect_unnamed_options(&args_info, &num_input);
  input_files = append_input_files(&args_info, input_files, &num_input);

  if (opt.aln)
    ggo_get_id_control(args_info, opt.id_control, "Alignment", "alignment", "_", 4, 1);
  else
    ggo_get_id_control(args_info, opt.id_control, "Sequence", "sequence", "_", 4, 1);

  /* free allocated memory of command line data structure */
  RNAeval_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */
  if (opt.md.circ && opt.md.gquad)
    vrna_message_error("G-Quadruplex support is currently not available for circular RNA structures");

  if (opt.keep_order)
    opt.output_queue = vrna_ostream_init(&flush_cstr_callback, NULL);

  int (*processing_func)(FILE           *stream,
                         const char     *filename,
                         struct options *opt);

  if (opt.aln)
    processing_func = &process_alignment_input;
  else
    processing_func = &process_input;

  /*
   ################################################
   # process input files or handle input from stdin
   ################################################
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
  vrna_ostream_free(opt.output_queue);


  free(input_files);
  free(opt.shape_file);
  free(opt.shape_method);
  free(opt.shape_conversion);

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

  /* print user help if we get input from tty */
  if (istty_in && istty_out) {
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
    vrna_message_input_seq("Use '&' to connect 2 sequences that shall form a complex.\n"
                           "Input sequence (upper or lower case) followed by structure");
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

    if (opt->output_queue)
      vrna_ostream_request(opt->output_queue, opt->next_record_number++);

    RUN_IN_PARALLEL(process_record, record);

    if (opt->shape) {
      ret = 0;
      break;
    }

    /* print user help for the next round if we get input from tty */
    if (istty_in && istty_out)
      vrna_message_input_seq("Use '&' to connect 2 sequences that shall form a complex.\n"
                             "Input sequence (upper or lower case) followed by structure");
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

    if (opt->output_queue)
      vrna_ostream_request(opt->output_queue, opt->next_record_number++);

    /* process the record we've just read */
    RUN_IN_PARALLEL(process_alignment_record, record);

    free(tmp_id);

    if (opt->shape) {
      ret = 0;
      break;
    }
  }

  return ret;
}


static void
process_record(struct record_data *record)
{
  struct options        *opt;
  struct output_stream  *o_stream;
  char                  *rec_sequence, *structure, *tmp;
  int                   n;
  float                 energy;
  vrna_fold_compound_t  *vc;

  opt           = record->options;
  o_stream      = (struct output_stream *)vrna_alloc(sizeof(struct output_stream));
  rec_sequence  = strdup(record->sequence);

  /* convert DNA alphabet to RNA if not explicitely switched off */
  if (!opt->noconv) {
    vrna_seq_toRNA(rec_sequence);
    vrna_seq_toRNA(record->sequence);
  }

  /* convert sequence to uppercase letters only */
  vrna_seq_toupper(rec_sequence);

  vc = vrna_fold_compound(rec_sequence,
                          &(opt->md),
                          VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY);

  if (!vc) {
    vrna_message_warning("Skipping computations for \"%s\"",
                         (record->id) ? record->id : "identifier unavailable");
    return;
  }

  n = (int)vc->length;

  if (opt->shape) {
    vrna_constraints_add_SHAPE(vc,
                               opt->shape_file,
                               opt->shape_method,
                               opt->shape_conversion,
                               opt->verbose,
                               VRNA_OPTION_DEFAULT);
  }

  /* retrieve string stream bound to stdout, 6*length should be enough memory to start with */
  o_stream->data = vrna_cstr(6 * n, stdout);
  /* retrieve string stream bound to stderr for any info messages */
  o_stream->err = vrna_cstr(n, stderr);

  tmp = vrna_extract_record_rest_structure((const char **)record->rest,
                                           0,
                                           (record->multiline_input) ? VRNA_OPTION_MULTILINE : 0);

  if (!tmp)
    vrna_message_error("structure missing for record %d\n", record->number);

  {
    int cp = -1;
    structure = vrna_cut_point_remove(tmp, &cp);
    if (cp != vc->cutpoint) {
      vrna_message_warning("cut_point = %d cut = %d", vc->cutpoint, cp);
      vrna_message_error("Sequence and Structure have different cut points.");
    }

    n = (int)strlen(structure);
    if (n != vc->length)
      vrna_message_error("structure and sequence differ in length!");

    free(tmp);
  }

  if (record->tty) {
    if (vc->cutpoint == -1) {
      vrna_message_info(stdout, "length = %d", n);
    } else {
      vrna_message_info(stdout,
                        "length1 = %d\nlength2 = %d",
                        vc->cutpoint - 1,
                        n - vc->cutpoint + 1);
    }
  }

  /*
   ########################################################
   # begin actual computations
   ########################################################
   */
  vrna_cstr_print_fasta_header(o_stream->data, record->id);

  energy = vrna_eval_structure_cstr(vc, structure, opt->verbose, o_stream->data);
  vrna_cstr_printf(o_stream->data, "%s\n", record->sequence);

  char *pstruct = vrna_cut_point_insert(structure, vc->cutpoint);
  vrna_cstr_printf_structure(o_stream->data,
                             pstruct,
                             record->tty ? "\n energy = %6.2f kcal/mol" : " (%6.2f)",
                             energy);
  free(pstruct);

  if (opt->output_queue)
    vrna_ostream_provide(opt->output_queue, record->number, (void *)o_stream);
  else
    flush_cstr_callback(NULL, 0, (void *)o_stream);

  /* clean up */
  vrna_fold_compound_free(vc);
  free(record->id);
  free(record->SEQ_ID);
  free(record->sequence);
  free(rec_sequence);
  free(structure);

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
  char                  **alignment, *consensus_sequence, *structure;
  unsigned int          n;
  double                real_en, cov_en;
  struct options        *opt;
  vrna_fold_compound_t  *vc;
  struct output_stream  *o_stream;

  o_stream = (struct output_stream *)vrna_alloc(sizeof(struct output_stream));

  if (!record->structure)
    vrna_message_error("structure missing for record %d\n", record->number);

  opt       = record->options;
  structure = vrna_db_from_WUSS(record->structure);

  /*
   *  create a new alignment for internal computations
   *  which require all-uppercase letters, and RNA/DNA alphabet
   */
  alignment = vrna_aln_copy((const char **)record->alignment,
                            VRNA_ALN_UPPERCASE |
                            (opt->noconv ? 0 : VRNA_ALN_RNA));

  vc = vrna_fold_compound_comparative((const char **)alignment,
                                      &(opt->md),
                                      VRNA_OPTION_DEFAULT | VRNA_OPTION_EVAL_ONLY);

  if (!vc) {
    vrna_message_warning("Skipping computations for \"%s\"",
                         (record->MSA_ID) ? record->MSA_ID : "identifier unavailable");
    return;
  }

  n = vc->length;

  /* retrieve string stream bound to stdout, 6*length should be enough memory to start with */
  o_stream->data = vrna_cstr(6 * n, stdout);
  /* retrieve string stream bound to stderr for any info messages */
  o_stream->err = vrna_cstr(n, stderr);

  /*
   ########################################################
   # begin actual calculations
   ########################################################
   */

  /* generate consensus sequence */
  consensus_sequence = (opt->mis) ?
                       vrna_aln_consensus_mis((const char **)alignment, &(opt->md)) :
                       vrna_aln_consensus_sequence((const char **)alignment, &(opt->md));

  /* put header + sequence into output string stream */
  vrna_cstr_print_fasta_header(o_stream->data, record->MSA_ID);

  real_en = vrna_eval_structure_cstr(vc, structure, opt->verbose, o_stream->data);
  cov_en  = vrna_eval_covar_structure(vc, structure);

  vrna_cstr_printf(o_stream->data, "%s\n", consensus_sequence);
  vrna_cstr_printf_structure(o_stream->data,
                             structure,
                             record->tty ?
                             "\n energy = %6.2f kcal/mol "
                             "(%6.2f + %6.2f)" :
                             " (%6.2f = %6.2f + %6.2f)",
                             DBL_ROUND(real_en - cov_en, 2),
                             DBL_ROUND(real_en, 2),
                             DBL_ROUND(-cov_en, 2));

  /* print what we've collected in output charstream */
  if (opt->output_queue)
    vrna_ostream_provide(opt->output_queue, record->number, (void *)o_stream);
  else
    flush_cstr_callback(NULL, record->number, (void *)o_stream);

  vrna_fold_compound_free(vc);
  vrna_aln_free(alignment);
  vrna_aln_free(record->alignment);
  vrna_aln_free(record->names);
  free(consensus_sequence);
  free(record->MSA_ID);
  free(record->structure);
  free(structure);

  free(record);
}
