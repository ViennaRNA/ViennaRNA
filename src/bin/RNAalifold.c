/*
 *                Access to alifold Routines
 *
 *                c Ivo L Hofacker
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
#include "ViennaRNA/plotting/probabilities.h"
#include "ViennaRNA/plotting/structures.h"
#include "ViennaRNA/plotting/alignments.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/alifold.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/io/file_formats_msa.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/params/io.h"
#include "ViennaRNA/MEA.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/constraints/SHAPE.h"
#include "ViennaRNA/plotting/utils.h"
#include "ViennaRNA/datastructures/char_stream.h"
#include "ViennaRNA/datastructures/stream_output.h"

#include "RNAalifold_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"
#include "parallel_helpers.h"

#include "ViennaRNA/color_output.inc"

#define DBL_ROUND(a, digits) (round((a) * pow(10., (double)(digits))) / pow(10., (double)(digits)))

struct options {
  unsigned int    input_format;
  int             filename_full;
  char            *filename_delim;
  int             pf;
  int             noPS;
  int             noDP;
  int             noconv;
  int             MEA;
  double          MEAgamma;
  double          bppmThreshold;
  int             verbose;
  int             quiet;
  vrna_md_t       md;

  dataset_id      id_control;
  int             continuous_names;

  int             n_back;
  int             eval_en;
  int             non_red;

  int             color;
  int             aln_PS;
  int             aln_PS_cols;
  int             mis;
  int             sci;
  int             endgaps;

  int             aln_out;
  char            *aln_out_prefix;

  char            *constraint;
  char            *constraint_file;
  int             constraint_SScons;
  int             constraint_batch;
  int             constraint_enforce;

  int             shape;
  char            **shape_files;
  char            *shape_method;
  int             *shape_file_association;

  int             jobs;
  int             keep_order;
  unsigned int    next_record_number;
  vrna_ostream_t  output_queue;
};


struct record_data {
  unsigned int    number;
  char            *MSA_ID;
  char            **alignment;
  char            **names;
  char            *consensus_structure;
  unsigned int    n_seq;

  int             tty;
  FILE            *output;

  struct options  *options;
};


struct output_stream {
  vrna_cstr_t data;
  vrna_cstr_t err;
};


struct nr_en_data {
  vrna_cstr_t           output;
  vrna_fold_compound_t  *fc;
  double                kT;
  double                ens_en;
};


PRIVATE void
print_nr_samples(const char *structure,
                 void       *data);


PRIVATE void
print_nr_samples_en(const char  *structure,
                    void        *data);


PRIVATE void
print_pi(const vrna_pinfo_t pi,
         FILE               *file);


PRIVATE void
print_aliout(vrna_fold_compound_t *vc,
             plist                *pl,
             double               threshold,
             const char           *mfe,
             FILE                 *aliout);


PRIVATE void
mark_endgaps(char *seq,
             char egap);


static void
apply_constraints(vrna_fold_compound_t  *fc,
                  const char            *SS_cons,
                  struct options        *opt);


static double
compute_sci(const char  **alignment,
            vrna_md_t   *md,
            double      mfe_comparative);


static void
postscript_layout(const char      *filename,
                  const char      **alignment,
                  const char      *consensus_sequence,
                  const char      *consensus_struct,
                  struct options  *opt);


static void
Boltzmann_sampling(vrna_fold_compound_t *fc,
                   double               dG,
                   struct options       *opt,
                   vrna_cstr_t          rec_output);


static void
compute_MEA(vrna_fold_compound_t  *fc,
            struct options        *opt,
            vrna_cstr_t           rec_output);


static void
compute_centroid(vrna_fold_compound_t *fc,
                 struct options       *opt,
                 vrna_cstr_t          rec_output);


static void
apply_SHAPE_data(vrna_fold_compound_t *fc,
                 struct options       *opt);


static char *
get_filename(const char     *id,
             const char     *suffix,
             const char     *filename_default,
             struct options *opt);


int
process_input(FILE            *input_stream,
              const char      *input_filename,
              struct options  *opt);


static void
process_record(struct record_data *record);


void
flush_cstr_callback(void          *auxdata,
                    unsigned int  i,
                    void          *data);


void
init_default_options(struct options *opt)
{
  opt->input_format   = VRNA_FILE_FORMAT_MSA_CLUSTAL; /* default to ClustalW format */
  opt->filename_full  = 0;
  opt->filename_delim = NULL;
  opt->pf             = 0;
  opt->noPS           = 0;
  opt->noDP           = 0;
  opt->noconv         = 0;
  opt->MEA            = 0;
  opt->MEAgamma       = 1.;
  opt->bppmThreshold  = 1e-6;
  opt->verbose        = 0;
  opt->quiet          = 0;
  set_model_details(&(opt->md));

  opt->continuous_names = 0;

  opt->n_back   = 0;
  opt->eval_en  = 0;
  opt->non_red  = 0;

  opt->color        = 0;
  opt->aln_PS       = 0;
  opt->aln_PS_cols  = 60;
  opt->mis          = 0;
  opt->sci          = 0;
  opt->endgaps      = 0;

  opt->aln_out        = 0;
  opt->aln_out_prefix = NULL;

  opt->constraint         = NULL;
  opt->constraint_file    = NULL;
  opt->constraint_SScons  = 0;
  opt->constraint_batch   = 0;
  opt->constraint_enforce = 0;

  opt->shape                  = 0;
  opt->shape_files            = NULL;
  opt->shape_file_association = NULL;
  opt->shape_method           = NULL;

  opt->jobs               = 1;
  opt->keep_order         = 1;
  opt->next_record_number = 0;
  opt->output_queue       = NULL;
}


static char **
collect_unnamed_options(struct RNAalifold_args_info *ggostruct,
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
  struct RNAalifold_args_info args_info;
  unsigned int                longest_string;
  char                        *tmp_string, **input_files;
  int                         s,
                              tmp_number;
  struct  options             opt;
  long int                    first_alignment_number;
  int                         num_input;

  num_input = 0;

  init_default_options(&opt);

  input_files = NULL;

  /*
   #############################################
   # check the command line prameters
   #############################################
   */
  if (RNAalifold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* get basic set of model details */
  ggo_get_md_eval(args_info, opt.md);
  ggo_get_md_fold(args_info, opt.md);
  ggo_get_md_part(args_info, opt.md);
  ggo_get_circ(args_info, opt.md.circ);

  /* temperature */
  ggo_get_temperature(args_info, opt.md.temperature);

  /* check dangle model */
  if (!((opt.md.dangles == 0) || (opt.md.dangles == 2))) {
    vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    opt.md.dangles = dangles = 2;
  }

  ggo_get_id_control(args_info, opt.id_control, "Alignment", "alignment", "_", 4, 1);

  /* do not treat first alignment special */
  if (args_info.continuous_ids_given || get_auto_id(opt.id_control))
    opt.continuous_names = 1;

  ggo_get_constraints_settings(args_info,
                               fold_constrained,
                               opt.constraint_file,
                               opt.constraint_enforce,
                               opt.constraint_batch);

  if (args_info.SS_cons_given) {
    fold_constrained      = 1;
    opt.constraint_SScons = 1;
  }

  /* do not produce postscript output */
  if (args_info.noPS_given)
    opt.noPS = 1;

  /* do not produce dot-plot output */
  if (args_info.noDP_given)
    opt.noDP = 1;

  /* partition function settings */
  if (args_info.partfunc_given) {
    opt.pf = 1;
    if (args_info.partfunc_arg != -1)
      opt.md.compute_bpp = args_info.partfunc_arg;
  }

  /* MEA (maximum expected accuracy) settings */
  if (args_info.MEA_given) {
    opt.pf = opt.MEA = 1;
    if (args_info.MEA_arg != -1)
      opt.MEAgamma = args_info.MEA_arg;
  }

  /* set the bppm threshold for the dotplot */
  if (args_info.bppmThreshold_given)
    opt.bppmThreshold = MIN2(1., MAX2(0., args_info.bppmThreshold_arg));

  /* set cfactor */
  if (args_info.cfactor_given)
    opt.md.cv_fact = args_info.cfactor_arg;

  /* set nfactor */
  if (args_info.nfactor_given)
    opt.md.nc_fact = args_info.nfactor_arg;

  if (args_info.endgaps_given)
    opt.endgaps = 1;

  if (args_info.mis_given)
    opt.mis = 1;

  if (args_info.color_given)
    opt.color = 1;

  if (args_info.aln_given)
    opt.aln_PS = 1;

  if (args_info.aln_EPS_cols_given)
    opt.aln_PS_cols = args_info.aln_EPS_cols_arg;

  if (args_info.aln_stk_given) {
    opt.aln_out = 1;
    if (args_info.aln_stk_arg)
      opt.aln_out_prefix = strdup(args_info.aln_stk_arg);
  }

  if (args_info.old_given)
    opt.md.oldAliEn = 1;

  if (args_info.stochBT_given) {
    opt.n_back          = args_info.stochBT_arg;
    opt.md.uniq_ML      = 1;
    opt.md.compute_bpp  = 0;
    opt.pf              = 1;
    vrna_init_rand();
  }

  if (args_info.stochBT_en_given) {
    opt.n_back          = args_info.stochBT_en_arg;
    opt.md.uniq_ML      = 1;
    opt.md.compute_bpp  = 0;
    opt.pf              = 1;
    opt.eval_en         = 1;
    vrna_init_rand();
  }

  /* non-redundant backtracing */
  if (args_info.nonRedundant_given)
    opt.non_red = 1;

  if (args_info.ribosum_file_given) {
    RibosumFile = strdup(args_info.ribosum_file_arg);
    opt.md.ribo = ribo = 1;
  }

  if (args_info.ribosum_scoring_given) {
    RibosumFile = NULL;
    opt.md.ribo = 1;
  }

  if (args_info.layout_type_given)
    rna_plot_type = args_info.layout_type_arg;

  if (args_info.verbose_given)
    opt.verbose = 1;

  if (args_info.quiet_given) {
    if (opt.verbose)
      vrna_message_warning(
        "Can not be verbose and quiet at the same time! I keep on being chatty...");
    else
      opt.quiet = 1;
  }

  /* SHAPE reactivity data */
  if (args_info.shape_given) {
    if (opt.verbose)
      vrna_message_info(stderr, "SHAPE reactivity data correction activated");

    opt.shape                   = 1;
    opt.shape_files             = (char **)vrna_alloc(sizeof(char *) * (args_info.shape_given + 1));
    opt.shape_file_association  = (int *)vrna_alloc(sizeof(int) * (args_info.shape_given + 1));

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
        opt.shape_files[s]            = strdup(tmp_string);
        opt.shape_file_association[s] = tmp_number - 1;
      } else {
        opt.shape_files[s]            = strdup(args_info.shape_arg[s]);
        opt.shape_file_association[s] = s;
      }

      if (opt.verbose) {
        vrna_message_info(stderr,
                          "Using SHAPE reactivity data provided in file %s for sequence %d",
                          opt.shape_files[s],
                          opt.shape_file_association[s] + 1);
      }
    }

    opt.shape_file_association[s] = -1;

    free(tmp_string);
  }

  if (opt.shape)
    opt.shape_method = strdup(args_info.shapeMethod_arg);

  /* sci computation */
  if (args_info.sci_given)
    opt.sci = 1;

  /* alignment file name(s) given as unnamed option? */
  input_files = collect_unnamed_options(&args_info, &num_input);

  if (num_input > 0)
    /*
     *  Use default alignment file formats.
     *  This may be overridden when we parse the
     *  --input-format parameter below
     */
    opt.input_format = VRNA_FILE_FORMAT_MSA_DEFAULT;

  if (args_info.input_format_given) {
    switch (args_info.input_format_arg[0]) {
      case 'C': /* ClustalW format */
        opt.input_format = VRNA_FILE_FORMAT_MSA_CLUSTAL;
        break;

      case 'S': /* Stockholm 1.0 format */
        opt.input_format = VRNA_FILE_FORMAT_MSA_STOCKHOLM;
        break;

      case 'F': /* FASTA format */
        opt.input_format = VRNA_FILE_FORMAT_MSA_FASTA;
        break;

      case 'M': /* MAF format */
        opt.input_format = VRNA_FILE_FORMAT_MSA_MAF;
        break;

      default:
        vrna_message_warning("Unknown input format specified");
        break;
    }
  }

  /* filename sanitize delimiter */
  if (args_info.filename_delim_given)
    opt.filename_delim = strdup(args_info.filename_delim_arg);
  else if (get_id_delim(opt.id_control))
    opt.filename_delim = strdup(get_id_delim(opt.id_control));

  if ((opt.filename_delim) && isspace(*opt.filename_delim)) {
    free(opt.filename_delim);
    opt.filename_delim = NULL;
  }

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (args_info.noconv_given)
    opt.noconv = 1;

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
      "This version of RNAfold has been built without parallel input processing capabilities");
#endif

    if (args_info.unordered_given)
      opt.keep_order = 0;
  }

  /* free allocated memory of command line data structure */
  RNAalifold_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */
  if (opt.md.circ && opt.md.gquad)
    vrna_message_error("G-Quadruplex support is currently not available for circular RNA structures");

  if (opt.md.circ && opt.md.noLP)
    vrna_message_warning("Depending on the origin of the circular sequence, "
                         "some structures may be missed when using --noLP\n"
                         "Try rotating your sequence a few times\n");

  first_alignment_number = get_current_id(opt.id_control);

  if (opt.keep_order)
    opt.output_queue = vrna_ostream_init(&flush_cstr_callback, NULL);

  /*
   ################################################
   # read constraint from stdin
   ################################################
   */
  if (fold_constrained && (!opt.constraint_file) && (!opt.constraint_SScons)) {
    if (isatty(fileno(stdin))) {
      vrna_message_constraint_options_all();
      vrna_message_input_msa("");
    }

    char          *input_string = NULL;
    unsigned int  input_type    = get_input_line(&input_string, VRNA_INPUT_NOSKIP_COMMENTS);

    if (input_type & VRNA_INPUT_QUIT)
      return EXIT_FAILURE;
    else if ((input_type & VRNA_INPUT_MISC) && (strlen(input_string) > 0))
      opt.constraint = strdup(input_string);

    free(input_string);
  }

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
          vrna_message_error("Unable to open %d. input file \"%s\" for reading",
                             i + 1,
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

  UNINIT_PARALLELIZATION
  /*
   ################################################
   # post processing
   ################################################
   */
  vrna_ostream_free(opt.output_queue);


  /* check whether we've actually processed any alignment so far */
  if (first_alignment_number == get_current_id(opt.id_control)) {
    char *msg = "Missing sequences in input file(s)! "
                "Either your file is empty, or not in %s format!";

    switch (opt.input_format) {
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

  free(opt.shape_files);
  free(input_files);
  free(opt.filename_delim);
  free_id_data(opt.id_control);

  return EXIT_SUCCESS;
}


int
process_input(FILE            *input_stream,
              const char      *input_filename,
              struct options  *opt)
{
  int           ret           = 1;
  unsigned int  input_format  = opt->input_format;
  int           istty_in      = isatty(fileno(input_stream));

  /* detect input file format if reading from file */
  if (input_filename) {
    unsigned int format_guess = vrna_file_msa_detect_format(input_filename, opt->input_format);

    if (format_guess == VRNA_FILE_FORMAT_MSA_UNKNOWN) {
      char *msg = "Your input file is missing sequences! "
                  "Either your file is empty, or not in %s format!";

      switch (opt->input_format) {
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
    struct record_data *record = (struct record_data *)vrna_alloc(sizeof(struct record_data));

    record->number  = opt->next_record_number;
    record->MSA_ID  = fileprefix_from_id_alifold(tmp_id,
                                                 opt->id_control,
                                                 opt->continuous_names);
    record->alignment           = alignment;
    record->names               = names;
    record->consensus_structure = tmp_structure;
    record->n_seq               = (unsigned int)n_seq;

    record->tty     = istty_in;
    record->output  = stdout;

    record->options = opt;

    if (opt->output_queue)
      vrna_ostream_request(opt->output_queue, opt->next_record_number++);

    /* process the record we've just read */
    RUN_IN_PARALLEL(process_record, record);

    free(tmp_id);

    /* break after first record if fold_constrained and not explicitly instructed otherwise */
    if (opt->shape || (fold_constrained && (!(opt->constraint_batch || opt->constraint_SScons)))) {
      ret = 0;
      break;
    }
  }

  return ret;
}


static void
process_record(struct record_data *record)
{
  char                  **alignment, *consensus_sequence, *mfe_structure;
  unsigned int          n, i, n_seq;
  double                min_en, real_en, cov_en;
  struct options        *opt;
  vrna_fold_compound_t  *vc;
  struct output_stream  *o_stream;

  o_stream = (struct output_stream *)vrna_alloc(sizeof(struct output_stream));

  opt   = record->options;
  n_seq = record->n_seq;

  /* construct output file names */
  char  *filename_plot  = get_filename(record->MSA_ID, "ss.ps", "alirna.ps", opt);
  char  *filename_dot   = get_filename(record->MSA_ID, "dp.ps", "alidot.ps", opt);
  char  *filename_aln   = get_filename(record->MSA_ID, "aln.ps", "aln.ps", opt);
  char  *filename_out   = get_filename(record->MSA_ID, "ali.out", "alifold.out", opt);

  /*
   *  create a new alignment for internal computations
   *  which require all-uppercase letters, and RNA/DNA alphabet
   */
  alignment = vrna_aln_copy((const char **)record->alignment,
                            VRNA_ALN_UPPERCASE |
                            (opt->noconv ? 0 : VRNA_ALN_RNA));

  if (opt->endgaps)
    for (i = 0; i < n_seq; i++)
      mark_endgaps(alignment[i], '~');

  vc = vrna_fold_compound_comparative((const char **)alignment,
                                      &(opt->md),
                                      VRNA_OPTION_DEFAULT);

  if (!vc) {
    vrna_message_warning("Skipping computations for \"%s\"",
                         (record->MSA_ID) ? record->MSA_ID : "identifier unavailable");
    return;
  }

  n = vc->length;

  if (fold_constrained)
    apply_constraints(vc, record->consensus_structure, opt);

  if (opt->shape)
    apply_SHAPE_data(vc, opt);

  /* retrieve string stream bound to stdout, 6*length should be enough memory to start with */
  o_stream->data = vrna_cstr(6 * n, stdout);
  /* retrieve string stream bound to stderr for any info messages */
  o_stream->err = vrna_cstr(n, stderr);

  if (!opt->quiet) {
    vrna_cstr_message_info(o_stream->err,
                           "%u sequences; length of alignment %u.",
                           n_seq,
                           n);
  }

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
  vrna_cstr_printf(o_stream->data, "%s\n", consensus_sequence);

  mfe_structure = (char *)vrna_alloc(sizeof(char) * (n + 1));
  min_en        = vrna_mfe(vc, mfe_structure);

  /* check whether the constraint allows for any solution */
  if (fold_constrained) {
    if (min_en == (double)(INF / 100.)) {
      vrna_message_error(
        "Supplied structure constraints create empty solution set for alignment Nr. %d",
        record->number);
      exit(EXIT_FAILURE);
    }
  }

  real_en = vrna_eval_structure(vc, mfe_structure);
  cov_en  = vrna_eval_covar_structure(vc, mfe_structure);

  if (opt->sci) {
    double sci = compute_sci((const char **)alignment, &(opt->md), min_en);

    if (opt->shape) {
      vrna_cstr_printf_structure(o_stream->data,
                                 mfe_structure,
                                 record->tty ?
                                 "\n minimum free energy = %6.2f kcal/mol "
                                 "(%6.2f + %6.2f + %6.2f)\n SCI = %2.4f" :
                                 " (%6.2f = %6.2f + %6.2f + %6.2f) [sci = %2.4f]",
                                 DBL_ROUND(min_en, 2),
                                 DBL_ROUND(real_en, 2),
                                 DBL_ROUND(-cov_en, 2),
                                 DBL_ROUND(min_en - real_en + cov_en, 2),
                                 DBL_ROUND(sci, 4));
    } else {
      vrna_cstr_printf_structure(o_stream->data,
                                 mfe_structure,
                                 record->tty ?
                                 "\n minimum free energy = %6.2f kcal/mol "
                                 "(%6.2f + %6.2f)\n SCI = %2.4f" :
                                 " (%6.2f = %6.2f + %6.2f) [sci = %2.4f]",
                                 DBL_ROUND(min_en, 2),
                                 DBL_ROUND(real_en, 2),
                                 DBL_ROUND(min_en - real_en, 2),
                                 DBL_ROUND(sci, 4));
    }
  } else {
    if (opt->shape) {
      vrna_cstr_printf_structure(o_stream->data,
                                 mfe_structure,
                                 record->tty ?
                                 "\n minimum free energy = %6.2f kcal/mol "
                                 "(%6.2f + %6.2f + %6.2f)" :
                                 " (%6.2f = %6.2f + %6.2f + %6.2f)",
                                 DBL_ROUND(min_en, 2),
                                 DBL_ROUND(real_en, 2),
                                 DBL_ROUND(-cov_en, 2),
                                 DBL_ROUND(min_en - real_en + cov_en, 2));
    } else {
      vrna_cstr_printf_structure(o_stream->data,
                                 mfe_structure,
                                 record->tty ?
                                 "\n minimum free energy = %6.2f kcal/mol "
                                 "(%6.2f + %6.2f)" :
                                 " (%6.2f = %6.2f + %6.2f)",
                                 DBL_ROUND(min_en, 2),
                                 DBL_ROUND(real_en, 2),
                                 DBL_ROUND(min_en - real_en, 2));
    }
  }

  if (!opt->noPS) {
    postscript_layout(filename_plot,
                      (const char **)alignment,
                      consensus_sequence,
                      mfe_structure,
                      opt);
  }

  if (opt->aln_PS) {
    THREADSAFE_FILE_OUTPUT(
      vrna_file_PS_aln(filename_aln,
                       (const char **)record->alignment,
                       (const char **)record->names,
                       mfe_structure,
                       opt->aln_PS_cols));
  }

  /* free mfe arrays */
  vrna_mx_mfe_free(vc);

  if (opt->pf) {
    double  energy;
    char    *pairing_propensity;

    pairing_propensity = (char *)vrna_alloc(sizeof(char) * (n + 1));

    vrna_exp_params_rescale(vc, &min_en);

    if ((n > 2000) && (!opt->quiet))
      vrna_cstr_message_info(o_stream->err,
                             "scaling factor %f\n",
                             vc->exp_params->pf_scale);

    energy = vrna_pf(vc, pairing_propensity);

    if (opt->n_back > 0) {
      Boltzmann_sampling(vc, energy, opt, o_stream->data);
    } else if (opt->md.compute_bpp) {
      FILE  *aliout;
      cpair *cp;
      plist *pl, *mfel;

      vrna_cstr_printf_structure(o_stream->data,
                                 pairing_propensity,
                                 record->tty ?
                                 "\n free energy of ensemble = %6.2f kcal/mol" :
                                 " [%6.2f]",
                                 DBL_ROUND(energy, 2));

      pl    = vrna_plist_from_probs(vc, opt->bppmThreshold);
      mfel  = vrna_plist(mfe_structure, 0.95 * 0.95);

      if (!opt->md.circ)
        compute_centroid(vc, opt, o_stream->data);

      if (opt->MEA)
        compute_MEA(vc, opt, o_stream->data);

      THREADSAFE_FILE_OUTPUT({
        aliout = fopen(filename_out, "w");
        if (!aliout)
          vrna_message_warning("can't open %s ... skipping output", filename_out);
        else
          print_aliout(vc, pl, opt->bppmThreshold, mfe_structure, aliout);

        fclose(aliout);
      });

      if (!opt->noDP) {
        cp = vrna_annotate_covar_pairs((const char **)alignment,
                                       pl,
                                       mfel,
                                       opt->bppmThreshold,
                                       &(opt->md));

        THREADSAFE_FILE_OUTPUT((void)PS_color_dot_plot(consensus_sequence,
                                                       cp,
                                                       filename_dot));

        free(cp);
      }

      free(pl);
      free(mfel);

      vrna_cstr_printf_structure(o_stream->data,
                                 NULL,
                                 " frequency of mfe structure in ensemble %g"
                                 "; ensemble diversity %-6.2f",
                                 vrna_pr_energy(vc, min_en),
                                 vrna_mean_bp_distance(vc));
    } else {
      vrna_cstr_printf_structure(o_stream->data,
                                 NULL,
                                 " free energy of ensemble = %6.2f kcal/mol\n"
                                 " frequency of mfe structure in ensemble %g;",
                                 DBL_ROUND(energy, 2),
                                 vrna_pr_energy(vc, min_en));
    }

    free(pairing_propensity);
  } /* end partition function block */

  /* write output multiple sequence alignment */
  if (opt->aln_out) {
    char *filename = NULL;
    if (!opt->aln_out_prefix) {
      filename = strdup("RNAalifold_results.stk");
    } else {
      char *tmp = vrna_strdup_printf("%s.stk", opt->aln_out_prefix);
      filename = vrna_filename_sanitize(tmp, opt->filename_delim);
      free(tmp);
    }

    THREADSAFE_FILE_OUTPUT(
      vrna_file_msa_write((const char *)filename,
                          (const char **)record->names,
                          (const char **)record->alignment,
                          record->MSA_ID,
                          (const char *)mfe_structure,
                          "RNAalifold prediction",
                          (opt->mis ? VRNA_FILE_FORMAT_MSA_MIS : 0) |
                          VRNA_FILE_FORMAT_MSA_STOCKHOLM |
                          VRNA_FILE_FORMAT_MSA_APPEND));

    free(filename);
  }

  /* print what we've collected in output charstream */
  if (opt->output_queue)
    vrna_ostream_provide(opt->output_queue, record->number, (void *)o_stream);
  else
    THREADSAFE_STREAM_OUTPUT(flush_cstr_callback(NULL, record->number, (void *)o_stream));

  free(consensus_sequence);
  free(mfe_structure);
  free(filename_plot);
  free(filename_dot);
  free(filename_aln);
  free(filename_out);
  vrna_fold_compound_free(vc);

  vrna_aln_free(alignment);

  free(record->MSA_ID);
  vrna_aln_free(record->alignment);
  vrna_aln_free(record->names);
  free(record->consensus_structure);

  free(record);
}


static void
apply_constraints(vrna_fold_compound_t  *fc,
                  const char            *SS_cons,
                  struct options        *opt)
{
  const char    *constraints_file;
  const char    *cstruc;
  unsigned int  constraint_options;

  constraints_file  = opt->constraint_file;
  cstruc            = opt->constraint;

  if (constraints_file) {
    /** [Adding hard constraints from file] */
    vrna_constraints_add(fc, constraints_file, VRNA_OPTION_DEFAULT);
    /** [Adding hard constraints from file] */
  } else if (cstruc) {
    unsigned int  length  = fc->length;
    unsigned int  cl      = strlen(cstruc);

    if (cl == 0)
      vrna_message_warning("structure constraint is missing");
    else if (cl < length)
      vrna_message_warning("structure constraint is shorter than sequence");
    else if (cl > length)
      vrna_message_error("structure constraint is too long");

    /** [Adding hard constraints from pseudo dot-bracket] */
    constraint_options = VRNA_CONSTRAINT_DB_DEFAULT;

    if (opt->constraint_enforce)
      constraint_options |= VRNA_CONSTRAINT_DB_ENFORCE_BP;

    vrna_constraints_add(fc, cstruc, constraint_options);
  } else if ((opt->constraint_SScons) && (SS_cons)) {
    constraint_options = VRNA_CONSTRAINT_DB_DEFAULT | VRNA_CONSTRAINT_DB_WUSS;
    if (opt->constraint_enforce)
      constraint_options |= VRNA_CONSTRAINT_DB_ENFORCE_BP;

    vrna_constraints_add(fc, SS_cons, constraint_options);
  } else {
    vrna_message_warning("Constraint missing");
  }
}


static double
compute_sci(const char  **alignment,
            vrna_md_t   *md,
            double      mfe_comparative)
{
  unsigned int  i;
  double        sci, e_mean;

  sci = e_mean = 0.;

  for (i = 0; alignment[i]; i++) {
    char *seq = vrna_seq_ungapped(alignment[i]);

    if (strlen(seq) > 0) {
      vrna_fold_compound_t *fc = vrna_fold_compound(seq, md, VRNA_OPTION_DEFAULT);
      e_mean += vrna_mfe(fc, NULL);
      vrna_fold_compound_free(fc);
    }

    free(seq);
  }

  if ((i > 0) && (e_mean != 0.)) {
    e_mean  /= i;
    sci     = mfe_comparative / e_mean;
  }

  return sci;
}


static void
postscript_layout(const char      *filename,
                  const char      **alignment,
                  const char      *consensus_sequence,
                  const char      *consensus_structure,
                  struct options  *opt)
{
  char **A;

  A = vrna_annotate_covar_db(alignment, consensus_structure, &(opt->md));

  THREADSAFE_FILE_OUTPUT(
    (void)vrna_file_PS_rnaplot_a(consensus_sequence,
                                 consensus_structure,
                                 filename,
                                 (opt->color) ? A[0] : NULL,
                                 A[1],
                                 &(opt->md)));

  free(A[0]);
  free(A[1]);
  free(A);
}


static void
Boltzmann_sampling(vrna_fold_compound_t *fc,
                   double               dG,
                   struct options       *opt,
                   vrna_cstr_t          rec_output)
{
  unsigned int i, options;

  options = (opt->non_red) ?
            VRNA_PBACKTRACK_NON_REDUNDANT :
            VRNA_PBACKTRACK_DEFAULT;

  /*stochastic sampling*/
  if (opt->eval_en) {
    struct nr_en_data dat;
    dat.output  = rec_output;
    dat.fc      = fc;
    dat.kT      = fc->exp_params->kT / 1000.;
    dat.ens_en  = dG;

    vrna_pbacktrack_cb(fc,
                       opt->n_back,
                       &print_nr_samples_en,
                       (void *)&dat,
                       options);
  } else {
    vrna_pbacktrack_cb(fc,
                       opt->n_back,
                       &print_nr_samples,
                       (void *)&rec_output,
                       options);
  }
}


static void
print_nr_samples(const char *structure,
                 void       *data)
{
  if (structure)
    vrna_cstr_printf_structure(*((vrna_cstr_t *)data), structure, NULL);
}


static void
print_nr_samples_en(const char  *structure,
                    void        *data)
{
  if (structure) {
    struct nr_en_data     *d      = (struct nr_en_data *)data;
    vrna_cstr_t           output  = d->output;
    vrna_fold_compound_t  *fc     = d->fc;
    double                kT      = d->kT;
    double                ens_en  = d->ens_en;

    double                e = vrna_eval_structure(fc, structure);
    e -= (double)vrna_eval_covar_structure(fc, structure);
    double                prob = exp((ens_en - e) / kT);

    vrna_cstr_printf_structure(output,
                               structure,
                               " %6.2f %6g",
                               e,
                               prob);
  }
}


static void
compute_MEA(vrna_fold_compound_t  *fc,
            struct options        *opt,
            vrna_cstr_t           rec_output)
{
  char  *MEA_structure;
  int   gq;
  float mea, *ens;

  /* compute MEA */
  MEA_structure = vrna_MEA(fc, opt->MEAgamma, &mea);

  /* evaluate MEA structure */
  ens     = (float *)vrna_alloc(2 * sizeof(float));
  ens[0]  = vrna_eval_structure(fc, MEA_structure);
  ens[1]  = vrna_eval_covar_structure(fc, MEA_structure);

  vrna_cstr_printf_structure(rec_output,
                             MEA_structure,
                             " {%6.2f = %6.2f + %6.2f MEA=%.2f}",
                             DBL_ROUND(ens[0] - ens[1], 2),
                             DBL_ROUND(ens[0], 2),
                             DBL_ROUND((-1) * ens[1], 2),
                             mea);

  /* cleanup */
  free(MEA_structure);
  free(ens);
}


static void
compute_centroid(vrna_fold_compound_t *fc,
                 struct options       *opt,
                 vrna_cstr_t          rec_output)
{
  char    *centroid_structure;
  float   *ens;
  double  dist;

  dist = 0;
  /* compute centroid structure */
  centroid_structure = vrna_centroid(fc, &dist);
  /* evaluate centroid structure */
  ens     = (float *)vrna_alloc(2 * sizeof(float));
  ens[0]  = vrna_eval_structure(fc, centroid_structure);
  ens[1]  = vrna_eval_covar_structure(fc, centroid_structure);

  vrna_cstr_printf_structure(rec_output,
                             centroid_structure,
                             " {%6.2f = %6.2f + %6.2f d=%.2f}",
                             DBL_ROUND(ens[0] - ens[1], 2),
                             DBL_ROUND(ens[0], 2),
                             DBL_ROUND((-1) * ens[1], 2),
                             dist);

  /* cleanup */
  free(centroid_structure);
  free(ens);
}


static void
apply_SHAPE_data(vrna_fold_compound_t *fc,
                 struct options       *opt)
{
  unsigned int s;

  /* count number of SHAPE reactivity data sets */
  for (s = 0; opt->shape_file_association[s] != -1; s++);

  if ((s != fc->n_seq) && (!opt->quiet)) {
    vrna_message_warning("Number of sequences in alignment (%d) does not match "
                         "number of provided SHAPE reactivity data files (%d)!",
                         fc->n_seq,
                         s);
  }

  /* shrink SHAPE data lists to required size */
  opt->shape_files = (char **)vrna_realloc(opt->shape_files,
                                           (s + 1) * sizeof(char *));
  opt->shape_file_association = (int *)vrna_realloc(opt->shape_file_association,
                                                    (s + 1) * sizeof(int));

  vrna_constraints_add_SHAPE_ali(fc,
                                 opt->shape_method,
                                 (const char **)opt->shape_files,
                                 opt->shape_file_association,
                                 opt->verbose,
                                 VRNA_OPTION_DEFAULT);
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


PRIVATE void
mark_endgaps(char *seq,
             char egap)
{
  int i, n;

  n = strlen(seq);
  for (i = 0; i < n && (seq[i] == '-'); i++)
    seq[i] = egap;
  for (i = n - 1; i > 0 && (seq[i] == '-'); i--)
    seq[i] = egap;
}


PRIVATE void
print_pi(const vrna_pinfo_t pi,
         FILE               *file)
{
  const char  *pname[8] = {
    "", "CG", "GC", "GU", "UG", "AU", "UA", "--"
  };
  int         i;

  /* numbering starts with 1 in output */
  fprintf(file, "%5d %5d %2d %5.1f%% %7.3f",
          pi.i, pi.j, pi.bp[0], 100. * pi.p, pi.ent);
  for (i = 1; i <= 7; i++)
    if (pi.bp[i])
      fprintf(file, " %s:%-4d", pname[i], pi.bp[i]);

  if (!pi.comp)
    fprintf(file, " +");

  fprintf(file, "\n");
}


/*-------------------------------------------------------------------------*/

PRIVATE void
print_aliout(vrna_fold_compound_t *vc,
             plist                *pl,
             double               threshold,
             const char           *mfe,
             FILE                 *aliout)
{
  int           k;
  vrna_pinfo_t  *pi;
  char          **AS  = vc->sequences;
  int           n_seq = vc->n_seq;

  pi = vrna_aln_pinfo(vc, (const char *)mfe, threshold);

  /* print it */
  fprintf(aliout, "%d sequence; length of alignment %d\n",
          n_seq, (int)strlen(AS[0]));
  fprintf(aliout, "alifold output\n");

  for (k = 0; pi[k].i > 0; k++)
    print_pi(pi[k], aliout);

  fprintf(aliout, "%s\n", mfe);
  free(pi);
}
