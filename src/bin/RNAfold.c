/*
 *                Ineractive Access to folding Routines
 *
 *                c Ivo L Hofacker
 *                Vienna RNA package
 */

/** \file
 *  \brief RNAfold program source code
 *
 *  This code provides an interface for MFE and Partition function folding
 *  of single linear or circular RNA molecules.
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
#include <string.h>

#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/plotting/probabilities.h"
#include "ViennaRNA/plotting/structures.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/params/io.h"
#include "ViennaRNA/centroid.h"
#include "ViennaRNA/MEA.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/basic.h"
#include "ViennaRNA/constraints/SHAPE.h"
#include "ViennaRNA/constraints/ligand.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/commands.h"
#include "ViennaRNA/equilibrium_probs.h"
#include "ViennaRNA/datastructures/char_stream.h"
#include "ViennaRNA/datastructures/stream_output.h"
#include "ViennaRNA/combinatorics.h"
#include "ViennaRNA/color_output.inc"

#include "RNAfold_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"
#include "parallel_helpers.h"


struct options {
  int             filename_full;
  char            *filename_delim;
  int             pf;
  int             noPS;
  int             noconv;
  int             lucky;
  int             MEA;
  double          MEAgamma;
  double          bppmThreshold;
  int             verbose;
  char            *ligandMotif;
  vrna_cmd_t      cmds;
  vrna_md_t       md;
  dataset_id      id_control;

  char            *constraint_file;
  int             constraint_batch;
  int             constraint_enforce;
  int             constraint_canonical;

  int             shape;
  char            *shape_file;
  char            *shape_method;
  char            *shape_conversion;

  int             jobs;
  int             tofile;
  char            *output_file;
  int             keep_order;
  FILE            *output_stream;
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


struct output_stream {
  vrna_cstr_t data;
  int         individual;
};


static char *
annotate_ligand_motif(vrna_fold_compound_t  *vc,
                      const char            *structure);


static void
print_ligand_motifs(vrna_fold_compound_t  *vc,
                    const char            *structure,
                    const char            *structure_name,
                    vrna_cstr_t           buf);


static void add_ligand_motif(vrna_fold_compound_t *vc,
                             char                 *motifstring,
                             int                  verbose,
                             unsigned int         options);


static char *
annotate_ud_motif(vrna_fold_compound_t  *vc,
                  vrna_ud_motif_t       *motifs);


static void
print_ud_motifs(vrna_fold_compound_t  *vc,
                vrna_ud_motif_t       *motifs,
                const char            *structure_name,
                vrna_cstr_t           buf);


static void
add_ligand_motifs_dot(vrna_fold_compound_t  *fc,
                      vrna_ep_t             **prob_list,
                      vrna_ep_t             **mfe_list,
                      const char            *structure);


static void
add_ligand_motifs_to_list(vrna_ep_t       **list,
                          vrna_sc_motif_t *motifs);


static void
compute_MEA(vrna_fold_compound_t  *fc,
            double                MEAgamma,
            const char            *ligandMotif,
            int                   verbose,
            vrna_cstr_t           buf);


static void
compute_centroid(vrna_fold_compound_t *fc,
                 const char           *ligandMotif,
                 int                  verbose,
                 vrna_cstr_t          buf);


static void
apply_constraints(vrna_fold_compound_t  *fc,
                  const char            *constraints_file,
                  const char            **rec_rest,
                  int                   maybe_multiline,
                  int                   enforceConstraints,
                  int                   canonicalBPonly);


static char *
generate_filename(const char  *pattern,
                  const char  *def_name,
                  const char  *id,
                  const char  *filename_delim);


int
process_input(FILE            *input_stream,
              const char      *input_filename,
              struct options  *opt);


static void
process_record(struct record_data *record);


/*--------------------------------------------------------------------------*/
void
flush_cstr_callback(void          *auxdata,
                    unsigned int  i,
                    void          *data)
{
  struct output_stream *s = (struct output_stream *)data;

  if (s) {
    /* flush data[k] */
    vrna_cstr_fflush(s->data);
    /* free/close data[k] */
    if (s->individual)
      vrna_cstr_close(s->data);
    else
      vrna_cstr_free(s->data);

    free(s);
  }
}


static void
postscript_layout(vrna_fold_compound_t  *fc,
                  const char            *orig_sequence,
                  const char            *structure,
                  const char            *SEQ_ID,
                  const char            *ligandMotif,
                  const char            *filename_delim,
                  int                   verbose)
{
  char      *filename_plot  = NULL;
  char      *annotation     = NULL;
  vrna_md_t *md             = &(fc->params->model_details);

  filename_plot = generate_filename("%s%sss.ps",
                                    "rna.ps",
                                    SEQ_ID,
                                    filename_delim);

  if (ligandMotif) {
    char *annote = annotate_ligand_motif(fc, structure);
    vrna_strcat_printf(&annotation, annote);
    free(annote);
  }

  if (fc->domains_up) {
    vrna_ud_motif_t *m  = vrna_ud_motifs_MFE(fc, structure);
    char            *a  = annotate_ud_motif(fc, m);
    vrna_strcat_printf(&annotation, a);
    free(a);
    free(m);
  }

  THREADSAFE_FILE_OUTPUT(
    vrna_file_PS_rnaplot_a(orig_sequence,
                           structure,
                           filename_plot,
                           annotation,
                           NULL,
                           md));
  free(annotation);
  free(filename_plot);
}


static void
ImFeelingLucky(vrna_fold_compound_t *fc,
               const char           *orig_sequence,
               const char           *SEQ_ID,
               int                  noPS,
               const char           *filename_delim,
               vrna_cstr_t          buf,
               int                  istty_in)
{
  vrna_md_t *md = &(fc->params->model_details);

  vrna_init_rand();

  char      *filename_plot  = NULL;
  char      *s              = vrna_pbacktrack(fc);
  float     e               = vrna_eval_structure(fc, (const char *)s);

  vrna_cstr_printf_structure(buf,
                             s,
                             (istty_in) ? "\n free energy = %6.2f kcal/mol" : " (%6.2f)",
                             e);

  if (!noPS) {
    filename_plot = generate_filename("%s%sss.ps",
                                      "rna.ps",
                                      SEQ_ID,
                                      filename_delim);

    THREADSAFE_FILE_OUTPUT(
      vrna_file_PS_rnaplot(orig_sequence,
                           s,
                           filename_plot,
                           md));
    free(filename_plot);
  }

  free(s);
}


static char *
generate_filename(const char  *pattern,
                  const char  *def_name,
                  const char  *id,
                  const char  *filename_delim)
{
  char *filename, *ptr;

  if (id) {
    filename  = vrna_strdup_printf(pattern, id, filename_delim);
    ptr       = vrna_filename_sanitize(filename, filename_delim);
    free(filename);
    filename = ptr;
  } else {
    filename = strdup(def_name);
  }

  return filename;
}


static char **
collect_unnamed_options(struct RNAfold_args_info  *ggostruct,
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
append_input_files(struct RNAfold_args_info *ggostruct,
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
init_default_options(struct options *opt)
{
  opt->filename_full  = 0;
  opt->filename_delim = NULL;
  opt->pf             = 0;
  opt->noPS           = 0;
  opt->noconv         = 0;
  opt->lucky          = 0;
  opt->MEA            = 0;
  opt->MEAgamma       = 1.;
  opt->bppmThreshold  = 1e-5;
  opt->verbose        = 0;
  opt->ligandMotif    = NULL;
  opt->cmds           = NULL;
  set_model_details(&(opt->md));

  opt->constraint_file      = NULL;
  opt->constraint_batch     = 0;
  opt->constraint_enforce   = 0;
  opt->constraint_canonical = 0;

  opt->shape            = 0;
  opt->shape_file       = NULL;
  opt->shape_method     = NULL;
  opt->shape_conversion = NULL;

  opt->jobs               = 1;
  opt->tofile             = 0;
  opt->output_file        = NULL;
  opt->keep_order         = 1;
  opt->output_stream      = NULL;
  opt->next_record_number = 0;
  opt->output_queue       = NULL;
}


int
main(int  argc,
     char *argv[])
{
  struct  RNAfold_args_info args_info;
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
  if (RNAfold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* get basic set of model details */
  ggo_get_md_eval(args_info, opt.md);
  ggo_get_md_fold(args_info, opt.md);
  ggo_get_md_part(args_info, opt.md);
  ggo_get_circ(args_info, opt.md.circ);

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

  /* always look on the bright side of life */
  if (args_info.ImFeelingLucky_given)
    opt.md.uniq_ML = opt.lucky = opt.pf = st_back = 1;

  /* set the bppm threshold for the dotplot */
  if (args_info.bppmThreshold_given)
    opt.bppmThreshold = MIN2(1., MAX2(0., args_info.bppmThreshold_arg));

  /* do not produce postscript output */
  if (args_info.noPS_given)
    opt.noPS = 1;

  /* partition function settings */
  if (args_info.partfunc_given) {
    opt.pf = 1;
    if (args_info.partfunc_arg != 1)
      opt.md.compute_bpp = do_backtrack = args_info.partfunc_arg;
    else
      opt.md.compute_bpp = do_backtrack = 1;
  }

  /* MEA (maximum expected accuracy) settings */
  if (args_info.MEA_given) {
    opt.pf = opt.MEA = 1;
    if (args_info.MEA_arg != -1)
      opt.MEAgamma = args_info.MEA_arg;
  }

  if (args_info.layout_type_given)
    rna_plot_type = args_info.layout_type_arg;

  if (args_info.verbose_given)
    opt.verbose = 1;

  if (args_info.outfile_given) {
    opt.tofile = 1;
    if (args_info.outfile_arg)
      opt.output_file = strdup(args_info.outfile_arg);
  }

  if (args_info.motif_given)
    opt.ligandMotif = strdup(args_info.motif_arg);

  if (args_info.commands_given)
    opt.cmds = vrna_file_commands_read(args_info.commands_arg, VRNA_CMD_PARSE_DEFAULTS);

  /* filename sanitize delimiter */
  if (args_info.filename_delim_given)
    opt.filename_delim = strdup(args_info.filename_delim_arg);
  else if (get_id_delim(opt.id_control))
    opt.filename_delim = strdup(get_id_delim(opt.id_control));

  if ((opt.filename_delim) && isspace(*opt.filename_delim)) {
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
      "This version of RNAfold has been built without parallel input processing capabilities");
#endif

    if (args_info.unordered_given)
      opt.keep_order = 0;
  }

  input_files = collect_unnamed_options(&args_info, &num_input);
  input_files = append_input_files(&args_info, input_files, &num_input);

  /* free allocated memory of command line data structure */
  RNAfold_cmdline_parser_free(&args_info);


  /*
   #############################################
   # begin initializing
   #############################################
   */
  if (opt.md.circ && opt.md.gquad) {
    vrna_message_error("G-Quadruplex support is currently not available for circular RNA structures");
    exit(EXIT_FAILURE);
  }

  if (opt.md.circ && opt.md.noLP)
    vrna_message_warning("depending on the origin of the circular sequence, some structures may be missed when using --noLP\n"
                         "Try rotating your sequence a few times");

  if ((opt.verbose) && (opt.jobs > 1))
    vrna_message_info(stderr, "Preparing %d parallel computation slots", opt.jobs);

  if (opt.keep_order)
    opt.output_queue = vrna_ostream_init(&flush_cstr_callback, NULL);

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

  /* close output stream if necessary */
  if ((opt.output_stream) && (opt.output_stream != stdout))
    fclose(opt.output_stream);

  vrna_ostream_free(opt.output_queue);

  free(input_files);
  free(opt.constraint_file);
  free(opt.ligandMotif);
  free(opt.shape_file);
  free(opt.shape_method);
  free(opt.shape_conversion);
  free(opt.filename_delim);
  vrna_commands_free(opt.cmds);

  free_id_data(opt.id_control);

  return EXIT_SUCCESS;
}


struct output_stream *
get_output_stream(unsigned int    init_size,
                  struct options  *opt,
                  const char      *SEQ_ID,
                  const char      *input_filename)
{
  struct output_stream  *o_stream;
  FILE                  *output;
  int                   individual_stream;

  individual_stream = 0; /* we default to using a single output sink */

  o_stream = (struct output_stream *)vrna_alloc(sizeof(struct output_stream));

  /* in case we do parallel processing of input, let's block access to the opt->output_stream pointer */
  ATOMIC_BLOCK(({
    /* default to stream that we've already opened */
    output = opt->output_stream;

    if ((!opt->tofile) && (!output)) {
      output = stdout;
      opt->output_stream = stdout;
    } else if (opt->tofile) {
      char *filename, *tmp;

      tmp = filename = NULL;

      if ((!opt->output_file) && (SEQ_ID)) {
        /* need to open new individual output file */
        tmp = vrna_strdup_printf("%s.fold", SEQ_ID);
        individual_stream = 1;

        filename = vrna_filename_sanitize(tmp, opt->filename_delim);

        if ((input_filename) && !strcmp(input_filename, filename))
          vrna_message_error("Input and output file names are identical");

        if (!(output = fopen(filename, "a")))
          vrna_message_error("Failed to open file for writing");
      } else if (!output) {
        /* we need to open global output file */
        tmp = (opt->output_file) ?
              vrna_strdup_printf("%s", opt->output_file) :
              vrna_strdup_printf("RNAfold_output.fold");

        filename = vrna_filename_sanitize(tmp, opt->filename_delim);

        if ((input_filename) && !strcmp(input_filename, filename))
          vrna_message_error("Input and output file names are identical");

        if (!(output = fopen(filename, "a")))
          vrna_message_error("Failed to open file for writing");

        opt->output_stream = output;
      }

      free(tmp);
      free(filename);
    }

    /* actually initialize vrna_cstr_t of the stream */
    o_stream->data = vrna_cstr(init_size, output);
    o_stream->individual = (individual_stream) ? 1 : 0;
  }));

  return o_stream;
}


/* main loop that processes an input stream */
int
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
    if (fold_constrained) {
      vrna_message_constraint_options_all();
      vrna_message_input_seq("Input sequence (upper or lower case) followed by structure constraint");
    } else {
      vrna_message_input_seq_simple();
    }
  }

  /* set options we wanna pass to vrna_file_fasta_read_record() */
  if (istty_in)
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;

  if (!fold_constrained)
    read_opt |= VRNA_INPUT_NO_REST;

  /* main loop that processes each record obtained from input stream */
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

    if (opt->shape || (opt->constraint_file && (!opt->constraint_batch))) {
      ret = 0;
      break;
    }

    /* print user help for the next round if we get input from tty */
    if (istty_in && istty_out) {
      if (fold_constrained) {
        vrna_message_constraint_options_all();
        vrna_message_input_seq(
          "Input sequence (upper or lower case) followed by structure constraint");
      } else {
        vrna_message_input_seq_simple();
      }
    }
  } while (1);

  return ret;
}


static void
process_record(struct record_data *record)
{
  unsigned int          length;
  struct options        *opt;
  char                  *rec_sequence, *mfe_structure;
  double                min_en, energy;
  vrna_fold_compound_t  *vc;
  struct output_stream  *o_stream;

  opt = record->options;

  rec_sequence = strdup(record->sequence);

  /* convert DNA alphabet to RNA if not explicitely switched off */
  if (!opt->noconv) {
    vrna_seq_toRNA(rec_sequence);
    vrna_seq_toRNA(record->sequence);
  }

  /* convert sequence to uppercase letters only */
  vrna_seq_toupper(rec_sequence);

  vc = vrna_fold_compound(rec_sequence, &(opt->md), VRNA_OPTION_DEFAULT);

  length = vc->length;

  if ((opt->md.circ) && (vrna_rotational_symmetry(rec_sequence) > 1))
    vrna_message_warning("Input sequence %ld is rotationally symmetric! "
                         "Symmetry correction might be required to compute actual MFE and equilibrium properties!",
                         record->number);

  /* retrieve string stream, 6*length should be enough memory to start with */
  o_stream = get_output_stream(6 * length,
                               opt,
                               record->SEQ_ID,
                               record->input_filename);

  if (record->tty)
    vrna_message_info(stdout, "length = %d\n", length);

  mfe_structure = (char *)vrna_alloc(sizeof(char) * (length + 1));

  /* parse the rest of the current dataset to obtain a structure constraint */
  if (fold_constrained) {
    apply_constraints(vc,
                      opt->constraint_file,
                      (const char **)record->rest,
                      record->multiline_input,
                      opt->constraint_enforce,
                      opt->constraint_canonical);
  }

  if (opt->shape) {
    vrna_constraints_add_SHAPE(vc,
                               opt->shape_file,
                               opt->shape_method,
                               opt->shape_conversion,
                               opt->verbose,
                               VRNA_OPTION_DEFAULT);
  }

  if (opt->ligandMotif) {
    add_ligand_motif(vc,
                     opt->ligandMotif,
                     opt->verbose,
                     VRNA_OPTION_MFE | ((opt->pf) ? VRNA_OPTION_PF : 0));
  }

  if (opt->cmds)
    vrna_commands_apply(vc,
                        opt->cmds,
                        VRNA_CMD_PARSE_DEFAULTS);

  /*
   ########################################################
   # begin actual computations
   ########################################################
   */


  /* put header + sequence into output string stream */
  vrna_cstr_print_fasta_header(o_stream->data, record->id);
  vrna_cstr_printf(o_stream->data, "%s\n", record->sequence);

  min_en = (double)vrna_mfe(vc, mfe_structure);

  /* check whether the constraint allows for any solution */
  if ((fold_constrained && opt->constraint_file) || (opt->cmds)) {
    if (min_en == (double)(INF / 100.)) {
      vrna_message_error(
        "Supplied structure constraints create empty solution set for sequence:\n%s",
        record->sequence);
      exit(EXIT_FAILURE);
    }
  }

  if (!opt->lucky) {
    vrna_cstr_printf_structure(o_stream->data,
                               mfe_structure,
                               record->tty ?  "\n minimum free energy = %6.2f kcal/mol" : " (%6.2f)",
                               min_en);

    if (opt->verbose) {
      if (opt->ligandMotif)
        print_ligand_motifs(vc, mfe_structure, "MFE", o_stream->data);

      if (vc->domains_up) {
        vrna_ud_motif_t *m = vrna_ud_motifs_MFE(vc, mfe_structure);
        print_ud_motifs(vc, m, "MFE", o_stream->data);
        free(m);
      }
    }

    if (!opt->noPS) {
      postscript_layout(vc,
                        record->sequence,
                        mfe_structure,
                        record->SEQ_ID,
                        opt->ligandMotif,
                        opt->filename_delim,
                        opt->verbose);
    }
  }

  if (length > 2000)
    vrna_mx_mfe_free(vc);

  if (opt->pf) {
    char *pf_struc = (char *)vrna_alloc(sizeof(char) * (length + 1));
    if (vc->params->model_details.dangles % 2) {
      int dang_bak = vc->params->model_details.dangles;
      vc->params->model_details.dangles = 2;   /* recompute with dangles as in pf_fold() */
      min_en                            = vrna_eval_structure(vc, mfe_structure);
      vc->params->model_details.dangles = dang_bak;
    }

    vrna_exp_params_rescale(vc, &min_en);

    if (length > 2000)
      vrna_message_info(stderr, "scaling factor %f", vc->exp_params->pf_scale);

    energy = (double)vrna_pf(vc, pf_struc);

    /* in case we abort because of floating point errors */
    if (length > 1600)
      vrna_message_info(stderr, "free energy = %8.2f", energy);

    if (opt->lucky) {
      ImFeelingLucky(vc,
                     record->sequence,
                     record->SEQ_ID,
                     opt->noPS,
                     opt->filename_delim,
                     o_stream->data,
                     record->tty);
    } else if (opt->md.compute_bpp) {
      vrna_cstr_printf_structure(o_stream->data,
                                 pf_struc,
                                 record->tty ? "\n free energy of ensemble = %6.2f kcal/mol" : " [%6.2f]",
                                 energy);

      char  *filename_dotplot = NULL;
      plist *pl1, *pl2;

      /* generate initial element probability lists for dot-plot */
      pl1 = vrna_plist_from_probs(vc, opt->bppmThreshold);
      pl2 = vrna_plist(mfe_structure, 0.95 * 0.95);

      /* add ligand motif annotation if necessary */
      if (opt->ligandMotif)
        add_ligand_motifs_dot(vc, &pl1, &pl2, mfe_structure);

      /* generate dot-plot file name */
      filename_dotplot = generate_filename("%s%sdp.ps",
                                           "dot.ps",
                                           record->SEQ_ID,
                                           opt->filename_delim);

      if (filename_dotplot) {
        THREADSAFE_FILE_OUTPUT(
          vrna_plot_dp_EPS(filename_dotplot,
                           record->sequence,
                           pl1,
                           pl2,
                           NULL,
                           VRNA_PLOT_PROBABILITIES_DEFAULT));
      }

      free(filename_dotplot);
      free(pl2);

      /* compute stack probabilities and generate dot-plot */
      if (opt->md.compute_bpp == 2) {
        char *filename_stackplot = generate_filename("%s%sdp2.ps",
                                                     "dot2.ps",
                                                     record->SEQ_ID,
                                                     opt->filename_delim);

        pl2 = vrna_stack_prob(vc, 1e-5);

        if (filename_stackplot) {
          THREADSAFE_FILE_OUTPUT(
            PS_dot_plot_list(record->sequence, filename_stackplot, pl1, pl2,
                             "Probabilities for stacked pairs (i,j)(i+1,j-1)"));
        }

        free(pl2);
        free(filename_stackplot);
      }

      free(pl1);

      /* compute centroid structure */
      compute_centroid(vc, opt->ligandMotif, opt->verbose, o_stream->data);

      /* compute MEA structure */
      if (opt->MEA) {
        compute_MEA(vc,
                    opt->MEAgamma,
                    opt->ligandMotif,
                    opt->verbose,
                    o_stream->data);
      }

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
                                 energy,
                                 vrna_pr_energy(vc, min_en));
    }

    free(pf_struc);
  }

  /* print what we've collected in output charstream */
  if (opt->output_queue) {
    if (o_stream->individual) {
      /* output immediately */
      ATOMIC_BLOCK(flush_cstr_callback(NULL, record->number, (void *)o_stream));

      /* use dummy element for insert into queue */
      o_stream = NULL;
    }

    vrna_ostream_provide(opt->output_queue, record->number, (void *)o_stream);
  } else {
    ATOMIC_BLOCK(flush_cstr_callback(NULL, record->number, (void *)o_stream));
  }

  /* clean up */
  vrna_fold_compound_free(vc);
  free(record->id);
  free(record->SEQ_ID);
  free(record->sequence);
  free(rec_sequence);
  free(mfe_structure);

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
apply_constraints(vrna_fold_compound_t  *fc,
                  const char            *constraints_file,
                  const char            **rec_rest,
                  int                   maybe_multiline,
                  int                   enforceConstraints,
                  int                   canonicalBPonly)
{
  if (constraints_file) {
    /** [Adding hard constraints from file] */
    vrna_constraints_add(fc, constraints_file, VRNA_OPTION_DEFAULT);
    /** [Adding hard constraints from file] */
  } else {
    char          *cstruc   = NULL;
    unsigned int  length    = fc->length;
    unsigned int  coptions  = (maybe_multiline) ? VRNA_OPTION_MULTILINE : 0;
    cstruc = vrna_extract_record_rest_structure((const char **)rec_rest, 0, coptions);
    unsigned int  cl = (cstruc) ? strlen(cstruc) : 0;

    if (cl == 0)
      vrna_message_warning("structure constraint is missing");
    else if (cl < length)
      vrna_message_warning("structure constraint is shorter than sequence");
    else if (cl > length)
      vrna_message_error("structure constraint is too long");

    if (cstruc) {
      /** [Adding hard constraints from pseudo dot-bracket] */
      unsigned int constraint_options = VRNA_CONSTRAINT_DB_DEFAULT;

      if (enforceConstraints)
        constraint_options |= VRNA_CONSTRAINT_DB_ENFORCE_BP;

      if (canonicalBPonly)
        constraint_options |= VRNA_CONSTRAINT_DB_CANONICAL_BP;

      vrna_constraints_add(fc, (const char *)cstruc, constraint_options);
      /** [Adding hard constraints from pseudo dot-bracket] */

      free(cstruc);
    }
  }
}


static void
compute_MEA(vrna_fold_compound_t  *fc,
            double                MEAgamma,
            const char            *ligandMotif,
            int                   verbose,
            vrna_cstr_t           rec_output)
{
  char  *structure;
  float mea, mea_en;
  /*  this is a hack since vrna_plist_from_probs() always resolves g-quad pairs,
   *  while MEA_seq() still expects unresolved gquads */
  int   gq = fc->exp_params->model_details.gquad;

  /* we need to create a string as long as the sequence for the MEA implementation :( */
  structure = strdup(fc->sequence);

  fc->exp_params->model_details.gquad = 0;
  plist *pl = vrna_plist_from_probs(fc, 1e-4 / (1 + MEAgamma));
  fc->exp_params->model_details.gquad = gq;

  if (gq)
    mea = MEA_seq(pl, fc->sequence, structure, MEAgamma, fc->exp_params);
  else
    mea = MEA(pl, structure, MEAgamma);

  mea_en = vrna_eval_structure(fc, (const char *)structure);

  vrna_cstr_printf_structure(rec_output, structure, " {%6.2f MEA=%.2f}", mea_en, mea);

  if ((ligandMotif) && (verbose))
    print_ligand_motifs(fc, structure, "MEA", rec_output);

  if ((fc->domains_up) && (verbose)) {
    vrna_ud_motif_t *m = vrna_ud_motifs_MEA(fc, structure, pl);
    print_ud_motifs(fc, m, "MEA", rec_output);
    free(m);
  }

  free(pl);
  free(structure);
}


static void
compute_centroid(vrna_fold_compound_t *fc,
                 const char           *ligandMotif,
                 int                  verbose,
                 vrna_cstr_t          rec_output)
{
  char    *cent;
  double  cent_en, dist;

  cent    = vrna_centroid(fc, &dist);
  cent_en = vrna_eval_structure(fc, (const char *)cent);

  vrna_cstr_printf_structure(rec_output, cent, " {%6.2f d=%.2f}", cent_en, dist);

  if ((ligandMotif) && (verbose))
    print_ligand_motifs(fc, cent, "centroid", rec_output);

  if ((fc->domains_up) && (verbose)) {
    vrna_ud_motif_t *m = vrna_ud_motifs_centroid(fc, cent);
    print_ud_motifs(fc, m, "centroid", rec_output);
    free(m);
  }

  free(cent);
}


static void
add_ligand_motif(vrna_fold_compound_t *vc,
                 char                 *motifstring,
                 int                  verbose,
                 unsigned int         options)
{
  int   r, l, error;
  char  *seq, *str, *ptr;
  float energy;

  l   = strlen(motifstring);
  seq = vrna_alloc(sizeof(char) * (l + 1));
  str = vrna_alloc(sizeof(char) * (l + 1));

  error = 1;

  if (motifstring) {
    error = 0;
    /* parse sequence */
    for (r = 0, ptr = motifstring; *ptr != '\0'; ptr++) {
      if (*ptr == ',')
        break;

      seq[r++] = toupper(*ptr);
    }
    seq[r]  = '\0';
    seq     = vrna_realloc(seq, sizeof(char) * (strlen(seq) + 1));

    for (ptr++, r = 0; *ptr != '\0'; ptr++) {
      if (*ptr == ',')
        break;

      str[r++] = *ptr;
    }
    str[r]  = '\0';
    str     = vrna_realloc(str, sizeof(char) * (strlen(seq) + 1));

    ptr++;
    if (!(sscanf(ptr, "%f", &energy) == 1)) {
      vrna_message_warning("Energy contribution in ligand motif missing!");
      error = 1;
    }

    if (strlen(seq) != strlen(str)) {
      vrna_message_warning("Sequence and structure length in ligand motif have unequal lengths!");
      error = 1;
    }

    if (strlen(seq) == 0) {
      vrna_message_warning("Sequence length in ligand motif is zero!");
      error = 1;
    }

    if (!error && verbose)
      vrna_message_info(stderr, "Read ligand motif: %s, %s, %f", seq, str, energy);
  }

  if (error || (!vrna_sc_add_hi_motif(vc, seq, str, energy, options)))
    vrna_message_warning("Malformatted ligand motif! Skipping stabilizing motif.");

  free(seq);
  free(str);
}


static char *
annotate_ligand_motif(vrna_fold_compound_t  *vc,
                      const char            *structure)
{
  char            *annote;
  vrna_sc_motif_t *motifs, *m_ptr;

  annote  = NULL;
  motifs  = vrna_sc_ligand_detect_motifs(vc, structure);

  if (motifs) {
    for (m_ptr = motifs; m_ptr->i; m_ptr++) {
      char *tmp_string, *annotation;
      annotation  = NULL;
      tmp_string  = annote;

      if (m_ptr->i != m_ptr->k) {
        annotation = vrna_strdup_printf(" %d %d %d %d 1. 0 0 BFmark",
                                        m_ptr->i,
                                        m_ptr->j,
                                        m_ptr->k,
                                        m_ptr->l);
      } else {
        annotation = vrna_strdup_printf(" %d %d 1. 0 0 Fomark",
                                        m_ptr->i,
                                        m_ptr->j);
      }

      if (tmp_string)
        annote = vrna_strdup_printf("%s %s", tmp_string, annotation);
      else
        annote = strdup(annotation);

      free(tmp_string);
      free(annotation);
    }
  }

  free(motifs);

  return annote;
}


static void
print_ligand_motifs(vrna_fold_compound_t  *vc,
                    const char            *structure,
                    const char            *structure_name,
                    vrna_cstr_t           buf)
{
  vrna_sc_motif_t *motifs, *m_ptr;

  motifs = vrna_sc_ligand_detect_motifs(vc, structure);

  if (motifs) {
    for (m_ptr = motifs; m_ptr->i; m_ptr++) {
      if (m_ptr->i != m_ptr->k) {
        /* put annotation into output vrna_cstr_t */
        vrna_cstr_message_info(buf,
                               "specified motif detected in %s structure: [%d:%d] & [%d:%d]",
                               structure_name,
                               m_ptr->i,
                               m_ptr->k,
                               m_ptr->l,
                               m_ptr->j);
      } else {
        /* put annotation into output vrna_cstr_t */
        vrna_cstr_message_info(buf,
                               "specified motif detected in %s structure: [%d:%d]",
                               structure_name,
                               m_ptr->i,
                               m_ptr->j);
      }
    }
  }

  free(motifs);
}


static char *
annotate_ud_motif(vrna_fold_compound_t  *vc,
                  vrna_ud_motif_t       *motifs)
{
  int   m, i, size;
  char  *annote;

  m       = 0;
  annote  = NULL;

  if (motifs) {
    while (motifs[m].start != 0) {
      char  *tmp_string = annote;
      i     = motifs[m].start;
      size  = vc->domains_up->motif_size[motifs[m].number];
      char  *annotation;

      annotation = vrna_strdup_printf(" %d %d 12 0.4 0.65 0.95 omark", i, i + size - 1);

      if (tmp_string)
        annote = vrna_strdup_printf("%s %s", tmp_string, annotation);
      else
        annote = strdup(annotation);

      free(tmp_string);
      free(annotation);
      m++;
    }
  }

  return annote;
}


static void
print_ud_motifs(vrna_fold_compound_t  *vc,
                vrna_ud_motif_t       *motifs,
                const char            *structure_name,
                vrna_cstr_t           buf)
{
  int m, i, size;

  m = 0;

  if (motifs) {
    while (motifs[m].start != 0) {
      i     = motifs[m].start;
      size  = vc->domains_up->motif_size[motifs[m].number];

      /* put annotation into output vrna_cstr_t */
      vrna_cstr_message_info(buf,
                             "ud motif %d detected in %s structure: [%d:%d]",
                             motifs[m].number,
                             structure_name,
                             i,
                             i + size - 1);
      m++;
    }
  }
}


static void
add_ligand_motifs_dot(vrna_fold_compound_t  *fc,
                      vrna_ep_t             **prob_list,
                      vrna_ep_t             **mfe_list,
                      const char            *structure)
{
  vrna_sc_motif_t *motifs;

  /* append motif positions to the plists of base pair probabilities */
  motifs = vrna_sc_ligand_get_all_motifs(fc);
  if (motifs) {
    add_ligand_motifs_to_list(prob_list, motifs);
    free(motifs);
  }

  /* now scan for the motif in MFE structure again */
  motifs = vrna_sc_ligand_detect_motifs(fc, structure);
  if (motifs) {
    add_ligand_motifs_to_list(mfe_list, motifs);
    free(motifs);
  }
}


static void
add_ligand_motifs_to_list(vrna_ep_t       **list,
                          vrna_sc_motif_t *motifs)
{
  unsigned int    cnt, add, size;
  vrna_ep_t       *ptr;
  vrna_sc_motif_t *m_ptr;

  cnt = 0;
  add = 10;

  /* get current size of list */
  for (size = 0, ptr = (*list); ptr->i; size++, ptr++);

  /* increase length of list */
  (*list) = vrna_realloc((*list), sizeof(vrna_ep_t) * (size + add + 1));

  for (m_ptr = motifs; m_ptr->i; m_ptr++) {
    if (m_ptr->i == m_ptr->k) {
      /* hairpin motif */
      (*list)[size + cnt].i     = m_ptr->i;
      (*list)[size + cnt].j     = m_ptr->j;
      (*list)[size + cnt].p     = 0.95 * 0.95;
      (*list)[size + cnt].type  = VRNA_PLIST_TYPE_H_MOTIF;
      cnt++;
      if (cnt == add) {
        add += 10;
        /* increase length of (*prob_list) */
        (*list) = vrna_realloc((*list), sizeof(vrna_ep_t) * (size + add + 1));
      }
    } else {
      /* interior loop motif */
      (*list)[size + cnt].i     = m_ptr->i;
      (*list)[size + cnt].j     = m_ptr->j;
      (*list)[size + cnt].p     = 0.95 * 0.95;
      (*list)[size + cnt].type  = VRNA_PLIST_TYPE_I_MOTIF;
      cnt++;
      (*list)[size + cnt].i     = m_ptr->k;
      (*list)[size + cnt].j     = m_ptr->l;
      (*list)[size + cnt].p     = 0.95 * 0.95;
      (*list)[size + cnt].type  = VRNA_PLIST_TYPE_I_MOTIF;
      cnt++;
      if (cnt == add) {
        add += 10;
        /* increase length of (*prob_list) */
        (*list) = vrna_realloc((*list), sizeof(vrna_ep_t) * (size + add + 1));
      }
    }
  }

  /* resize pl1 to actual needs */
  (*list)               = vrna_realloc((*list), sizeof(vrna_ep_t) * (size + cnt + 1));
  (*list)[size + cnt].i = 0;
  (*list)[size + cnt].j = 0;
}
