/*
 *                    Heat Capacity of RNA molecule
 *
 *                  c Ivo Hofacker and Peter Stadler
 *                        Vienna RNA package
 *
 *
 *          calculates specific heat using C = - T d^2/dT^2 G(T)
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/params/io.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/datastructures/char_stream.h"
#include "ViennaRNA/datastructures/stream_output.h"
#include "ViennaRNA/color_output.inc"

#include "RNAheat_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"
#include "parallel_helpers.h"


#define MAXWIDTH  201

struct options {
  int             filename_full;
  int             noconv;
  float           T_min;
  float           T_max;
  float           h;
  int             mpoints;
  vrna_md_t       md;
  dataset_id      id_control;

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
  char            *input_filename;
  int             multiline_input;
  struct options  *options;
  int             tty;
};


struct output_stream {
  vrna_cstr_t data;
  vrna_cstr_t err;
};


PRIVATE float
ddiff(float f[],
      float h,
      int   m);


static int
process_input(FILE            *input_stream,
              const char      *input_filename,
              struct options  *opt);


static void
process_record(struct record_data *record);


void
init_default_options(struct options *opt)
{
  opt->filename_full  = 0;
  opt->noconv         = 0;

  vrna_md_set_default(&(opt->md));
  opt->md.backtrack   = 0;
  opt->md.compute_bpp = 0;

  opt->T_min    = 0.;
  opt->T_max    = 100.;
  opt->h        = 1;
  opt->mpoints  = 2;

  opt->jobs               = 1;
  opt->keep_order         = 1;
  opt->next_record_number = 0;
  opt->output_queue       = NULL;
}


static char **
collect_unnamed_options(struct RNAheat_args_info  *ggostruct,
                        int                       *num_files)
{
  char          **input_files = NULL;
  unsigned int  i;

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
append_input_files(struct RNAheat_args_info *ggostruct,
                   char                     **files,
                   int                      *numfiles)
{
  unsigned int i;

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
  struct RNAheat_args_info  args_info;
  char                      **input_files;
  int                       num_input;
  struct options            opt;

  num_input = 0;

  init_default_options(&opt);

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNAheat_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /*
   *  - dangles
   *  - special_hp
   *  - gquad
   *  - energy_set
   *  - ns_bases
   *  - parameter file
   */
  ggo_get_md_eval(args_info, opt.md);

  /* check dangle model */
  if (!((opt.md.dangles == 0) || (opt.md.dangles == 2))) {
    vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    opt.md.dangles = 2;
  }

  /*
   *  - noLP
   *  - noGU
   *  - noGUclosure
   *  - maxBPspan
   */
  ggo_get_md_fold(args_info, opt.md);
  ggo_get_circ(args_info, opt.md.circ);

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (args_info.noconv_given)
    opt.noconv = 1;

  /* Tmin */
  if (args_info.Tmin_given)
    opt.T_min = args_info.Tmin_arg;

  /* Tmax */
  if (args_info.Tmax_given)
    opt.T_max = args_info.Tmax_arg;

  /* step size */
  if (args_info.stepsize_given)
    opt.h = args_info.stepsize_arg;

  /* ipoints */
  if (args_info.ipoints_given) {
    opt.mpoints = args_info.ipoints_arg;
    if (opt.mpoints < 1)
      opt.mpoints = 1;

    if (opt.mpoints > 100)
      opt.mpoints = 100;
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
      "This version of RNAheat has been built without parallel input processing capabilities");
#endif

    if (args_info.unordered_given)
      opt.keep_order = 0;
  }

  input_files = collect_unnamed_options(&args_info, &num_input);
  input_files = append_input_files(&args_info, input_files, &num_input);

  /* parse options for ID manipulation */
  ggo_get_id_control(args_info, opt.id_control, "Sequence", "sequence", "_", 4, 1);

  /* free allocated memory of command line data structure */
  RNAheat_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */
  if (opt.md.circ && opt.md.gquad)
    vrna_message_error("G-Quadruplex support is currently not available for circular RNA structures");

  if (opt.keep_order)
    opt.output_queue = vrna_ostream_init(&flush_cstr_callback, NULL);

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


  free(input_files);

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

  unsigned int  read_opt = VRNA_INPUT_NO_REST;

  /* print user help if we get input from tty */
  if (istty_in && istty_out) {
    vrna_message_input_seq_simple();
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
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
    record->multiline_input = maybe_multiline;
    record->options         = opt;
    record->tty             = istty_in && istty_out;
    record->input_filename  = (input_filename) ? strdup(input_filename) : NULL;

    if (opt->output_queue)
      vrna_ostream_request(opt->output_queue, opt->next_record_number++);

    RUN_IN_PARALLEL(process_record, record);

    /* print user help for the next round if we get input from tty */
    if (istty_in && istty_out)
      vrna_message_input_seq_simple();
  } while (1);

  return ret;
}


static void
process_record(struct record_data *record)
{
  char                  *rec_sequence;
  int                   i, n, m;
  float                 hc, F[MAXWIDTH], T_min, T_max, h;
  double                min_en;
  vrna_fold_compound_t  *fc;
  vrna_md_t             md;

  struct options        *opt;
  struct output_stream  *o_stream;

  opt           = record->options;
  o_stream      = (struct output_stream *)vrna_alloc(sizeof(struct output_stream));
  rec_sequence  = strdup(record->sequence);

  T_min = opt->T_min;
  T_max = opt->T_max;
  h     = opt->h;
  m     = opt->mpoints;

  /* convert DNA alphabet to RNA if not explicitely switched off */
  if (!opt->noconv) {
    vrna_seq_toRNA(rec_sequence);
    vrna_seq_toRNA(record->sequence);
  }

  /* convert sequence to uppercase letters only */
  vrna_seq_toupper(rec_sequence);

  fc = vrna_fold_compound(rec_sequence,
                          &(opt->md),
                          VRNA_OPTION_DEFAULT);

  n = (int)fc->length;

  /* retrieve string stream bound to stdout, 6*length should be enough memory to start with */
  o_stream->data = vrna_cstr(6 * n, stdout);
  /* retrieve string stream bound to stderr for any info messages */
  o_stream->err = vrna_cstr(n, stderr);

  if (record->tty)
    vrna_message_info(stdout, "length = %d", n);

  /*
   ########################################################
   # begin actual computations
   ########################################################
   */
  vrna_cstr_print_fasta_header(o_stream->data, record->id);

  md = fc->params->model_details;

  /* required for vrna_exp_param_rescale() in subsequent calls */
  md.sfact = 1.;

  md.temperature = T_min - m * h;
  vrna_params_reset(fc, &md);

  /* initialize_fold(length); <- obsolete */
  min_en  = (double)vrna_mfe(fc, NULL);
  min_en  *= md.sfact;

  vrna_exp_params_rescale(fc, &min_en);

  for (i = 0; i < 2 * m + 1; i++) {
    F[i] = vrna_pf(fc, NULL);
    /* increase temperature */
    md.temperature += h;
    /* reset all energy parameters according to temperature changes */
    vrna_params_reset(fc, &md);

    min_en = F[i] + h * 0.00727 * n;

    vrna_exp_params_rescale(fc, &min_en);
  }

  while (md.temperature <= (T_max + m * h + h)) {
    hc = -ddiff(F, h, m) * (md.temperature + K0 - m * h - h);

    vrna_cstr_printf_tbody(o_stream->data,
                           "%g\t%g",
                           (md.temperature - m * h - h),
                           hc);

    for (i = 0; i < 2 * m; i++)
      F[i] = F[i + 1];

    F[2 * m] = vrna_pf(fc, NULL);

    /*       printf("%g\n", F[2*m]);*/
    md.temperature += h;

    vrna_params_reset(fc, &md);

    min_en = F[i] + h * 0.00727 * n;

    vrna_exp_params_rescale(fc, &min_en);
  }

  if (opt->output_queue)
    vrna_ostream_provide(opt->output_queue, record->number, (void *)o_stream);
  else
    flush_cstr_callback(NULL, 0, (void *)o_stream);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(record->id);
  free(record->SEQ_ID);
  free(record->sequence);
  free(rec_sequence);
  free(record->input_filename);
  free(record);
}


/* ------------------------------------------------------------------------- */

PRIVATE float
ddiff(float f[],
      float h,
      int   m)
{
  int   i;
  float fp, A, B;

  A = (float)(m * (m + 1) * (2 * m + 1) / 3);                                     /* 2*sum(x^2) */
  B = (float)(m * (m + 1) * (2 * m + 1)) * (float)(3 * m * m + 3 * m - 1) / 15.;  /* 2*sum(x^4) */

  fp = 0.;
  for (i = 0; i < 2 * m + 1; i++)
    fp += f[i] * (A - (float)((2 * m + 1) * (i - m) * (i - m)));

  fp /= ((A * A - B * ((float)(2 * m + 1))) * h * h / 2.);
  return (float)fp;
}
