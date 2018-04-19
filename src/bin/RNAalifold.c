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
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/alifold.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/MEA.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/constraints_SHAPE.h"
#include "ViennaRNA/plot_utils.h"
#include "RNAalifold_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"

#include "ViennaRNA/color_output.inc"

#define DBL_ROUND(a, digits) (round((a) * pow(10., (double)(digits))) / pow(10., (double)(digits)))

struct options {
  int         filename_full;
  char        *filename_delim;
  int         pf;
  int         noPS;
  int         noconv;
  int         MEA;
  double      MEAgamma;
  double      bppmThreshold;
  int         verbose;
  int         quiet;
  vrna_md_t   md;

  dataset_id  id_control;
  int         continuous_names;

  int         n_back;
  int         eval_en;

  int         color;
  int         aln_PS;
  int         aln_PS_cols;
  int         mis;
  int         sci;
  int         endgaps;

  int         aln_out;
  char        *aln_out_prefix;

  char        *constraint_file;
  int         constraint_SScons;
  int         constraint_batch;
  int         constraint_enforce;

  int         shape;
  char        **shape_files;
  char        *shape_method;
  int         *shape_file_association;
};


PRIVATE void  print_pi(const vrna_pinfo_t pi,
                       FILE               *file);


PRIVATE void  print_aliout(vrna_fold_compound_t *vc,
                           plist                *pl,
                           double               threshold,
                           char                 *mfe,
                           FILE                 *aliout);


PRIVATE void  mark_endgaps(char *seq,
                           char egap);


static void
apply_constraints(vrna_fold_compound_t  *fc,
                  const char            *constraints_file,
                  const char            *cstruc,
                  const char            *SS_cons,
                  struct options        *opt);


void
init_default_options(struct options *opt)
{
  opt->filename_full  = 0;
  opt->filename_delim = NULL;
  opt->pf             = 0;
  opt->noPS           = 0;
  opt->noconv         = 0;
  opt->MEA            = 0;
  opt->MEAgamma       = 1.;
  opt->bppmThreshold  = 1e-6;
  opt->verbose        = 0;
  opt->quiet          = 0;
  set_model_details(&(opt->md));

  opt->continuous_names = 0;

  opt->n_back           = 0;
  opt->eval_en          = 0;

  opt->color            = 0;
  opt->aln_PS           = 0;
  opt->aln_PS_cols      = 60;
  opt->mis              = 0;
  opt->sci              = 0;
  opt->endgaps          = 0;

  opt->aln_out          = 0;
  opt->aln_out_prefix   = NULL;

  opt->constraint_file = NULL;
  opt->constraint_SScons  = 0;
  opt->constraint_batch = 0;
  opt->constraint_enforce = 0;

  opt->shape            = 0;
  opt->shape_files      = NULL;
  opt->shape_file_association = NULL;
  opt->shape_method     = NULL;
}


int
main(int  argc,
     char *argv[])
{
  struct RNAalifold_args_info args_info;
  FILE                        *clust_file;
  unsigned int                input_type, longest_string, input_format_options,
                              aln_options;
  char                        *input_string, *string, *structure, *cstruc, **AS, **names,
                              *filename_plot,
                              *filename_dot, *filename_aln, *filename_out, *filename_in,
                              *tmp_id, *tmp_structure, *tmp_string, **input_files;
  int                         s, n_seq, i, length, istty_out, istty_in, istty, 
                              tmp_number, input_file_num;
  double                      min_en, real_en, cov_en;
  vrna_fold_compound_t        *vc;
  struct  options             opt;

  input_file_num = 0;

  init_default_options(&opt);

  string      = NULL;
  structure   = NULL;
  cstruc      = NULL;

  filename_plot           = NULL;
  filename_dot            = NULL;
  filename_aln            = NULL;
  filename_out            = NULL;
  filename_in             = NULL;
  input_files             = NULL;
  tmp_id                  = NULL;
  tmp_structure           = NULL;
  input_format_options    = VRNA_FILE_FORMAT_MSA_CLUSTAL; /* default to ClustalW format */
  aln_options             = VRNA_ALN_UPPERCASE;           /* we always require uppercase sequence letters internally */
  vc                      = NULL;
  clust_file              = stdin;

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
    opt.constraint_SScons  = 1;
  }

  /* do not produce postscript output */
  if (args_info.noPS_given)
    opt.noPS = 1;

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
    opt.md.cv_fact = cv_fact = args_info.cfactor_arg;

  /* set nfactor */
  if (args_info.nfactor_given)
    opt.md.nc_fact = nc_fact = args_info.nfactor_arg;

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
    opt.eval_en     = 1;
    vrna_init_rand();
  }

  if (args_info.ribosum_file_given) {
    RibosumFile = strdup(args_info.ribosum_file_arg);
    opt.md.ribo     = ribo = 1;
  }

  if (args_info.ribosum_scoring_given) {
    RibosumFile = NULL;
    opt.md.ribo     = 1;
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

    opt.shape             = 1;
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

  /* sci computation */
  if (args_info.sci_given)
    opt.sci = 1;

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
  else
    aln_options |= VRNA_ALN_RNA;

  /* free allocated memory of command line data structure */
  RNAalifold_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */
  istty     = isatty(fileno(stdout)) && isatty(fileno(stdin));
  istty_out = isatty(fileno(stdout));
  istty_in  = isatty(fileno(stdin)) && (!filename_in);

  if (opt.md.circ && opt.md.gquad)
    vrna_message_error("G-Quadruplex support is currently not available for circular RNA structures");

  if (opt.md.circ && opt.md.noLP)
    vrna_message_warning("Depending on the origin of the circular sequence, "
                         "some structures may be missed when using --noLP\n"
                         "Try rotating your sequence a few times\n");

  /*
   ########################################################
   # handle user input from 'stdin' if necessary
   ########################################################
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

  if (fold_constrained && (!opt.constraint_file) && (!opt.constraint_SScons)) {
    if (isatty(fileno(stdin))) {
      vrna_message_constraint_options_all();
      vrna_message_input_seq("");
    }

    input_string  = NULL;
    input_type    = get_input_line(&input_string, VRNA_INPUT_NOSKIP_COMMENTS);
    if (input_type & VRNA_INPUT_QUIT) {
      return 0;
    } else if ((input_type & VRNA_INPUT_MISC) && (strlen(input_string) > 0)) {
      cstruc = strdup(input_string);
      free(input_string);
    }
  }

  long int first_alignment_number = get_current_id(opt.id_control);

  if (opt.aln_out) {
    if (!opt.aln_out_prefix) {
      opt.aln_out_prefix = strdup("RNAalifold_results.stk");
    } else {
      char *tmp = vrna_strdup_printf("%s.stk", opt.aln_out_prefix);
      free(opt.aln_out_prefix);
      opt.aln_out_prefix = vrna_filename_sanitize(tmp, opt.filename_delim);
      free(tmp);
    }
  }

  while (!feof(clust_file)) {
    char  *MSA_ID     = NULL;
    char  **MSA_orig  = NULL;

    fflush(stdout);

    if (istty && (clust_file == stdin)) {
      switch (input_format_options) {
        case VRNA_FILE_FORMAT_MSA_CLUSTAL:
          vrna_message_input_seq("Input aligned sequences in ClustalW format\n"
                                 "press Ctrl+d when finished to indicate the end of your input)");
          break;

        case VRNA_FILE_FORMAT_MSA_STOCKHOLM:
          vrna_message_input_seq("Input aligned sequences in Stockholm format (Insert one alignment at a time!)\n"
                                 "press Ctrl+d when finished to indicate the end of your input)");
          break;

        case VRNA_FILE_FORMAT_MSA_FASTA:
          vrna_message_input_seq("Input aligned sequences in FASTA format\n"
                                 "press Ctrl+d when finished to indicate the end of your input)");
          break;

        case VRNA_FILE_FORMAT_MSA_MAF:
          vrna_message_input_seq("Input aligned sequences in MAF format (Insert one alignment at a time!)\n"
                                 "press Ctrl+d when finished to indicate the end of your input)");
          break;

        default:
          vrna_message_error("Which input format are you using?");
          break;
      }
    }

    if (opt.quiet)
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
    MSA_ID = fileprefix_from_id_alifold(tmp_id, opt.id_control, opt.continuous_names);

    /* construct output file names */
    if (MSA_ID) {
      /* construct file names */
      filename_plot = vrna_strdup_printf("%s%sss.ps", MSA_ID, opt.filename_delim);
      filename_dot  = vrna_strdup_printf("%s%sdp.ps", MSA_ID, opt.filename_delim);
      filename_aln  = vrna_strdup_printf("%s%saln.ps", MSA_ID, opt.filename_delim);
      filename_out  = vrna_strdup_printf("%s%sali.out", MSA_ID, opt.filename_delim);

      /* sanitize file names */
      tmp_string = vrna_filename_sanitize(filename_plot, opt.filename_delim);
      free(filename_plot);
      filename_plot = tmp_string;
      tmp_string    = vrna_filename_sanitize(filename_dot, opt.filename_delim);
      free(filename_dot);
      filename_dot  = tmp_string;
      tmp_string    = vrna_filename_sanitize(filename_aln, opt.filename_delim);
      free(filename_aln);
      filename_aln  = tmp_string;
      tmp_string    = vrna_filename_sanitize(filename_out, opt.filename_delim);
      free(filename_out);
      filename_out = tmp_string;
    } else {
      filename_plot = vrna_strdup_printf("alirna.ps");
      filename_dot  = vrna_strdup_printf("alidot.ps");
      filename_aln  = vrna_strdup_printf("aln.ps");
      filename_out  = vrna_strdup_printf("alifold.out");
    }

    /*
     ##############################################################
     # done with retrieving alignment, now init everything properly
     ##############################################################
     */

    /*
     *  store original alignment and create a new one for internal computations
     *  which require all-uppercase letters, and RNA/DNA alphabet
     */
    MSA_orig  = AS;
    AS        = vrna_aln_copy((const char **)MSA_orig, aln_options);
    length    = (int)strlen(AS[0]);
    structure = (char *)vrna_alloc(sizeof(char) * (length + 1));

    if (opt.endgaps)
      for (i = 0; i < n_seq; i++)
        mark_endgaps(AS[i], '~');

    vc = vrna_fold_compound_comparative((const char **)AS, &opt.md, VRNA_OPTION_DEFAULT);

    if (fold_constrained)
      apply_constraints(vc, opt.constraint_file, cstruc, tmp_structure, &opt);

    if (opt.shape) {
      for (s = 0; opt.shape_file_association[s] != -1; s++);

      if ((s != n_seq) && (!opt.quiet))
        vrna_message_warning("Number of sequences in alignment (%d) does not match "
                             "number of provided SHAPE reactivity data files (%d)!",
                             n_seq,
                             s);

      opt.shape_files             = (char **)vrna_realloc(opt.shape_files,
                                                          (n_seq + 1) * sizeof(char *));
      opt.shape_file_association  = (int *)vrna_realloc(opt.shape_file_association,
                                                    (n_seq + 1) * sizeof(int));

      vrna_constraints_add_SHAPE_ali(vc,
                                     opt.shape_method,
                                     (const char **)opt.shape_files,
                                     opt.shape_file_association,
                                     opt.verbose,
                                     VRNA_OPTION_DEFAULT);
    }

    /*
     ########################################################
     # begin actual calculations
     ########################################################
     */

    /* generate consensus sequence */
    string = (opt.mis) ? consens_mis((const char **)AS) : consensus((const char **)AS);

    /* put header + sequence into output string stream */

    print_fasta_header(stdout, MSA_ID);
    fprintf(stdout, "%s\n", string);

    min_en  = vrna_mfe(vc, structure);
    real_en = vrna_eval_structure(vc, structure);
    cov_en  = vrna_eval_covar_structure(vc, structure);


    float sci     = min_en;
    float e_mean  = 0;

    if (opt.sci) {
      for (i = 0; AS[i] != NULL; i++) {
        char *seq = get_ungapped_sequence(AS[i]);
        if (strlen(seq) > 0) {
          vrna_fold_compound_t *fc = vrna_fold_compound(seq, &opt.md, VRNA_OPTION_DEFAULT);
          e_mean += vrna_mfe(fc, NULL);
          vrna_fold_compound_free(fc);
        }

        free(seq);
      }
      e_mean /= i;
      if (e_mean == 0.)
        sci = 0.;
      else
        sci /= e_mean;
    }

    char  *energy_string  = NULL;
    char  *e_individual   = NULL;

    if (opt.shape) {
      e_individual = vrna_strdup_printf("%6.2f + %6.2f + %6.2f",
                                        DBL_ROUND(real_en, 2),
                                        DBL_ROUND(-cov_en, 2),
                                        DBL_ROUND(min_en - real_en + cov_en, 2));
    } else {
      e_individual = vrna_strdup_printf("%6.2f + %6.2f",
                                        DBL_ROUND(real_en, 2),
                                        DBL_ROUND(min_en - real_en, 2));
    }

    if (istty_in) {
      if (opt.sci) {
        energy_string = vrna_strdup_printf(
          "\n minimum free energy = %6.2f kcal/mol (%s)\n SCI = %2.4f",
          DBL_ROUND(min_en, 2),
          e_individual,
          DBL_ROUND(sci, 4));
      } else {
        energy_string = vrna_strdup_printf("\n minimum free energy = %6.2f kcal/mol (%s)",
                                           DBL_ROUND(min_en, 2),
                                           e_individual);
      }
    } else {
      if (opt.sci)
        energy_string = vrna_strdup_printf(" (%6.2f = %s) [sci = %2.4f]",
                                           DBL_ROUND(min_en, 2),
                                           e_individual,
                                           sci);
      else
        energy_string = vrna_strdup_printf(" (%6.2f = %s)",
                                           DBL_ROUND(min_en, 2),
                                           e_individual);
    }

    print_structure(stdout, structure, energy_string);

    free(energy_string);
    free(e_individual);

    if (!opt.noPS) {
      char **A;
      A = vrna_annotate_covar_struct((const char **)AS, structure, &opt.md);

      if (opt.color)
        (void)vrna_file_PS_rnaplot_a(string, structure, filename_plot, A[0], A[1], &opt.md);
      else
        (void)vrna_file_PS_rnaplot_a(string, structure, filename_plot, NULL, A[1], &opt.md);

      free(A[0]);
      free(A[1]);
      free(A);
    }

    if (opt.aln_PS)
      vrna_file_PS_aln(filename_aln, (const char **)MSA_orig, (const char **)names, structure,
                       opt.aln_PS_cols);

    if (opt.aln_out) {
      unsigned int options = VRNA_FILE_FORMAT_MSA_STOCKHOLM
                             | VRNA_FILE_FORMAT_MSA_APPEND;
      if (opt.mis)
        options |= VRNA_FILE_FORMAT_MSA_MIS;

      vrna_file_msa_write((const char *)opt.aln_out_prefix,
                          (const char **)names,
                          (const char **)MSA_orig,
                          MSA_ID,
                          (const char *)structure,
                          "RNAalifold prediction",
                          options);
    }

    /* free mfe arrays */
    vrna_mx_mfe_free(vc);

    if (opt.pf) {
      double energy, kT;
      char  *mfe_struc;

      mfe_struc = strdup(structure);

      vrna_exp_params_rescale(vc, &min_en);
      pf_scale = vc->exp_params->pf_scale;


      kT = vc->exp_params->kT / 1000.;

      if ((length > 2000) && (!opt.quiet))
        vrna_message_info(stderr, "scaling factor %f\n", vc->exp_params->pf_scale);

      fflush(stdout);

      energy = vrna_pf(vc, structure);


      if (opt.n_back > 0) {
        /*stochastic sampling*/
        for (i = 0; i < opt.n_back; i++) {
          char    *s, *e_string = NULL;

#if CHECK_PROBABILITIES
          double  prob2, prob = 1.;
          s = alipbacktrack(&prob);
          double  e = (double)vrna_eval_structure(vc, s);
          e     -= (double)vrna_eval_covar_structure(vc, s);
          prob2 = exp((energy - e) / kT);
          if (opt.eval_en)
            e_string =
              vrna_strdup_printf(" %6g (%6g) %.2f (%.2f)",
                                 prob,
                                 prob2,
                                 -1 * (kT * log(prob) - energy),
                                 e);

#else
          double prob = 1.;
          s = vrna_pbacktrack(vc);

          if (opt.eval_en) {
            double e = (double)vrna_eval_structure(vc, s);
            e         -= (double)vrna_eval_covar_structure(vc, s);
            prob      = exp((energy - e) / kT);
            e_string  = vrna_strdup_printf(" %6g %.2f", prob, -1 * (kT * log(prob) - energy));
          }

#endif
          print_structure(stdout, s, e_string);
          free(s);
          free(e_string);
        }
      }

      if (opt.md.compute_bpp) {
        char *msg = NULL;
        if (istty_in)
          msg = vrna_strdup_printf("\n free energy of ensemble = %6.2f kcal/mol",
                                   DBL_ROUND(energy, 2));
        else
          msg = vrna_strdup_printf(" [%6.2f]",
                                   DBL_ROUND(energy, 2));

        print_structure(stdout, structure, msg);
        free(msg);

        FILE    *aliout;
        cpair   *cp;
        char    *cent;
        double  dist;
        plist   *pl, *mfel;

        pl    = vrna_plist_from_probs(vc, opt.bppmThreshold);
        mfel  = vrna_plist(mfe_struc, 0.95 * 0.95);

        if (!opt.md.circ) {
          float *ens;
          cent    = vrna_centroid(vc, &dist);
          ens     = (float *)vrna_alloc(2 * sizeof(float));
          ens[0]  = vrna_eval_structure(vc, cent);
          ens[1]  = vrna_eval_covar_structure(vc, cent);

          char *energy_string = vrna_strdup_printf(" {%6.2f = %6.2f + %6.2f d=%.2f}",
                                                   DBL_ROUND(ens[0] - ens[1], 2),
                                                   DBL_ROUND(ens[0], 2),
                                                   DBL_ROUND((-1) * ens[1], 2),
                                                   dist);
          print_structure(stdout, cent, energy_string);
          free(energy_string);
          free(cent);
          free(ens);
        }

        if (opt.MEA) {
          float mea, *ens;
          plist *pl2;
          pl2     = vrna_plist_from_probs(vc, 1e-4 / (1 + opt.MEAgamma));
          mea     = MEA(pl2, structure, opt.MEAgamma);
          ens     = (float *)vrna_alloc(2 * sizeof(float));
          ens[0]  = vrna_eval_structure(vc, structure);
          ens[1]  = vrna_eval_covar_structure(vc, structure);

          char *energy_string = vrna_strdup_printf(" {%6.2f = %6.2f + %6.2f MEA=%.2f}",
                                                   DBL_ROUND(ens[0] - ens[1], 2),
                                                   DBL_ROUND(ens[0], 2),
                                                   DBL_ROUND((-1) * ens[1], 2),
                                                   mea);
          print_structure(stdout, structure, energy_string);
          free(energy_string);
          free(ens);
          free(pl2);
        }

        aliout = fopen(filename_out, "w");
        if (!aliout)
          vrna_message_warning("can't open %s    skipping output", filename_out);
        else
          print_aliout(vc, pl, opt.bppmThreshold, mfe_struc, aliout);

        fclose(aliout);
        cp = vrna_annotate_covar_pairs((const char **)AS, pl, mfel, opt.bppmThreshold, &opt.md);
        (void)PS_color_dot_plot(string, cp, filename_dot);
        free(cp);
        free(pl);
        free(mfel);
      } else {
        char *msg = vrna_strdup_printf(" free energy of ensemble = %6.2f kcal/mol",
                                       DBL_ROUND(energy, 2));
        print_structure(stdout, NULL, msg);
        free(msg);
      }

      {
        char *msg = NULL;
        if (opt.md.compute_bpp) {
          msg = vrna_strdup_printf(" frequency of mfe structure in ensemble %g"
                                   "; ensemble diversity %-6.2f",
                                   exp((energy - min_en) / kT),
                                   vrna_mean_bp_distance(vc));
        } else {
          msg = vrna_strdup_printf(" frequency of mfe structure in ensemble %g;",
                                   exp((energy - min_en) / kT));
        }

        print_structure(stdout, NULL, msg);
        free(msg);
      }


      free(mfe_struc);
    } /* end partition function block */

    (void)fflush(stdout);

    free(string);
    free(structure);
    free(filename_plot);
    free(filename_dot);
    free(filename_aln);
    free(filename_out);
    vrna_fold_compound_free(vc);

    vrna_aln_free(AS);
    vrna_aln_free(MSA_orig);
    vrna_aln_free(names);

    free(tmp_id);
    free(tmp_structure);

    free(MSA_ID);

    /* break after first record if fold_constrained and not explicitly instructed otherwise */
    if (opt.shape || (fold_constrained && (!(opt.constraint_batch || opt.constraint_SScons))))
      break;
  } /* end of input */

  /* check whether we've actually processed any alignment so far */
  if (first_alignment_number == get_current_id(opt.id_control)) {
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

  if (cstruc != NULL)
    free(cstruc);

  (void)fflush(stdout);
  if (opt.shape_files)
    free(opt.shape_files);

  free(filename_in);

  for (i = 0; i < input_file_num; i++)
    free(input_files[i]);
  free(input_files);

  free(opt.filename_delim);

  free_id_data(opt.id_control);

  return EXIT_SUCCESS;
}


void
process_input(FILE            *input_stream,
              const char      *input_filename,
              struct options  *opt)
{
  int           istty_in  = isatty(fileno(input_stream));
  int           istty_out = isatty(fileno(stdout));

}


static void
apply_constraints(vrna_fold_compound_t  *fc,
                  const char            *constraints_file,
                  const char            *cstruc,
                  const char            *SS_cons,
                  struct options        *opt)
{
  unsigned int constraint_options;
  
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
             char                 *mfe,
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
