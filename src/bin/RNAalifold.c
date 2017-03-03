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
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/alifold.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/MEA.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/constraints_SHAPE.h"
#include "RNAalifold_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helper.h"

#include "ViennaRNA/color_output.inc"

PRIVATE char **annote(const char  *structure,
                      const char  *AS[]);


PRIVATE void  print_pi(const vrna_pinfo_t pi,
                       FILE               *file);


PRIVATE void  print_aliout(vrna_fold_compound_t *vc,
                           plist                *pl,
                           double               threshold,
                           char                 *mfe,
                           FILE                 *aliout);


PRIVATE void  mark_endgaps(char *seq,
                           char egap);


PRIVATE cpair *make_color_pinfo(char    **sequences,
                                plist   *pl,
                                double  threshold,
                                int     n_seq,
                                plist   *mfel);


PRIVATE void  dot_bracketify(char *string);


/*--------------------------------------------------------------------------*/
int
main(int  argc,
     char *argv[])
{
  struct RNAalifold_args_info args_info;
  FILE                        *clust_file;
  unsigned int                input_type, options, longest_string,
                              input_format_options;
  char                        *input_string, *string, *structure, *cstruc, **AS, **names,
                              *constraints_file, **shape_files, *shape_method, *filename_plot,
                              *filename_dot, *filename_aln, *filename_out, *filename_in,
                              *tmp_id, *tmp_structure, *tmp_string, **input_files, *id_prefix,
                              *aln_prefix, *id_delim, *filename_delim;
  int                         s, n_seq, i, length, noPS, with_shapes, verbose, with_sci, endgaps,
                              aln_columns, mis, circular, doAlnPS, doColor, doMEA, n_back, istty_out,
                              istty_in, eval_energy, pf, istty, *shape_file_association, quiet,
                              tmp_number, batch, continuous_names, id_digits, auto_id, aln_out,
                              input_file_num, consensus_constraint, enforceConstraints;
  long int                    alignment_number;
  double                      min_en, real_en, cov_en, MEAgamma, bppmThreshold;
  vrna_md_t                   md;
  vrna_fold_compound_t        *vc;

  string                  = structure = cstruc = NULL;
  endgaps                 = mis = pf = circular = doAlnPS = doColor = n_back = eval_energy = oldAliEn = doMEA = ribo = noPS = 0;
  aln_columns             = 60;
  aln_out                 = 0;
  aln_prefix              = NULL;
  do_backtrack            = 1;
  bppmThreshold           = 1e-6;
  MEAgamma                = 1.0;
  shape_files             = NULL;
  shape_file_association  = 0;
  shape_method            = NULL;
  with_shapes             = 0;
  verbose                 = 0;
  quiet                   = 0;
  with_sci                = 0;
  filename_plot           = NULL;
  filename_dot            = NULL;
  filename_aln            = NULL;
  filename_out            = NULL;
  filename_in             = NULL;
  input_files             = NULL;
  tmp_id                  = NULL;
  tmp_structure           = NULL;
  input_format_options    = VRNA_FILE_FORMAT_MSA_CLUSTAL; /* default to ClustalW format */
  vc                      = NULL;
  clust_file              = stdin;
  continuous_names        = 0;
  auto_id                 = 0;
  input_file_num          = 0;
  consensus_constraint    = 0;
  set_model_details(&md);

  /*
   #############################################
   # check the command line prameters
   #############################################
   */
  if (RNAalifold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* get basic set of model details */
  ggo_get_md_eval(args_info, md);
  ggo_get_md_fold(args_info, md);
  ggo_get_md_part(args_info, md);
  ggo_get_circ(args_info, md.circ);

  /* check dangle model */
  if (!((md.dangles == 0) || (md.dangles == 2))) {
    vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    md.dangles = dangles = 2;
  }

  ggo_get_ID_manipulation(args_info,
                          auto_id,
                          id_prefix, "alignment",
                          id_delim, "_",
                          id_digits, 4,
                          alignment_number, 1);

  /* do not treat first alignment special */
  if (args_info.continuous_ids_given || auto_id)
    continuous_names = 1;

  ggo_get_constraints_settings(args_info,
                               fold_constrained,
                               constraints_file,
                               enforceConstraints,
                               batch);

  if (args_info.SS_cons_given) {
    fold_constrained      = 1;
    consensus_constraint  = 1;
  }

  /* do not produce postscript output */
  if (args_info.noPS_given)
    noPS = 1;

  /* partition function settings */
  if (args_info.partfunc_given) {
    pf = 1;
    if (args_info.partfunc_arg != -1)
      md.compute_bpp = do_backtrack = args_info.partfunc_arg;
  }

  /* MEA (maximum expected accuracy) settings */
  if (args_info.MEA_given) {
    pf = doMEA = 1;
    if (args_info.MEA_arg != -1)
      MEAgamma = args_info.MEA_arg;
  }

  /* set the bppm threshold for the dotplot */
  if (args_info.bppmThreshold_given)
    bppmThreshold = MIN2(1., MAX2(0., args_info.bppmThreshold_arg));

  /* set cfactor */
  if (args_info.cfactor_given)
    md.cv_fact = cv_fact = args_info.cfactor_arg;

  /* set nfactor */
  if (args_info.nfactor_given)
    md.nc_fact = nc_fact = args_info.nfactor_arg;

  if (args_info.endgaps_given)
    endgaps = 1;

  if (args_info.mis_given)
    mis = 1;

  if (args_info.color_given)
    doColor = 1;

  if (args_info.aln_given)
    doAlnPS = 1;

  if (args_info.aln_EPS_cols_given)
    aln_columns = args_info.aln_EPS_cols_arg;

  if (args_info.aln_stk_given) {
    aln_out = 1;
    if (args_info.aln_stk_arg)
      aln_prefix = strdup(args_info.aln_stk_arg);
  }

  if (args_info.old_given)
    md.oldAliEn = oldAliEn = 1;

  if (args_info.stochBT_given) {
    n_back          = args_info.stochBT_arg;
    md.uniq_ML      = 1;
    md.compute_bpp  = do_backtrack = 0;
    pf              = 1;
    vrna_init_rand();
  }

  if (args_info.stochBT_en_given) {
    n_back          = args_info.stochBT_en_arg;
    md.uniq_ML      = 1;
    md.compute_bpp  = do_backtrack = 0;
    pf              = 1;
    eval_energy     = 1;
    vrna_init_rand();
  }

  if (args_info.ribosum_file_given) {
    RibosumFile = strdup(args_info.ribosum_file_arg);
    md.ribo     = ribo = 1;
  }

  if (args_info.ribosum_scoring_given) {
    RibosumFile = NULL;
    md.ribo     = ribo = 1;
  }

  if (args_info.layout_type_given)
    rna_plot_type = args_info.layout_type_arg;

  if (args_info.verbose_given)
    verbose = 1;

  if (args_info.quiet_given) {
    if (verbose)
      vrna_message_warning("Can not be verbose and quiet at the same time! I keep on being chatty...");
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

  if (args_info.shapeMethod_given)
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

  /* sci computation */
  if (args_info.sci_given)
    with_sci = 1;

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
  else
    filename_delim = strdup(id_delim);

  if (isspace(*filename_delim)) {
    free(filename_delim);
    filename_delim = NULL;
  }

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

  if (circular && gquad)
    vrna_message_error("G-Quadruplex support is currently not available for circular RNA structures");

  make_pair_matrix(); /* for make_color_pinfo */

  if (circular && noLonelyPairs)
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
      vrna_message_error("Your input file is missing sequences! Either your file is empty, or not in %s format!",
                         format);
      free(format);
    }

    input_format_options = format_guess;
  }

  if (fold_constrained && (!constraints_file) && (!consensus_constraint)) {
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

  long int first_alignment_number = alignment_number;

  if (aln_out) {
    if (!aln_prefix) {
      aln_prefix = strdup("RNAalifold_results.stk");
    } else {
      char *tmp = vrna_strdup_printf("%s.stk", aln_prefix);
      free(aln_prefix);
      aln_prefix = vrna_filename_sanitize(tmp, filename_delim);
      free(tmp);
    }
  }

  while (!feof(clust_file)) {
    char *MSA_ID = NULL;
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

    if (quiet)
      input_format_options |= VRNA_FILE_FORMAT_MSA_QUIET;

    /* read the first record from input file */
    n_seq = vrna_file_msa_read_record(clust_file, &names, &AS, &tmp_id, &tmp_structure, input_format_options);
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
    if (tmp_id && (!auto_id))                                       /* we've read an ID from file, so we use it */
      MSA_ID = strdup(tmp_id);
    else if (auto_id || (alignment_number > 1) || continuous_names) /* we have nuffin', Jon Snow (...so we simply generate an ID) */
      MSA_ID = vrna_strdup_printf("%s%s%0*ld", id_prefix, id_delim, id_digits, alignment_number);

    /* construct output file names */
    if (MSA_ID) {
      /* construct file names */
      filename_plot = vrna_strdup_printf("%s%sss.ps", MSA_ID, id_delim);
      filename_dot  = vrna_strdup_printf("%s%sdp.ps", MSA_ID, id_delim);
      filename_aln  = vrna_strdup_printf("%s%saln.ps", MSA_ID, id_delim);
      filename_out  = vrna_strdup_printf("%s%sali.out", MSA_ID, id_delim);

      /* sanitize file names */
      tmp_string = vrna_filename_sanitize(filename_plot, filename_delim);
      free(filename_plot);
      filename_plot = tmp_string;
      tmp_string    = vrna_filename_sanitize(filename_dot, filename_delim);
      free(filename_dot);
      filename_dot  = tmp_string;
      tmp_string    = vrna_filename_sanitize(filename_aln, filename_delim);
      free(filename_aln);
      filename_aln  = tmp_string;
      tmp_string    = vrna_filename_sanitize(filename_out, filename_delim);
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

    length    = (int)strlen(AS[0]);
    structure = (char *)vrna_alloc((unsigned)length + 1);

    if (endgaps)
      for (i = 0; i < n_seq; i++)
        mark_endgaps(AS[i], '~');

    /*
     ########################################################
     # begin actual calculations
     ########################################################
     */
    options = VRNA_OPTION_MFE;

    if (pf)
      options |= VRNA_OPTION_PF;

    vc = vrna_fold_compound_comparative((const char **)AS, &md, VRNA_OPTION_DEFAULT);

    if (fold_constrained) {
      if (cstruc) {
        strncpy(structure, cstruc, length);
        if (enforceConstraints)
          vrna_constraints_add(vc, (const char *)structure, VRNA_CONSTRAINT_DB_DEFAULT | VRNA_CONSTRAINT_DB_ENFORCE_BP);
        else
          vrna_constraints_add(vc, (const char *)structure, VRNA_CONSTRAINT_DB_DEFAULT);
      } else if (constraints_file) {
        vrna_constraints_add(vc, constraints_file, 0);
      } else if (tmp_structure) {
        dot_bracketify(tmp_structure);
        if (enforceConstraints)
          vrna_constraints_add(vc, (const char *)tmp_structure, VRNA_CONSTRAINT_DB_DEFAULT | VRNA_CONSTRAINT_DB_ENFORCE_BP);
        else
          vrna_constraints_add(vc, (const char *)tmp_structure, VRNA_CONSTRAINT_DB_DEFAULT);
      } else {
        vrna_message_warning("Constraint missing");
      }
    }

    if (with_shapes) {
      for (s = 0; shape_file_association[s] != -1; s++);

      if ((s != n_seq) && (!quiet))
        vrna_message_warning("Number of sequences in alignment does not match number of provided SHAPE reactivity data files! ");

      shape_files             = (char **)vrna_realloc(shape_files, (n_seq + 1) * sizeof(char *));
      shape_file_association  = (int *)vrna_realloc(shape_file_association, (n_seq + 1) * sizeof(int));
    }

    if (with_shapes) {
      vrna_constraints_add_SHAPE_ali(vc, \
                                     shape_method, \
                                     (const char **)shape_files, \
                                     shape_file_association, \
                                     verbose, \
                                     VRNA_OPTION_MFE | ((pf) ? VRNA_OPTION_PF : 0));
    }

    min_en  = vrna_mfe(vc, structure);
    real_en = vrna_eval_structure(vc, structure);
    cov_en  = vrna_eval_covar_structure(vc, structure);

    string = (mis) ? consens_mis((const char **)AS) : consensus((const char **)AS);

    float sci     = min_en;
    float e_mean  = 0;

    if (with_sci) {
      for (i = 0; AS[i] != NULL; i++) {
        char *seq = get_ungapped_sequence(AS[i]);
        e_mean += vrna_fold(seq, NULL);
        free(seq);
      }
      e_mean /= i;
      if (e_mean == 0.)
        sci = 0.;
      else
        sci /= e_mean;
    }

    print_fasta_header(stdout, MSA_ID);
    fprintf(stdout, "%s\n", string);
    char *energy_string = NULL;
    char *e_individual  = NULL;

    if (with_shapes)
      e_individual = vrna_strdup_printf("%6.2f + %6.2f + %6.2f", real_en, -cov_en, min_en - real_en + cov_en);
    else
      e_individual = vrna_strdup_printf("%6.2f + %6.2f", real_en, min_en - real_en);

    if (istty_in) {
      if (with_sci)
        energy_string = vrna_strdup_printf("\n minimum free energy = %6.2f kcal/mol (%s)\n SCI = %2.4f",
                                           min_en, e_individual, sci);
      else
        energy_string = vrna_strdup_printf("\n minimum free energy = %6.2f kcal/mol (%s)",
                                           min_en, e_individual);
    } else {
      if (with_sci)
        energy_string = vrna_strdup_printf(" (%6.2f = %s) [sci = %2.4f]",
                                           min_en, e_individual, sci);
      else
        energy_string = vrna_strdup_printf(" (%6.2f = %s)",
                                           min_en, e_individual);
    }

    print_structure(stdout, structure, energy_string);

    free(energy_string);
    free(e_individual);

    if (!noPS) {
      char **A;
      A = annote(structure, (const char **)AS);

      if (doColor)
        (void)vrna_file_PS_rnaplot_a(string, structure, filename_plot, A[0], A[1], &md);
      else
        (void)vrna_file_PS_rnaplot_a(string, structure, filename_plot, NULL, A[1], &md);

      free(A[0]);
      free(A[1]);
      free(A);
    }

    if (doAlnPS)
      vrna_file_PS_aln(filename_aln, (const char **)AS, (const char **)names, structure, aln_columns);

    if (aln_out) {
      unsigned int options = VRNA_FILE_FORMAT_MSA_STOCKHOLM
                             | VRNA_FILE_FORMAT_MSA_APPEND;
      if (mis)
        options |= VRNA_FILE_FORMAT_MSA_MIS;

      vrna_file_msa_write((const char *)aln_prefix,
                          (const char **)names,
                          (const char **)AS,
                          MSA_ID,
                          (const char *)structure,
                          "RNAalifold prediction",
                          options);
    }

    /* free mfe arrays */
    vrna_mx_mfe_free(vc);

    if (pf) {
      float energy, kT;
      char  *mfe_struc;

      mfe_struc = strdup(structure);

      vrna_exp_params_rescale(vc, &min_en);
      pf_scale = vc->exp_params->pf_scale;


      kT = vc->exp_params->kT / 1000.;

      if ((length > 2000) && (!quiet))
        vrna_message_info(stderr, "scaling factor %f\n", vc->exp_params->pf_scale);

      fflush(stdout);

      if (cstruc != NULL)
        strncpy(structure, cstruc, length + 1);

      energy = vrna_pf(vc, structure);


      if (n_back > 0) {
        /*stochastic sampling*/
        for (i = 0; i < n_back; i++) {
          char    *s, *e_string = NULL;

#if CHECK_PROBABILITIES
          double  prob2, prob = 1.;
          s = alipbacktrack(&prob);
          double  e = (double)vrna_eval_structure(vc, s);
          e     -= (double)vrna_eval_covar_structure(vc, s);
          prob2 = exp((energy - e) / kT);
          if (eval_energy)
            e_string = vrna_strdup_printf(" %6g (%6g) %.2f (%.2f)", prob, prob2, -1 * (kT * log(prob) - energy), e);

#else
          double prob = 1.;
          s = vrna_pbacktrack(vc);

          if (eval_energy) {
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

      if (do_backtrack) {
        char *msg = NULL;
        if (istty_in)
          msg = vrna_strdup_printf("\n free energy of ensemble = %6.2f kcal/mol", energy);
        else
          msg = vrna_strdup_printf(" [%6.2f]", energy);

        print_structure(stdout, structure, msg);
        free(msg);
      } else {
        char *msg = vrna_strdup_printf(" free energy of ensemble = %6.2f kcal/mol", energy);
        print_structure(stdout, NULL, msg);
        free(msg);
      }

      if (do_backtrack) {
        FILE    *aliout;
        cpair   *cp;
        char    *cent;
        double  dist;
        plist   *pl, *mfel;

        pl    = vrna_plist_from_probs(vc, bppmThreshold);
        mfel  = vrna_plist(mfe_struc, 0.95 * 0.95);

        if (!circular) {
          float *ens;
          cent    = vrna_centroid(vc, &dist);
          ens     = (float *)vrna_alloc(2 * sizeof(float));
          ens[0]  = vrna_eval_structure(vc, cent);
          ens[1]  = vrna_eval_covar_structure(vc, cent);

          char *energy_string = vrna_strdup_printf(" {%6.2f = %6.2f + %6.2f d=%.2f}",
                                                   ens[0] - ens[1], ens[0], (-1) * ens[1], dist);
          print_structure(stdout, cent, energy_string);
          free(energy_string);
          free(cent);
          free(ens);
        }

        if (doMEA) {
          float mea, *ens;
          plist *pl2;
          pl2     = vrna_plist_from_probs(vc, 1e-4 / (1 + MEAgamma));
          mea     = MEA(pl2, structure, MEAgamma);
          ens     = (float *)vrna_alloc(2 * sizeof(float));
          ens[0]  = vrna_eval_structure(vc, structure);
          ens[1]  = vrna_eval_covar_structure(vc, structure);

          char *energy_string = vrna_strdup_printf(" {%6.2f = %6.2f + %6.2f MEA=%.2f}",
                                                   ens[0] - ens[1], ens[0], (-1) * ens[1], mea);
          print_structure(stdout, structure, energy_string);
          free(energy_string);
          free(ens);
          free(pl2);
        }

        aliout = fopen(filename_out, "w");
        if (!aliout)
          vrna_message_warning("can't open %s    skipping output", filename_out);
        else
          print_aliout(vc, pl, bppmThreshold, mfe_struc, aliout);

        fclose(aliout);
        cp = make_color_pinfo(AS, pl, bppmThreshold, n_seq, mfel);
        (void)PS_color_dot_plot(string, cp, filename_dot);
        free(cp);
        free(pl);
        free(mfel);
      }

      {
        char *msg = NULL;
        if (do_backtrack) {
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

    ID_number_increase(alignment_number, "Alignment");

    (void)fflush(stdout);

    free(string);
    free(structure);
    free(filename_plot);
    free(filename_dot);
    free(filename_aln);
    free(filename_out);
    vrna_fold_compound_free(vc);

    for (i = 0; AS[i]; i++) {
      free(AS[i]);
      free(names[i]);
    }
    free(AS);
    free(names);

    free(tmp_id);
    free(tmp_structure);

    free(MSA_ID);

    /* break after first record if fold_constrained and not explicitly instructed otherwise */
    if (with_shapes || (fold_constrained && (!(batch || consensus_constraint))))
      break;
  } /* end of input */

  if (first_alignment_number == alignment_number) {
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
    vrna_message_error("Your input file is missing sequences! Either your file is empty, or not in %s format!",
                       format);
    free(format);
  }

  if (clust_file != stdin)
    fclose(clust_file);

  if (cstruc != NULL)
    free(cstruc);

  (void)fflush(stdout);
  if (shape_files)
    free(shape_files);

  free(filename_in);

  for (i = 0; i < input_file_num; i++)
    free(input_files[i]);
  free(input_files);

  free(id_prefix);
  free(id_delim);
  free(filename_delim);

  return EXIT_SUCCESS;
}


PRIVATE void
dot_bracketify(char *string)
{
  int i;

  if (string) {
    for (i = 0; string[i] != '\0'; i++) {
      switch (string[i]) {
        case '<':
        case '(':
          string[i] = '(';
          break;
        case '>':
        case ')':
          string[i] = ')';
          break;
        default:
          string[i] = '.';
          break;
      }
    }
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

PRIVATE char **
annote(const char *structure,
       const char *AS[])
{
  /* produce annotation for colored drawings from vrna_file_PS_rnaplot_a() */
  char  *ps, *colorps, **A;
  int   i, n, s, pairings, maxl;
  short *ptable;
  char  *colorMatrix[6][3] = {
    { "0.0 1",  "0.0 0.6",  "0.0 0.2"  }, /* red    */
    { "0.16 1", "0.16 0.6", "0.16 0.2" }, /* ochre  */
    { "0.32 1", "0.32 0.6", "0.32 0.2" }, /* turquoise */
    { "0.48 1", "0.48 0.6", "0.48 0.2" }, /* green  */
    { "0.65 1", "0.65 0.6", "0.65 0.2" }, /* blue   */
    { "0.81 1", "0.81 0.6", "0.81 0.2" }  /* violet */
  };

  n     = strlen(AS[0]);
  maxl  = 1024;

  A       = (char **)vrna_alloc(sizeof(char *) * 2);
  ps      = (char *)vrna_alloc(maxl);
  colorps = (char *)vrna_alloc(maxl);
  ptable  = vrna_ptable(structure);
  for (i = 1; i <= n; i++) {
    char  pps[64], ci = '\0', cj = '\0';
    int   j, type, pfreq[8] = {
      0, 0, 0, 0, 0, 0, 0, 0
    }, vi = 0, vj = 0;
    if ((j = ptable[i]) < i)
      continue;

    for (s = 0; AS[s] != NULL; s++) {
      type = pair[encode_char(AS[s][i - 1])][encode_char(AS[s][j - 1])];
      pfreq[type]++;
      if (type) {
        if (AS[s][i - 1] != ci) {
          ci = AS[s][i - 1];
          vi++;
        }

        if (AS[s][j - 1] != cj) {
          cj = AS[s][j - 1];
          vj++;
        }
      }
    }
    for (pairings = 0, s = 1; s <= 7; s++)
      if (pfreq[s])
        pairings++;

    if ((maxl - strlen(ps) < 192) || ((maxl - strlen(colorps)) < 64)) {
      maxl    *= 2;
      ps      = (char *)vrna_realloc(ps, sizeof(char) * maxl);
      colorps = (char *)vrna_realloc(colorps, sizeof(char) * maxl);
      if ((ps == NULL) || (colorps == NULL))
        vrna_message_error("out of memory in realloc");
    }

    if (pfreq[0] <= 2 && pairings > 0) {
      snprintf(pps, 64, "%d %d %s colorpair\n",
               i, j, colorMatrix[pairings - 1][pfreq[0]]);
      strcat(colorps, pps);
    }

    if (pfreq[0] > 0) {
      snprintf(pps, 64, "%d %d %d gmark\n", i, j, pfreq[0]);
      strcat(ps, pps);
    }

    if (vi > 1) {
      snprintf(pps, 64, "%d cmark\n", i);
      strcat(ps, pps);
    }

    if (vj > 1) {
      snprintf(pps, 64, "%d cmark\n", j);
      strcat(ps, pps);
    }
  }
  free(ptable);
  A[0]  = colorps;
  A[1]  = ps;
  return A;
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


PRIVATE cpair *
make_color_pinfo(char   **sequences,
                 plist  *pl,
                 double threshold,
                 int    n_seq,
                 plist  *mfel)
{
  /* produce info for PS_color_dot_plot */
  cpair *cp;
  int   i, n, s, a, b, z, t, j, c;
  int   pfreq[7];

  for (n = 0; pl[n].i > 0; n++);
  c   = 0;
  cp  = (cpair *)vrna_alloc(sizeof(cpair) * (n + 1));
  for (i = 0; i < n; i++) {
    int ncomp = 0;
    if (pl[i].p > threshold) {
      cp[c].i = pl[i].i;
      cp[c].j = pl[i].j;
      cp[c].p = pl[i].p;
      for (z = 0; z < 7; z++)
        pfreq[z] = 0;
      for (s = 0; s < n_seq; s++) {
        a = encode_char(toupper(sequences[s][cp[c].i - 1]));
        b = encode_char(toupper(sequences[s][cp[c].j - 1]));
        if ((sequences[s][cp[c].j - 1] == '~') || (sequences[s][cp[c].i - 1] == '~'))
          continue;

        pfreq[pair[a][b]]++;
      }
      for (z = 1; z < 7; z++)
        if (pfreq[z] > 0)
          ncomp++;

      cp[c].hue = (ncomp - 1.0) / 6.2;   /* hue<6/6.9 (hue=1 ==  hue=0) */
      cp[c].sat = 1 - MIN2(1.0, (float)(pfreq[0] * 2. /*pi[i].bp[0]*/ / (n_seq)));
      c++;
    }
  }
  for (t = 0; mfel[t].i > 0; t++) {
    int nofound = 1;
    for (j = 0; j < c; j++) {
      if ((cp[j].i == mfel[t].i) && (cp[j].j == mfel[t].j)) {
        cp[j].mfe = 1;
        nofound   = 0;
        break;
      }
    }
    if (nofound) {
      vrna_message_warning("mfe base pair with very low prob in pf: %d %d",
                           mfel[t].i,
                           mfel[t].j);

      cp        = (cpair *)vrna_realloc(cp, sizeof(cpair) * (c + 2));
      cp[c].i   = mfel[t].i;
      cp[c].j   = mfel[t].j;
      cp[c].p   = 0.;
      cp[c].hue = 0;
      cp[c].sat = 0;
      cp[c].mfe = 1;
      c++;
      cp[c].i = cp[c].j = 0;
    }
  }
  return cp;
}
