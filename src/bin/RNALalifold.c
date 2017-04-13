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
#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/plot_aln.h"
#include "ViennaRNA/plot_structure.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/alifold.h"
#include "ViennaRNA/Lfold.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/file_formats_msa.h"
#include "ViennaRNA/file_utils.h"
#include "RNALalifold_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helper.h"

#include "ViennaRNA/color_output.inc"

#define MAX_NUM_NAMES    500

typedef struct {
  char      **names;
  char      **strings;
  char      *prefix;
  int       columns;
  vrna_md_t *md;
  int       ss_eps;
  int       msa_eps;
  int       msa_stk;
  int       csv;
  int       mis;
} hit_data;


PRIVATE char **annote(const char  *structure,
                      const char  *AS[],
                      vrna_md_t   *md);


PRIVATE char **get_subalignment(const char  *AS[],
                                int         i,
                                int         j);


PRIVATE void  delete_alignment(char *AS[]);


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
  char                          *string, *structure, *ParamFile, *ns_bases, *prefix,
                                **AS, **names, *filename_in, *id_prefix, *id_delim, *filename_delim,
                                **input_files, *tmp_id, *tmp_structure;
  unsigned int                  input_format_options;
  int                           n_seq, i, maxdist, unchangednc, unchangedcv, quiet, auto_id, gquad,
                                mis, istty, alnPS, aln_columns, aln_out, ssPS, input_file_num, id_digits;
  long int                      alignment_number, first_alignment_number;
  vrna_md_t                     md;
  vrna_fold_compound_t          *fc;

  clust_file            = stdin;
  string                = structure = ParamFile = ns_bases = prefix = NULL;
  mis                   = 0;
  maxdist               = 70;
  do_backtrack          = unchangednc = unchangedcv = 1;
  dangles               = 2;
  ribo                  = 0;
  alnPS                 = 0;
  ssPS                  = 0;
  aln_out               = 0;
  aln_columns           = 60;
  input_format_options  = VRNA_FILE_FORMAT_MSA_CLUSTAL;   /* default to ClustalW format */
  filename_in           = NULL;
  input_files           = NULL;
  input_file_num        = 0;
  quiet                 = 0;
  auto_id               = 0;
  gquad                 = 0;

  vrna_md_set_default(&md);

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNALalifold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  ggo_get_ID_manipulation(args_info,
                          auto_id,
                          id_prefix, "alignment",
                          id_delim, "_",
                          id_digits, 4,
                          alignment_number, 1);

  /* temperature */
  if (args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;

  /* structure constraint */
  if (args_info.noTetra_given)
    md.special_hp = tetra_loop = 0;

  /* set dangle model */
  if (args_info.dangles_given) {
    if ((args_info.dangles_arg < 0) || (args_info.dangles_arg > 3))
      vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    else
      md.dangles = dangles = args_info.dangles_arg;
  }

  /* do not allow weak pairs */
  if (args_info.noLP_given)
    md.noLP = noLonelyPairs = 1;

  /* do not allow wobble pairs (GU) */
  if (args_info.noGU_given)
    md.noGU = noGU = 1;

  /* do not allow weak closing pairs (AU,GU) */
  if (args_info.noClosingGU_given)
    md.noGUclosure = no_closingGU = 1;

  /* set energy model */
  if (args_info.energyModel_given)
    md.energy_set = energy_set = args_info.energyModel_arg;

  /* take another energy parameter set */
  if (args_info.paramFile_given)
    ParamFile = strdup(args_info.paramFile_arg);

  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if (args_info.nsp_given)
    ns_bases = strdup(args_info.nsp_arg);

  /* set cfactor */
  if (args_info.cfactor_given) {
    md.cv_fact  = cv_fact = args_info.cfactor_arg;
    unchangedcv = 0;
  }

  /* set nfactor */
  if (args_info.nfactor_given) {
    md.nc_fact  = nc_fact = args_info.nfactor_arg;
    unchangednc = 0;
  }

  /* set the maximum base pair span */
  if (args_info.span_given)
    md.max_bp_span = maxdist = args_info.span_arg;

  /* calculate most informative sequence */
  if (args_info.mis_given)
    mis = 1;

  if (args_info.csv_given)
    csv = 1;

  if (args_info.ribosum_file_given) {
    md.ribo     = ribo = 1;
    RibosumFile = strdup(args_info.ribosum_file_arg);
  }

  if (args_info.ribosum_scoring_given) {
    RibosumFile = NULL;
    md.ribo     = ribo = 1;
  }

  /* gquadruplex support */
  if (args_info.gquad_given)
    md.gquad = gquad = 1;

  if (args_info.aln_given) {
    alnPS = aln_out = ssPS = 1;
    if (args_info.aln_arg)
      prefix = strdup(args_info.aln_arg);
  }

  if (args_info.aln_stk_given) {
    aln_out = 1;
    if (args_info.aln_stk_arg) {
      if (prefix) {
        vrna_message_info(stdout, "multiple output prefixes detected, using \"%s\"", args_info.aln_stk_arg);
        free(prefix);
      }

      prefix = strdup(args_info.aln_stk_arg);
    }
  }

  if (args_info.aln_EPS_ss_given)
    ssPS = 1;

  if (args_info.aln_EPS_given)
    alnPS = 1;

  if (args_info.aln_EPS_cols_given)
    aln_columns = args_info.aln_EPS_cols_arg;

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
  else
    filename_delim = strdup(id_delim);

  if (isspace(*filename_delim)) {
    free(filename_delim);
    filename_delim = NULL;
  }

  if (args_info.quiet_given)
    quiet = 1;

  /* free allocated memory of command line data structure */
  RNALalifold_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */
  if ((ribo == 1) && (unchangednc))
    md.nc_fact = nc_fact = 0.5;

  if ((ribo == 1) && (unchangedcv))
    md.cv_fact = cv_fact = 0.6;

  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (ns_bases != NULL)
    vrna_md_set_nonstandards(&md, ns_bases);

  istty = isatty(fileno(stdout)) && isatty(fileno(stdin));

  if (istty && (clust_file == stdin))
    vrna_message_input_seq("Input aligned sequences in clustalw format");

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
      vrna_message_error("Your input file is missing sequences! Either your file is empty, or not in %s format!",
                         format);
      free(format);
    }

    input_format_options = format_guess;
  }

  first_alignment_number = alignment_number;

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
    if (tmp_id && (!auto_id)) /* we've read an ID from file, so we use it */
      MSA_ID = strdup(tmp_id);
    else if (auto_id)         /* we have nuffin', Jon Snow (...so we simply generate an ID) */
      MSA_ID = vrna_strdup_printf("%s%s%0*ld", id_prefix, id_delim, id_digits, alignment_number);

    print_fasta_header(stdout, MSA_ID);

    fc = vrna_fold_compound_comparative((const char **)AS, &md, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

    hit_data data;

    data.names    = names;
    data.strings  = AS;
    data.prefix   = MSA_ID;
    data.columns  = aln_columns;
    data.md       = &md;
    data.ss_eps   = ssPS;
    data.msa_eps  = alnPS;
    data.msa_stk  = aln_out;
    data.csv      = csv;
    data.mis      = mis;

    (void)vrna_mfe_window_cb(fc, &print_hit_cb, (void *)&data);

    string = (mis) ? consens_mis((const char **)AS) : consensus((const char **)AS);
    printf("%s\n", string);

    (void)fflush(stdout);
    free(string);
    free(structure);

    for (i = 0; AS[i]; i++) {
      free(AS[i]);
      free(names[i]);
    }
    free(AS);
    free(names);

    vrna_fold_compound_free(fc);
    free(MSA_ID);
    free(tmp_id);
    free(tmp_structure);

    ID_number_increase(alignment_number, "Alignment");
  }

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
  char      *prefix;
  char      *msg;
  int       columns;
  vrna_md_t *md;
  int       with_ss, with_msa, with_stk, with_csv, with_mis;

  names     = ((hit_data *)data)->names;
  strings   = ((hit_data *)data)->strings;
  prefix    = ((hit_data *)data)->prefix;
  columns   = ((hit_data *)data)->columns;
  md        = ((hit_data *)data)->md;
  with_ss   = ((hit_data *)data)->ss_eps;
  with_msa  = ((hit_data *)data)->msa_eps;
  with_stk  = ((hit_data *)data)->msa_stk;
  with_csv  = ((hit_data *)data)->csv;
  with_mis  = ((hit_data *)data)->mis;

  sub   = get_subalignment((const char **)strings, start, end);
  cons  = (with_mis) ? consens_mis((const char **)sub) : consensus((const char **)sub);
  A     = annote(structure, (const char **)sub, md);

  if (with_csv == 1)
    msg = vrna_strdup_printf(",%6.2f,%d,%d", en, start, end);
  else
    msg = vrna_strdup_printf(" (%6.2f) %4d - %4d", en, start, end);

  print_structure(stdout, structure, msg);
  free(msg);

  if (with_ss) {
    if (prefix)
      fname = vrna_strdup_printf("%s_ss_%d_%d.eps", prefix, start, end);
    else
      fname = vrna_strdup_printf("ss_%d_%d.eps", start, end);

    tmp_string = vrna_filename_sanitize(fname, "_");
    free(fname);
    fname = tmp_string;

    (void)vrna_file_PS_rnaplot_a(cons, structure, fname, A[0], A[1], md);
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
    free(fname);
    fname = tmp_string;

    vrna_file_PS_aln_sub(fname, (const char **)sub, (const char **)names, structure, start, -1, columns);
    free(fname);
  }

  if (with_stk) {
    unsigned int options = VRNA_FILE_FORMAT_MSA_STOCKHOLM
                           | VRNA_FILE_FORMAT_MSA_APPEND;
    if (with_mis)
      options |= VRNA_FILE_FORMAT_MSA_MIS;

    if (prefix) {
      id    = vrna_strdup_printf("%s_aln_%d_%d", prefix, start, end);
      fname = vrna_strdup_printf("%s.stk", prefix);

      tmp_string = vrna_filename_sanitize(fname, "_");
      free(fname);
      fname = tmp_string;

      vrna_file_msa_write(fname, (const char **)names, (const char **)sub, id, structure, "RNALalifold prediction", options);
      free(fname);
    } else {
      id = vrna_strdup_printf("aln_%d_%d", start, end);
      vrna_file_msa_write("RNALalifold_results.stk", (const char **)names, (const char **)sub, id, structure, "RNALalifold prediction", options);
    }

    free(id);
  }

  delete_alignment(sub);
}


PRIVATE char **
annote(const char *structure,
       const char *AS[],
       vrna_md_t  *md)
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
      type = md->pair[vrna_nucleotide_encode(AS[s][i - 1], md)][vrna_nucleotide_encode(AS[s][j - 1], md)];
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
      ps      = realloc(ps, maxl);
      colorps = realloc(colorps, maxl);
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


PRIVATE char **
get_subalignment(const char **AS,
                 int        i,
                 int        j)
{
  char  **sub;
  int   n_seq, s;

  /* get number of sequences in alignment */
  for (n_seq = 0; AS[n_seq] != NULL; n_seq++);

  sub = (char **)vrna_alloc(sizeof(char *) * (n_seq + 1));
  for (s = 0; s < n_seq; s++)
    sub[s] = vrna_alloc(sizeof(char) * (j - i + 2));
  sub[s] = NULL;

  /* copy subalignment */
  for (s = 0; s < n_seq; s++) {
    sub[s]              = memcpy(sub[s], AS[s] + i - 1, sizeof(char) * (j - i + 1));
    sub[s][(j - i + 1)] = '\0';
  }

  return sub;
}


PRIVATE void
delete_alignment(char **AS)
{
  int n_seq;

  /* get number of sequences in alignment */
  for (n_seq = 0; AS[n_seq] != NULL; n_seq++)
    free(AS[n_seq]);
  free(AS);
}
