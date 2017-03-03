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


#define MAX_NUM_NAMES    500

typedef struct {
  char      **names;
  char      **strings;
  char      *prefix;
  int        columns;
  vrna_md_t  *md;
  int        ss_eps;
  int        msa_eps;
  int        msa_stk;
  int        csv;
} hit_data;


PRIVATE char **annote(const char  *structure,
                      const char  *AS[],
                      vrna_md_t *md);


PRIVATE void  get_subalignment(const char *AS[],
                               char       *sub[],
                               int        i,
                               int        j);

PRIVATE void  delete_alignment(char *AS[]);

PRIVATE void
print_hit_cb(int start, int end, const char *structure, float en, void *data);


int
main(int  argc,
     char *argv[])
{
  FILE                          *clust_file;
  struct RNALalifold_args_info  args_info;
  char                          *string, *structure, *ParamFile, *ns_bases, *prefix,
                                *AS[MAX_NUM_NAMES], *names[MAX_NUM_NAMES];
  int                           n_seq, i, length, maxdist, unchangednc, unchangedcv,
                                mis, pf, istty, alnPS, aln_columns, aln_out, ssPS;
  double                        min_en;
  vrna_md_t                     md;

  clust_file    = stdin;
  string        = structure = ParamFile = ns_bases = prefix = NULL;
  mis           = pf = 0;
  maxdist       = 70;
  do_backtrack  = unchangednc = unchangedcv = 1;
  dangles       = 2;
  ribo          = 0;
  alnPS         = 0;
  ssPS          = 0;
  aln_out       = 0;
  aln_columns   = 60;

  vrna_md_set_default(&md);

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNALalifold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

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

  if (args_info.aln_EPS_ss_given) {
    ssPS = 1;
    if (args_info.aln_EPS_ss_arg) {
      if (prefix) {
        vrna_message_info(stdout, "multiple output prefixes detected, using \"%s\"", args_info.aln_EPS_ss_arg);
        free(prefix);
      }

      prefix = strdup(args_info.aln_EPS_ss_arg);
    }
  }

  if (args_info.aln_EPS_given) {
    alnPS = 1;
    if (args_info.aln_EPS_arg) {
      if (prefix) {
        vrna_message_info(stdout, "multiple output prefixes detected, using \"%s\"", args_info.aln_EPS_arg);
        free(prefix);
      }

      prefix = strdup(args_info.aln_EPS_arg);
    }
  }

  if (args_info.aln_EPS_cols_given)
    aln_columns = args_info.aln_EPS_cols_arg;

  /* check unnamed options a.k.a. filename of input alignment */
  if (args_info.inputs_num == 1) {
    clust_file = fopen(args_info.inputs[0], "r");
    if (clust_file == NULL)
      vrna_message_warning("can't open %s", args_info.inputs[0]);
  } else {
    RNALalifold_cmdline_parser_print_help();
    exit(1);
  }

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

  n_seq = read_clustal(clust_file, AS, names);
  if (clust_file != stdin)
    fclose(clust_file);

  if (n_seq == 0)
    vrna_message_error("no sequences found");

  length = (int)strlen(AS[0]);
  if (length < maxdist) {
    vrna_message_warning("Alignment length < window size: setting L=%d", length);
    maxdist = length;
  }

  structure = (char *)vrna_alloc((unsigned)length + 1);

  /*
   #############################################
   # begin calculations
   #############################################
   */
  update_fold_params();

  if (!pf) {
    hit_data  data;

    data.names    = names;
    data.strings  = AS;
    data.prefix   = prefix;
    data.columns  = aln_columns;
    data.md       = &md;
    data.ss_eps   = ssPS;
    data.msa_eps  = alnPS;
    data.msa_stk  = aln_out;
    data.csv      = csv;

    min_en = aliLfold_cb((const char **)AS, maxdist, &print_hit_cb, (void *)&data);
  }

  {
    eos_debug = -1; /* shut off warnings about nonstandard pairs */
    /*   for (i=0; AS[i]!=NULL; i++)
     * s += energy_of_struct(AS[i], structure);
     * real_en = s/i;*/
  }
  string = (mis) ? consens_mis((const char **)AS) : consensus((const char **)AS);
  printf("%s\n%s\n", string, structure);

  free(base_pair);
  (void)fflush(stdout);
  free(string);
  free(structure);
  free(prefix);

  for (i = 0; AS[i]; i++) {
    free(AS[i]);
    free(names[i]);
  }

  return EXIT_SUCCESS;
}


PRIVATE void
print_hit_cb(int start, int end, const char *structure, float en, void *data)
{
  char  *sub[500];
  char  *fname, *tmp_string;
  char  *id;
  char  **A, *cons;
  char  **names;
  char  **strings;
  char  *prefix;
  int        columns;
  vrna_md_t  *md;
  int         with_ss, with_msa, with_stk, with_csv;

  names   = ((hit_data *)data)->names;
  strings = ((hit_data *)data)->strings;
  prefix  = ((hit_data *)data)->prefix;
  columns = ((hit_data *)data)->columns;
  md      = ((hit_data *)data)->md;
  with_ss   = ((hit_data *)data)->ss_eps;
  with_msa  = ((hit_data *)data)->msa_eps;
  with_stk  = ((hit_data *)data)->msa_stk;
  with_csv  = ((hit_data *)data)->csv;

  get_subalignment((const char **)strings, sub, start, end);
  cons  = consensus((const char **)sub);
  A     = annote(structure, (const char **)sub, md);

  if (with_csv == 1)
    printf("%s ,%6.2f, %4d, %4d\n", structure, en, start, end);
  else
    printf("%s (%6.2f) %4d - %4d\n", structure, en, start, end);

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
    if (prefix) {
      id    = vrna_strdup_printf("%s_aln_%d_%d", prefix, start, end);
      fname = vrna_strdup_printf("%s.stk", prefix);

      tmp_string = vrna_filename_sanitize(fname, "_");
      free(fname);
      fname = tmp_string;

      vrna_file_msa_write(fname, (const char **)names, (const char **)sub, id, structure, "RNALalifold prediction", VRNA_FILE_FORMAT_MSA_STOCKHOLM | VRNA_FILE_FORMAT_MSA_APPEND);
      free(fname);
    } else {
      id = vrna_strdup_printf("aln_%d_%d", start, end);
      vrna_file_msa_write("RNALalifold_results.stk", (const char **)names, (const char **)sub, id, structure, "RNALalifold prediction", VRNA_FILE_FORMAT_MSA_STOCKHOLM | VRNA_FILE_FORMAT_MSA_APPEND);
    }

    free(id);
  }

  delete_alignment(sub);
}


PRIVATE char **
annote(const char *structure,
       const char *AS[],
       vrna_md_t *md)
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

PRIVATE void
get_subalignment(const char *AS[],
                 char       *sub[],
                 int        i,
                 int        j)
{
  int n_seq, s;

  /* get number of sequences in alignment */
  for (n_seq = 0; AS[n_seq] != NULL; n_seq++) ;

  for (s = 0; s < n_seq; s++)
    sub[s] = vrna_alloc(sizeof(char) * (j - i + 2));
  sub[s] = NULL;

  /* copy subalignment */
  for (s = 0; s < n_seq; s++) {
    sub[s]              = memcpy(sub[s], AS[s] + i - 1, sizeof(char) * (j - i + 1));
    sub[s][(j - i + 1)] = '\0';
  }
}


PRIVATE void
delete_alignment(char *AS[])
{
  int n_seq;

  /* get number of sequences in alignment */
  for (n_seq = 0; AS[n_seq] != NULL; n_seq++)
    free(AS[n_seq]);
}


