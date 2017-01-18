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
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/alifold.h"
#include "ViennaRNA/Lfold.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/read_epars.h"
#include "RNALalifold_cmdl.h"


#define MAX_NUM_NAMES    500

int
main(int  argc,
     char *argv[])
{
  FILE                          *clust_file;
  struct RNALalifold_args_info  args_info;
  char                          *string, *structure, *ParamFile, *ns_bases,
                                ffname[FILENAME_MAX_LENGTH], gfname[FILENAME_MAX_LENGTH],
                                *AS[MAX_NUM_NAMES], *names[MAX_NUM_NAMES];
  int                           n_seq, i, length, maxdist, unchangednc, unchangedcv,
                                mis, pf, istty;
  double                        min_en;
  vrna_md_t                     md;

  clust_file    = stdin;
  string        = structure = ParamFile = ns_bases = NULL;
  mis           = pf = 0;
  maxdist       = 70;
  do_backtrack  = unchangednc = unchangedcv = 1;
  dangles       = 2;
  ribo          = 0;

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
  if (!pf)
    min_en = aliLfold((const char **)AS, structure, maxdist);

  {
    eos_debug = -1; /* shut off warnings about nonstandard pairs */
    /*   for (i=0; AS[i]!=NULL; i++)
     * s += energy_of_struct(AS[i], structure);
     * real_en = s/i;*/
  }
  string = (mis) ? consens_mis((const char **)AS) : consensus((const char **)AS);
  printf("%s\n%s\n", string, structure);
  strcpy(ffname, "alirna.ps");
  strcpy(gfname, "alirna.g");

  free(base_pair);
  (void)fflush(stdout);
  free(string);
  free(structure);
  for (i = 0; AS[i]; i++) {
    free(AS[i]);
    free(names[i]);
  }
  return 0;
}
