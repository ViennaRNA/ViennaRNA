#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/2Dfold.h"
#include "ViennaRNA/2Dpfold.h"
#include "RNA2Dfold_cmdl.h"

#include "ViennaRNA/color_output.inc"

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct nbhoods {
  int             k;
  int             l;
  struct nbhoods  *next;
} nbhoods;

int
main(int  argc,
     char *argv[])
{
  struct        RNA2Dfold_args_info args_info;
  struct        nbhoods             *neighborhoods, *neighborhoods_cur;
  unsigned int                      input_type;
  char                              *string, *input_string, *orig_sequence;
  char                              *mfe_structure, *structure1, *structure2, *reference_struc1,
                                    *reference_struc2, *ParamFile;
  int                               i, length, l, pf, istty, noconv, circ, maxDistance1, maxDistance2,
                                    do_backtrack, stBT, nstBT;
  double                            min_en;
  vrna_md_t                         md;

  string            = input_string = orig_sequence = ParamFile = NULL;
  mfe_structure     = structure1 = structure2 = reference_struc1 = reference_struc2 = NULL;
  dangles           = 2;
  pf                = 0;
  noconv            = 0;
  circ              = 0;
  maxDistance1      = -1;
  maxDistance2      = -1;
  do_backtrack      = 1;
  stBT              = 0;
  nstBT             = 0;
  neighborhoods     = NULL;
  neighborhoods_cur = NULL;

  vrna_md_set_default(&md);


  /*
   #############################################
   # check the command line prameters
   #############################################
   */
  if (RNA2Dfold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* temperature */
  if (args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;

  /* max distance to 1st reference structure */
  if (args_info.maxDist1_given)
    maxDistance1 = args_info.maxDist1_arg;

  /* max distance to 2nd reference structure */
  if (args_info.maxDist2_given)
    maxDistance2 = args_info.maxDist2_arg;

  /* compute partition function and boltzmann probabilities */
  if (args_info.partfunc_given)
    pf = 1;

  /* do stachastic backtracking */
  if (args_info.stochBT_given) {
    pf    = 1;
    stBT  = 1;
    nstBT = args_info.stochBT_arg;
  }

  if (args_info.noTetra_given)
    md.special_hp = tetra_loop = 0;

  /* assume RNA sequence to be circular */
  if (args_info.circ_given)
    md.circ = circ = 1;

  /* dangle options */
  if (args_info.dangles_given) {
    if ((args_info.dangles_arg != 0) && (args_info.dangles_arg != 2))
      vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    else
      md.dangles = dangles = args_info.dangles_arg;
  }

  /* set number of threads for parallel computation */
  if (args_info.numThreads_given)
#ifdef _OPENMP
    omp_set_num_threads(args_info.numThreads_arg);

#else
    vrna_message_error("\'j\' option is available only if compiled with OpenMP support!");
#endif

  /* get energy parameter file name */
  if (args_info.parameterFile_given)
    ParamFile = strdup(args_info.parameterFile_arg);

  /* do not allow GU pairs ? */
  if (args_info.noGU_given)
    md.noGU = noGU = 1;

  /* do not allow GU pairs at the end of helices? */
  if (args_info.noClosingGU_given)
    md.noGUclosure = no_closingGU = 1;

  /* pf scaling factor */
  if (args_info.pfScale_given)
    md.sfact = args_info.pfScale_arg;

  /* do not backtrack structures ? */
  if (args_info.noBT_given)
    md.backtrack = do_backtrack = 0;

  for (i = 0; i < args_info.neighborhood_given; i++) {
    int kappa, lambda;
    kappa = lambda = 0;
    if (sscanf(args_info.neighborhood_arg[i], "%d:%d", &kappa, &lambda) == 2) {
      if ((kappa > -2) && (lambda > -2)) {
        if (neighborhoods_cur != NULL) {
          neighborhoods_cur->next = (nbhoods *)vrna_alloc(sizeof(nbhoods));
          neighborhoods_cur       = neighborhoods_cur->next;
        } else {
          neighborhoods     = (nbhoods *)vrna_alloc(sizeof(nbhoods));
          neighborhoods_cur = neighborhoods;
        }

        neighborhoods_cur->k    = kappa;
        neighborhoods_cur->l    = lambda;
        neighborhoods_cur->next = NULL;
      }
    }
  }
  /* free allocated memory of command line data structure */
  RNA2Dfold_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin actual program code
   #############################################
   */
  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  istty = isatty(fileno(stdout)) && isatty(fileno(stdin));

  /*
   #############################################
   # main loop, continue until end of file
   #############################################
   */
  do {
    if (istty)
      vrna_message_input_seq("Input strings\n1st line: sequence (upper or lower case)\n2nd + 3rd line: reference structures (dot bracket notation)\n@ to quit\n");

    char *rec_id = NULL;

    while ((input_type = get_input_line(&input_string, 0)) & VRNA_INPUT_FASTA_HEADER)
      rec_id = input_string;

    /* break on any error, EOF or quit request */
    if (input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) {
      break;
    }
    /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else {
      length  = (int)strlen(input_string);
      string  = strdup(input_string);
      free(input_string);
    }

    mfe_structure = (char *)vrna_alloc((unsigned)length + 1);
    structure1    = (char *)vrna_alloc((unsigned)length + 1);
    structure2    = (char *)vrna_alloc((unsigned)length + 1);

    input_type = get_input_line(&input_string, VRNA_INPUT_NOSKIP_COMMENTS);
    if (input_type & VRNA_INPUT_QUIT) {
      break;
    } else if ((input_type & VRNA_INPUT_MISC) && (strlen(input_string) > 0)) {
      reference_struc1 = strdup(input_string);
      free(input_string);
      if (strlen(reference_struc1) != length)
        vrna_message_error("sequence and 1st reference structure have unequal length");
    } else {
      vrna_message_error("1st reference structure missing\n");
    }

    strncpy(structure1, reference_struc1, length);

    input_type = get_input_line(&input_string, VRNA_INPUT_NOSKIP_COMMENTS);
    if (input_type & VRNA_INPUT_QUIT) {
      break;
    } else if ((input_type & VRNA_INPUT_MISC) && (strlen(input_string) > 0)) {
      reference_struc2 = strdup(input_string);
      free(input_string);
      if (strlen(reference_struc2) != length)
        vrna_message_error("sequence and 2nd reference structure have unequal length");
    } else {
      vrna_message_error("2nd reference structure missing\n");
    }

    strncpy(structure2, reference_struc2, length);

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if (!noconv)
      vrna_seq_toRNA(string);

    /* store case-unmodified sequence */
    orig_sequence = strdup(string);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(string);

    if (istty)
      printf("length = %d\n", length);

    vrna_fold_compound_t *vc_global = vrna_fold_compound(string, &md, VRNA_OPTION_MFE);

    min_en = vrna_mfe(vc_global, mfe_structure);

    print_fasta_header(stdout, rec_id);
    fprintf(stdout, "%s\n", orig_sequence);

    char *msg = NULL;
    if (istty)
      msg = vrna_strdup_printf("\n minimum free energy = %6.2f kcal/mol", min_en);
    else
      msg = vrna_strdup_printf(" (%6.2f)", min_en);

    print_structure(stdout, mfe_structure, msg);
    free(msg);
    (void)fflush(stdout);

    msg = vrna_strdup_printf(" (%6.2f) <ref 1>", vrna_eval_structure(vc_global, structure1));
    print_structure(stdout, structure1, msg);
    free(msg);

    msg = vrna_strdup_printf(" (%6.2f) <ref 2>", vrna_eval_structure(vc_global, structure2));
    print_structure(stdout, structure2, msg);
    free(msg);

    vrna_fold_compound_free(vc_global);

    /* get all variables need for the folding process (some memory will be preallocated here too) */
    vrna_fold_compound_t  *vc     = vrna_fold_compound_TwoD(string, structure1, structure2, &md, VRNA_OPTION_MFE | (pf ? VRNA_OPTION_PF : 0));
    vrna_sol_TwoD_t       *mfe_s  = vrna_mfe_TwoD(vc, maxDistance1, maxDistance2);

    if (!pf) {
#ifdef COUNT_STATES
      print_table(stdout, "k\tl\tn\tMFE\tMFE-structure", NULL);
      for (i = 0; mfe_s[i].k != INF; i++) {
        char *tline = vrna_strdup_printf("%d\t%d\t%lu\t%6.2f\t%s",
                                         mfe_s[i].k, mfe_s[i].l,
                                         vc->N_F5[length][mfe_s[i].k][mfe_s[i].l / 2],
                                         mfe_s[i].en, mfe_s[i].s);
        print_table(stdout, NULL, tline);
        if (mfe_s[i].s)
          free(mfe_s[i].s);

        free(tline);
      }
      free(mfe_s);
#else
      print_table(stdout, "k\tl\tMFE\tMFE-structure", NULL);
      for (i = 0; mfe_s[i].k != INF; i++) {
        char *tline = vrna_strdup_printf("%d\t%d\t%6.2f\t%s",
                                         mfe_s[i].k, mfe_s[i].l,
                                         mfe_s[i].en, mfe_s[i].s);
        print_table(stdout, NULL, tline);
        if (mfe_s[i].s)
          free(mfe_s[i].s);

        free(tline);
      }
      free(mfe_s);
#endif
    }

    if (pf) {
      int     maxD1 = (int)vc->maxD1;
      int     maxD2 = (int)vc->maxD2;
      double  mmfe  = INF;
      double  Q;
      for (i = 0; mfe_s[i].k != INF; i++)
        if (mmfe > mfe_s[i].en)
          mmfe = (double)mfe_s[i].en;

      vrna_exp_params_rescale(vc, &mmfe);

      /* we dont need the mfe DP arrays anymore, so we can savely free their occupying memory */
      vrna_mx_mfe_free(vc);

      vrna_sol_TwoD_pf_t *pf_s = vrna_pf_TwoD(vc, maxD1, maxD2);

      Q = 0.;

      for (i = 0; pf_s[i].k != INF; i++)
        Q += pf_s[i].q;

      double fee = (-log(Q) - length * log(vc->exp_params->pf_scale)) * (vc->exp_params->kT / 1000.);

      if (!stBT) {
        char *msg = NULL;
        msg = vrna_strdup_printf("free energy of ensemble = %6.2f kcal/mol", fee);
        print_structure(stdout, NULL, msg);
        free(msg);
        print_table(stdout, "k\tl\tP(neighborhood)\tP(MFE in neighborhood)\tP(MFE in ensemble)\tMFE\tE_gibbs\tMFE-structure", NULL);
        for (i = 0; pf_s[i].k != INF; i++) {
          float free_energy = (-log((float)pf_s[i].q) - length * log(vc->exp_params->pf_scale)) * (vc->exp_params->kT / 1000.);
          if ((pf_s[i].k != mfe_s[i].k) || (pf_s[i].l != mfe_s[i].l))
            vrna_message_error("This should never happen!");

          char  *tline = vrna_strdup_printf("%d\t%d\t%2.8f\t%2.8f\t%2.8f\t%6.2f\t%6.2f\t%s",
                                            pf_s[i].k,
                                            pf_s[i].l,
                                            (float)(pf_s[i].q) / (float)Q,
                                            exp((free_energy - mfe_s[i].en) / (vc->exp_params->kT / 1000.)),
                                            exp((fee - mfe_s[i].en) / (vc->exp_params->kT / 1000.)),
                                            mfe_s[i].en,
                                            free_energy,
                                            mfe_s[i].s);
          print_table(stdout, NULL, tline);
          free(tline);
        }
      } else {
        vrna_init_rand();
        print_table(stdout, "k\tl\ten\tstructure", NULL);
        if (neighborhoods != NULL) {
          nbhoods *tmp;
          for (tmp = neighborhoods; tmp != NULL; tmp = tmp->next) {
            int k, l;
            k = tmp->k;
            l = tmp->l;
            for (i = 0; i < nstBT; i++) {
              char  *s      = vrna_pbacktrack_TwoD(vc, k, l);
              char  *tline  = vrna_strdup_printf("%d\t%d\t%6.2f\t%s", k, l, vrna_eval_structure(vc, s), s);
              print_table(stdout, NULL, tline);
              free(tline);
              free(s);
            }
          }
        } else {
          for (i = 0; pf_s[i].k != INF; i++) {
            for (l = 0; l < nstBT; l++) {
              char  *s      = vrna_pbacktrack_TwoD(vc, pf_s[i].k, pf_s[i].l);
              char  *tline  = vrna_strdup_printf("%d\t%d\t%6.2f\t%s", pf_s[i].k, pf_s[i].l, vrna_eval_structure(vc, s), s);
              print_table(stdout, NULL, tline);
              free(tline);
              free(s);
            }
          }
        }
      }

      for (i = 0; mfe_s[i].k != INF; i++)
        if (mfe_s[i].s)
          free(mfe_s[i].s);

      free(pf_s);
      free(mfe_s);
    }

    vrna_fold_compound_free(vc);

    free(string);
    free(orig_sequence);
    free(mfe_structure);
    free(structure1);
    free(structure2);
    free(reference_struc1);
    free(reference_struc2);
    free(rec_id);
    string = orig_sequence = mfe_structure = rec_id = NULL;
  } while (1);
  return 0;
}
