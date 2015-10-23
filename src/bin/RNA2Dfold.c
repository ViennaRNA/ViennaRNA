
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

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct nbhoods{
  int k;
  int l;
  struct nbhoods *next;
} nbhoods;

int main(int argc, char *argv[]){
  struct        RNA2Dfold_args_info args_info;
  unsigned int  input_type;
  char *string, *input_string, *orig_sequence;
  char *mfe_structure=NULL, *structure1=NULL, *structure2=NULL, *reference_struc1=NULL, *reference_struc2=NULL;
  char  *ParamFile=NULL;
  int   i, j, length, l;
  double min_en;
  int   pf=0,istty;
  int noconv=0;
  int circ=0;
  int maxDistance1 = -1;
  int maxDistance2 = -1;
  int do_backtrack = 1;
  int stBT = 0;
  int nstBT = 0;
  string=NULL;
  dangles = 2;
  struct nbhoods *neighborhoods = NULL;
  struct nbhoods *neighborhoods_cur = NULL;
  vrna_md_t   md;

  vrna_md_set_default(&md);

  string = input_string = orig_sequence = NULL;
  /*
  #############################################
  # check the command line prameters
  #############################################
  */
  if(RNA2Dfold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);

  /* temperature */
  if(args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;

  /* max distance to 1st reference structure */
  if(args_info.maxDist1_given)
    maxDistance1 = args_info.maxDist1_arg;

  /* max distance to 2nd reference structure */
  if(args_info.maxDist2_given)
    maxDistance2 = args_info.maxDist2_arg;

  /* compute partition function and boltzmann probabilities */
  if(args_info.partfunc_given)
    pf = 1;

  /* do stachastic backtracking */
  if(args_info.stochBT_given){
    pf = 1;
    stBT = 1;
    nstBT = args_info.stochBT_arg;
  }

  if(args_info.noTetra_given)
    md.special_hp = tetra_loop = 0;

  /* assume RNA sequence to be circular */
  if(args_info.circ_given)
    md.circ = circ = 1;

  /* dangle options */
  if(args_info.dangles_given){
    if((args_info.dangles_arg != 0) && (args_info.dangles_arg != 2))
      vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    else
      md.dangles = dangles = args_info.dangles_arg;
  }
  /* set number of threads for parallel computation */
  if(args_info.numThreads_given)
#ifdef _OPENMP
  omp_set_num_threads(args_info.numThreads_arg);
#else
  vrna_message_error("\'j\' option is available only if compiled with OpenMP support!");
#endif

  /* get energy parameter file name */
  if(args_info.parameterFile_given)
    ParamFile = strdup(args_info.parameterFile_arg);

  /* do not allow GU pairs ? */
  if(args_info.noGU_given)
    md.noGU = noGU = 1;

  /* do not allow GU pairs at the end of helices? */
  if(args_info.noClosingGU_given)
    md.noGUclosure = no_closingGU = 1;

  /* pf scaling factor */
  if(args_info.pfScale_given)
    md.sfact = args_info.pfScale_arg;

  /* do not backtrack structures ? */
  if(args_info.noBT_given)
    md.backtrack = do_backtrack = 0;

  for (i = 0; i < args_info.neighborhood_given; i++){
    int kappa, lambda;
    kappa = lambda = 0;
    if(sscanf(args_info.neighborhood_arg[i], "%d:%d", &kappa, &lambda) == 2);
    if ((kappa>-2) && (lambda>-2)){
      if(neighborhoods_cur != NULL){
        neighborhoods_cur->next = (nbhoods *)vrna_alloc(sizeof(nbhoods));
        neighborhoods_cur = neighborhoods_cur->next;
      }
      else{
        neighborhoods = (nbhoods *)vrna_alloc(sizeof(nbhoods));
        neighborhoods_cur = neighborhoods;
      }
      neighborhoods_cur->k = kappa;
      neighborhoods_cur->l = lambda;
      neighborhoods_cur->next = NULL;
    }
  }
  /* free allocated memory of command line data structure */
  RNA2Dfold_cmdline_parser_free (&args_info);

  /*
  #############################################
  # begin actual program code
  #############################################
  */
  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  istty = isatty(fileno(stdout))&&isatty(fileno(stdin));

 /*
  #############################################
  # main loop, continue until end of file
  #############################################
  */
  do {
    if (istty)
      vrna_message_input_seq("Input strings\n1st line: sequence (upper or lower case)\n2nd + 3rd line: reference structures (dot bracket notation)\n@ to quit\n");

    while((input_type = get_input_line(&input_string, 0)) & VRNA_INPUT_FASTA_HEADER){
      printf(">%s\n", input_string); /* print fasta header if available */
      free(input_string);
    }

    /* break on any error, EOF or quit request */
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)){ break;}
    /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else{
      length = (int)    strlen(input_string);
      string = strdup(input_string);
      free(input_string);
    }

    mfe_structure = (char *) vrna_alloc((unsigned) length+1);
    structure1 = (char *) vrna_alloc((unsigned) length+1);
    structure2 = (char *) vrna_alloc((unsigned) length+1);

    input_type = get_input_line(&input_string, VRNA_INPUT_NOSKIP_COMMENTS);
    if(input_type & VRNA_INPUT_QUIT){ break;}
    else if((input_type & VRNA_INPUT_MISC) && (strlen(input_string) > 0)){
      reference_struc1 = strdup(input_string);
      free(input_string);
      if(strlen(reference_struc1) != length)
        vrna_message_error("sequence and 1st reference structure have unequal length");
    }
    else vrna_message_error("1st reference structure missing\n");
    strncpy(structure1, reference_struc1, length);

    input_type = get_input_line(&input_string, VRNA_INPUT_NOSKIP_COMMENTS);
    if(input_type & VRNA_INPUT_QUIT){ break;}
    else if((input_type & VRNA_INPUT_MISC) && (strlen(input_string) > 0)){
      reference_struc2 = strdup(input_string);
      free(input_string);
      if(strlen(reference_struc2) != length)
        vrna_message_error("sequence and 2nd reference structure have unequal length");
    }
    else vrna_message_error("2nd reference structure missing\n");
    strncpy(structure2, reference_struc2, length);

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv) vrna_seq_toRNA(string);
    /* store case-unmodified sequence */
    orig_sequence = strdup(string);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(string);

    if (istty)  printf("length = %d\n", length);

    vrna_fold_compound  *vc_global = vrna_get_fold_compound(string,&md, VRNA_OPTION_MFE);

    min_en = vrna_fold(vc_global, mfe_structure);

    printf("%s\n%s", orig_sequence, mfe_structure);

    if (istty)
      printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
    else
      printf(" (%6.2f)\n", min_en);

    printf("%s (%6.2f) <ref 1>\n", structure1, vrna_eval_structure(vc_global, structure1));
    printf("%s (%6.2f) <ref 2>\n", structure2, vrna_eval_structure(vc_global, structure2));

    vrna_free_fold_compound(vc_global);

    /* get all variables need for the folding process (some memory will be preallocated here too) */
    vrna_fold_compound *vc  = vrna_get_fold_compound_2D(string, structure1, structure2, &md, VRNA_OPTION_MFE | (pf ? VRNA_OPTION_PF : 0));
    vrna_sol_TwoD_t *mfe_s  = vrna_TwoD_fold(vc, maxDistance1, maxDistance2);

    if(!pf){
#ifdef COUNT_STATES
      printf("k\tl\tn\tMFE\tMFE-structure\n");
      for(i = 0; mfe_s[i].k != INF; i++){
        printf("%d\t%d\t%lu\t%6.2f\t%s\n", mfe_s[i].k, mfe_s[i].l, vc->N_F5[length][mfe_s[i].k][mfe_s[i].l/2], mfe_s[i].en, mfe_s[i].s);
        if(mfe_s[i].s) free(mfe_s[i].s);
      }
      free(mfe_s);
#else
      printf("k\tl\tMFE\tMFE-structure\n");
      for(i = 0; mfe_s[i].k != INF; i++){
        printf("%d\t%d\t%6.2f\t%s\n", mfe_s[i].k, mfe_s[i].l, mfe_s[i].en, mfe_s[i].s);
        if(mfe_s[i].s) free(mfe_s[i].s);
      }
      free(mfe_s);
#endif
    }

    if(pf){
      int maxD1 = (int) vc->maxD1;
      int maxD2 = (int) vc->maxD2;
      double mmfe = INF;
      double Q;
      for(i = 0; mfe_s[i].k != INF; i++){
        if(mmfe > mfe_s[i].en) mmfe = (double)mfe_s[i].en;
      }
      vrna_exp_params_rescale(vc, &mmfe);

      /* we dont need the mfe DP arrays anymore, so we can savely free their occupying memory */
      vrna_mx_mfe_free(vc);

      vrna_sol_TwoD_pf_t *pf_s = vrna_TwoD_pf_fold(vc, maxD1, maxD2);

      Q = 0.;
      
      for(i = 0; pf_s[i].k != INF; i++){
        Q += pf_s[i].q;
      }

      double fee = (-log(Q)-length*log(vc->exp_params->pf_scale))*(vc->exp_params->kT/1000.);

      if(!stBT){
        printf("free energy of ensemble = %6.2f kcal/mol\n",fee);
        printf("k\tl\tP(neighborhood)\tP(MFE in neighborhood)\tP(MFE in ensemble)\tMFE\tE_gibbs\tMFE-structure\n");
        for(i=0; pf_s[i].k != INF;i++){
          float free_energy = (-log((float)pf_s[i].q)-length*log(vc->exp_params->pf_scale))*(vc->exp_params->kT/1000.);
          if((pf_s[i].k != mfe_s[i].k) || (pf_s[i].l != mfe_s[i].l))
            vrna_message_error("This should never happen!");
          fprintf(stdout,
                  "%d\t%d\t%2.8f\t%2.8f\t%2.8f\t%6.2f\t%6.2f\t%s\n",
                  pf_s[i].k,
                  pf_s[i].l,
                  (float)(pf_s[i].q)/(float)Q,
                  exp((free_energy-mfe_s[i].en)/(vc->exp_params->kT/1000.)),
                  exp((fee-mfe_s[i].en)/(vc->exp_params->kT/1000.)),
                  mfe_s[i].en,
                  free_energy,
                  mfe_s[i].s);
        }
      }
      else{
        vrna_init_rand();
        printf("k\tl\ten\tstructure\n");
        if(neighborhoods != NULL){
          nbhoods *tmp, *tmp2;
          for(tmp = neighborhoods; tmp != NULL; tmp = tmp->next){
            int k,l;
            k = tmp->k;
            l = tmp->l;
            for(i = 0; i < nstBT; i++){
              char *s = vrna_TwoD_pbacktrack(vc, k, l);
              printf("%d\t%d\t%6.2f\t%s\n", k, l, vrna_eval_structure(vc, s), s);
            }
          }
        }
        else{
          for(i=0; pf_s[i].k != INF;i++){
            for(l = 0; l < nstBT; l++){
              char *s = vrna_TwoD_pbacktrack(vc, pf_s[i].k, pf_s[i].l);
              printf("%d\t%d\t%6.2f\t%s\n", pf_s[i].k, pf_s[i].l, vrna_eval_structure(vc, s), s);
            }
          }
        }
      }

      for(i=0; mfe_s[i].k != INF;i++){
        if(mfe_s[i].s) free(mfe_s[i].s);
      }
      free(pf_s);
      free(mfe_s);
    }

    vrna_free_fold_compound(vc);

    free(string);
    free(orig_sequence);
    free(mfe_structure);
    free(structure1);
    free(structure2);
    free(reference_struc1);
    free(reference_struc2);
    string = orig_sequence = mfe_structure = NULL;
  } while (1);
  return 0;
}
