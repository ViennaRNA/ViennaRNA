#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>


#include <math.h>
#include <unistd.h>
#include <string.h>
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/constraints_SHAPE.h"
#include "ViennaRNA/perturbation_fold.h"
#include "ViennaRNA/file_formats.h"
#include "RNApvmin_cmdl.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static size_t     g_length    = 0;
static const char *g_statpath = 0;
static const char *g_sequence = 0;

static void
print_perturbation_vector(FILE    *f,
                          double  *epsilon)
{
  size_t i;

  for (i = 1; i <= g_length; ++i)
    fprintf(f, "%zu %c %f\n", i, g_sequence[i - 1], epsilon[i]);
}


static void
print_progress(int    iteration,
               double score,
               double *epsilon)
{
  FILE  *f;
  char  *path;

  vrna_message_info(stderr, "Iteration: %d\t Score: %f", iteration, score);

  if (!g_statpath)
    return;

  path  = vrna_strdup_printf("%s_%04d", g_statpath, iteration);
  f     = fopen(path, "w");
  if (!f) {
    vrna_message_warning("Couldn't open file '%s'", path);
    return;
  }

  fprintf(f, "#iteration %d\n#score %f\n", iteration, score);
  print_perturbation_vector(f, epsilon);

  fclose(f);
  free(path);
}


static void
init_perturbation_vector(double *epsilon,
                         int    length,
                         double max_energy)
{
  int i;

  if (max_energy == 0)
    return;

  for (i = 1; i <= length; ++i)
    epsilon[i] = max_energy * (vrna_urn() * 2 - 1);
}


int
main(int  argc,
     char *argv[])
{
  struct RNApvmin_args_info args_info;
  vrna_md_t                 md;
  char                      *rec_id, *rec_sequence, **rec_rest, *shape_sequence;
  size_t                    length;
  unsigned int              read_opt, rec_type;
  int                       istty, algorithm, i;
  double                    *shape_data, initialStepSize, minStepSize, minImprovement,
                            minimizerTolerance;

  istty               = 0;
  read_opt            = 0;
  initialStepSize     = 0.01;
  minStepSize         = 1e-15;
  minImprovement      = 1e-3;
  minimizerTolerance  = 1e-3;
  algorithm           = VRNA_MINIMIZER_DEFAULT;

  if (RNApvmin_cmdline_parser(argc, argv, &args_info))
    return 1;

  if (args_info.inputs_num != 1) {
    RNApvmin_cmdline_parser_print_help();
    return 1;
  }

  if (args_info.tauSigmaRatio_arg <= 0) {
    vrna_message_warning("invalid value for tauSigmaRatio");
    return 1;
  }

  vrna_init_rand();
  set_model_details(&md);

  /* set number of threads for parallel computation */
  if (args_info.numThreads_given)
#ifdef _OPENMP
    omp_set_num_threads(args_info.numThreads_arg);

#else
    vrna_message_error("\'j\' option is available only if compiled with OpenMP support!");
#endif

  if (args_info.paramFile_given)
    read_parameter_file(args_info.paramFile_arg);

  if (args_info.temp_given)
    md.temperature = args_info.temp_arg;

  if (args_info.noTetra_given)
    md.special_hp = 0;

  if (args_info.dangles_given) {
    if (args_info.dangles_given > 3)
      vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    else
      md.dangles = args_info.dangles_arg;
  }

  if (args_info.noLP_given)
    md.noLP = 1;

  if (args_info.noGU_given)
    md.noGU = 1;

  if (args_info.noClosingGU_given)
    md.noGUclosure = 1;

  if (args_info.energyModel_given)
    md.energy_set = args_info.energyModel_arg;

  if (args_info.intermediatePath_given)
    g_statpath = args_info.intermediatePath_arg;

  if (args_info.objectiveFunction_arg < 0 || args_info.objectiveFunction_arg > 1) {
    vrna_message_warning("required objective function mode not implemented, falling back to default");
    args_info.objectiveFunction_arg = 0;
  }

  if (args_info.minimizer_given) {
    struct {
      int algorithm;
      int arg;
    } mapper[] = { { VRNA_MINIMIZER_CONJUGATE_FR,     minimizer_arg_conjugate_fr     },
                   { VRNA_MINIMIZER_CONJUGATE_PR,     minimizer_arg_conjugate_pr     },
                   { VRNA_MINIMIZER_VECTOR_BFGS,      minimizer_arg_vector_bfgs      },
                   { VRNA_MINIMIZER_VECTOR_BFGS2,     minimizer_arg_vector_bfgs2     },
                   { VRNA_MINIMIZER_STEEPEST_DESCENT, minimizer_arg_steepest_descent },
                   { 0,                               0                              } };
    for (i = 0; mapper[i].algorithm; ++i)
      if (args_info.minimizer_arg == mapper[i].arg) {
        algorithm = mapper[i].algorithm;
        break;
      }
  }

  if (args_info.initialStepSize_given)
    initialStepSize = args_info.initialStepSize_arg;

  if (args_info.minStepSize_given)
    minStepSize = args_info.minStepSize_arg;

  if (args_info.minImprovement_given)
    minImprovement = args_info.minImprovement_arg;

  if (args_info.minimizerTolerance_given)
    minimizerTolerance = args_info.minimizerTolerance_arg;

  if (args_info.pfScale_given)
    md.sfact = args_info.pfScale_arg;

  istty = isatty(fileno(stdout)) && isatty(fileno(stdin));
  if (istty) {
    vrna_message_input_seq_simple();
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  }

  rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, NULL, read_opt);
  if (rec_type & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))
    return 0;

  md.uniq_ML = 1;

  g_sequence  = rec_sequence;
  g_length    = length = strlen(rec_sequence);

  shape_sequence  = vrna_alloc(sizeof(char) * (length + 1));
  shape_data      = vrna_alloc(sizeof(double) * (length + 1));

  if (vrna_file_SHAPE_read(args_info.inputs[0], length, -1, shape_sequence, shape_data)) {
    double                *epsilon;
    vrna_fold_compound_t  *vc;

    double                mfe;
    double                tau = 0.01;
    if (args_info.tauSigmaRatio_arg >= 10.)
      tau *= 10;

    if (args_info.tauSigmaRatio_arg >= 100.)
      tau *= 10;

    if (args_info.tauSigmaRatio_arg >= 1000.)
      tau *= 10;

    double sigma = tau / args_info.tauSigmaRatio_arg;

    vrna_sc_SHAPE_to_pr(args_info.shapeConversion_arg, shape_data, length, -1);

    vc  = vrna_fold_compound(rec_sequence, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    mfe = (double)vrna_mfe(vc, NULL);
    vrna_exp_params_rescale(vc, &mfe);

    epsilon = vrna_alloc(sizeof(double) * (length + 1));
    init_perturbation_vector(epsilon, length, args_info.initialVector_arg);
    vrna_sc_minimize_pertubation(vc,
                                 shape_data,
                                 args_info.objectiveFunction_arg,
                                 sigma, tau,
                                 algorithm,
                                 args_info.sampleSize_arg,
                                 epsilon,
                                 initialStepSize,
                                 minStepSize,
                                 minImprovement,
                                 minimizerTolerance,
                                 print_progress);

    vrna_fold_compound_free(vc);

    print_perturbation_vector(stdout, epsilon);

    free(epsilon);
  }

  free(shape_data);
  free(shape_sequence);

  free(rec_sequence);
  free(rec_id);
  for (i = 0; rec_rest[i]; ++i)
    free(rec_rest[i]);
  free(rec_rest);

  RNApvmin_cmdline_parser_free(&args_info);

  return 0;
}
