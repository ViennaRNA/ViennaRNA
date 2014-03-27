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
#include "ViennaRNA/perturbation_fold.h"
#include "ViennaRNA/file_formats.h"
#include "RNApvmin_cmdl.h"

static size_t g_length = 0;
static const char *g_statpath = 0;
static const char *g_sequence = 0;

static void print_perturbation_vector(FILE *f, double *epsilon)
{
  size_t i;
  for(i = 1; i <= g_length; ++i)
    fprintf(f, "%zu %c %f\n", i, g_sequence[i-1], epsilon[i]);
}

static void print_progress(int iteration, double score, double *epsilon)
{
  FILE *f;
  char path[256];

  fprintf(stderr, "Iteration: %d\t Score: %f\n", iteration, score);

  if(!g_statpath)
    return;

  sprintf(path, "%s_%04d", g_statpath, iteration);
  f = fopen(path, "w");
  if(!f) {
    fprintf(stderr, "Couldn't open file '%s'\n", path);
    return;
  }

  fprintf(f, "#iteration %d\n#score %f\n", iteration, score);
  print_perturbation_vector(f, epsilon);

  fclose(f);
}

static void init_perturbation_vector(double *epsilon, int length, double max_energy)
{
  int i;

  if (max_energy == 0)
    return;

  for (i = 1; i <= length; ++i)
    epsilon[i] = max_energy * (urn() * 2 - 1);
}

int main(int argc, char *argv[]){
  struct RNApvmin_args_info args_info;
  model_detailsT md;
  int istty = 0;
  unsigned int read_opt = 0;
  unsigned int rec_type;
  size_t length;
  char *rec_id;
  char *rec_sequence;
  char **rec_rest;
  char *shape_sequence;
  double *shape_data;
  int algorithm = VRNA_MINIMIZER_DEFAULT;
  int i;

  if (RNApvmin_cmdline_parser(argc, argv, &args_info))
    return 1;

  if(args_info.inputs_num != 1){
    RNApvmin_cmdline_parser_print_help();
    return 1;
  }

  init_rand();
  set_model_details(&md);

  if(args_info.paramFile_given)
    read_parameter_file(args_info.paramFile_arg);

  if(args_info.temp_given)
    md.temperature = args_info.temp_arg;
  if(args_info.noTetra_given)
    md.special_hp = 0;
  if(args_info.dangles_given){
    if(args_info.dangles_given > 3)
      warn_user("required dangle model not implemented, falling back to default dangles=2");
    else
      md.dangles = args_info.dangles_arg;
  }
  if(args_info.noLP_given)
    md.noLP =  1;
  if(args_info.noGU_given)
    md.noGU = 1;
  if(args_info.noClosingGU_given)
    md.noGUclosure = 1;
  if(args_info.energyModel_given)
    md.energy_set = args_info.energyModel_arg;
  if(args_info.intermediatePath_given)
    g_statpath = args_info.intermediatePath_arg;
  if(args_info.objectiveFunction_arg < 0 || args_info.objectiveFunction_arg > 1)
  {
    warn_user("required objective function mode not implemented, falling back to default");
    args_info.objectiveFunction_arg = 0;
  }
  if(args_info.minimizer_given)
  {
    struct {int algorithm; int arg;} mapper[] = {{VRNA_MINIMIZER_CONJUGATE_FR, minimizer_arg_conjugate_fr},
                                                 {VRNA_MINIMIZER_CONJUGATE_PR, minimizer_arg_conjugate_pr},
                                                 {VRNA_MINIMIZER_VECTOR_BFGS, minimizer_arg_vector_bfgs},
                                                 {VRNA_MINIMIZER_VECTOR_BFGS2, minimizer_arg_vector_bfgs2},
                                                 {VRNA_MINIMIZER_STEEPEST_DESCENT, minimizer_arg_steepest_descent},
                                                 {0, 0}};
    for(i = 0; mapper[i].algorithm; ++i)
      if(args_info.minimizer_arg == mapper[i].arg)
      {
        algorithm = mapper[i].algorithm;
        break;
      }
  }


  istty = isatty(fileno(stdout)) && isatty(fileno(stdin));
  if(istty)
  {
    print_tty_input_seq();
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  }

  rec_type = vrna_read_fasta_record(&rec_id, &rec_sequence, &rec_rest, NULL, read_opt);
  if (rec_type & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))
    return 0;

  md.uniq_ML = 1;

  g_sequence = rec_sequence;
  g_length = length  = strlen(rec_sequence);

  shape_sequence = space(sizeof(char) * (length + 1));
  shape_data = space(sizeof(double) * (length + 1));

  if(parse_soft_constraints_file(args_info.inputs[0], length, -1, shape_sequence, shape_data))
  {
    double *epsilon;
    vrna_fold_compound *vc;
    pf_paramT *pf_parameters;
    size_t i;
    float mfe;
    const double kT = (md.temperature + K0) * GASCONST / 1000.;

    convert_shape_reactivities_to_probabilities(args_info.shapeConversion_arg, shape_data, length, -1);

    vc = vrna_get_fold_compound(rec_sequence, &md, VRNA_OPTION_MFE);
    mfe = vrna_fold(vc, NULL);
    vrna_free_fold_compound(vc);

    vc = vrna_get_fold_compound(rec_sequence, &md, VRNA_OPTION_PF);
    pf_scale = exp(-(args_info.pfScale_arg * mfe) / kT / length);
    pf_parameters = get_boltzmann_factors(md.temperature, md.betaScale, md, pf_scale);
    vrna_update_pf_params(vc, pf_parameters);

    epsilon = space(sizeof(double) * (length + 1));
    init_perturbation_vector(epsilon, length, args_info.initialVector_arg);
    vrna_find_perturbation_vector(vc, shape_data, args_info.objectiveFunction_arg, args_info.sigma_arg, args_info.tau_arg, algorithm, args_info.sampleSize_arg, epsilon, print_progress);

    vrna_free_fold_compound(vc);
    free(pf_parameters);

    print_perturbation_vector(stdout, epsilon);

    free(epsilon);
  }

  free(shape_data);
  free(shape_sequence);

  free(rec_sequence);
  free(rec_id);
  for(i = 0; rec_rest[i]; ++i)
    free(rec_rest[i]);
  free(rec_rest);

  RNApvmin_cmdline_parser_free (&args_info);

  return 0;
}
