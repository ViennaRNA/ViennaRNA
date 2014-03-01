#include <stdio.h>
#include <stdlib.h>


#include <unistd.h>
#include <string.h>
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/perturbation_fold.h"
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

  if (RNApvmin_cmdline_parser(argc, argv, &args_info))
    return 1;

  if(args_info.inputs_num != 1){
    RNApvmin_cmdline_parser_print_help();
    return 1;
  }

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


  istty = isatty(fileno(stdout)) && isatty(fileno(stdin));
  if(istty)
  {
    print_tty_input_seq();
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  }

  rec_type = read_record(&rec_id, &rec_sequence, &rec_rest, read_opt);
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
    size_t i;

    //use a cuttoff approach to divide into paired/unpaired
    for (i = 1; i <= length; ++i)
      shape_data[i] = shape_data[i] < args_info.cutoff_arg ? 0 : 1;

    epsilon = space(sizeof(double) * (length + 1));
    vc = vrna_get_fold_compound(rec_sequence, &md, VRNA_OPTION_PF);

    vrna_find_perturbation_vector(vc, shape_data, args_info.sigma_arg, args_info.tau_arg, args_info.sampleSize_arg, epsilon, print_progress);

    destroy_fold_compound(vc);

    print_perturbation_vector(stdout, epsilon);

    free(epsilon);
  }

  free(shape_data);
  free(shape_sequence);

  free(rec_sequence);
  free(rec_id);

  RNApvmin_cmdline_parser_free (&args_info);

  return 0;
}
