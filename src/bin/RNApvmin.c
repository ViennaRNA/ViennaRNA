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
#include "ViennaRNA/pertubation_fold.h"
#include "RNApvmin_cmdl.h"

static void print_progress(int iteration, double score)
{
  fprintf(stderr, "Iteration: %d\t Score: %f\n", iteration, score);
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
  double *prob_unpaired;

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

  length  = strlen(rec_sequence);

  shape_sequence = space(sizeof(char) * (length + 1));
  prob_unpaired = space(sizeof(double) * (length + 1));

  if(parse_soft_constraints_file(args_info.inputs[0], length, 0.5, shape_sequence, prob_unpaired))
  {
    double *epsilon;
    pf_paramT *pf_parameters;
    vrna_fold_compound *vc;
    size_t i;

    epsilon = space(sizeof(double) * (length + 1));
    pf_parameters = vrna_get_boltzmann_factors(md);
    vc = get_fold_compound_pf_constrained(rec_sequence, NULL, NULL, pf_parameters);

    vrna_find_pertubation_vector(vc, prob_unpaired, args_info.sigma_arg, args_info.tau_arg, args_info.sampleSize_arg, epsilon, print_progress);
    free(pf_parameters);
    destroy_fold_compound(vc);

    for (i = 1; i <= length; ++i)
      printf("%zu %c %f\n", i, rec_sequence[i-1], epsilon[i]);

    free(epsilon);
  }

  free(prob_unpaired);
  free(shape_sequence);

  free(rec_sequence);
  free(rec_id);

  RNApvmin_cmdline_parser_free (&args_info);

  return 0;
}
