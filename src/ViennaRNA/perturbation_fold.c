#include "perturbation_fold.h"

#include "constraints.h"
#include "eval.h"
#include "fold_vars.h"
#include "part_func.h"
#include "utils.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void calculate_probability_unpaired(vrna_fold_compound *vc, double *probability)
{
  int length = vc->length;
  double *probs = vc->exp_matrices->probs;
  int *iidx = vc->iindx;
  int i, j;

  for (i = 0; i <= length; ++i)
    probability[i] = 1;

  for (i = 1; i <= length; ++i)
    for (j = i + 1; j <= length; ++j)
    {
      probability[i] -= probs[iidx[i]-j];
      probability[j] -= probs[iidx[i]-j];
    }
}

static double calculate_norm(double *vector, int length)
{
  double sum = 0;
  int i;

  for (i = 1; i <= length; ++i)
    sum += vector[i] * vector[i];

  return sqrt(sum);
}

static soft_constraintT* addSoftConstraint(vrna_fold_compound *vc, const double *epsilon, int length)
{
  soft_constraintT *sc;
  int i, j;
  double kT = vc->exp_params->kT / 1000;

  sc = space(sizeof(soft_constraintT));

  sc->boltzmann_factors = space(sizeof(double*) * (length + 2));
  sc->boltzmann_factors[0] = space(1);
  for (i = 1; i <= length; ++i)
    sc->boltzmann_factors[i] = space(sizeof(double) * (length - i + 2));

  for (i = 1; i <= length; ++i)
  {
    sc->boltzmann_factors[i][0] = 1;
    for (j = 1; j <= length - i + 1; ++j)
      sc->boltzmann_factors[i][j] = sc->boltzmann_factors[i][j-1] * exp(-epsilon[i + j - 1] / kT);
  }

  vc->sc = sc;
}

double vrna_evaluate_perturbation_vector_score(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared)
{
  double ret = 0;
  double *p_prob_unpaired;
  int i;
  soft_constraintT *sc;
  int length = vc->length;

  //calculate pairing probabilty in the pertubated energy model
  p_prob_unpaired = space(sizeof(double) * (length + 1));

  addSoftConstraint(vc, epsilon, length);
  vc->exp_params->model_details.do_backtrack = 1;
  vrna_pf_fold(vc, NULL);

  calculate_probability_unpaired(vc, p_prob_unpaired);

  remove_soft_constraints(vc);

  for (i = 1; i <= length; ++i)
  {
    double diff;

    //add penalty for pertubation energies
    ret += epsilon[i] * epsilon[i] / tau_squared;

    //add penalty for mismatches between observed and predicted probabilities
    diff = p_prob_unpaired[i] - q_prob_unpaired[i];
    ret += diff * diff / sigma_squared;
  }

  free(p_prob_unpaired);

  return ret;
}

void pairing_probabilities_from_restricted_pf(vrna_fold_compound *vc, const double *epsilon, double *prob_unpaired, double **conditional_prob_unpaired)
{
  int length = vc->length;
  int i;
  char *hc_string;
  hard_constraintT *hc_backup;

  addSoftConstraint(vc, epsilon, length);
  vc->exp_params->model_details.do_backtrack = 1;

  vrna_pf_fold(vc, NULL);
  calculate_probability_unpaired(vc, prob_unpaired);

  hc_string = space(sizeof(char) * (length + 1));
  hc_backup = vc->hc;

  for (i = 1; i <= length; ++i)
  {
    unsigned int constraint_options = VRNA_CONSTRAINT_IINDX
                                      | VRNA_CONSTRAINT_DB
                                      | VRNA_CONSTRAINT_PIPE
                                      | VRNA_CONSTRAINT_DOT
                                      | VRNA_CONSTRAINT_X
                                      | VRNA_CONSTRAINT_ANG_BRACK
                                      | VRNA_CONSTRAINT_RND_BRACK;

    memset(hc_string, '.', length);
    hc_string[i - 1] = 'x';

    vc->hc = get_hard_constraints(vc->sequence, hc_string, &(vc->exp_params->model_details), TURN, constraint_options);
    vrna_pf_fold(vc, NULL);
    calculate_probability_unpaired(vc, conditional_prob_unpaired[i]);

    destroy_hard_constraints(vc->hc);
    vc->hc = NULL;
  }

  remove_soft_constraints(vc);
  vc->hc = hc_backup;

  free(hc_string);
}

void pairing_probabilities_from_sampling(vrna_fold_compound *vc, const double *epsilon, int sample_size, double *prob_unpaired, double **conditional_prob_unpaired)
{
  int length = vc->length;
  int i, j, s;
  st_back = 1;

  addSoftConstraint(vc, epsilon, length);
  vc->exp_params->model_details.do_backtrack = 0;

  vrna_pf_fold(vc, NULL);

  for (s = 0; s < sample_size; ++s)
  {
    char *sample = vrna_pbacktrack5(vc->length, vc);

    for (i = 1; i <= length; ++i)
    {
      if (sample[i-1] != '.')
        continue;

      ++prob_unpaired[i];

      for (j = 1; j <= length; ++j)
        if (sample[j-1] == '.')
          ++conditional_prob_unpaired[i][j];
    }

    free(sample);
  }

  for (i = 1; i <= length; ++i)
  {
    if (prob_unpaired[i])
      for (j = 1; j <= length; ++j)
        conditional_prob_unpaired[i][j] /= prob_unpaired[i];

    prob_unpaired[i] /= sample_size;

    assert(prob_unpaired[i] >= 0 && prob_unpaired[i] <= 1);
  }

  remove_soft_constraints(vc);
}

static void allocateProbabilityArrays(double **unpaired, double ***conditional_unpaired, int length)
{
  int i;

  *unpaired = space(sizeof(double) * (length + 1));
  *conditional_unpaired = space(sizeof(double*) * (length + 1));

  for (i = 1; i <= length; ++i)
    (*conditional_unpaired)[i] = space(sizeof(double) * (length + 1));
}

static void freeProbabilityArrays(double *unpaired, double **conditional_unpaired, int length)
{
  int i;

  free(unpaired);
  for (i = 1; i <= length; ++i)
    free(conditional_unpaired[i]);
  free(conditional_unpaired);
}

void vrna_evaluate_perturbation_vector_gradient(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int sample_size, double *gradient)
{
  double *p_prob_unpaired;
  double **p_conditional_prob_unpaired;
  int i, mu;
  int length = vc->length;
  double kT = vc->exp_params->kT / 1000;

  allocateProbabilityArrays(&p_prob_unpaired, &p_conditional_prob_unpaired, length);

  if (sample_size > 0)
    pairing_probabilities_from_sampling(vc, epsilon, sample_size, p_prob_unpaired, p_conditional_prob_unpaired);
  else
    pairing_probabilities_from_restricted_pf(vc, epsilon, p_prob_unpaired, p_conditional_prob_unpaired);

  for (mu = 1; mu <= length; ++mu)
  {
    double sum = 0;

    for (i = 1; i <= length; ++i)
    {
      sum +=  (p_prob_unpaired[i] - q_prob_unpaired[i]) *
              p_prob_unpaired[i] *
              (p_prob_unpaired[mu] - p_conditional_prob_unpaired[i][mu]) /
              sigma_squared;
    }

    gradient[mu] = 2 * (epsilon[mu] / tau_squared + sum/kT);
  }

  freeProbabilityArrays(p_prob_unpaired, p_conditional_prob_unpaired, length);
}

void vrna_find_perturbation_vector(vrna_fold_compound *vc, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int sample_size, double *epsilon, progress_callback callback)
{
  int iteration_count = 0;
  int length = vc->length;
  double improvement;

  double *new_epsilon = space(sizeof(double) * (length + 1));
  double *gradient = space(sizeof(double) * (length + 1));

  double score = vrna_evaluate_perturbation_vector_score(vc, epsilon, q_prob_unpaired, sigma_squared, tau_squared);

  if (callback)
    callback(0, score);

  do
  {
    double new_score;
    double step_size;

    ++iteration_count;

    vrna_evaluate_perturbation_vector_gradient(vc, epsilon, q_prob_unpaired, sigma_squared, tau_squared, sample_size, gradient);

    step_size = 0.5 / calculate_norm(gradient, length);

    do
    {
      int i;
      for (i = 1; i <= length; ++i)
        new_epsilon[i] = epsilon[i] - step_size * gradient[i];

      new_score = vrna_evaluate_perturbation_vector_score(vc, new_epsilon, q_prob_unpaired, sigma_squared, tau_squared);
      step_size /= 2;
    } while (new_score > score && step_size >= 1e-15);

    if (new_score > score)
      break;

    if (callback)
      callback(iteration_count, new_score);

    improvement = 1 - new_score / score;
    score = new_score;
    memcpy(epsilon, new_epsilon, sizeof(double) * (length+1));
  } while (improvement >= 0.0001);

  free(gradient);
  free(new_epsilon);
}

