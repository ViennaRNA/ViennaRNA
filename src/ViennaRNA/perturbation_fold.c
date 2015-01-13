#include "perturbation_fold.h"

#include <config.h>
#include "eval.h"
#include "fold_vars.h"
#include "constraints.h"
#include "fold.h"
#include "part_func.h"
#include "utils.h"
#include "params.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WITH_GSL
#include <gsl/gsl_multimin.h>
#endif

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

static void addSoftConstraint(vrna_fold_compound *vc, const double *epsilon, int length)
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
      sc->boltzmann_factors[i][j] = sc->boltzmann_factors[i][j-1] * exp(-(epsilon[i + j - 1]) / kT);
  }

  /* also add sc for MFE computation */
  sc->free_energies = space(sizeof(int*) * (length + 2));
  sc->free_energies[0] = space(sizeof(int));
  for (i = 1; i <= length; ++i)
    sc->free_energies[i] = space(sizeof(int) * (length - i + 2));

  for (i = 1; i <= length; ++i){
    sc->free_energies[i][0] = 0;
    for (j = 1; j <= length - i + 1; ++j)
      sc->free_energies[i][j] = sc->free_energies[i][j-1] + (epsilon[i + j - 1]*100.);
  }

  vc->sc = sc;
}

static double evaluate_objective_function_contribution(double value, int objective_function)
{
  if (objective_function == VRNA_OBJECTIVE_FUNCTION_QUADRATIC)
    return value * value;
  if (objective_function == VRNA_OBJECTIVE_FUNCTION_ABSOLUTE)
    return fabs(value);

  assert(0);
  return 0;
}

static double evaluate_perturbation_vector_score(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int objective_function)
{
  double kT, ret = 0;
  double ret2 = 0.;
  double *p_prob_unpaired;
  int i;
  int length = vc->length;

  /* calculate pairing probabilty in the pertubated energy model */
  p_prob_unpaired = space(sizeof(double) * (length + 1));

  addSoftConstraint(vc, epsilon, length);

  vc->exp_params->model_details.compute_bpp = 1;

  /* get new (constrained) MFE to scale pf computations properly */
  double mfe = (double)vrna_fold(vc, NULL);
  vrna_rescale_pf_params(vc, &mfe);

  vrna_pf_fold(vc, NULL);

  calculate_probability_unpaired(vc, p_prob_unpaired);

  vrna_sc_remove(vc);

  
  for (i = 1; i <= length; ++i)
  {
    /* add penalty for pertubation energies */
    ret += evaluate_objective_function_contribution(epsilon[i], objective_function) / tau_squared;

    /* add penalty for mismatches between observed and predicted probabilities */
    if (q_prob_unpaired[i] >= 0) /* ignore positions with missing data */
      ret2 += evaluate_objective_function_contribution(p_prob_unpaired[i] - q_prob_unpaired[i], objective_function) / sigma_squared;
  }

  fprintf(stderr, "Score: pertubation: %g\tdiscrepancy: %g\n", ret, ret2);
  free(p_prob_unpaired);

  return ret + ret2;
}

static void pairing_probabilities_from_restricted_pf(vrna_fold_compound *vc, const double *epsilon, double *prob_unpaired, double **conditional_prob_unpaired)
{
  int length = vc->length;
  int i;

  addSoftConstraint(vc, epsilon, length);
  vc->exp_params->model_details.compute_bpp = 1;

  /* get new (constrained) MFE to scale pf computations properly */
  double mfe = (double)vrna_fold(vc, NULL);
  vrna_rescale_pf_params(vc, &mfe);

  vrna_pf_fold(vc, NULL);

  calculate_probability_unpaired(vc, prob_unpaired);

  #pragma omp parallel for private(i)
  for (i = 1; i <= length; ++i)
  {
    vrna_fold_compound *restricted_vc;
    char *hc_string;
    unsigned int constraint_options = VRNA_CONSTRAINT_DB
                                      | VRNA_CONSTRAINT_PIPE
                                      | VRNA_CONSTRAINT_DOT
                                      | VRNA_CONSTRAINT_X
                                      | VRNA_CONSTRAINT_ANG_BRACK
                                      | VRNA_CONSTRAINT_RND_BRACK;

    hc_string = space(sizeof(char) * (length + 1));
    memset(hc_string, '.', length);
    hc_string[i - 1] = 'x';

    restricted_vc = vrna_get_fold_compound(vc->sequence, &(vc->exp_params->model_details), VRNA_OPTION_PF);
    vrna_hc_add(restricted_vc, hc_string, constraint_options);
    free(hc_string);

    vrna_update_pf_params(restricted_vc, vc->exp_params);

    vrna_pf_fold(restricted_vc, NULL);
    calculate_probability_unpaired(restricted_vc, conditional_prob_unpaired[i]);

    restricted_vc->sc = NULL;
    vrna_free_fold_compound(restricted_vc);
  }

  vrna_sc_remove(vc);
}

static void pairing_probabilities_from_sampling(vrna_fold_compound *vc, const double *epsilon, int sample_size, double *prob_unpaired, double **conditional_prob_unpaired)
{
  double kT;
  int length = vc->length;
  int i, j, s;
  st_back = 1; /* is this really required? */

  addSoftConstraint(vc, epsilon, length);

  vc->exp_params->model_details.compute_bpp = 0;

  /* get new (constrained) MFE to scale pf computations properly */
  double mfe = (double)vrna_fold(vc, NULL);
  vrna_rescale_pf_params(vc, &mfe);

  vrna_pf_fold(vc, NULL);


  #pragma omp parallel for private(s)
  for (s = 0; s < sample_size; ++s)
  {
    char *sample = vrna_pbacktrack(vc);

    #pragma omp critical
    {
      for (i = 1; i <= length; ++i)
      {
        if (sample[i-1] != '.')
          continue;

        ++prob_unpaired[i];

        for (j = 1; j <= length; ++j)
          if (sample[j-1] == '.')
            ++conditional_prob_unpaired[i][j];
      }
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

  vrna_sc_remove(vc);
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

static void evaluate_perturbation_vector_gradient(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int objective_function, int sample_size, double *gradient)
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

    if (objective_function == VRNA_OBJECTIVE_FUNCTION_QUADRATIC)
    {
      for (i = 1; i <= length; ++i)
      {
        if (q_prob_unpaired[i] < 0) /* ignore positions with missing data */
          continue;

        sum += (p_prob_unpaired[i] - q_prob_unpaired[i])
               * p_prob_unpaired[i] * (p_prob_unpaired[mu] - p_conditional_prob_unpaired[i][mu])
               / sigma_squared;
      }

      gradient[mu] = 2 * (epsilon[mu] / tau_squared + sum/kT);
    }
    else if (objective_function == VRNA_OBJECTIVE_FUNCTION_ABSOLUTE)
    {
      for (i = 1; i <= length; ++i)
        if (q_prob_unpaired[i] >= 0 && p_prob_unpaired[i] != q_prob_unpaired[i])
          sum += (p_prob_unpaired[i] * (p_prob_unpaired[mu] - p_conditional_prob_unpaired[i][mu])) / kT
                 / sigma_squared
                 * (p_prob_unpaired[i] > q_prob_unpaired[i] ? 1. : -1.);

      if (epsilon[mu])
        sum += (epsilon[mu] > 0 ? 1. : -1.) / tau_squared;

      gradient[mu] = sum;
    }
  }

  freeProbabilityArrays(p_prob_unpaired, p_conditional_prob_unpaired, length);
}

#ifdef WITH_GSL
typedef struct parameters_gsl {
  vrna_fold_compound *vc;
  const double *q_prob_unpaired;
  double sigma_squared;
  double tau_squared;
  int objective_function;
  int sample_size;
} parameters_gsl;

static double f_gsl(const gsl_vector *x, void *params)
{
  parameters_gsl *p = params;

  return evaluate_perturbation_vector_score(p->vc, x->data, p->q_prob_unpaired, p->sigma_squared, p->tau_squared, p->objective_function);
}

static void df_gsl(const gsl_vector *x, void *params, gsl_vector *df)
{
  parameters_gsl *p = params;

  gsl_vector_set(df, 0, 0);
  evaluate_perturbation_vector_gradient(p->vc, x->data, p->q_prob_unpaired, p->sigma_squared, p->tau_squared, p->objective_function, p->sample_size, df->data);
}

static void fdf_gsl(const gsl_vector *x, void *params, double *f, gsl_vector *g)
{
  *f = f_gsl(x, params);
  df_gsl(x, params, g);
}
#endif /* WITH_GSL */

PUBLIC void
vrna_find_perturbation_vector(vrna_fold_compound *vc,
                              const double *q_prob_unpaired,
                              int objective_function,
                              double sigma_squared,
                              double tau_squared,
                              int algorithm,
                              int sample_size,
                              double *epsilon,
                              double initialStepSize,
                              double minStepSize,
                              double minImprovement,
                              double minimizerTolerance,
                              progress_callback callback){

  int iteration_count = 0;
  const int max_iterations = 100;
  int length = vc->length;

#ifdef WITH_GSL
  const gsl_multimin_fdfminimizer_type *minimizer_type = 0;

  struct {int type; const gsl_multimin_fdfminimizer_type *gsl_type;} algorithms[] = {{VRNA_MINIMIZER_CONJUGATE_FR, gsl_multimin_fdfminimizer_conjugate_fr},
                                                                                     {VRNA_MINIMIZER_CONJUGATE_PR, gsl_multimin_fdfminimizer_conjugate_pr},
                                                                                     {VRNA_MINIMIZER_VECTOR_BFGS, gsl_multimin_fdfminimizer_vector_bfgs},
                                                                                     {VRNA_MINIMIZER_VECTOR_BFGS2, gsl_multimin_fdfminimizer_vector_bfgs2},
                                                                                     {VRNA_MINIMIZER_STEEPEST_DESCENT, gsl_multimin_fdfminimizer_steepest_descent},
                                                                                     {0, NULL}};
  int i;
  for (i = 0; algorithms[i].type; ++i)
    if (algorithms[i].type == algorithm)
    {
      minimizer_type = algorithms[i].gsl_type;
      break;
    }

  if (minimizer_type)
  {
    parameters_gsl parameters;
    gsl_multimin_function_fdf fdf;
    gsl_multimin_fdfminimizer *minimizer;
    gsl_vector *vector;

    int status;

    parameters.vc = vc;
    parameters.q_prob_unpaired = q_prob_unpaired;
    parameters.sigma_squared = sigma_squared;
    parameters.tau_squared = tau_squared;
    parameters.objective_function = objective_function;
    parameters.sample_size = sample_size;

    fdf.n = length + 1;
    fdf.f = &f_gsl;
    fdf.df = &df_gsl;
    fdf.fdf = &fdf_gsl;
    fdf.params = (void*)&parameters;

    minimizer = gsl_multimin_fdfminimizer_alloc(minimizer_type, length + 1);
    vector = gsl_vector_calloc(length + 1);

    /* gsl_multimin_fdfminimizer_set(minimizer, &fdf, vector, 0.01, 1e-4); */
    gsl_multimin_fdfminimizer_set(minimizer, &fdf, vector, initialStepSize, minimizerTolerance);

    if (callback)
      callback(0, minimizer->f, minimizer->x->data);

    do
    {
      ++iteration_count;
      status = gsl_multimin_fdfminimizer_iterate(minimizer);

      if (callback)
        callback(iteration_count, minimizer->f, minimizer->x->data);

      if (status)
        break;

      status = gsl_multimin_test_gradient(minimizer->gradient, minimizerTolerance);
    }
    while (status == GSL_CONTINUE && iteration_count < max_iterations);

    memcpy(epsilon, minimizer->x->data, sizeof(double) * (length + 1));

    gsl_multimin_fdfminimizer_free(minimizer);
    gsl_vector_free(vector);

    return;
  }
#endif /* WITH_GSL */

  double improvement;
  const double min_improvement = minImprovement;

  double *new_epsilon = space(sizeof(double) * (length + 1));
  double *gradient = space(sizeof(double) * (length + 1));

  double score = evaluate_perturbation_vector_score(vc, epsilon, q_prob_unpaired, sigma_squared, tau_squared, objective_function);

  if (callback)
    callback(0, score, epsilon);

  do
  {
    double new_score;
    double step_size;

    ++iteration_count;

    evaluate_perturbation_vector_gradient(vc, epsilon, q_prob_unpaired, sigma_squared, tau_squared, objective_function, sample_size, gradient);

    /*    step_size = 0.5 / calculate_norm(gradient, length);*/
    step_size = initialStepSize;

    do
    {
      int i;
      for (i = 1; i <= length; ++i)
        new_epsilon[i] = epsilon[i] - step_size * gradient[i];

      new_score = evaluate_perturbation_vector_score(vc, new_epsilon, q_prob_unpaired, sigma_squared, tau_squared, objective_function);
      improvement = 1 - new_score / score;
      step_size /= 2;
    } while ((improvement < min_improvement) && (step_size >= minStepSize));

    if (new_score > score)
      break;

    if (callback)
      callback(iteration_count, new_score, new_epsilon);

    score = new_score;
    memcpy(epsilon, new_epsilon, sizeof(double) * (length+1));
  } while (improvement >= min_improvement && iteration_count < max_iterations);

  free(gradient);
  free(new_epsilon);
}

