#ifndef __VIENNA_RNA_PACKAGE_PERTURBATION_FOLD_H__
#define __VIENNA_RNA_PACKAGE_PERTURBATION_FOLD_H__

#include "data_structures.h"

#define VRNA_OBJECTIVE_FUNCTION_QUADRATIC 0
#define VRNA_OBJECTIVE_FUNCTION_ABSOLUTE 1

#define VRNA_MINIMIZER_DEFAULT 0
#define VRNA_MINIMIZER_CONJUGATE_FR 1
#define VRNA_MINIMIZER_CONJUGATE_PR 2
#define VRNA_MINIMIZER_VECTOR_BFGS 3
#define VRNA_MINIMIZER_VECTOR_BFGS2 4
#define VRNA_MINIMIZER_STEEPEST_DESCENT 5

double vrna_evaluate_perturbation_vector_score(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int objective_function);
void vrna_evaluate_perturbation_vector_gradient(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int objective_function, int sample_size, double *gradient);

typedef void (*progress_callback)(int iteration, double score, double *epsilon);
void vrna_find_perturbation_vector(vrna_fold_compound *vc, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int objective_function, int algorithm, int sample_size, double *epsilon, progress_callback callback);

#endif
