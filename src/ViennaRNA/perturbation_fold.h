#ifndef __VIENNA_RNA_PACKAGE_PERTURBATION_FOLD_H__
#define __VIENNA_RNA_PACKAGE_PERTURBATION_FOLD_H__

#include "data_structures.h"

#define VRNA_OBJECTIVE_FUNCTION_QUADRATIC 0
#define VRNA_OBJECTIVE_FUNCTION_ABSOLUTE 1

double vrna_evaluate_perturbation_vector_score(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int method);
void vrna_evaluate_perturbation_vector_gradient(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int method, int sample_size, double *gradient);

typedef void (*progress_callback)(int iteration, double score, double *epsilon);
void vrna_find_perturbation_vector(vrna_fold_compound *vc, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int method, int sample_size, double *epsilon, progress_callback callback);

#endif
