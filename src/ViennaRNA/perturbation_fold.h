#ifndef __VIENNA_RNA_PACKAGE_PERTURBATION_FOLD_H__
#define __VIENNA_RNA_PACKAGE_PERTURBATION_FOLD_H__

#include "data_structures.h"

double vrna_evaluate_perturbation_vector_score(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared);
void vrna_evaluate_perturbation_vector_gradient(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int sample_size, double *gradient);

typedef void (*progress_callback)(int iteration, double score);
void vrna_find_perturbation_vector(vrna_fold_compound *vc, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int sample_size, double *epsilon, progress_callback callback);

#endif
