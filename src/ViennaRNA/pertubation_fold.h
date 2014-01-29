#ifndef __VIENNA_RNA_PACKAGE_PERTUBATION_FOLD_H__
#define __VIENNA_RNA_PACKAGE_PERTUBATION_FOLD_H__

#include "data_structures.h"

double vrna_evaluate_pertubation_vector_score(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared);
void vrna_evaluate_pertubation_vector_gradient(vrna_fold_compound *vc, const double *epsilon, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int sample_size, double *gradient);

typedef void (*progress_callback)(int iteration, double score);
void vrna_find_pertubation_vector(vrna_fold_compound *vc, const double *q_prob_unpaired, double sigma_squared, double tau_squared, int sample_size, double *epsilon, progress_callback callback);

#endif
