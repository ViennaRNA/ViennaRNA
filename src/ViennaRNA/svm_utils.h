#ifndef VIENNA_RNA_PACKAGE_SUBOPT_H
#define VIENNA_RNA_PACKAGE_SUBOPT_H

#include "svm.h"

extern  char *avg_model_string;
extern  char *sd_model_string;

float     get_z(char *sequence,
                double energy);
double    avg_regression (int N,
                          int A,
                          int C,
                          int G,
                          int T,
                          struct svm_model *avg_model,
                          int *info );
double    sd_regression  (int N,
                          int A,
                          int C,
                          int G,
                          int T,
                          struct svm_model  *sd_model);
double    minimal_sd     (int N,
                          int A,
                          int C,
                          int G,
                          int T);
struct svm_model *svm_load_model_string(char *modelString);
int       *get_seq_composition( short *S,
                                unsigned int start,
                                unsigned int stop,
                                unsigned int length);

#endif
