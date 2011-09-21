#ifndef __VIENNA_RNA_PACKAGE_SUBOPT_H__
#define __VIENNA_RNA_PACKAGE_SUBOPT_H__

typedef struct svm_model{
  struct svm_parameter param;
  int nr_class;
  int l;
  struct svm_node **SV;
  double **sv_coef;
  double *rho;
  double *probA;
  double *probB;
  int *label;
  int *nSV;
  int free_sv;
} svm_model;

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
svm_model *svm_load_model_string(char *modelString);
int       *get_seq_composition( short *S,
                                unsigned int start,
                                unsigned int stop);

#endif
