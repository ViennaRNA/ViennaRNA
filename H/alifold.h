#ifndef __VIENNA_RNA_PACKAGE_ALIFOLD_H__
#define __VIENNA_RNA_PACKAGE_ALIFOLD_H__

#include "data_structures.h"

extern  double  cv_fact /* =1 */;
extern  double  nc_fact /* =1 */;

float   alifold(char **strings, char *structure);
void    free_alifold_arrays(void);


/* extern float aliLfold(char **strings, char *structure, int maxdist); */
/* extern float alipfW_fold(char **sequences, char *structure, struct plist **pl,int winSize, float cutoff, int pairsize); */
/* extern struct plist *get_plistW(struct plist *pl, int length, double cutoff, int start, FLT_OR_DBL **Tpr, int winSize); */
float   **readribosum(char *name);
void    energy_of_alistruct(char **sequences, const char *structure, int n_seq, float *energy);
double  circalifold(const char **strings, char *structure);

/* partition function variants */
float   alipf_fold(char **sequences, char *structure, struct plist **pl);
void    free_alipf_arrays(void);
char    *centroid_ali(int length, double *dist,struct plist *pl);
char    *alipbacktrack(double *prob) ;
float   alipf_circ_fold(char **sequences, char *structure, struct plist **pl);

#endif
