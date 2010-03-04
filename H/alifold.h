#ifndef __VIENNA_RNA_PACKAGE_ALIFOLD_H__
#define __VIENNA_RNA_PACKAGE_ALIFOLD_H__

#include "data_structures.h"

/** \file alifold.h
***
**/

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

/* some helper functions that might be useful in the library */

/** get arrays with encoded sequence of the alignment
***
*** this function assumes that in S, S5, s3, ss and as enough
*** space is already allocated (size must be at least sequence length+2)
*** \param sequence The gapped sequence from the alignment
*** \param S        pointer to an array that holds encoded sequence
*** \params s5      pointer to an array that holds the next base 5' of alignment position i
*** \params s3      pointer to an array that holds the next base 3' of alignment position i
*** \params ss      
*** \params as      
**/
void encode_ali_sequence(const char *sequence, short *S, short *s5, short *s3, char *ss, unsigned short *as);

#endif
