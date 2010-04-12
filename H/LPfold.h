#ifndef __VIENNA_RNA_PACKAGE_LPFOLD_H__
#define __VIENNA_RNA_PACKAGE_LPFOLD_H__

void    update_pf_paramsLP(int length);
struct  plist *pfl_fold(char *sequence, int winSize, int pairSize, float cutoffb, double **pU, struct plist **dpp2, FILE *pUfp, FILE *spup);
void    putoutpU_prob(double **pU,int length, int ulength, FILE *fp, int energies);


#endif
