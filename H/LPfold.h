#ifndef LPFOLD_H
#define LPFOLD_H

extern void  update_pf_paramsLP(int length);
extern struct plist *pfl_fold(char *sequence, int winSize, int pairSize,
			      float cutoff, double *pup);

#endif
