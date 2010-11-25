#ifndef __VIENNA_RNA_PACKAGE_LPFOLD_H__
#define __VIENNA_RNA_PACKAGE_LPFOLD_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
*** \file LPfold.h
*** \brief
***
**/

/**
*** \brief
***
*** \param  length
**/
void    update_pf_paramsLP(int length);

/**
*** \brief Compute partition functions for locally stable secondary structures <b>(berni! update me)</b>
***
*** pfl_fold computes partition functions for every window of size
*** winSize possible in a RNA molecule, allowing only pairs with a span
*** smaller than pairSize. It returns the mean pair probabilities averaged
*** over all windows containing the pair in 'pl'. 'winSize' should
*** always be >= 'pairSize'. Note that in contrast to Lfold(),
*** bases outside of the window do not influence the structure at all. Only
*** probabilities higher than 'cutoff' are kept.
***
*** If 'pup' is supplied (i.e is not the NULL pointer), pfl_fold()
*** will also compute the mean probability that regions of length 'u' are
*** unpaired. The parameter 'u' is supplied in 'pup[0]'. On return
*** the 'pup' array will contain these probabilities, with the entry on
*** 'pup[x]' containing the mean probability that x and the u-1
*** preceding bases are unpaired. The 'pup' array needs to be large
*** enough to hold n+1 float entries, where n is the sequence length.
***
*** \param  sequence  
*** \param  winSize   
*** \param  pairSize  
*** \param  cutoffb   
*** \param  pU        
*** \param  dpp2      
*** \param  pUfp      
*** \param  spup      
*** \return           
**/
plist *pfl_fold(char *sequence, int winSize, int pairSize, float cutoffb, double **pU, struct plist **dpp2, FILE *pUfp, FILE *spup);

/**
*** \brief
***
*** \param  pU        
*** \param  length    
*** \param  ulength   
*** \param  fp        
*** \param  energies  
**/
void    putoutpU_prob(double **pU,int length, int ulength, FILE *fp, int energies);

/**
*** Dunno if this function was ever used by external programs linking to RNAlib, but it
*** was declared PUBLIC before.
*** Anyway, never use this function as it will be removed soon and does nothing at all
**/
DEPRECATED(void init_pf_foldLP(int length));

#endif
