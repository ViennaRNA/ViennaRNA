#ifndef __VIENNA_RNA_PACKAGE_PART_FUNC_H__
#define __VIENNA_RNA_PACKAGE_PART_FUNC_H__

#include "data_structures.h"

#define FLT_OR_DBL double

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif


/**
*** \file part_func.h
**/

/* functions from part_func.c */
float   pf_fold(const char *sequence, char *structure);
float   pf_circ_fold(const char *sequence, char *structure);
/* calculate partition function and base pair probabilities */

/** allocate space for pf_fold()
*** \deprecated {This function is obsolete and will be removed soon!}
**/
void    DEPRECATED(init_pf_fold(int length));

void    free_pf_arrays(void);        /* free arrays from pf_fold() */
void    update_pf_params(int length); /*recalculate energy parameters */
char    bppm_symbol(float *x);
double  mean_bp_dist(int length); /* mean pair distance of ensemble */
int     get_pf_arrays(short **S_p, short **S1_p, char **ptype_p, FLT_OR_DBL **qb_p, FLT_OR_DBL **qm_p, FLT_OR_DBL **q1k_p, FLT_OR_DBL **qln_p);

/**
*** <H2>Sample a secondary structure from the Boltzmann ensemble according its probability</H2>
*** \param  sequence  The RNA sequence
*** \return           A sampled secondary structure in dot-bracket notation
**/
char    *pbacktrack(char *sequence);

/**
*** <H2>Sample a secondary structure of a circular RNA from the Boltzmann ensemble according its probability</H2>
*** This function does the same as \func pbacktrack() but assumes the RNA molecule to be circular
*** \param  sequence  The RNA sequence
*** \return           A sampled secondary structure in dot-bracket notation
**/
char    *pbacktrack_circ(char *seq);

/*
#################################################
# OTHER PARTITION FUNCTION RELATED DECLARATIONS #
#################################################
*/
char    DEPRECATED(*centroid(int length, double *dist));     /* mean pair distance of ensemble */
/* this function is a threadsafe replacement for centroid() with a 'struct plist' input */
char    *get_centroid_struct_pl(int length, double *dist, plist *pl);
/* this function is a threadsafe replacement for centroid() with a probability array input */
char    *get_centroid_struct_pr(int length, double *dist, double *pr);

/**
*** Create a dot-bracket like structure string from base pair probability matrix
**/
void    bppm_to_structure(char *structure, FLT_OR_DBL *pr, unsigned int length);

/**
*** Create a plist from a probability matrix
*** The probability matrix given is parsed and all pair probabilities above
*** the given threshold are used to create an entry in the plist
***
*** The end of the plist is marked by sequence positions i as well as j
*** equal to 0. This condition should be used to stop looping over its
*** entries
***
*** \note This function is threadsafe
***
*** \param pl     A pointer to the plist that is to be created
*** \param probs  The probability matrix used for creting the plist
*** \param length The length of the RNA sequence
*** \param cutoff The cutoff value
**/
void    assign_plist_from_pr(plist **pl, double *probs, int length, double cut_off);

extern  int st_back;

plist   *stackProb(double cutoff);

/* deprecated, use exp_E_IntLoop() from loop_energies.h instead */
double  DEPRECATED(expLoopEnergy(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1));
/* deprecated, use exp_E_Hairpin() from loop_energies.h instead */
double  DEPRECATED(expHairpinEnergy(int u, int type, short si1, short sj1, const char *string));

#endif
