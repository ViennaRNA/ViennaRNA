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
void    init_pf_fold(int length);    /* allocate space for pf_fold() */
void    free_pf_arrays(void);        /* free arrays from pf_fold() */
void    update_pf_params(int length); /*recalculate energy parameters */
char    bppm_symbol(float *x);  /* string representation of structure */
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
char    *get_centroid_struct_pr(int length, double *dist, double *pl);

/**
*** Create a plist from a probability matrix
*** The probability matrix given is parsed and all pair probabilities above
*** the given threshold are used to create an entry in the plist
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
