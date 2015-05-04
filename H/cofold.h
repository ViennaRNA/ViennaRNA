#ifndef __VIENNA_RNA_PACKAGE_COFOLD_H__
#define __VIENNA_RNA_PACKAGE_COFOLD_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  \file cofold.h
 * 
 *  \brief MFE version of cofolding routines
 * 
 *  This file includes (almost) all function declarations within the <b>RNAlib</b> that are related to
 *  MFE Cofolding...
 *  This also includes the Zuker suboptimals calculations, since they are implemented using the cofold
 *  routines.
 */

/**
 *  \brief Compute the minimum free energy of two interacting RNA molecules
 * 
 *  The code is analog to the fold() function. If #cut_point ==-1 results
 *  should be the same as with fold().
 * 
 *  \param    sequence  The two sequences concatenated
 *  \param    structure Will hold the barcket dot structure of the dimer molecule
 *  \return   minimum free energy of the structure
 */
float cofold( const char *sequence,
              char *structure);

float cofold_par( const char *string,
                  char *structure,
                  paramT *parameters,
                  int is_constrained);

/**
 *  \brief Free memory occupied by cofold()
 */
void      free_co_arrays(void);

/**
 *  \brief Recalculate parameters
 */
void      update_cofold_params(void);

void      update_cofold_params_par(paramT *parameters);

/**
 *  \brief Compute Zuker type suboptimal structures
 *
 *  Compute Suboptimal structures according to M. Zuker, i.e. for every 
 *  possible base pair the minimum energy structure containing the resp. base pair. 
 *  Returns a list of these structures and their energies.
 *
 *  \param  string  RNA sequence
 *  \return         List of zuker suboptimal structures
 */
SOLUTION  *zukersubopt(const char *string);

SOLUTION  *zukersubopt_par( const char *string,
                            paramT *parameters);

/**
 *  \brief get_monomer_free_energies
 *
 *  Export monomer free energies out of cofold arrays
 * 
 *  \param e1 A pointer to a variable where the energy of molecule A will be written to
 *  \param e2 A pointer to a variable where the energy of molecule B will be written to
 */
void get_monomere_mfes( float *e1,
                        float *e2);

/**
 *  \brief Export the arrays of partition function cofold
 * 
 *  Export the cofold arrays for use e.g. in the concentration
 *  Computations or suboptimal secondary structure backtracking
 *
 *  \param  f5_p    A pointer to the 'f5' array, i.e. array conatining best free energy in interval [1,j]
 *  \param  c_p     A pointer to the 'c' array, i.e. array containing best free energy in interval [i,j] given that i pairs with j
 *  \param  fML_p   A pointer to the 'M' array, i.e. array containing best free energy in interval [i,j] for any multiloop segment with at least one stem
 *  \param  fM1_p   A pointer to the 'M1' array, i.e. array containing best free energy in interval [i,j] for multiloop segment with exactly one stem
 *  \param  fc_p    A pointer to the 'fc' array, i.e. array ...
 *  \param  indx_p  A pointer to the indexing array used for accessing the energy matrices
 *  \param  ptype_p A pointer to the ptype array containing the base pair types for each possibility (i,j)
 */
void export_cofold_arrays(int **f5_p,
                          int **c_p,
                          int **fML_p,
                          int **fM1_p,
                          int **fc_p,
                          int **indx_p,
                          char **ptype_p);

/**
 *  allocate arrays for folding
 *  \deprecated{This function is obsolete and will be removed soon!}
 */
DEPRECATED(void initialize_cofold(int length));

#endif
