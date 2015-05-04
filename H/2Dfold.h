/*
      minimum free energy
      RNA secondary structure with
      basepair distance d to reference structure prediction

*/
#ifndef __VIENNA_RNA_PACKAGE_TWO_D_FOLD_H__
#define __VIENNA_RNA_PACKAGE_TWO_D_FOLD_H__

/**
 *  \addtogroup kl_neighborhood
 *  \brief Compute Thermodynamic properties for a Distance Class Partitioning of the Secondary Structure Space
 *
 *  All functions related to this group implement the basic recursions for MFE folding, partition function
 *  computation and stochastic backtracking with a \e classified \e dynamic \e programming approach.
 *  The secondary structure space is divided into partitions according to the base pair distance to two
 *  given reference structures and all relevant properties are calculated for each of the resulting partitions
 *  \see  For further details have a look into \cite lorenz:2009
 */

/**
 *  \addtogroup kl_neighborhood_mfe
 *  \brief Compute the minimum free energy (MFE) and secondary structures for a partitioning of
 *  the secondary structure space according to the base pair distance to two fixed reference structures
 *  basepair distance to two fixed reference structures
 *  @{
 *
 *  \file 2Dfold.h
 *
 */

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  \brief Get a structure of type TwoDfold_vars prefilled with current global settings
 * 
 *  This function returns a datastructure of type TwoDfold_vars.
 *  The data fields inside the TwoDfold_vars are prefilled by global settings and all memory
 *  allocations necessary to start a computation are already done for the convenience of the user
 * 
 *  \note Make sure that the reference structures are compatible with the sequence according to Watson-Crick- and Wobble-base pairing
 * 
 *  \see destroy_TwoDfold_variables(), TwoDfold(), TwoDfold_circ
 * 
 *  \param seq          The RNA sequence
 *  \param structure1   The first reference structure in dot-bracket notation
 *  \param structure2   The second reference structure in dot-bracket notation
 *  \param circ         A switch to indicate the assumption to fold a circular instead of linear RNA (0=OFF, 1=ON)
 *  \returns            A datastructure prefilled with folding options and allocated memory
 */
TwoDfold_vars *get_TwoDfold_variables(const char *seq,
                                      const char *structure1,
                                      const char *structure2,
                                      int circ);

/**
 *  \brief Destroy a TwoDfold_vars datastructure without memory loss
 * 
 *  This function free's all allocated memory that depends on the datastructure given.
 * 
 *  \see get_TwoDfold_variables()
 * 
 *  \param our_variables  A pointer to the datastructure to be destroyed
 */
void          destroy_TwoDfold_variables(TwoDfold_vars *our_variables);

/**
 * 
 */
DEPRECATED(TwoDfold_solution **TwoDfold(TwoDfold_vars *our_variables,
                                        int distance1,
                                        int distance2));

/**
 * \brief Compute MFE's and representative for distance partitioning
 *
 * This function computes the minimum free energies and a representative
 * secondary structure for each distance class according to the two references
 * specified in the datastructure 'vars'.
 * The maximum basepair distance to each of both references may be set
 * by the arguments 'distance1' and 'distance2', respectively.
 * If both distance arguments are set to '-1', no restriction is assumed and
 * the calculation is performed for each distance class possible.
 *
 * The returned list contains an entry for each distance class. If a maximum
 * basepair distance to either of the references was passed, an entry with
 * k=l=-1 will be appended in the list, denoting the class where all structures
 * exceeding the maximum will be thrown into.
 * The end of the list is denoted by an attribute value of #INF in
 * the k-attribute of the list entry.
 *
 * \see get_TwoDfold_variables(), destroy_TwoDfold_variables(), #TwoDfold_solution
 *
 * \param vars      the datastructure containing all predefined folding attributes
 * \param distance1 maximum distance to reference1 (-1 means no restriction)
 * \param distance2 maximum distance to reference2 (-1 means no restriction)
 */
TwoDfold_solution *TwoDfoldList(TwoDfold_vars *vars,
                                int distance1,
                                int distance2);

/**
 * \brief Backtrack a minimum free energy structure from a 5' section of specified length
 *
 * This function allows to backtrack a secondary structure beginning at the 5' end, a specified
 * length and residing in a specific distance class.
 * If the argument 'k' gets a value of -1, the structure that is backtracked is assumed to
 * reside in the distance class where all structures exceeding the maximum basepair distance
 * specified in TwoDfoldList() belong to.
 * \note The argument 'vars' must contain precalculated energy values in the energy matrices,
 * i.e. a call to TwoDfoldList() preceding this function is mandatory!
 *
 * \see TwoDfoldList(), get_TwoDfold_variables(), destroy_TwoDfold_variables()
 *
 * \param j     The length in nucleotides beginning from the 5' end
 * \param k     distance to reference1 (may be -1)
 * \param l     distance to reference2
 * \param vars  the datastructure containing all predefined folding attributes
 */
char *TwoDfold_backtrack_f5(unsigned int j,
                            int k,
                            int l,
                            TwoDfold_vars *vars);

/**
 *  @}
 */
#endif
