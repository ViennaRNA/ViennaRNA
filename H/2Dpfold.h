/*
      minimum free energy
      RNA secondary structure with
      basepair distance d to reference structure prediction

*/
#ifndef __VIENNA_RNA_PACKAGE_TWO_D_PF_FOLD_H__
#define __VIENNA_RNA_PACKAGE_TWO_D_PF_FOLD_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  \file 2Dpfold.h
 *  \brief Compute the partition function and stochastically sample secondary structures for a partitioning of
 *  the secondary structure space according to the base pair distance to two fixed reference structures
 */

/**
 * \brief Get a datastructure containing all necessary attributes and global folding switches
 *
 * This function prepares all necessary attributes and matrices etc which are needed for a call
 * of TwoDpfoldList.
 * A snapshot of all current global model switches (dangles, temperature and so on) is done and
 * stored in the returned datastructure. Additionally, all matrices that will hold the partition
 * function values are prepared.
 *
 * \param seq         the RNA sequence in uppercase format with letters from the alphabet {AUCG}
 * \param structure1  the first reference structure in dot-bracket notation
 * \param structure2  the second reference structure in dot-bracket notation
 * \param circ        a switch indicating if the sequence is linear (0) or circular (1)
 * \returns           the datastructure containing all necessary partition function attributes
 */
TwoDpfold_vars  *get_TwoDpfold_variables( const char *seq,
                                          const char *structure1,
                                          char *structure2,
                                          int circ);

/**
 * \brief Get the datastructure containing all necessary attributes and global folding switches from 
 * a pre-filled mfe-datastructure
 *
 * This function actually does the same as get_TwoDpfold_variables but takes its switches and
 * settings from a pre-filled MFE equivalent datastructure
 *
 * \see get_TwoDfold_variables(), get_TwoDpfold_variables()
 *
 * \param mfe_vars    the pre-filled mfe datastructure
 * \returns           the datastructure containing all necessary partition function attributes
 */
TwoDpfold_vars  *get_TwoDpfold_variables_from_MFE(TwoDfold_vars *mfe_vars);

/**
 * \brief Free all memory occupied by a TwoDpfold_vars datastructure
 *
 * This function free's all memory occupied by a datastructure obtained from from
 * get_TwoDpfold_variables() or get_TwoDpfold_variables_from_MFE()
 *
 * \see get_TwoDpfold_variables(), get_TwoDpfold_variables_from_MFE()
 *
 * \param vars   the datastructure to be free'd
 */
void            destroy_TwoDpfold_variables(TwoDpfold_vars *vars);

/**
 * \brief
 *
 *
 */
DEPRECATED(FLT_OR_DBL          **TwoDpfold(TwoDpfold_vars *our_variables,
                                int maxDistance1,
                                int maxDistance2));

/**
 * \brief
 *
 *
 */
DEPRECATED(FLT_OR_DBL          **TwoDpfold_circ(
                                TwoDpfold_vars *our_variables,
                                int maxDistance1,
                                int maxDistance2));

/**
 * \brief Compute the partition function for all distance classes
 *
 * This function computes the partition functions for all distance classes
 * according the two reference structures specified in the datastructure 'vars'.
 * Similar to TwoDfoldList() the arguments maxDistance1 and maxDistance2 specify
 * the maximum distance to both reference structures. A value of '-1' in either of
 * them makes the appropriate distance restrictionless, i.e. all basepair distancies
 * to the reference are taken into account during computation.
 * In case there is a restriction, the returned solution contains an entry where
 * the attribute k=l=-1 contains the partition function for all structures exceeding
 * the restriction.
 * A values of #INF in the attribute 'k' of the returned list denotes the end of the list
 *
 * \see get_TwoDpfold_variables(), destroy_TwoDpfold_variables(), #TwoDpfold_solution
 *
 * \param vars          the datastructure containing all necessary folding attributes and matrices
 * \param maxDistance1  the maximum basepair distance to reference1 (may be -1)
 * \param maxDistance2  the maximum basepair distance to reference2 (may be -1)
 * \returns             a list of partition funtions for the appropriate distance classes
 */
TwoDpfold_solution  *TwoDpfoldList( TwoDpfold_vars *vars,
                                    int maxDistance1,
                                    int maxDistance2);

/**
 *  \brief Sample secondary structure representatives from a set of distance classes according to their 
 *  Boltzmann probability
 *
 * If the argument 'd1' is set to '-1', the structure will be backtracked in the distance class
 * where all structures exceeding the maximum basepair distance to either of the references reside.
 *
 * \note The argument 'vars' must contain precalculated partition function matrices,
 * i.e. a call to TwoDpfoldList() preceding this function is mandatory!
 *
 * \see TwoDpfoldList()
 *
 *  \param vars   the datastructure containing all necessary folding attributes and matrices
 *  \param d1     the distance to reference1 (may be -1)
 *  \param d2     the distance to reference2
 *  \returns      a sampled secondary structure in dot-bracket notation
 */
char            *TwoDpfold_pbacktrack(TwoDpfold_vars *vars,
                                      int d1,
                                      int d2);

/**
 * \brief Sample secondary structure representatives with a specified length from a set of distance classes according to their 
 *  Boltzmann probability
 *
 * This function does essentially the same as TwoDpfold_pbacktrack with the only difference that partial structures,
 * i.e. structures beginning from the 5' end with a specified length of the sequence, are backtracked
 *
 * \note The argument 'vars' must contain precalculated partition function matrices,
 * i.e. a call to TwoDpfoldList() preceding this function is mandatory!
 * \note This function does not work (since it makes no sense) for circular RNA sequences!
 *
 * \see TwoDpfold_pbacktrack(), TwoDpfoldList()
 *
 *  \param vars   the datastructure containing all necessary folding attributes and matrices
 *  \param d1     the distance to reference1 (may be -1)
 *  \param d2     the distance to reference2
 *  \param length the length of the structure beginning from the 5' end
 *  \returns      a sampled secondary structure in dot-bracket notation
 */
char            *TwoDpfold_pbacktrack5( TwoDpfold_vars *vars,
                                          int d1,
                                          int d2,
                                          unsigned int length);

#endif
