/* subopt_zuker.h */
#ifndef VIENNA_RNA_PACKAGE_SUBOPT_ZUKER_H
#define VIENNA_RNA_PACKAGE_SUBOPT_ZUKER_H

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/subopt.h>

/**
 *  @brief Compute Zuker type suboptimal structures
 *
 *  Compute Suboptimal structures according to M. Zuker @cite zuker:1989 , i.e. for every
 *  possible base pair the minimum energy structure containing the resp. base pair.
 *  Returns a list of these structures and their energies.
 *
 *  @note This function internally uses the cofold implementation to compute
 *        the suboptimal structures. For that purpose, the function doubles
 *        the sequence and enlarges the DP matrices, which in fact will grow
 *        by a factor of 4 during the computation!
 *        At the end of the structure prediction, everything will be re-set
 *        to its original requriements, i.e. normal sequence, normal (empty)
 *        DP matrices.
 *
 *  @bug  Due to resizing, any pre-existing constraints will be lost!
 *
 *  @ingroup subopt_zuker
 *
 *  @see vrna_subopt(), zukersubopt(), zukersubopt_par()
 *
 *  @param  vc  fold compound
 *  @return     List of zuker suboptimal structures
 */
vrna_subopt_solution_t *
vrna_subopt_zuker(vrna_fold_compound_t *fc);


#endif
