/* subopt_zuker.h */
#ifndef VIENNA_RNA_PACKAGE_SUBOPT_ZUKER_H
#define VIENNA_RNA_PACKAGE_SUBOPT_ZUKER_H

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/subopt/basic.h>

/**
 *  @brief Compute Zuker type suboptimal structures
 *
 *  Compute Suboptimal structures according to @rstinline :cite:t:`zuker:1989` @endrst , i.e. for every
 *  possible base pair the minimum energy structure containing the resp. base pair.
 *  Returns a list of these structures and their energies.
 *
 *  @ingroup subopt_zuker
 *
 *  @see vrna_subopt(), zukersubopt(), zukersubopt_par()
 *
 *  @param  fc  fold compound
 *  @return     List of zuker suboptimal structures
 */
vrna_subopt_solution_t *
vrna_subopt_zuker(vrna_fold_compound_t *fc);


#endif
