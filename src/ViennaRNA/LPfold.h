#ifndef VIENNA_RNA_PACKAGE_LPFOLD_H
#define VIENNA_RNA_PACKAGE_LPFOLD_H

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @file     LPfold.h
 *  @ingroup  part_func_window_deprecated
 *  @brief    Partition function and equilibrium probability implementation for the sliding window algorithm
 *
 *  This file contains the implementation for sliding window partition function and equilibrium
 *  probabilities. It also provides the unpaired probability implementation from
*   @rstinline :cite:t:`bernhart:2011` @endrst
 */

#include <stdio.h>

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/partfunc/local.h>

#ifdef VRNA_WARN_DEPRECATED
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

/**
 *  \brief
 *
 *  \ingroup part_func_window_deprecated
 *
 *  \param  length
 */
DEPRECATED(void update_pf_paramsLP(int length),
"This function is obsolete");


/**
 *  \brief
 *
 *  \ingroup part_func_window_deprecated
 *
 */
DEPRECATED(void update_pf_paramsLP_par(int              length,
                                       vrna_exp_param_t *parameters),
"Use the new API with vrna_folc_compound_t instead");


/**
 *  \brief Compute partition functions for locally stable secondary structures
 *
 *  pfl_fold computes partition functions for every window of size
 *  'winSize' possible in a RNA molecule, allowing only pairs with a span
 *  smaller than 'pairSize'. It returns the mean pair probabilities averaged
 *  over all windows containing the pair in 'pl'. 'winSize' should
 *  always be >= 'pairSize'. Note that in contrast to Lfold(),
 *  bases outside of the window do not influence the structure at all. Only
 *  probabilities higher than 'cutoffb' are kept.
 *
 *  If 'pU' is supplied (i.e is not the NULL pointer), pfl_fold()
 *  will also compute the mean probability that regions of length 'u' and smaller are
 *  unpaired. The parameter 'u' is supplied in 'pup[0][0]'. On return
 *  the 'pup' array will contain these probabilities, with the entry on
 *  'pup[x][y]' containing the mean probability that x and the y-1
 *  preceding bases are unpaired. The 'pU' array needs to be large
 *  enough to hold n+1 float* entries, where n is the sequence length.
 *
 *  If an array dpp2 is supplied, the probability of base pair (i,j)
 *  given that there already exists a base pair (i+1,j-1) is also
 *  computed and saved in this array. If pUfp is given (i.e. not NULL), pU
 *  is not saved but put out imediately. If spup is given (i.e. is not NULL),
 *  the pair probabilities in pl are not saved but put out imediately.
 *
 *  \ingroup part_func_window_deprecated
 *
 *  \param  sequence  RNA sequence
 *  \param  winSize   size of the window
 *  \param  pairSize  maximum size of base pair
 *  \param  cutoffb   cutoffb for base pairs
 *  \param  pU        array holding all unpaired probabilities
 *  \param  dpp2  array of dependent pair probabilities
 *  \param  pUfp     file pointer for pU
 *  \param  spup     file pointer for pair probabilities
 *  \return list of pair probabilities
 */
DEPRECATED(vrna_ep_t *pfl_fold(char          *sequence,
                               int           winSize,
                               int           pairSize,
                               float         cutoffb,
                               double        **pU,
                               vrna_ep_t     **dpp2,
                               FILE          *pUfp,
                               FILE          *spup),
"Use vrna_pfl_fold(), vrna_pfl_fold_cb(), vrna_pfl_fold_up(), or vrna_pfl_fold_up_cb() instead");


/**
 *  \brief Compute partition functions for locally stable secondary structures
 *
 *  \ingroup part_func_window_deprecated
 *
 */
DEPRECATED(vrna_ep_t *pfl_fold_par(char              *sequence,
                                      int               winSize,
                                      int               pairSize,
                                      float             cutoffb,
                                      double            **pU,
                                      vrna_ep_t         **dpp2,
                                      FILE              *pUfp,
                                      FILE              *spup,
                                      vrna_exp_param_t  *parameters),
"Use the new API and vrna_probs_window() instead");


DEPRECATED(void putoutpU_prob_par(double            **pU,
                                  int               length,
                                  int               ulength,
                                  FILE              *fp,
                                  int               energies,
                                  vrna_exp_param_t  *parameters),
"");


/**
 *  \brief Writes the unpaired probabilities (pU) or opening energies into a file
 *
 *  Can write either the unpaired probabilities (accessibilities) pU or
 *  the opening energies -log(pU)kT into a file
 *
 *  \ingroup part_func_window_deprecated
 *
 *  \param  pU       pair probabilities
 *  \param  length   length of RNA sequence
 *  \param  ulength  maximum length of unpaired stretch
 *  \param  fp file pointer of destination file
 *  \param  energies  switch to put out as  opening energies
 */
DEPRECATED(void    putoutpU_prob(double **pU,
                                 int    length,
                                 int    ulength,
                                 FILE   *fp,
                                 int    energies),
"");


DEPRECATED(void putoutpU_prob_bin_par(double            **pU,
                                      int               length,
                                      int               ulength,
                                      FILE              *fp,
                                      int               energies,
                                      vrna_exp_param_t  *parameters),
"");


/**
 *  \brief Writes the unpaired probabilities (pU) or opening energies into a binary file
 *
 *  Can write either the unpaired probabilities (accessibilities) pU or
 *  the opening energies -log(pU)kT into a file
 *
 *  \ingroup part_func_window_deprecated
 *
 *  \param  pU       pair probabilities
 *  \param  length   length of RNA sequence
 *  \param  ulength  maximum length of unpaired stretch
 *  \param  fp file pointer of destination file
 *  \param  energies  switch to put out as  opening energies
 */
DEPRECATED(void    putoutpU_prob_bin(double **pU,
                                     int    length,
                                     int    ulength,
                                     FILE   *fp,
                                     int    energies),
"");


/**
 *  Dunno if this function was ever used by external programs linking to RNAlib, but it
 *  was declared PUBLIC before.
 *  Anyway, never use this function as it will be removed soon and does nothing at all
 */
DEPRECATED(void init_pf_foldLP(int length),
"This function is obsolete");

#endif

#endif
