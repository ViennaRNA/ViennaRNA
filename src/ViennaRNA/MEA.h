#ifndef VIENNA_RNA_PACKAGE_MEA_H
#define VIENNA_RNA_PACKAGE_MEA_H

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>

/**
 *  @file MEA.h
 *  @ingroup subopt_and_representatives
 *  @brief Computes a MEA (maximum expected accuracy) structure.
 */

/**
 *  @brief Computes a MEA (maximum expected accuracy) structure.
 *
 *  @ingroup  mea_fold
 *
 *  The algorithm maximizes the expected accuracy
 *  @f[ A(S) = \sum_{(i,j) \in S} 2 \gamma p_{ij} + \sum_{i \notin S} p^u_i @f]
 *  Higher values of @f$\gamma@f$ result in more base pairs of lower
 *  probability and thus higher sensitivity. Low values of @f$\gamma@f$ result in structures
 *  containing only highly likely pairs (high specificity).
 *  The code of the MEA function also demonstrates the use of sparse dynamic
 *  programming scheme to reduce the time and memory complexity of folding.
 */
float MEA(plist *p,
          char *structure,
          double gamma);

float MEA_seq(plist *p,
              const char *sequence,
              char *structure,
              double gamma,
              vrna_exp_param_t *pf);

#endif
