#ifndef VIENNA_RNA_PACKAGE_STRUCTURES_BENCHMARK_H
#define VIENNA_RNA_PACKAGE_STRUCTURES_BENCHMARK_H

/**
 *  @file     ViennaRNA/structures/benchmark.h
 *  @ingroup  struct_utils
 *  @brief    Various utility- and helper-functions for secondary structure parsing, converting, etc.
 */

#include <stdio.h>

#include <ViennaRNA/datastructures/basic.h>

/**
 *  @addtogroup struct_utils_benchmark
 *  @{
 */

/**
 *  @brief Typename for the score data structure #vrna_score_s
 *  @ingroup  struct_utils_benchmark
 */
typedef struct vrna_score_s vrna_score_t;


/**
 *  @brief The data structure that contains statistic result of two structures comparaison
 *
 */
struct vrna_score_s {
  double  TP;             /**< @brief True Positive count */
  double  TN;             /**< @brief True Negative count */
  double  FP;             /**< @brief False Positive count */
  double  FN;             /**< @brief False Negative count */
  double  TPR;            /**< @brief True Positive Rate */
  double  PPV;            /**< @brief Positive Predictive Value */
  double  FPR;            /**< @brief False Positive Rate */
  double  FOR;            /**< @brief False Omission Rate */
  double  TNR;            /**< @brief True Negative Rate */
  double  FDR;            /**< @brief False Discovery Rate */
  double  FNR;            /**< @brief False Negative Rate  */
  double  NPV;            /**< @brief Negative Predictive Value */
  double  F1;             /**< @brief F1 Score */
  double  MCC;            /**< @brief Matthews Correlation Coefficient */
};


/**
 *  @brief Construct score data structure from given confusion matrix
 *
 *  @param TP          True positive count
 *  @param TN          True negative count
 *  @param FP          False positive count
 *  @param FN          False negative count
 *  @return  The score data structure to write
 */
vrna_score_t
vrna_score_from_confusion_matrix(double TP,
                                 double TN,
                                 double FP,
                                 double FN);


/**
 *  @brief Return statistic of two structure (in pair table) comparaison
 *
 *  @see   vrna_compare_structure
 *
 *  @param pt_gold   Gold standard structure in pair table
 *  @param pt_other  Structure to compare in pair table
 *  @param fuzzy     Allows for base pair slippage. Hence, for any base pair (i,j) in the gold standard, a base pair (p, q) in the second structure is considered a true positive, if i - fuzzy <= p <= i + fuzzy, and j - fuzzy <= q <= j + fuzzy. Each gold standard base pair may be claimed as a true positive by at most one predicted base pair: once matched (within the fuzzy window) it is consumed and can't be matched again by another predicted pair. Matching is greedy in the order predicted pairs are encountered, not a globally optimal assignment, but for typical slippage values (fuzzy of a few nucleotides) this does not affect the result.
 *  @return          The #vrna_score_s data structure
 */
vrna_score_t
vrna_compare_structure_pt(const short *pt_gold,
                          const short *pt_other,
                          int         fuzzy);


/**
 *  @brief Return statistic of two structure (in dbn) comparaison
 *
 *  @param structure_gold   Gold standard structure
 *  @param structure_other  Structure to compare
 *  @param fuzzy            Allows for base pair slippage. Hence, for any base pair (i,j) in the gold standard, a base pair (p, q) in the second structure is considered a true positive, if i - fuzzy <= p <= i + fuzzy, and j - fuzzy <= q <= j + fuzzy. Each gold standard base pair may be claimed as a true positive by at most one predicted base pair: once matched (within the fuzzy window) it is consumed and can't be matched again by another predicted pair. Matching is greedy in the order predicted pairs are encountered, not a globally optimal assignment, but for typical slippage values (fuzzy of a few nucleotides) this does not affect the result.
 *  @param  options         A bitmask to specify which brackets are recognized during conversion to pair table
 *  @return                 The #vrna_score_s data structure
 */
vrna_score_t
vrna_compare_structure(const char   *structure_gold,
                       const char   *structure_other,
                       int          fuzzy,
                       unsigned int options);


/* End benchmark */
/** @} */

#endif
