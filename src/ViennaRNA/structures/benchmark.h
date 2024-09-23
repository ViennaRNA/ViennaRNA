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
  int   TP;             /**< @brief True Positive count */
  int   TN;             /**< @brief True Negative count */
  int   FP;             /**< @brief False Positive count */
  int   FN;             /**< @brief False Negative count */
  float TPR;            /**< @brief True Positive Rate */
  float PPV;            /**< @brief Positive Predictive Value */
  float FPR;            /**< @brief False Positive Rate */
  float FOR;            /**< @brief False Omission Rate */
  float TNR;            /**< @brief True Negative Rate */
  float FDR;            /**< @brief False Discovery Rate */
  float FNR;            /**< @brief False Negative Rate  */
  float NPV;            /**< @brief Negative Predictive Value */
  float F1;             /**< @brief F1 Score */
  float MCC;            /**< @brief Matthews Correlation Coefficient */
};


/**
 *  @brief Construct score data structure from given confusion matrix
 *
 *  @param TP     True positive count
 *  @param TN     True negative count
 *  @param FP     False positive count
 *  @param FN     False negative count
 *  @return  The score data structure to write
 */
vrna_score_t
vrna_score_from_confusion_matrix(int  TP,
                                 int  TN,
                                 int  FP,
                                 int  FN);


/**
 *  @brief Return statistic of two structure (in pair table) comparaison
 *
 *  @see   vrna_compare_structure
 *
 *  @param pt_gold   Gold standard structure in pair table
 *  @param pt_other  Structure to compare in pair table
 *  @param fuzzy     Allows for base pair slippage. Hence, for any base pair (i,j) in the gold standard, a base pair (p, q) in the second structure is considered a true positive, if i - fuzzy <= p <= i + fuzzy, and j - fuzzy <= q <= j + fuzzy.
 *  @return          The #vrna_score_s data structure
 */
vrna_score_t
vrna_compare_structure_pt(const short *pt_gold,
                          const short *pt_other,
                          int         fuzzy);


/**
 *  @brief Return statistic of two structure (in dbn) comparaison
 *
 *  @param pt_gold   Gold standard structure
 *  @param pt_other  Structure to compare
 *  @param fuzzy     Allows for base pair slippage. Hence, for any base pair (i,j) in the gold standard, a base pair (p, q) in the second structure is considered a true positive, if i - fuzzy <= p <= i + fuzzy, and j - fuzzy <= q <= j + fuzzy.
 *  @param  options  A bitmask to specify which brackets are recognized during conversion to pair table
 *  @return          The #vrna_score_s data structure
 */
vrna_score_t
vrna_compare_structure(const char   *structure_gold,
                       const char   *structure_other,
                       int          fuzzy,
                       unsigned int options);


/* End benchmark */
/** @} */

#endif
