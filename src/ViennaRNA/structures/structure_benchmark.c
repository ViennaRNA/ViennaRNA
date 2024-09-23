#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <math.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/structures/benchmark.h"

PUBLIC vrna_score_t
vrna_score_from_confusion_matrix(int  TP,
                                 int  TN,
                                 int  FP,
                                 int  FN)
{
  vrna_score_t score;

  score.TP  = TP;
  score.TN  = TN;
  score.FP  = FP;
  score.FN  = FN;

  score.TPR = TP + FN > 0 ? (float)TP / (TP + FN) : 0.0;
  score.PPV = TP + FP > 0 ? (float)TP / (TP + FP) : 0.0;
  score.FPR = TN + FP > 0 ? (float)FP / (TN + FP) : 0.0;
  score.FOR = TN + FN > 0 ? (float)FN / (TN + FN) : 0.0;
  score.FDR = 1 - score.PPV;
  score.FNR = 1 - score.TPR;
  score.TNR = 1 - score.FPR;
  score.NPV = 1 - score.FOR;
  score.F1  = score.TPR + score.PPV >
              0 ? (float)2 * score.TPR * score.PPV / (score.TPR + score.PPV) : 0.0;

  score.MCC = sqrt(score.PPV * score.TPR * score.TNR * score.NPV) -
              sqrt(score.FDR * score.FNR * score.FPR * score.FOR);

  return score;
}


PUBLIC vrna_score_t
vrna_compare_structure_pt(const short *pt_gold,
                          const short *pt_other,
                          int         fuzzy)
{
  int TP = 0, TN = 0, FP = 0, FN = 0;
  /* number of false positive, but compatible pairs */
  int compatible = 0;
  /* count the number of base pairs in pt_gold, and pt_other */
  int bps_gold = 0, bps_other = 0;

  /* Hard coded fuzzy as 0 */
  fuzzy = 0;

  for (int i = 1; i <= pt_gold[0]; i++)
    if (pt_gold[i] > i)
      bps_gold++;

  for (int i = 1; i <= pt_other[0]; i++)
    if (pt_other[i] > i)
      bps_other++;

  /* got through pt_other and compare its pairs to pt_gold */
  for (int i = 1; i <= pt_other[0]; i++) {
    int j = pt_other[i];
    if (j < i)
      continue;

    int is_true_positive = 0, is_inconsistent = 0, is_contradicting = 0;

    /*
     * True Positives
     * let's see if position i matches in gold position j +/- fuzzy
     */
    for (int add = 0; add <= fuzzy; add++) {
      if ((pt_gold[i] == j + add) && (j + add <= pt_other[0])) {
        is_true_positive = 1;
        bps_gold--;
        break;
      }

      if ((pt_gold[i] == j - add) && (j - add >= 1)) {
        is_true_positive = 1;
        bps_gold--;
        break;
      }
    }
    /* let's see if position j matches in gold position i +/- fuzzy */
    if (is_true_positive == 0) {
      for (int add = 0; add <= fuzzy; add++) {
        if ((pt_gold[j] == i + add) && (i + add <= pt_other[0])) {
          is_true_positive = 1;
          bps_gold--;
          break;
        }

        if ((pt_gold[j] == i - add) && (i - add >= 1)) {
          is_true_positive = 1;
          bps_gold--;
          break;
        }
      }
    }

    /*
     * False Positives
     * let's check if base-pair i, j is inconsistent with reference structure
     * base i or base j in reference is paired to something else
     */
    if (is_true_positive == 0)
      is_inconsistent = 1;

    /*
     * let's check if a base-pair is contradicting
     * there is a pair k * l in ther reference that k < i < l < j
     */
    for (int k = 1; k <= pt_gold[0]; k++) {
      if (pt_gold[k] < k)
        continue;

      int l = pt_gold[k];
      if ((k < i) && (i < l) && (l < j)) {
        is_contradicting = 1;
        break;
      }
    }

    if (is_true_positive == 0) {
      if ((is_inconsistent == 0) && (is_contradicting == 0))
        compatible++;

      FP++;
    } else {
      TP++;
    }
  }

  FN  = bps_gold;
  TN  = (int)(0.5 * (pt_gold[0] * (pt_gold[0] - 1))) - TP - FP - FN;

  return vrna_score_from_confusion_matrix(TP, TN, FP, FN);
}


PUBLIC vrna_score_t
vrna_compare_structure(const char   *structure_gold,
                       const char   *structure_other,
                       int          fuzzy,
                       unsigned int options)
{
  short         *pt_gold  = NULL;
  short         *pt_other = NULL;
  vrna_score_t  score;

  pt_gold   = vrna_ptable_from_string(structure_gold, options);
  pt_other  = vrna_ptable_from_string(structure_other, options);
  score     = vrna_compare_structure_pt(pt_gold, pt_other, fuzzy);

  free(pt_gold);
  free(pt_other);
  return score;
}
