#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* malloc, free, rand */

#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/constraints.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>

#suite  MFE_Prediction

#tcase  Backward_Compatibility

#test test_fold
  /* unit test code */
  const char *seq1 = "CGCAGGGAUACCCGCG";
  const char *str1 = "(((.(((...))))))";
  char *structure = (char*)malloc(sizeof(char) * 17);
  float en = fold(seq1, structure);
  ck_assert(strcmp(str1, structure) == 0);
  free(structure);

#suite  Partition_Function

#tcase Stochastic_Backtracking

#test test_sample_structure
  vrna_md_t md;
  vrna_fold_compound_t *vc;
  const char sequence[] = "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU";
  char *sample;

  vrna_md_set_default(&md);
  md.uniq_ML = 1;
  md.compute_bpp = 0;

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_PF);

  vrna_pf(vc, NULL);

  sample = vrna_pbacktrack(vc);
  ck_assert_int_eq(strlen(sample), sizeof(sequence) - 1);
  free(sample);

  sample = vrna_pbacktrack5(vc, 16);
  ck_assert_int_eq(strlen(sample), 16);
  free(sample);

  vrna_fold_compound_free(vc);

#suite  Constraints_Implementation

#tcase  Soft_Constraints

#test test_sc_sanity_check
  vrna_md_t md;
  vrna_fold_compound_t *vc;
  const char sequence[] = "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU";
  FLT_OR_DBL *sc_up;
  FLT_OR_DBL **sc_bp;
  int i, j;
  const int length = sizeof(sequence) - 1;

  char mfe_structure_unconstrained[length + 1];
  char mfe_structure_constrained[length + 1];
  char pf_structure_unconstrained[length + 1];
  char pf_structure_constrained[length + 1];

  double mfe_energy_unconstrained;
  double mfe_energy_constrained;
  double pf_energy_unconstrained;
  double pf_energy_constrained;

  plist *plist_unconstrained;
  plist *plist_constrained;


  vrna_md_set_default(&md);
  md.compute_bpp = 1;

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);

  mfe_energy_unconstrained = (double)vrna_mfe(vc, mfe_structure_unconstrained);
  vrna_exp_params_rescale(vc, &mfe_energy_unconstrained);
  pf_energy_unconstrained = (double)vrna_pf(vc, pf_structure_unconstrained);
  plist_unconstrained = vrna_plist_from_probs(vc, 0);

  sc_up = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (length + 1));
  sc_bp = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (length + 1));

  for(i = 1; i <= length; ++i)
  {
    sc_up[i] = -1;
    sc_bp[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (length + 1));
    for(j = i + 1; j <= length; ++j)
      sc_bp[i][j] = -2;
  }

  vrna_sc_set_up(vc, (const FLT_OR_DBL *)sc_up, VRNA_OPTION_MFE | VRNA_OPTION_PF);
  vrna_sc_set_bp(vc, (const FLT_OR_DBL **)sc_bp, VRNA_OPTION_MFE | VRNA_OPTION_PF);

  mfe_energy_constrained = (double)vrna_mfe(vc, mfe_structure_constrained);
  vrna_exp_params_rescale(vc, &mfe_energy_constrained);
  pf_energy_constrained = (double)vrna_pf(vc, pf_structure_constrained);
  plist_constrained = vrna_plist_from_probs(vc, 0);

  ck_assert_int_eq(strlen(mfe_structure_constrained), sizeof(sequence) - 1);
  ck_assert_int_eq(strlen(mfe_structure_unconstrained), sizeof(sequence) - 1);
  ck_assert_int_eq(strlen(pf_structure_constrained), sizeof(sequence) - 1);
  ck_assert_int_eq(strlen(pf_structure_unconstrained), sizeof(sequence) - 1);

  ck_assert(mfe_energy_constrained != mfe_energy_unconstrained);
  ck_assert(pf_energy_constrained != pf_energy_unconstrained);

  ck_assert_str_eq(mfe_structure_constrained, mfe_structure_unconstrained);
  ck_assert_str_eq(pf_structure_constrained, pf_structure_unconstrained);

  for (i = 0; plist_constrained[i].i || plist_constrained[i].j || plist_unconstrained[i].i || plist_unconstrained[i].j; ++i)
  {
    ck_assert_int_eq(plist_constrained[i].i, plist_unconstrained[i].i);
    ck_assert_int_eq(plist_constrained[i].j, plist_unconstrained[i].j);
    ck_assert_int_eq(plist_constrained[i].type, plist_unconstrained[i].type);
    ck_assert(plist_constrained[i].p == plist_unconstrained[i].p);
  }

  vrna_fold_compound_free(vc);

  for(i = 1; i <= length; ++i)
    free(sc_bp[i]);
  free(sc_bp);
  free(sc_up);

  free(plist_constrained);
  free(plist_unconstrained);
