#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* malloc, free, rand */

#include <check.h>
#include <ViennaRNA/constraints.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/utils.h>


START_TEST(test_fold)
{
  /* unit test code */
  const char *seq1 = "CGCAGGGAUACCCGCG";
  const char *str1 = "(((.(((...))))))";
  char *structure = (char*)malloc(sizeof(char) * 17);
  float en = fold(seq1, structure);
  ck_assert(strcmp(str1, structure) == 0);
  free(structure);
}
END_TEST

START_TEST(test_sample_structure)
{
  vrna_md_t md;
  vrna_fold_compound *vc;
  const char sequence[] = "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU";
  char *sample;

  set_model_details(&md);
  md.uniq_ML = 1;
  md.compute_bpp = 0;

  vc = vrna_get_fold_compound(sequence, &md, VRNA_OPTION_PF);

  vrna_pf_fold(vc, NULL);

  sample = vrna_pbacktrack(vc);
  ck_assert_int_eq(strlen(sample), sizeof(sequence) - 1);
  free(sample);

  sample = vrna_pbacktrack5(vc, 16);
  ck_assert_int_eq(strlen(sample), 16);
  free(sample);

  vrna_free_fold_compound(vc);
}
END_TEST

START_TEST(test_sc_sanity_check)
{
  vrna_md_t md;
  vrna_fold_compound *vc;
  const char sequence[] = "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU";
  double *sc_up;
  double **sc_bp;
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


  set_model_details(&md);
  md.compute_bpp = 1;

  vc = vrna_get_fold_compound(sequence, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);

  mfe_energy_unconstrained = vrna_fold(vc, mfe_structure_unconstrained);
  pf_energy_unconstrained = vrna_pf_fold(vc, pf_structure_unconstrained);
  plist_unconstrained = vrna_get_plist_from_pr(vc, 0);

  sc_up = (double *)space(sizeof(double) * (length + 1));
  sc_bp = (double **)space(sizeof(double *) * (length + 1));

  for(i = 1; i <= length; ++i)
  {
    sc_up[i] = -1.;
    sc_bp[i] = (double *)space(sizeof(double) * (length + 1));
    for(j = i + 1; j <= length; ++j)
      sc_bp[i][j] = -2.;
  }

  vrna_sc_add_up(vc, sc_up, VRNA_CONSTRAINT_SOFT_MFE | VRNA_CONSTRAINT_SOFT_PF);
  vrna_sc_add_bp(vc, sc_bp, VRNA_CONSTRAINT_SOFT_MFE | VRNA_CONSTRAINT_SOFT_PF);

  mfe_energy_constrained = vrna_fold(vc, mfe_structure_constrained);
  pf_energy_constrained = vrna_pf_fold(vc, pf_structure_constrained);
  plist_constrained = vrna_get_plist_from_pr(vc, 0);

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

  vrna_free_fold_compound(vc);

  for(i = 1; i <= length; ++i)
    free(sc_bp[i]);
  free(sc_bp);
  free(sc_up);

  free(plist_constrained);
  free(plist_unconstrained);
}
END_TEST

TCase* constraints_testcase();
TCase* utils_testcase();
TCase* loop_energies_testcase();

Suite *
fold_suite(void)
{
  Suite *s = suite_create("Fold");

  /* Core test case */
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_fold);
  tcase_add_test(tc_core, test_sample_structure);
  tcase_add_test(tc_core, test_sc_sanity_check);
  suite_add_tcase(s, tc_core);

  suite_add_tcase(s, constraints_testcase());
  suite_add_tcase(s, utils_testcase());
  suite_add_tcase(s, loop_energies_testcase());

  return s;
}


int
main(void)
{
  int number_failed;
  Suite *s      = fold_suite();
  SRunner *sr   = srunner_create(s);
  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

