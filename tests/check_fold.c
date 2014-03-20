#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* malloc, free, rand */

#include <check.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>


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
  model_detailsT md;
  vrna_fold_compound *vc;
  const char sequence[] = "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU";
  char *sample;

  set_model_details(&md);
  md.uniq_ML = 1;
  md.do_backtrack = 0;

  vc = vrna_get_fold_compound(sequence, &md, VRNA_OPTION_PF);

  vrna_pf_fold(vc, NULL);

  sample = vrna_pbacktrack(vc);
  ck_assert_int_eq(strlen(sample), sizeof(sequence) - 1);
  free(sample);

  sample = vrna_pbacktrack5(vc, 16);
  ck_assert_int_eq(strlen(sample), 16);
  free(sample);

  destroy_fold_compound(vc);
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

