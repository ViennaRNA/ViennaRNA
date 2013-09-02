#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* malloc, free, rand */

#include <check.h>
#include <ViennaRNA/fold.h>


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

Suite *
fold_suite(void)
{
  Suite *s = suite_create("Fold");

  /* Core test case */
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_fold);
  suite_add_tcase(s, tc_core);

  return s;
}


int
main(void)
{
  int number_failed;
  Suite *s      = fold_suite();
  SRunner *sr   = srunner_create(s);
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

