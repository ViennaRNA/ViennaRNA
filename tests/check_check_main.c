#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* malloc, free, rand */

#include <check.h>

Suite *make_fold_suite();
Suite *make_constraints_suite();
Suite *make_utils_suite();
Suite *make_energies_suite();

int
main(void)
{
  int number_failed;
  SRunner *sr   = srunner_create(make_fold_suite());
  srunner_add_suite(sr, make_constraints_suite());
  srunner_add_suite(sr, make_utils_suite());
  srunner_add_suite(sr, make_energies_suite());

  srunner_set_fork_status(sr, CK_NOFORK);

  srunner_set_log (sr, "test.log");

  srunner_run_all(sr, CK_VERBOSE);

  number_failed = srunner_ntests_failed(sr);

  srunner_free(sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

