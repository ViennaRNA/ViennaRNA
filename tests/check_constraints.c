#include <check.h>

#include <ViennaRNA/constraints.h>

#include <math.h>
#include <stdio.h>
#include <unistd.h>

static int deltaCompare(double a, double b)
{
  if (fabs(a - b) < 1e-5)
    return 1;

  printf("%f != %f\n", a, b);
  return 0;
}

static void writeTempFile(char *tempfile, const char *data)
{
  FILE *f;

  ck_assert(tmpnam(tempfile) != NULL);
  f = fopen(tempfile, "w");
  ck_assert(f != NULL);
  fputs(data, f);
  fclose(f);
}

START_TEST(test_parse_soft_constraints_file)
{
  char tempfile[L_tmpnam + 1];
  const size_t len = 5;
  char sequence[len];
  double values[len];

  //1 entry
  writeTempFile(tempfile, "1 A 0.5\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 1, 0, sequence, values), 1);
  ck_assert_str_eq(sequence, "A");
  ck_assert(deltaCompare(values[1], 0.5));
  unlink(tempfile);

  //multiple entries
  writeTempFile(tempfile, "1 A 0.1\n2 T 0.2\n3 G 0.3\n4 C 0.4\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 4, 0, sequence, values), 1);
  ck_assert_str_eq(sequence, "ATGC");
  ck_assert(deltaCompare(values[1], 0.1));
  ck_assert(deltaCompare(values[2], 0.2));
  ck_assert(deltaCompare(values[3], 0.3));
  ck_assert(deltaCompare(values[4], 0.4));
  unlink(tempfile);

  //value formats
  writeTempFile(tempfile, "1 A 1\n2 T 2.\n3 G .3\n4 C 1e-1\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 4, 0, sequence, values), 1);
  ck_assert_str_eq(sequence, "ATGC");
  ck_assert(deltaCompare(values[1], 1));
  ck_assert(deltaCompare(values[2], 2));
  ck_assert(deltaCompare(values[3], 0.3));
  ck_assert(deltaCompare(values[4], 0.1));
  unlink(tempfile);

  //missing value
  writeTempFile(tempfile, "1 A\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 1, 123, sequence, values), 1);
  ck_assert_str_eq(sequence, "A");
  ck_assert(deltaCompare(values[1], 123));
  unlink(tempfile);

  //missing nucleotide & value
  writeTempFile(tempfile, "1\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 1, 123, sequence, values), 1);
  ck_assert_str_eq(sequence, "N");
  ck_assert(deltaCompare(values[1], 123));
  unlink(tempfile);

  //upper limit
  writeTempFile(tempfile, "1 A 0.1\n2 T 0.2\n3 G 0.3\n4 C 0.4\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 3, 0, sequence, values), 0);
  unlink(tempfile);

  //lower limit
  writeTempFile(tempfile, "0 A 0.1\n1 T 0.2\n2 G 0.3\n3 C 0.4\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 4, 0, sequence, values), 0);
  unlink(tempfile);

  //unordered
  writeTempFile(tempfile, "3 G 0.3\n2 T 0.2\n4 C 0.4\n1 A 0.1\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 4, 0, sequence, values), 1);
  ck_assert_str_eq(sequence, "ATGC");
  ck_assert(deltaCompare(values[1], 0.1));
  ck_assert(deltaCompare(values[2], 0.2));
  ck_assert(deltaCompare(values[3], 0.3));
  ck_assert(deltaCompare(values[4], 0.4));
  unlink(tempfile);

  //missing indices middle
  writeTempFile(tempfile, "1 A 0.1\n4 C 0.4\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 4, 123, sequence, values), 1);
  ck_assert_str_eq(sequence, "ANNC");
  ck_assert(deltaCompare(values[1], 0.1));
  ck_assert(deltaCompare(values[2], 123));
  ck_assert(deltaCompare(values[3], 123));
  ck_assert(deltaCompare(values[4], 0.4));
  unlink(tempfile);

  //missing indices start end
  writeTempFile(tempfile, "2 T 0.2\n3 G 0.3\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 4, 123, sequence, values), 1);
  ck_assert_str_eq(sequence, "NTGN");
  ck_assert(deltaCompare(values[1], 123));
  ck_assert(deltaCompare(values[2], 0.2));
  ck_assert(deltaCompare(values[3], 0.3));
  ck_assert(deltaCompare(values[4], 123));
  unlink(tempfile);

  //invalid file
  ck_assert_int_eq(parse_soft_constraints_file(NULL, 0, 0, sequence, values), 0);

  //missing linebreak
  writeTempFile(tempfile, "1 A 0.5");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 1, 0, sequence, values), 1);
  ck_assert_str_eq(sequence, "A");
  ck_assert(deltaCompare(values[1], 0.5));
  unlink(tempfile);

  //garbage + entry
  writeTempFile(tempfile, "\nblablabla\n1 A 0.5\n\ngarbage\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 1, 0, sequence, values), 1);
  ck_assert_str_eq(sequence, "A");
  ck_assert(deltaCompare(values[1], 0.5));
  unlink(tempfile);

  //garbage only
  writeTempFile(tempfile, "\nblablabla\n\ngarbage\n");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 1, 0, sequence, values), 0);
  unlink(tempfile);

  //empty file
  writeTempFile(tempfile, "");
  ck_assert_int_eq(parse_soft_constraints_file(tempfile, 1, 0, sequence, values), 0);
  unlink(tempfile);
}
END_TEST

TCase* constraints_testcase()
{
  TCase *tc = tcase_create("constraints");
  tcase_add_test(tc, test_parse_soft_constraints_file);

  return tc;
}
