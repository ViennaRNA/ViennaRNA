#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include <ViennaRNA/file_formats.h>
#include <ViennaRNA/constraints.h>

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

/* end of prologue */

#suite Constraints

#tcase  SoftConstraints

#test test_vrna_sc_SHAPE_to_pr
  int ret;
  double negative_values[] = {0, -100, -1, -1e-10};
  double hardcoded_range_values[] = { 0, 0.125, 0.25, 0.275, 0.3, 0.5, 0.7};
  double upper_range_values[] = { 0, 0.8, 0.9, 1};
  double upper_range_values2[] = { 0, 1.2, 1.7};
  double cutoff_values_default[] = { 0, -1, 0.24, 0.25, 0.26 };
  double cutoff_values[] = { 0, -1, 0.49, 0.50, 0.51 };
  double skip_values[] = { 0, -1, 0.5, 2 };
  double linear_values[] = { 0, -1, 0, 0.25, 0.5, 0.75, 1, 2 };
  double linear_custom_values[] = { 0, -1, 0, 0.25, 0.5, 0.75, 1, 2 };
  double log_values[] = { 0, -1, 0, 0.25, 0.5, 0.75, 1, 2 };
  double log_values_custom[] = { 0, -1, 0, 0.25, 0.5, 0.75, 1, 2 };

  ret = vrna_sc_SHAPE_to_pr(NULL, NULL, 0, 0);
  ck_assert_int_eq(ret, 0);

  ret = vrna_sc_SHAPE_to_pr("", NULL, 0, 0);
  ck_assert_int_eq(ret, 0);

  ret = vrna_sc_SHAPE_to_pr("X", NULL, 0, 0);
  ck_assert_int_eq(ret, 0);

  ret = vrna_sc_SHAPE_to_pr("M", NULL, 0, 0);
  ck_assert_int_eq(ret, 0);

  ret = vrna_sc_SHAPE_to_pr("M", negative_values, 3, 0.123);
  ck_assert_int_eq(ret, 1);
  ck_assert(deltaCompare(negative_values[1], 0.123));
  ck_assert(deltaCompare(negative_values[2], 0.123));
  ck_assert(deltaCompare(negative_values[3], 0.123));

  ret = vrna_sc_SHAPE_to_pr("M", hardcoded_range_values, 6, 0);
  ck_assert_int_eq(ret, 1);
  ck_assert(deltaCompare(hardcoded_range_values[1], 0.175));
  ck_assert(deltaCompare(hardcoded_range_values[2], 0.35));
  ck_assert(deltaCompare(hardcoded_range_values[3], 0.45));
  ck_assert(deltaCompare(hardcoded_range_values[4], 0.55));
  ck_assert(deltaCompare(hardcoded_range_values[5], 0.7));
  ck_assert(deltaCompare(hardcoded_range_values[6], 0.85));

  ret = vrna_sc_SHAPE_to_pr("M", upper_range_values, 3, 0);
  ck_assert_int_eq(ret, 1);
  ck_assert(deltaCompare(upper_range_values[1], 0.9));
  ck_assert(deltaCompare(upper_range_values[2], 0.95));
  ck_assert(deltaCompare(upper_range_values[3], 1));

  ret = vrna_sc_SHAPE_to_pr("M", upper_range_values2, 2, 0);
  ck_assert_int_eq(ret, 1);
  ck_assert(deltaCompare(upper_range_values2[1], 0.925));
  ck_assert(deltaCompare(upper_range_values2[2], 1));

  ret = vrna_sc_SHAPE_to_pr("C", cutoff_values_default, 4, 0.5);
  ck_assert_int_eq(ret, 1);
  ck_assert(deltaCompare(cutoff_values_default[1], 0.5));
  ck_assert(deltaCompare(cutoff_values_default[2], 0));
  ck_assert(deltaCompare(cutoff_values_default[3], 1));
  ck_assert(deltaCompare(cutoff_values_default[4], 1));

  ret = vrna_sc_SHAPE_to_pr("C0.5", cutoff_values, 4, 0.5);
  ck_assert_int_eq(ret, 1);
  ck_assert(deltaCompare(cutoff_values[1], 0.5));
  ck_assert(deltaCompare(cutoff_values[2], 0));
  ck_assert(deltaCompare(cutoff_values[3], 1));
  ck_assert(deltaCompare(cutoff_values[4], 1));

  ret = vrna_sc_SHAPE_to_pr("S", skip_values, 3, 0.5);
  ck_assert_int_eq(ret, 1);
  ck_assert(deltaCompare(skip_values[1], -1));
  ck_assert(deltaCompare(skip_values[2], 0.5));
  ck_assert(deltaCompare(skip_values[3], 2));

  ret = vrna_sc_SHAPE_to_pr("L", linear_values, 7, 0.5);
  ck_assert_int_eq(ret, 1);
  ck_assert(deltaCompare(linear_values[1], 0.5));
  ck_assert(deltaCompare(linear_values[2], 0));
  ck_assert(deltaCompare(linear_values[3], 0.073529));
  ck_assert(deltaCompare(linear_values[4], 0.441176));
  ck_assert(deltaCompare(linear_values[5], 0.808824));
  ck_assert(deltaCompare(linear_values[6], 1));
  ck_assert(deltaCompare(linear_values[7], 1));

  ret = vrna_sc_SHAPE_to_pr("Ls0.5i0.1", linear_custom_values, 7, 0.5);
  ck_assert_int_eq(ret, 1);
  ck_assert(deltaCompare(linear_custom_values[1], 0.5));
  ck_assert(deltaCompare(linear_custom_values[2], 0));
  ck_assert(deltaCompare(linear_custom_values[3], 0.3));
  ck_assert(deltaCompare(linear_custom_values[4], 0.8));
  ck_assert(deltaCompare(linear_custom_values[5], 1));
  ck_assert(deltaCompare(linear_custom_values[6], 1));
  ck_assert(deltaCompare(linear_custom_values[7], 1));

  ret = vrna_sc_SHAPE_to_pr("O", log_values, 7, 0.5);
  ck_assert_int_eq(ret, 1);
  ck_assert(deltaCompare(log_values[1], 0.5));
  ck_assert(deltaCompare(log_values[2], 0));
  ck_assert(deltaCompare(log_values[3], 0.564816));
  ck_assert(deltaCompare(log_values[4], 0.998033));
  ck_assert(deltaCompare(log_values[5], 1));
  ck_assert(deltaCompare(log_values[6], 1));
  ck_assert(deltaCompare(log_values[7], 1));

  ret = vrna_sc_SHAPE_to_pr("Os1.5i-2", log_values_custom, 7, 0.5);
  ck_assert_int_eq(ret, 1);
  ck_assert(deltaCompare(log_values_custom[1], 0.5));
  ck_assert(deltaCompare(log_values_custom[2], 0));
  ck_assert(deltaCompare(log_values_custom[3], 0.409137));
  ck_assert(deltaCompare(log_values_custom[4], 0.871235));
  ck_assert(deltaCompare(log_values_custom[5], 1));
  ck_assert(deltaCompare(log_values_custom[6], 1));
  ck_assert(deltaCompare(log_values_custom[7], 1));

#test test_vrna_file_SHAPE_read
  char tempfile[L_tmpnam + 1];
  const size_t len = 5;
  char sequence[len];
  double values[len];
  int ret;

  //1 entry
  writeTempFile(tempfile, "1 A 0.5\n");
  ret = vrna_file_SHAPE_read(tempfile, 1, 0, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "A");
  ck_assert(deltaCompare(values[1], 0.5));
  unlink(tempfile);

  //1 entry 2 columns
  writeTempFile(tempfile, "1 0.5\n");
  ret = vrna_file_SHAPE_read(tempfile, 1, 0, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "N");
  ck_assert(deltaCompare(values[1], 0.5));
  unlink(tempfile);

  //multiple entries
  writeTempFile(tempfile, "1 A 0.1\n2 T 0.2\n3 G 0.3\n4 C 0.4\n");
  ret = vrna_file_SHAPE_read(tempfile, 4, 0, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "ATGC");
  ck_assert(deltaCompare(values[1], 0.1));
  ck_assert(deltaCompare(values[2], 0.2));
  ck_assert(deltaCompare(values[3], 0.3));
  ck_assert(deltaCompare(values[4], 0.4));
  unlink(tempfile);

  //value formats
  writeTempFile(tempfile, "1 A 1\n2 T 2.\n3 G .3\n4 C 1e-1\n");
  ret = vrna_file_SHAPE_read(tempfile, 4, 0, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "ATGC");
  ck_assert(deltaCompare(values[1], 1));
  ck_assert(deltaCompare(values[2], 2));
  ck_assert(deltaCompare(values[3], 0.3));
  ck_assert(deltaCompare(values[4], 0.1));
  unlink(tempfile);

  //whitespaces
  writeTempFile(tempfile, "1 \t 0.5\n2    A\t1\n");
  ret = vrna_file_SHAPE_read(tempfile, 2, 0, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "NA");
  ck_assert(deltaCompare(values[1], 0.5));
  ck_assert(deltaCompare(values[2], 1));
  unlink(tempfile);

  //missing value
  writeTempFile(tempfile, "1 A\n");
  ret = vrna_file_SHAPE_read(tempfile, 1, 123, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "A");
  ck_assert(deltaCompare(values[1], 123));
  unlink(tempfile);

  //missing nucleotide & value
  writeTempFile(tempfile, "1\n");
  ret = vrna_file_SHAPE_read(tempfile, 1, 123, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "N");
  ck_assert(deltaCompare(values[1], 123));
  unlink(tempfile);

  //upper limit
  writeTempFile(tempfile, "1 A 0.1\n2 T 0.2\n3 G 0.3\n4 C 0.4\n");
  ret = vrna_file_SHAPE_read(tempfile, 3, 0, sequence, values);
  ck_assert_int_eq(ret, 0);
  unlink(tempfile);

  //lower limit
  writeTempFile(tempfile, "0 A 0.1\n1 T 0.2\n2 G 0.3\n3 C 0.4\n");
  ret = vrna_file_SHAPE_read(tempfile, 4, 0, sequence, values);
  ck_assert_int_eq(ret, 0);
  unlink(tempfile);

  //unordered
  writeTempFile(tempfile, "3 G 0.3\n2 T 0.2\n4 C 0.4\n1 A 0.1\n");
  ret = vrna_file_SHAPE_read(tempfile, 4, 0, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "ATGC");
  ck_assert(deltaCompare(values[1], 0.1));
  ck_assert(deltaCompare(values[2], 0.2));
  ck_assert(deltaCompare(values[3], 0.3));
  ck_assert(deltaCompare(values[4], 0.4));
  unlink(tempfile);

  //missing indices middle
  writeTempFile(tempfile, "1 A 0.1\n4 C 0.4\n");
  ret = vrna_file_SHAPE_read(tempfile, 4, 123, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "ANNC");
  ck_assert(deltaCompare(values[1], 0.1));
  ck_assert(deltaCompare(values[2], 123));
  ck_assert(deltaCompare(values[3], 123));
  ck_assert(deltaCompare(values[4], 0.4));
  unlink(tempfile);

  //missing indices start end
  writeTempFile(tempfile, "2 T 0.2\n3 G 0.3\n");
  ret = vrna_file_SHAPE_read(tempfile, 4, 123, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "NTGN");
  ck_assert(deltaCompare(values[1], 123));
  ck_assert(deltaCompare(values[2], 0.2));
  ck_assert(deltaCompare(values[3], 0.3));
  ck_assert(deltaCompare(values[4], 123));
  unlink(tempfile);

  //invalid file
  ret = vrna_file_SHAPE_read(NULL, 0, 0, sequence, values);
  ck_assert_int_eq(ret, 0);

  //missing linebreak
  writeTempFile(tempfile, "1 A 0.5");
  ret = vrna_file_SHAPE_read(tempfile, 1, 0, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "A");
  ck_assert(deltaCompare(values[1], 0.5));
  unlink(tempfile);

  //garbage + entry
  writeTempFile(tempfile, "\nblablabla\n#evil_comment 123\n1 A 0.5\n\ngarbage\n");
  ret = vrna_file_SHAPE_read(tempfile, 1, 0, sequence, values);
  ck_assert_int_eq(ret, 1);
  ck_assert_str_eq(sequence, "A");
  ck_assert(deltaCompare(values[1], 0.5));
  unlink(tempfile);

  //garbage only
  writeTempFile(tempfile, "\nblablabla\n\ngarbage\n");
  ret = vrna_file_SHAPE_read(tempfile, 1, 0, sequence, values);
  ck_assert_int_eq(ret, 0);
  unlink(tempfile);

  //empty file
  writeTempFile(tempfile, "");
  ret = vrna_file_SHAPE_read(tempfile, 1, 0, sequence, values);
  ck_assert_int_eq(ret, 0);
  unlink(tempfile);

#test test_vrna_sc_SHAPE_parse_method
  float p1, p2;
  char method;
  int ret;

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method(NULL, &method, &p1, &p2);
  ck_assert_int_eq(ret, 0);
  ck_assert_int_eq(method, 0);
  ck_assert(deltaCompare(p1, 0));
  ck_assert(deltaCompare(p2, 0));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("", &method, &p1, &p2);
  ck_assert_int_eq(ret, 0);
  ck_assert_int_eq(method, 0);
  ck_assert(deltaCompare(p1, 0));
  ck_assert(deltaCompare(p2, 0));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("X", &method, &p1, &p2);
  ck_assert_int_eq(ret, 0);
  ck_assert_int_eq(method, 0);
  ck_assert(deltaCompare(p1, 0));
  ck_assert(deltaCompare(p2, 0));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("D", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'D');
  ck_assert(deltaCompare(p1, 1.8));
  ck_assert(deltaCompare(p2, -0.6));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Dm", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'D');
  ck_assert(deltaCompare(p1, 1.8));
  ck_assert(deltaCompare(p2, -0.6));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Db", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'D');
  ck_assert(deltaCompare(p1, 1.8));
  ck_assert(deltaCompare(p2, -0.6));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Dmb", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'D');
  ck_assert(deltaCompare(p1, 1.8));
  ck_assert(deltaCompare(p2, -0.6));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Dm3b4", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'D');
  ck_assert(deltaCompare(p1, 3));
  ck_assert(deltaCompare(p2, 4));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Dm3.4b4.5", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'D');
  ck_assert(deltaCompare(p1, 3.4));
  ck_assert(deltaCompare(p2, 4.5));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Dm3.4", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'D');
  ck_assert(deltaCompare(p1, 3.4));
  ck_assert(deltaCompare(p2, -0.6));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Db4.5", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'D');
  ck_assert(deltaCompare(p1, 1.8));
  ck_assert(deltaCompare(p2, 4.5));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Z", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'Z');
  ck_assert(deltaCompare(p1, 0.89));
  ck_assert(deltaCompare(p2, 0));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Zb4.5", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'Z');
  ck_assert(deltaCompare(p1, 4.5));
  ck_assert(deltaCompare(p2, 0));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Zx", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'Z');
  ck_assert(deltaCompare(p1, 0.89));
  ck_assert(deltaCompare(p2, 0));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("W", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'W');
  ck_assert(deltaCompare(p1, 0));
  ck_assert(deltaCompare(p2, 0));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Wb4.5", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'W');
  ck_assert(deltaCompare(p1, 0));
  ck_assert(deltaCompare(p2, 0));

  p1 = p2 = method = 0;
  ret = vrna_sc_SHAPE_parse_method("Wx", &method, &p1, &p2);
  ck_assert_int_eq(ret, 1);
  ck_assert_int_eq(method, 'W');
  ck_assert(deltaCompare(p1, 0));
  ck_assert(deltaCompare(p2, 0));
