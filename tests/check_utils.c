#include <check.h>

#include <ViennaRNA/utils.h>
#include <ViennaRNA/utils.c>

START_TEST(test_get_char_encoding)
{
  model_detailsT details = {0};

  ck_assert_int_eq(get_char_encoding('\0', &details), 0);
  ck_assert_int_eq(get_char_encoding('_', &details), 0);
  ck_assert_int_eq(get_char_encoding('A', &details), 1);
  ck_assert_int_eq(get_char_encoding('C', &details), 2);
  ck_assert_int_eq(get_char_encoding('G', &details), 3);
  ck_assert_int_eq(get_char_encoding('U', &details), 4);
  ck_assert_int_eq(get_char_encoding('T', &details), 4);
  ck_assert_int_eq(get_char_encoding('X', &details), 0);
  ck_assert_int_eq(get_char_encoding('K', &details), 0);
  ck_assert_int_eq(get_char_encoding('I', &details), 0);

  details.energy_set = 1;

  ck_assert_int_eq(get_char_encoding('\0', &details), -64);
  ck_assert_int_eq(get_char_encoding('_', &details), 31);
  ck_assert_int_eq(get_char_encoding('A', &details), 1);
  ck_assert_int_eq(get_char_encoding('C', &details), 3);
  ck_assert_int_eq(get_char_encoding('G', &details), 7);
  ck_assert_int_eq(get_char_encoding('U', &details), 21);
  ck_assert_int_eq(get_char_encoding('T', &details), 20);
  ck_assert_int_eq(get_char_encoding('X', &details), 24);
  ck_assert_int_eq(get_char_encoding('K', &details), 11);
  ck_assert_int_eq(get_char_encoding('I', &details), 9);
}
END_TEST

START_TEST(test_get_encoded_char)
{
  const char *characters = "_ACGUTXKI";
  const char *p;
  model_detailsT details = {0};

  ck_assert_int_eq(get_encoded_char(0, &details), '_');
  ck_assert_int_eq(get_encoded_char(1, &details), 'A');
  ck_assert_int_eq(get_encoded_char(2, &details), 'C');
  ck_assert_int_eq(get_encoded_char(3, &details), 'G');
  ck_assert_int_eq(get_encoded_char(4, &details), 'U');
  ck_assert_int_eq(get_encoded_char(5, &details), 'T');
  ck_assert_int_eq(get_encoded_char(6, &details), 'X');
  ck_assert_int_eq(get_encoded_char(7, &details), 'K');
  ck_assert_int_eq(get_encoded_char(8, &details), 'I');

  details.energy_set = 1;

  ck_assert_int_eq(get_encoded_char(-64, &details), '\0');
  ck_assert_int_eq(get_encoded_char(31, &details), '_');
  ck_assert_int_eq(get_encoded_char(1, &details), 'A');
  ck_assert_int_eq(get_encoded_char(3, &details), 'C');
  ck_assert_int_eq(get_encoded_char(7, &details), 'G');
  ck_assert_int_eq(get_encoded_char(21, &details), 'U');
  ck_assert_int_eq(get_encoded_char(20, &details), 'T');
  ck_assert_int_eq(get_encoded_char(24, &details), 'X');
  ck_assert_int_eq(get_encoded_char(11, &details), 'K');
  ck_assert_int_eq(get_encoded_char(9, &details), 'I');
}
END_TEST

START_TEST(test_get_sequence_encoding)
{
  model_detailsT details = {0};
  short *data;

  data = get_sequence_encoding("_AUGC", 0, &details);
  ck_assert_int_eq(data[0], 5); //sequence length
  ck_assert_int_eq(data[1], 0);
  ck_assert_int_eq(data[2], 1);
  ck_assert_int_eq(data[3], 4);
  ck_assert_int_eq(data[4], 3);
  ck_assert_int_eq(data[5], 2);
  ck_assert_int_eq(data[6], 0); //value stored in data[1]
  free(data);

  data = get_sequence_encoding("augc", 0, &details);
  ck_assert_int_eq(data[0], 4); //sequence length
  ck_assert_int_eq(data[1], 1);
  ck_assert_int_eq(data[2], 4);
  ck_assert_int_eq(data[3], 3);
  ck_assert_int_eq(data[4], 2);
  ck_assert_int_eq(data[5], 1); //value stored in data[1]
  free(data);

  //@TODO: type = 1
  //@TODO: details.energy_set = 1
}
END_TEST

START_TEST(test_fill_pair_matrices)
{
  model_detailsT details = {0};
  int i, j, count;

  fill_pair_matrices(&details);

  for (i = 0; i <= 4; ++i)
    ck_assert_int_eq(details.alias[i], i);

  ck_assert_int_eq(details.alias[5], 3);
  ck_assert_int_eq(details.alias[6], 2);

  for (i = 7; i <= MAXALPHA; ++i)
    ck_assert_int_eq(details.alias[i], 0);

  ck_assert_int_eq(details.pair[1][4], 5);
  ck_assert_int_eq(details.pair[1][7], 5);
  ck_assert_int_eq(details.pair[2][3], 1);
  ck_assert_int_eq(details.pair[3][2], 2);
  ck_assert_int_eq(details.pair[3][4], 3);
  ck_assert_int_eq(details.pair[4][1], 6);
  ck_assert_int_eq(details.pair[4][3], 4);
  ck_assert_int_eq(details.pair[4][7], 6);
  ck_assert_int_eq(details.pair[5][6], 2);
  ck_assert_int_eq(details.pair[6][5], 1);
  ck_assert_int_eq(details.pair[7][1], 6);
  ck_assert_int_eq(details.pair[7][4], 5);

  for (i = 0, count = 0; i < NBASES; ++i)
    for (j = 0; j < NBASES; ++j)
      if (details.pair[i][j] == 0)
        ++count;

  ck_assert_int_eq(count, 52);

  ck_assert_int_eq(details.rtype[0], 0);
  ck_assert_int_eq(details.rtype[1], 2);
  ck_assert_int_eq(details.rtype[2], 1);
  ck_assert_int_eq(details.rtype[3], 4);
  ck_assert_int_eq(details.rtype[4], 3);
  ck_assert_int_eq(details.rtype[5], 6);
  ck_assert_int_eq(details.rtype[6], 5);
  ck_assert_int_eq(details.rtype[7], 0);

  //@TODO: details.noGU = 1
  //@TODO: details.nonstandards
  //@TODO: details.energyset = [1, 2, 3]
}
END_TEST

START_TEST(test_get_ptypes)
{
  model_detailsT details = {0};
  const int len = 12;
  short sequence[] = {len, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1}; //ACGUACGUACGU
  int i, j, count;
  int *idx;
  char *ptype;

  int results[][3] = {{7, 2, 1},
                      {8, 1, 5},
                      {8, 3, 3},
                      {9, 4, 6},
                      {10, 3, 2},
                      {11, 2, 1},
                      {11, 4, 4},
                      {11, 6, 1},
                      {12, 1, 5},
                      {12, 3, 3},
                      {12, 5, 5},
                      {12, 7, 3},
                      {0, 0, 0}};

  fill_pair_matrices(&details);

  ptype = get_ptypes(sequence, &details, 0);
  idx = get_indx(len);

  for (i = 0; results[i][0]; ++i)
    ck_assert_int_eq(ptype[idx[results[i][0]] + results[i][1]], results[i][2]);

  for (j = 1, count = 0; j <= len; ++j)
    for (i = 1; i <= j; ++i)
      if (ptype[idx[j] + i] == 0)
        ++count;

  ck_assert_int_eq(count, 66);

  free(idx);
  free(ptype);

  //@TODO: extend alphabeth
  //@TODO: details.noLP = 1
  //@TODO: idx_type = 1
}
END_TEST

TCase* utils_testcase()
{
  TCase *tc = tcase_create("utils");
  tcase_add_test(tc, test_get_char_encoding);
  tcase_add_test(tc, test_get_encoded_char);
  tcase_add_test(tc, test_get_sequence_encoding);
  tcase_add_test(tc, test_fill_pair_matrices);
  tcase_add_test(tc, test_get_ptypes);

  return tc;
}
