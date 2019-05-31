#include <stdio.h>
#include <stdlib.h>
#include <ViennaRNA/datastructures/hash_tables.h>
#include <stdarg.h>


static unsigned
hash_function_test(void           *hash_entry,
                   unsigned long  hashtable_size)
{
  return 0; //*((unsigned int*)hash_entry);
}


static int
hash_comparison_test(void *x,
                     void *y)
{
  unsigned int  *hem_x  = ((unsigned int *)x);
  unsigned int  *hem_y  = ((unsigned int *)y);

  if ((x == NULL) ^ (y == NULL))
    return 1;

  return !(*hem_x == *hem_y);
}


static int
free_dummy(void *x)
{
  return 0;
}


#suite Hash Table

#test test_vrna_hash_table
{
  //test hash table.
  vrna_callback_ht_free_entry       *my_free          = free_dummy;
  vrna_callback_ht_compare_entries  *my_comparison    = hash_comparison_test;
  vrna_callback_ht_hash_function    *my_hash_function = hash_function_test;
  vrna_hash_table_t                 ht                = vrna_ht_init(27,
                                                                     my_comparison,
                                                                     my_hash_function,
                                                                     my_free);

  unsigned int                      val = 5;
  int                               res = vrna_ht_insert(ht, (void *)&val);
  //printf("insert res: %d\n",res);
  ck_assert_int_eq(res, 0);

  unsigned int                      val2 = 3;
  res = vrna_ht_insert(ht, (void *)&val2);
  ck_assert_int_eq(res, 0);

  unsigned int                      *res_p = vrna_ht_get(ht, (void *)&val);
  //printf("get res5: %d\n",*res_p);
  ck_assert_int_eq(*res_p, 5);

  res_p = vrna_ht_get(ht, (void *)&val2);
  //printf("get res3: %d\n",*res_p);
  ck_assert_int_eq(*res_p, 3);

  vrna_ht_remove(ht, (void *)&val);
  res_p = vrna_ht_get(ht, (void *)&val);
  ck_assert_ptr_eq(res_p, NULL);

  res_p = vrna_ht_get(ht, (void *)&val2);
  ck_assert_int_eq(*res_p, 3);

  vrna_ht_remove(ht, (void *)&val2);
  res_p = vrna_ht_get(ht, (void *)&val2);
  ck_assert_ptr_eq(res_p, NULL);

  //vrna_ht_clear(ht);
  vrna_ht_free(ht);
}
