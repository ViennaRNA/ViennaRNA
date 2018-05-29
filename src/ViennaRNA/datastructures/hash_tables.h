#ifndef VIENNA_RNA_PACKAGE_HASH_UTIL_H
#define VIENNA_RNA_PACKAGE_HASH_UTIL_H

/* Taken from the barriers tool and modified by GE. */

/**
 *  @file ViennaRNA/datastructures/hash_tables.h
 *  @ingroup  utils
 *  @brief  Implementations of hash table functions
 */


/**
 *  @addtogroup utils
 *  @{
 */

typedef int (vrna_hash_entry_comparison)(void *hash_entry_a,
                                         void *hash_entry_b);


typedef unsigned (vrna_hash_function)(void          *hash_entry,
                                      unsigned long hashtable_size);


typedef int (vrna_free_hash_entry)(void *hash_entry);


typedef struct vrna_hash_table_s *vrna_hash_table_t;

/**
 *  @brief  Get an initialized hash table
 *
 *  This function returns a ready-to-use hash table with pre-allocated
 *  memory for a particular number of entries.
 *
 *  @warning  If hash_bits is larger than 27 you have to compile it with
 *            the flag gcc -mcmodel=large.
 *
 *  @param  hash_bits   Number of bits for the hash table. This determines the size (2^hash_bits -1).
 *  @param  compare_function  A function pointer to compare any two entries in the hash table
 *  @param  hash_function     A function pointer to retrieve the hash value of any entry
 *  @param  free_hash_entry   A function pointer to free the memory occupied by any entry
 *  @return                   An initialized, empty hash table
 */
vrna_hash_table_t
vrna_hash_init(unsigned int               hash_bits,
               vrna_hash_entry_comparison *compare_function,
               vrna_hash_function         *hash_function,
               vrna_free_hash_entry       *free_hash_entry);


/**
 *  @brief  Get an element from the hash table
 *
 *  This function takes an object @p x and performs a look-up whether
 *  the object is stored within the hash table @p ht. If the object is
 *  already stored in @p ht, the function simply returns the entry,
 *  otherwise it returns @p NULL.
 *
 *  @see vrna_hash_insert(), vrna_hash_delete(), vrna_hash_init()
 *
 *  @param  ht  The hash table
 *  @param  x   The hash entry to look-up
 *  @return     The entry @p x if it is stored in @p ht, @p NULL otherwise
 */
void *
vrna_hash_get(vrna_hash_table_t ht,
              void              *x);


/**
 *  @brief  Insert an object into a hash table
 *
 *  Writes the pointer to your hash entry into the table.
 *
 *  @see vrna_hash_init(), vrna_hash_delete(), vrna_hash_clear()
 *
 *  @param  ht  The hash table
 *  @param  x   The hash entry
 *  @return     0 on success, 1 if the value is already in the hash table, -1 on error.
 */
int
vrna_hash_insert(vrna_hash_table_t  ht,
                 void               *x);


/**
 *  @brief  Remove an object from the hash table
 *
 *  Deletes the pointer to your hash entry from the table.
 *
 *  @note This function doesn't free any memory occupied by
 *        the hash entry.
 *
 *  @param  ht  The hash table
 *  @param  x   The hash entry
 */
void
vrna_hash_remove(vrna_hash_table_t  ht,
                 void               *x);


/**
 *  @brief  Clear the hash table
 *
 *  This function removes all entries from the hash table and
 *  automatically free's the memory occupied by each entry using
 *  the bound @vrna_free_hash_entry() function.
 *
 *  @see vrna_hash_free(), vrna_hash_init()
 *
 *  @param  ht  The hash table
 */
void
vrna_hash_clear(vrna_hash_table_t ht);


/**
 *  @brief  Free all memory occupied by the hash table
 *
 *  This function removes all entries from the hash table
 *  by calling the #vrna_free_hash_entry() function for each
 *  entry. Finally, the memory occupied by the hash table itself
 *  is free'd as well.
 *
 *  @param  ht  The hash table
 */
void
vrna_hash_free(vrna_hash_table_t ht);


/* modify hash_f(), hash_comp() and the typedef of hash_entry to suit your application */
int vrna_standard_hash_comparison(void  *x,
                                  void  *y);


unsigned vrna_standard_hash_function(void           *x,
                                     unsigned long  hashtable_size);


int vrna_standard_free_hash_entry(void *hash_entry);


typedef struct _hash_entry {
  char  *structure;
  float energy;
} vrna_standard_hash_entry;


/**
 *  @}
 */

#endif
