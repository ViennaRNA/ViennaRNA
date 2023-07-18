#ifndef VIENNA_RNA_PACKAGE_HASH_UTIL_H
#define VIENNA_RNA_PACKAGE_HASH_UTIL_H

/* Taken from the barriers tool and modified by GE. */

#ifdef VRNA_WARN_DEPRECATED
# if defined(DEPRECATED)
#   undef DEPRECATED
# endif
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

/**
 *  @file ViennaRNA/datastructures/hash_tables.h
 *  @ingroup  hash_table_utils
 *  @brief  Implementations of hash table functions
 */


/**
 *  @addtogroup hash_table_utils
 *  @{
 */

/**
 *  @name Abstract interface
 *  @{
 */

/**
 *  @brief  A hash table object
 *
 *  @see  vrna_ht_init(), vrna_ht_free()
 */
typedef struct vrna_hash_table_s *vrna_hash_table_t;


/**
 *  @brief  Callback function to compare two hash table entries
 *
 *  @see    vrna_ht_init(), vrna_ht_db_comp()
 *
 *  @param  x   A hash table entry
 *  @param  y   A hash table entry
 *  @return     -1 if x is smaller, +1 if x is larger than y. 0 if @f$x == y @f$
 */
typedef int (*vrna_ht_cmp_f)(void *x,
                                               void *y);

DEPRECATED(typedef int (vrna_callback_ht_compare_entries)(void *x,
                                               void *y),
           "Use vrna_ht_cmp_f instead!");



/**
 *  @brief  Callback function to generate a hash key, i.e. hash function
 *
 *  @see    vrna_ht_init(), vrna_ht_db_hash_func()
 *
 *  @param  x               A hash table entry
 *  @param  hashtable_size  The size of the hash table
 *  @return                 The hash table key for entry @p x
 */
typedef unsigned int (*vrna_ht_hashfunc_f)(void          *x,
                                                      unsigned long hashtable_size);

DEPRECATED(typedef unsigned int (vrna_callback_ht_hash_function)(void          *x,
                                                      unsigned long hashtable_size),
          "Use vrna_ht_hashfunc_f instead!");


/**
 *  @brief  Callback function to free a hash table entry
 *
 *  @see    vrna_ht_init(), vrna_ht_db_free_entry()
 *
 *  @param  x   A hash table entry
 *  @return     0 on success
 */
typedef int (*vrna_ht_free_f)(void *x);

DEPRECATED(typedef int (vrna_callback_ht_free_entry)(void *x),
           "Use vrna_ht_free_f instead!");


/**
 *  @brief  Get an initialized hash table
 *
 *  This function returns a ready-to-use hash table with pre-allocated
 *  memory for a particular number of entries.
 *
 *  @note
 *  @parblock
 *  If all function pointers are @p NULL, this function initializes
 *  the hash table with <em>default functions</em>, i.e.
 *  - vrna_ht_db_comp() for the @p compare_function,
 *  - vrna_ht_db_hash_func() for the @p hash_function, and
 *  - vrna_ht_db_free_entry() for the @p free_hash_entry
 *
 *  arguments.
 *  @endparblock
 *
 *  @warning  If @p hash_bits is larger than 27 you have to compile it with
 *            the flag gcc -mcmodel=large.
 *
 *  @param  b                 Number of bits for the hash table. This determines the size (@f$2^b -1@f$).
 *  @param  compare_function  A function pointer to compare any two entries in the hash table (may be @p NULL)
 *  @param  hash_function     A function pointer to retrieve the hash value of any entry (may be @p NULL)
 *  @param  free_hash_entry   A function pointer to free the memory occupied by any entry (may be @p NULL)
 *  @return                   An initialized, empty hash table, or @p NULL on any error
 */
vrna_hash_table_t
vrna_ht_init(unsigned int                     b,
             vrna_ht_cmp_f compare_function,
             vrna_ht_hashfunc_f   hash_function,
             vrna_ht_free_f      free_hash_entry);


/**
 *  @brief  Get the size of the hash table
 *
 *  @param  ht  The hash table
 *  @return     The size of the hash table, i.e. the maximum number of entries
 */
unsigned long
vrna_ht_size(vrna_hash_table_t ht);


/**
 *  @brief  Get the number of collisions in the hash table
 *
 *  @param  ht  The hash table
 *  @return     The number of collisions in the hash table
 */
unsigned long
vrna_ht_collisions(struct vrna_hash_table_s *ht);


/**
 *  @brief  Get an element from the hash table
 *
 *  This function takes an object @p x and performs a look-up whether
 *  the object is stored within the hash table @p ht. If the object is
 *  already stored in @p ht, the function simply returns the entry,
 *  otherwise it returns @p NULL.
 *
 *  @see vrna_ht_insert(), vrna_hash_delete(), vrna_ht_init()
 *
 *  @param  ht  The hash table
 *  @param  x   The hash entry to look-up
 *  @return     The entry @p x if it is stored in @p ht, @p NULL otherwise
 */
void *
vrna_ht_get(vrna_hash_table_t ht,
            void              *x);


/**
 *  @brief  Insert an object into a hash table
 *
 *  Writes the pointer to your hash entry into the table.
 *
 *  @warning  In case of collisions, this function simply
 *            increments the hash key until a free entry in
 *            the hash table is found.
 *
 *  @see vrna_ht_init(), vrna_hash_delete(), vrna_ht_clear()
 *
 *  @param  ht  The hash table
 *  @param  x   The hash entry
 *  @return     0 on success, 1 if the value is already in the hash table, -1 on error.
 */
int
vrna_ht_insert(vrna_hash_table_t  ht,
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
vrna_ht_remove(vrna_hash_table_t  ht,
               void               *x);


/**
 *  @brief  Clear the hash table
 *
 *  This function removes all entries from the hash table and
 *  automatically free's the memory occupied by each entry using
 *  the bound #vrna_ht_free_f() function.
 *
 *  @see vrna_ht_free(), vrna_ht_init()
 *
 *  @param  ht  The hash table
 */
void
vrna_ht_clear(vrna_hash_table_t ht);


/**
 *  @brief  Free all memory occupied by the hash table
 *
 *  This function removes all entries from the hash table
 *  by calling the #vrna_ht_free_f() function for each
 *  entry. Finally, the memory occupied by the hash table itself
 *  is free'd as well.
 *
 *  @param  ht  The hash table
 */
void
vrna_ht_free(vrna_hash_table_t ht);


/* End of abstract interface */
/**@}*/

/**
 *  @name Dot-Bracket / Free Energy entries
 *  @{
 */

/**
 *  @brief  Default hash table entry
 *
 *  @see  vrna_ht_init(), vrna_ht_db_comp(), vrna_ht_db_hash_func(), vrna_ht_db_free_entry()
 */
typedef struct {
  char  *structure; /**< A secondary structure in dot-bracket notation */
  float energy;     /**< The free energy of @p structure */
} vrna_ht_entry_db_t;


/**
 *  @brief  Default hash table entry comparison
 *
 *  This is the default comparison function for hash table entries.
 *  It assumes the both entries @p x and @p y are of type #vrna_ht_entry_db_t
 *  and compares the @p structure attribute of both entries
 *
 *  @see #vrna_ht_entry_db_t, vrna_ht_init(), vrna_ht_db_hash_func(), vrna_ht_db_free_entry()
 *
 *  @param  x   A hash table entry of type #vrna_ht_entry_db_t
 *  @param  y   A hash table entry of type #vrna_ht_entry_db_t
 *  @return     -1 if x is smaller, +1 if x is larger than y. 0 if both are equal.
 */
int
vrna_ht_db_comp(void  *x,
                void  *y);


/**
 *  @brief  Default hash function
 *
 *  This is the default hash function for hash table insertion/lookup. It
 *  assumes that entries are of type #vrna_ht_entry_db_t and uses
 *  the Bob Jenkins 1996 mix function to create a hash key from the
 *  @p structure attribute of the hash entry.
 *
 *  @see  #vrna_ht_entry_db_t, vrna_ht_init(), vrna_ht_db_comp(), vrna_ht_db_free_entry()
 *
 *  @param  x               A hash table entry to compute the key for
 *  @param  hashtable_size  The size of the hash table
 *  @return                 The hash key for entry @p x
 */
unsigned int
vrna_ht_db_hash_func(void           *x,
                     unsigned long  hashtable_size);


/**
 *  @brief  Default function to free memory occupied by a hash entry
 *
 *  This function assumes that hash entries are of type #vrna_ht_entry_db_t
 *  and free's the memory occupied by that entry.
 *
 *  @see  #vrna_ht_entry_db_t, vrna_ht_init(), vrna_ht_db_comp(), vrna_ht_db_hash_func()
 *
 *  @param  hash_entry  The hash entry to remove from memory
 *  @return             0 on success
 */
int vrna_ht_db_free_entry(void *hash_entry);


/* End of dot-bracket interface */
/**@}*/

/**
 *  @}
 */

#endif
