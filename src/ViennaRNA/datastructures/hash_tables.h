/* Taken from the barriers tool and modified by GE. */


#ifndef VIENNA_RNA_PACKAGE_HASH_UTIL_H
#define VIENNA_RNA_PACKAGE_HASH_UTIL_H


typedef int (vrna_hash_entry_comparison)(void *hash_entry_a,
                                         void *hash_entry_b);
typedef unsigned (vrna_hash_function)(void          *hash_entry,
                                      unsigned long hashtable_size);
typedef int (vrna_free_hash_entry)(void *hash_entry);

typedef struct _hash_table {
  unsigned int                hash_bits;
  unsigned long               Hash_size;
  /* must be power of 2^hash_bits -1 (example: HASHSIZE 67108864 -1 = 2^26 -1 )*/
  void                        **Hash_table;
  unsigned long               Collisions;
  vrna_hash_entry_comparison  *Compare_function;
  vrna_hash_function          *Hash_function;
  vrna_free_hash_entry        *Free_hash_entry;
} vrna_hash_table;


/**
 * ! Initialize the hash table and allocate memory.
 * Warning: If hash_bits is larger than 27 you have to compile it with the flag gcc -mcmodel=large.
 * @param hash_bits - Number of bits for the hash table. This determines the size (2^hash_bits -1).
 */
vrna_hash_table *vrna_initialize_hash_table(unsigned int                hash_bits,
                                            vrna_hash_entry_comparison  *compare_function,
                                            vrna_hash_function          *hash_function,
                                            vrna_free_hash_entry        *free_hash_entry);


/**
 * ! returns the hash entry if it is in the table, else NULL.
 * @param x - the hash entry
 * @param ht - the hash table
 */
void *vrna_lookup_hash(void             *x,
                       vrna_hash_table  *ht);


/**
 * ! Writes the pointer to your hash entry into the table.
 * @param x - the hash entry
 * @param ht - the hash table
 * @return 0 if everything is ok, 1 if the value is already in the hash table.
 */
int vrna_write_hash(void            *x,
                    vrna_hash_table *ht);


/**
 * ! Deletes the pointer to your hash entry from the table.
 * @param x - the hash entry
 * @param ht - the hash table
 */
void vrna_delete_hash(void            *x,
                      vrna_hash_table *ht);


/**
 * ! Applies the free_hash_entry function to all hash entries within the table
 *   and frees the memory that has been allocated for the table.
 * @param ht - the hash table
 */
void vrna_free_hash_table(vrna_hash_table *ht);


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


#endif

/* End of file */
