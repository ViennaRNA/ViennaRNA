/* Taken from the barriers tool and modified by GE. */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/hash_tables.h"


struct vrna_hash_table_s {
  unsigned int                      hash_bits;
  unsigned long                     Hash_size;
  /* must be power of 2^hash_bits -1 (example: HASHSIZE 67108864 -1 = 2^26 -1 )*/
  void                              **Hash_table;
  unsigned long                     Collisions;
  vrna_ht_cmp_f  Compare_function;
  vrna_ht_hashfunc_f    Hash_function;
  vrna_ht_free_f       Free_hash_entry;
};


typedef struct vrna_hash_entry_list_s {
  unsigned long length;
  unsigned long allocated_length;
  void          **hash_entries;
} vrna_hash_entry_list_t;

/* ----------------------------------------------------------------- */

PUBLIC struct vrna_hash_table_s *
vrna_ht_init(unsigned int                     hash_bits,
             vrna_ht_cmp_f compare_function,
             vrna_ht_hashfunc_f   hash_function,
             vrna_ht_free_f      free_hash_entry)
{
  struct vrna_hash_table_s *ht = NULL;

  if (hash_bits > 0) {
    ht = (struct vrna_hash_table_s *)vrna_alloc(sizeof(struct vrna_hash_table_s));

    ht->hash_bits = hash_bits;
    /* must be power of 2^hash_bits -1 (example: HASHSIZE 67108864 -1 = 2^26 -1 )
    * if we want to use '&' instead of modulo '%' for limiting the hash value. */
    ht->Hash_size   = (((unsigned long)1 << hash_bits) - 1);
    ht->Hash_table  = calloc((ht->Hash_size + 1), sizeof(void *));
    if (!ht->Hash_table) {
      fprintf(stderr, "Error: could not allocate space for the hash table!\n");
      free(ht);
      return NULL;
    }

    ht->Collisions = 0;

    if ((!compare_function) &&
        (!hash_function) &&
        (!free_hash_entry)) {
      /*
       *  Fall-back to expect dot-bracket structure string and
       *  free energy value as entries in hash table, i.e. pointers
       *  to vrna_ht_entry_db_t
       */
      ht->Compare_function  = &vrna_ht_db_comp;
      ht->Hash_function     = &vrna_ht_db_hash_func;
      ht->Free_hash_entry   = &vrna_ht_db_free_entry;
    } else if ((compare_function) &&
               (hash_function) &&
               (free_hash_entry)) {
      /* Bind user-defined compare, free, and hash functions */
      ht->Compare_function  = compare_function;
      ht->Hash_function     = hash_function;
      ht->Free_hash_entry   = free_hash_entry;
    } else {
      /*
       *  One of the function pointers is missing, so we don't initialize
       *  anything!
       */
      free(ht);
      ht = NULL;
    }
  }

  return ht;
}


unsigned long
vrna_ht_size(struct vrna_hash_table_s *ht)
{
  if (ht)
    return ht->Hash_size;

  return 0L;
}


unsigned long
vrna_ht_collisions(struct vrna_hash_table_s *ht)
{
  if (ht)
    return ht->Collisions;

  return 0L;
}


PUBLIC void *
vrna_ht_get(struct vrna_hash_table_s  *ht,
            void                      *x)             /* returns NULL unless x is in the hash */
{
  unsigned int hashval;

  if ((ht) && (x)) {
    hashval = ht->Hash_function(x, ht->Hash_size);
    if (hashval >= ht->Hash_size) {
      fprintf(stderr,
              "Error: hash function returns a value that is larger than the size of the hash map!\n");
      return NULL;
    }

    vrna_hash_entry_list_t *entries = ht->Hash_table[hashval];
    if (entries) {
      size_t i;
      for (i = 0; i < entries->length; i++) {
        if (ht->Compare_function(x, entries->hash_entries[i]) == 0)
          //found the same entry
          return entries->hash_entries[i]; /* success */
      }
    }
  }

  return NULL;
}


/* ----------------------------------------------------------------- */

PUBLIC int
vrna_ht_insert(struct vrna_hash_table_s *ht,
               void                     *x)         /* returns 1 if x already was in the hash */
{
  unsigned int hashval;

  if ((ht) && (x)) {
    hashval = ht->Hash_function(x, ht->Hash_size);
    if (hashval >= ht->Hash_size) {
      fprintf(stderr,
              "Error: hash function returns a value that is larger than the size of the hash map!\n");
      return -1; /* failure */;
    } else {
      vrna_hash_entry_list_t *entries = (vrna_hash_entry_list_t *)ht->Hash_table[hashval];
      if (entries) {
        //at first test if the entry is already in the list.
        size_t i;
        for (i = 0; i < entries->length; i++) {
          if (ht->Compare_function(x, entries->hash_entries[i]) == 0)
            //found the same entry
            return 0; /* success */
        }

        ht->Collisions++;
        //append the entry to the list with same hash values.
        if (i >= entries->length) {
          if (i >= entries->allocated_length) {
            //we have to extend the list.
            entries->allocated_length += 100;
            entries->hash_entries     = vrna_realloc(entries->hash_entries,
                                                     entries->allocated_length * sizeof(void *));
          }

          entries->hash_entries[entries->length] = x;
          entries->length++;
        }
      } else {
        //allocate new list and insert the value
        entries                   = malloc(sizeof(vrna_hash_entry_list_t));
        entries->allocated_length = 2;
        entries->hash_entries     = vrna_alloc(entries->allocated_length * sizeof(void *));
        entries->hash_entries[0]  = x;
        entries->length           = 1;
        ht->Hash_table[hashval]   = entries;
      }

      return 0; /* success */
    }
  }

  return -1; /* failure */
}


PUBLIC void
vrna_ht_clear(struct vrna_hash_table_s *ht)
{
  unsigned int k;

  if (ht) {
    for (k = 0; k < ht->Hash_size + 1; k++) {
      vrna_hash_entry_list_t *entries = ht->Hash_table[k];
      if (entries) {
        size_t i;
        for (i = 0; i < entries->length; i++) {
          ht->Free_hash_entry(entries->hash_entries[i]);
          entries->hash_entries[i] = NULL;
        }
        free(entries->hash_entries);
        free(entries);
      }
    }

    ht->Collisions = 0;
  }
}


PUBLIC void
vrna_ht_free(struct vrna_hash_table_s *ht)
{
  if (ht) {
    vrna_ht_clear(ht);
    free(ht->Hash_table);
    free(ht);
  }
}


/* ----------------------------------------------------------------- */

PUBLIC void
vrna_ht_remove(struct vrna_hash_table_s *ht,
               void                     *x)
{
  /* doesn't free anything ! */
  unsigned int hashval;

  if ((ht) && (x)) {
    hashval = ht->Hash_function(x, ht->Hash_size);
    if (hashval >= ht->Hash_size) {
      fprintf(stderr,
              "Error: hash function returns a value that is larger than the size of the hash map!\n");
      return;
    }

    vrna_hash_entry_list_t *entries = ht->Hash_table[hashval];
    if (entries) {
      size_t i;
      for (i = 0; i < entries->length; i++) {
        if (ht->Compare_function(x, entries->hash_entries[i]) == 0) {
          //found the same entry --> shift the list to the left in order to delete the value.
          size_t size_rest = entries->length - i - 1;
          if (size_rest == 0) {
            entries->hash_entries[i] = 0;
          } else {
            void  *offset     = entries->hash_entries + i;
            void  *next_entry = entries->hash_entries + i + 1;
            memcpy(offset, next_entry, size_rest * sizeof(void *));
          }

          entries->hash_entries[entries->length - 1] = NULL;
          entries->length--;
          return; /* success */
        }
      }
    }
  }
}


/* ----------------------------------------------------------------- */

/*
 * --------------------------------------------------------------------
 * mix -- mix 3 32-bit values reversibly.
 * For every delta with one or two bits set, and the deltas of all three
 * high bits or all three low bits, whether the original value of a,b,c
 * is almost all zero or is uniformly distributed,
 * If mix() is run forward or backward, at least 32 bits in a,b,c
 * have at least 1/4 probability of changing.
 * If mix() is run forward, every bit of c will change between 1/3 and
 * 2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
 * mix() takes 36 machine instructions, but only 18 cycles on a superscalar
 * machine (like a Pentium or a Sparc).  No faster mixer seems to work,
 * that's the result of my brute-force search.  There were about 2^^68
 * hashes to choose from.  I only tested about a billion of those.
 * --------------------------------------------------------------------
 */
#define mix(a, b, c) \
  { \
    a -= b; a -= c; a ^= (c >> 13); \
    b -= c; b -= a; b ^= (a << 8); \
    c -= a; c -= b; c ^= (b >> 13); \
    a -= b; a -= c; a ^= (c >> 12);  \
    b -= c; b -= a; b ^= (a << 16); \
    c -= a; c -= b; c ^= (b >> 5); \
    a -= b; a -= c; a ^= (c >> 3);  \
    b -= c; b -= a; b ^= (a << 10); \
    c -= a; c -= b; c ^= (b >> 15); \
  }

/*
 * --------------------------------------------------------------------
 * hash() -- hash a variable-length key into a 32-bit value
 * k       : the key (the unaligned variable-length array of bytes)
 * len     : the length of the key, counting by bytes
 * initval : can be any 4-byte value
 * Returns a 32-bit value.  Every bit of the key affects every bit of
 * the return value.  Every 1-bit and 2-bit delta achieves avalanche.
 * About 6*len+35 instructions.
 *
 * The best hash table sizes are powers of 2.  There is no need to do
 * mod a prime (mod is sooo slow!).  If you need less than 32 bits,
 * use a bitmask.  For example, if you need only 10 bits, do
 * h = (h & hashmask(10));
 * In which case, the hash table should have hashsize(10) elements.
 *
 * If you are hashing n strings (char **)k, do it like this:
 * for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);
 *
 * By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
 * code any way you wish, private, educational, or commercial.  It's free.
 *
 * See http://burtleburtle.net/bob/hash/evahash.html
 * Use for hash table lookup, or anything where one collision in 2^^32 is
 * acceptable.  Do NOT use for cryptographic purposes.
 * --------------------------------------------------------------------
 */
PUBLIC unsigned int
vrna_ht_db_hash_func(void           *x,
                     unsigned long  hashtable_size)
{
  register char           *k;           /* the key */
  register unsigned int   length;       /* the length of the key */
  register unsigned int   initval = 0;  /* the previous hash, or an arbitrary value */
  register unsigned int   a, b, c, len;

  /* Set up the internal state */
  k   = ((vrna_ht_entry_db_t *)x)->structure;
  len = length = strlen(k);
  a   = b = 0x9e3779b9; /* the golden ratio; an arbitrary value */
  c   = initval;        /* the previous hash value */

  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a +=
      (k[0] + ((unsigned int)k[1] << 8) + ((unsigned int)k[2] << 16) + ((unsigned int)k[3] << 24));
    b +=
      (k[4] + ((unsigned int)k[5] << 8) + ((unsigned int)k[6] << 16) + ((unsigned int)k[7] << 24));
    c +=
      (k[8] + ((unsigned int)k[9] << 8) + ((unsigned int)k[10] << 16) +
       ((unsigned int)k[11] << 24));
    mix(a, b, c);
    k   += 12;
    len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch (len) {
    /* all the case statements fall through */
    case 11:
      c += ((unsigned int)k[10] << 24);
      /* fall through */
    case 10:
      c += ((unsigned int)k[9] << 16);
      /* fall through */
    case 9:
      c += ((unsigned int)k[8] << 8);
    /* the first byte of c is reserved for the length */
      /* fall through */
    case 8:
      b += ((unsigned int)k[7] << 24);
      /* fall through */
    case 7:
      b += ((unsigned int)k[6] << 16);
      /* fall through */
    case 6:
      b += ((unsigned int)k[5] << 8);
      /* fall through */
    case 5:
      b += k[4];
      /* fall through */
    case 4:
      a += ((unsigned int)k[3] << 24);
      /* fall through */
    case 3:
      a += ((unsigned int)k[2] << 16);
      /* fall through */
    case 2:
      a += ((unsigned int)k[1] << 8);
      /* fall through */
    case 1:
      a += k[0];
      /* case 0: nothing left to add */
      break;
  }
  mix(a, b, c);
  /*-------------------------------------------- report the result */
  return c % hashtable_size;
}


/* ----------------------------------------------------------------- */
PUBLIC int
vrna_ht_db_comp(void  *x,
                void  *y)
{
  return strcmp(((vrna_ht_entry_db_t *)x)->structure,
                ((vrna_ht_entry_db_t *)y)->structure);
}


PUBLIC int
vrna_ht_db_free_entry(void *hash_entry)
{
  free(((vrna_ht_entry_db_t *)hash_entry)->structure);
  return 0;
}
