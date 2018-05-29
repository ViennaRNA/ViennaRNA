/* Taken from the barriers tool and modified by GE. */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/hash_tables.h"


struct vrna_hash_table_s {
  unsigned int                hash_bits;
  unsigned long               Hash_size;
  /* must be power of 2^hash_bits -1 (example: HASHSIZE 67108864 -1 = 2^26 -1 )*/
  void                        **Hash_table;
  unsigned long               Collisions;
  vrna_hash_entry_comparison  *Compare_function;
  vrna_hash_function          *Hash_function;
  vrna_free_hash_entry        *Free_hash_entry;
};

/* ----------------------------------------------------------------- */

PUBLIC struct vrna_hash_table_s *
vrna_hash_init(unsigned int               hash_bits,
               vrna_hash_entry_comparison *compare_function,
               vrna_hash_function         *hash_function,
               vrna_free_hash_entry       *free_hash_entry)
{
  struct vrna_hash_table_s *ht = NULL;

  if (hash_bits > 0) {
    ht = (struct vrna_hash_table_s *)vrna_alloc(sizeof(struct vrna_hash_table_s));

    ht->hash_bits = hash_bits;
    /* must be power of 2^hash_bits -1 (example: HASHSIZE 67108864 -1 = 2^26 -1 )*/
    ht->Hash_size   = (((unsigned long)1 << hash_bits) - 1);
    ht->Hash_table  = calloc(ht->Hash_size + 1, sizeof(void *));
    ht->Collisions  = 0;

    if ((!compare_function) &&
        (!hash_function) &&
        (!free_hash_entry)) {
      /*
       *  Fall-back to expect dot-bracket structure string and
       *  free energy value as entries in hash table, i.e. pointers
       *  to vrna_hash_entry_db_t
       */
      ht->Compare_function  = &vrna_hash_db_comp;
      ht->Hash_function     = &vrna_hash_db_hash_func;
      ht->Free_hash_entry   = &vrna_hash_db_free_entry;
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


/* ----------------------------------------------------------------- */

/* ----------------------------------------------------------------- */

PUBLIC void *
vrna_hash_get(struct vrna_hash_table_s  *ht,
              void                      *x)           /* returns NULL unless x is in the hash */
{
  unsigned int hashval;

  if ((ht) && (x)) {
    hashval = ht->Hash_function(x, ht->Hash_size);
    if (ht->Hash_table[hashval] == NULL)
      return NULL;

    while (ht->Hash_table[hashval]) {
      if (ht->Compare_function(x, ht->Hash_table[hashval]) == 0)
        return ht->Hash_table[hashval];

      hashval = ((hashval + 1) & (ht->Hash_size));
    }
  }

  return NULL;
}


/* ----------------------------------------------------------------- */

PUBLIC int
vrna_hash_insert(struct vrna_hash_table_s *ht,
                 void                     *x)         /* returns 1 if x already was in the hash */
{
  unsigned int hashval;

  if ((ht) && (x)) {
    hashval = ht->Hash_function(x, ht->Hash_size);
    if (hashval < ht->Hash_size) {
      void *entry = ht->Hash_table[hashval];
      while (entry != NULL) {
        if (ht->Compare_function(x, ht->Hash_table[hashval]) == 0)
          return 1;

        hashval = ((hashval + 1) & (ht->Hash_size));
        ht->Collisions++;
      }
      ht->Hash_table[hashval] = x;
      return 0; /* success */
    } else {
      vrna_message_warning("vrna_hash_insert: "
                           "The hash table (size %d) is too small for entry with key %d",
                           ht->Hash_size,
                           hashval);
    }
  }

  return -1; /* failure */
}


PUBLIC void
vrna_hash_clear(struct vrna_hash_table_s *ht)
{
  unsigned int i;

  if (ht) {
    for (i = 0; i < ht->Hash_size + 1; i++) {
      if (ht->Hash_table[i]) {
        ht->Free_hash_entry(ht->Hash_table[i]);
        ht->Hash_table[i] = NULL;
      }
    }

    ht->Collisions = 0;
  }
}


PUBLIC void
vrna_hash_free(struct vrna_hash_table_s *ht)
{
  if (ht) {
    vrna_hash_clear(ht);
    free(ht->Hash_table);
    free(ht);
  }
}


/* ----------------------------------------------------------------- */

PUBLIC void
vrna_hash_remove(struct vrna_hash_table_s *ht,
                 void                     *x)         /* doesn't work in case of collisions */
{
  /* doesn't free anything ! */
  unsigned int hashval;

  if ((ht) && (x)) {
    hashval = ht->Hash_function(x, ht->Hash_size);
    while (ht->Hash_table[hashval]) {
      if (ht->Compare_function(x, ht->Hash_table[hashval]) == 0) {
        ht->Hash_table[hashval] = NULL;
        return;
      }

      hashval = ((hashval + 1) & (ht->Hash_size));
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
vrna_hash_db_hash_func(void           *x,
                       unsigned long  hashtable_size)
{
  register unsigned char  *k;           /* the key */
  register unsigned int   length;       /* the length of the key */
  register unsigned int   initval = 0;  /* the previous hash, or an arbitrary value */
  register unsigned int   a, b, c, len;

  /* Set up the internal state */
  k   = ((vrna_hash_entry_db_t *)x)->structure;
  len = length = (unsigned int)strlen(k);
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
    case 10:
      c += ((unsigned int)k[9] << 16);
    case 9:
      c += ((unsigned int)k[8] << 8);
    /* the first byte of c is reserved for the length */
    case 8:
      b += ((unsigned int)k[7] << 24);
    case 7:
      b += ((unsigned int)k[6] << 16);
    case 6:
      b += ((unsigned int)k[5] << 8);
    case 5:
      b += k[4];
    case 4:
      a += ((unsigned int)k[3] << 24);
    case 3:
      a += ((unsigned int)k[2] << 16);
    case 2:
      a += ((unsigned int)k[1] << 8);
    case 1:
      a += k[0];
      /* case 0: nothing left to add */
  }
  mix(a, b, c);
  /*-------------------------------------------- report the result */
  return c & hashtable_size;
}


/* ----------------------------------------------------------------- */
PUBLIC int
vrna_hash_db_comp(void  *x,
                  void  *y)
{
  return strcmp(((vrna_hash_entry_db_t *)x)->structure,
                ((vrna_hash_entry_db_t *)y)->structure);
}


PUBLIC int
vrna_hash_db_free_entry(void *hash_entry)
{
  free(((vrna_hash_entry_db_t *)hash_entry)->structure);
  return 0;
}
