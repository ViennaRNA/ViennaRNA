/* Taken from the barriers tool and modified by GE. */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ViennaRNA/datastructures/hash_tables.h"


/* ----------------------------------------------------------------- */

vrna_hash_table *
vrna_initialize_hash_table(unsigned int               hash_bits,
                           vrna_hash_entry_comparison *compare_function,
                           vrna_hash_function         *hash_function,
                           vrna_free_hash_entry       *free_hash_entry)
{
  vrna_hash_table *ht = malloc(sizeof(vrna_hash_table));

  ht->Compare_function  = compare_function;
  ht->Hash_function     = hash_function;
  ht->Free_hash_entry   = free_hash_entry;
  ht->hash_bits         = hash_bits;
  /* must be power of 2^hash_bits -1 (example: HASHSIZE 67108864 -1 = 2^26 -1 )*/
  ht->Hash_size   = (((unsigned long)1 << hash_bits) - 1);
  ht->Hash_table  = calloc(ht->Hash_size + 1, sizeof(void *));
  ht->Collisions  = 0;
  return ht;
}


/* ----------------------------------------------------------------- */

/* ----------------------------------------------------------------- */
int
vrna_standard_hash_comparison(void  *x,
                              void  *y)
{
  return strcmp(((vrna_standard_hash_entry *)x)->structure,
                ((vrna_standard_hash_entry *)y)->structure);
}


/* ----------------------------------------------------------------- */

void *
vrna_lookup_hash(void             *x,
                 vrna_hash_table  *ht)         /* returns NULL unless x is in the hash */
{
  unsigned int hashval;

  hashval = ht->Hash_function(x, ht->Hash_size);
  if (ht->Hash_table[hashval] == NULL)
    return NULL;

  while (ht->Hash_table[hashval]) {
    if (ht->Compare_function(x, ht->Hash_table[hashval]) == 0)
      return ht->Hash_table[hashval];

    hashval = ((hashval + 1) & (ht->Hash_size));
  }
  return NULL;
}


/* ----------------------------------------------------------------- */

int
vrna_write_hash(void            *x,
                vrna_hash_table *ht)          /* returns 1 if x already was in the hash */
{
  unsigned int hashval;

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
  } else {
    printf(stderr, "Error: The hash table is too small: %d %d\n", hashval, ht->Hash_size);
  }

  return 0;
}


void
vrna_free_hash_table(vrna_hash_table *ht)
{
  unsigned int i;

  for (i = 0; i < ht->Hash_size + 1; i++) {
    if (ht->Hash_table[i]) {
      ht->Free_hash_entry(ht->Hash_table[i]);
      ht->Hash_table[i] = NULL;
    }
  }
  free(ht->Hash_table);
}


/* ----------------------------------------------------------------- */

void
vrna_delete_hash(void             *x,
                 vrna_hash_table  *ht)         /* doesn't work in case of collisions */
{
  /* doesn't free anything ! */
  unsigned int hashval;

  hashval = ht->Hash_function(x, ht->Hash_size);
  while (ht->Hash_table[hashval]) {
    if (ht->Compare_function(x, ht->Hash_table[hashval]) == 0) {
      ht->Hash_table[hashval] = NULL;
      return;
    }

    hashval = ((hashval + 1) & (ht->Hash_size));
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
unsigned
vrna_standard_hash_function(void          *x,
                            unsigned long hashtable_size)
{
  register unsigned char  *k;           /* the key */
  register unsigned       length;       /* the length of the key */
  register unsigned       initval = 0;  /* the previous hash, or an arbitrary value */
  register unsigned       a, b, c, len;

  /* Set up the internal state */
  k   = ((vrna_standard_hash_entry *)x)->structure;
  len = length = (unsigned)strlen(k);
  a   = b = 0x9e3779b9; /* the golden ratio; an arbitrary value */
  c   = initval;        /* the previous hash value */

  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a += (k[0] + ((unsigned)k[1] << 8) + ((unsigned)k[2] << 16) + ((unsigned)k[3] << 24));
    b += (k[4] + ((unsigned)k[5] << 8) + ((unsigned)k[6] << 16) + ((unsigned)k[7] << 24));
    c += (k[8] + ((unsigned)k[9] << 8) + ((unsigned)k[10] << 16) + ((unsigned)k[11] << 24));
    mix(a, b, c);
    k   += 12;
    len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch (len) {
    /* all the case statements fall through */
    case 11:
      c += ((unsigned)k[10] << 24);
    case 10:
      c += ((unsigned)k[9] << 16);
    case 9:
      c += ((unsigned)k[8] << 8);
    /* the first byte of c is reserved for the length */
    case 8:
      b += ((unsigned)k[7] << 24);
    case 7:
      b += ((unsigned)k[6] << 16);
    case 6:
      b += ((unsigned)k[5] << 8);
    case 5:
      b += k[4];
    case 4:
      a += ((unsigned)k[3] << 24);
    case 3:
      a += ((unsigned)k[2] << 16);
    case 2:
      a += ((unsigned)k[1] << 8);
    case 1:
      a += k[0];
      /* case 0: nothing left to add */
  }
  mix(a, b, c);
  /*-------------------------------------------- report the result */
  return c & hashtable_size;
}


int
vrna_standard_free_hash_entry(void *hash_entry)
{
  free(((vrna_standard_hash_entry *)hash_entry)->structure);
  return 0;
}
