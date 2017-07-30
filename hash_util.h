#ifndef _hash_util_h
#define _hash_util_h

#include <unordered_map>
#include <unordered_set>
#include <map>

extern "C" {
  #include "utils.h"
  #include "move_set_inside.h"
}

#include "pknots.h"

// help struct for hash
struct gw_struct {
  int count;
  struct_en he; // does not contain memory
  gw_struct(){
    he.structure = NULL;
    count = 0;
  }
};

#ifndef HASHBITS
#define HASHBITS 24
#endif
#define HASHSIZE (((unsigned) 1<<HASHBITS)-1)

#define mix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

struct hash_eq {
  bool operator()(const struct_en &lhs, const struct_en &rhs) const{
    int i=1;
    while (i<=lhs.structure[0] && lhs.structure[i]==rhs.structure[i]) {
      i++;
    }
    if (i>lhs.structure[0]) return true;
    else return false;
  }

  bool operator()(const short *lhs, const short *rhs) const {
    int i=1;
    while (i<=lhs[0] && lhs[i]==rhs[i]) {
      i++;
    }
    if (i>lhs[0]) return true;
    else return false;
  }

  bool operator()(const struct_en *lhs, const struct_en *rhs) const{
    int i=1;
    while (i<=lhs->structure[0] && lhs->structure[i]==rhs->structure[i]) {
      i++;
    }
    if (i>lhs->structure[0]) return true;
    else return false;
  }
  bool operator()(const Structure &lhs, const Structure &rhs) const{
    return lhs == rhs;
  }

  bool operator()(const Structure *lhs, const Structure *rhs) const{
    return *lhs == *rhs;
  }
};

struct hash_fncts{
  size_t operator()(const Structure &x) const {

  register short *k;        /* the key */
  register unsigned  length;   /* the length of the key */
  register unsigned  initval=0;  /* the previous hash, or an arbitrary value */
  register unsigned a,b,c,len;

  /* Set up the internal state */
  k = x.str;
  len = length = (unsigned) k[0];
  a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
  c = initval;         /* the previous hash value */

  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a += (k[0] +((unsigned)k[1]<<8) +((unsigned)k[2]<<16) +((unsigned)k[3]<<24));
    b += (k[4] +((unsigned)k[5]<<8) +((unsigned)k[6]<<16) +((unsigned)k[7]<<24));
    c += (k[8] +((unsigned)k[9]<<8) +((unsigned)k[10]<<16)+((unsigned)k[11]<<24));
    mix(a,b,c);
    k += 12; len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch(len) {             /* all the case statements fall through */
    case 11: c+=((unsigned)k[10]<<24);
    case 10: c+=((unsigned)k[9]<<16);
    case 9 : c+=((unsigned)k[8]<<8);
      /* the first byte of c is reserved for the length */
    case 8 : b+=((unsigned)k[7]<<24);
    case 7 : b+=((unsigned)k[6]<<16);
    case 6 : b+=((unsigned)k[5]<<8);
    case 5 : b+=k[4];
    case 4 : a+=((unsigned)k[3]<<24);
    case 3 : a+=((unsigned)k[2]<<16);
    case 2 : a+=((unsigned)k[1]<<8);
    case 1 : a+=k[0];
     /* case 0: nothing left to add */
  }
  mix(a,b,c);
   /*-------------------------------------------- report the result */
  return (c & HASHSIZE);
  }

  size_t operator()(const struct_en &x) const {

  register short *k;        /* the key */
  register unsigned  length;   /* the length of the key */
  register unsigned  initval=0;  /* the previous hash, or an arbitrary value */
  register unsigned a,b,c,len;

  /* Set up the internal state */
  k = x.structure;
  len = length = (unsigned) k[0];
  a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
  c = initval;         /* the previous hash value */

  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a += (k[0] +((unsigned)k[1]<<8) +((unsigned)k[2]<<16) +((unsigned)k[3]<<24));
    b += (k[4] +((unsigned)k[5]<<8) +((unsigned)k[6]<<16) +((unsigned)k[7]<<24));
    c += (k[8] +((unsigned)k[9]<<8) +((unsigned)k[10]<<16)+((unsigned)k[11]<<24));
    mix(a,b,c);
    k += 12; len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch(len) {             /* all the case statements fall through */
    case 11: c+=((unsigned)k[10]<<24);
    case 10: c+=((unsigned)k[9]<<16);
    case 9 : c+=((unsigned)k[8]<<8);
      /* the first byte of c is reserved for the length */
    case 8 : b+=((unsigned)k[7]<<24);
    case 7 : b+=((unsigned)k[6]<<16);
    case 6 : b+=((unsigned)k[5]<<8);
    case 5 : b+=k[4];
    case 4 : a+=((unsigned)k[3]<<24);
    case 3 : a+=((unsigned)k[2]<<16);
    case 2 : a+=((unsigned)k[1]<<8);
    case 1 : a+=k[0];
     /* case 0: nothing left to add */
  }
  mix(a,b,c);
   /*-------------------------------------------- report the result */
  return (c & HASHSIZE);
  }
  size_t operator()(const Structure *x) const {
  register short *k;        /* the key */
  register unsigned  length;   /* the length of the key */
  register unsigned  initval=0;  /* the previous hash, or an arbitrary value */
  register unsigned a,b,c,len;

  /* Set up the internal state */
  k = x->str;
  len = length = (unsigned) k[0];
  a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
  c = initval;         /* the previous hash value */

  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a += (k[0] +((unsigned)k[1]<<8) +((unsigned)k[2]<<16) +((unsigned)k[3]<<24));
    b += (k[4] +((unsigned)k[5]<<8) +((unsigned)k[6]<<16) +((unsigned)k[7]<<24));
    c += (k[8] +((unsigned)k[9]<<8) +((unsigned)k[10]<<16)+((unsigned)k[11]<<24));
    mix(a,b,c);
    k += 12; len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch(len) {             /* all the case statements fall through */
    case 11: c+=((unsigned)k[10]<<24);
    case 10: c+=((unsigned)k[9]<<16);
    case 9 : c+=((unsigned)k[8]<<8);
      /* the first byte of c is reserved for the length */
    case 8 : b+=((unsigned)k[7]<<24);
    case 7 : b+=((unsigned)k[6]<<16);
    case 6 : b+=((unsigned)k[5]<<8);
    case 5 : b+=k[4];
    case 4 : a+=((unsigned)k[3]<<24);
    case 3 : a+=((unsigned)k[2]<<16);
    case 2 : a+=((unsigned)k[1]<<8);
    case 1 : a+=k[0];
     /* case 0: nothing left to add */
  }
  mix(a,b,c);
   /*-------------------------------------------- report the result */
  return (c & HASHSIZE);
  }

  size_t operator()(const struct_en *x) const {

  register short *k;        /* the key */
  register unsigned  length;   /* the length of the key */
  register unsigned  initval=0;  /* the previous hash, or an arbitrary value */
  register unsigned a,b,c,len;

  /* Set up the internal state */
  k = x->structure;
  len = length = (unsigned) k[0];
  a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
  c = initval;         /* the previous hash value */

  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a += (k[0] +((unsigned)k[1]<<8) +((unsigned)k[2]<<16) +((unsigned)k[3]<<24));
    b += (k[4] +((unsigned)k[5]<<8) +((unsigned)k[6]<<16) +((unsigned)k[7]<<24));
    c += (k[8] +((unsigned)k[9]<<8) +((unsigned)k[10]<<16)+((unsigned)k[11]<<24));
    mix(a,b,c);
    k += 12; len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch(len) {             /* all the case statements fall through */
    case 11: c+=((unsigned)k[10]<<24);
    case 10: c+=((unsigned)k[9]<<16);
    case 9 : c+=((unsigned)k[8]<<8);
      /* the first byte of c is reserved for the length */
    case 8 : b+=((unsigned)k[7]<<24);
    case 7 : b+=((unsigned)k[6]<<16);
    case 6 : b+=((unsigned)k[5]<<8);
    case 5 : b+=k[4];
    case 4 : a+=((unsigned)k[3]<<24);
    case 3 : a+=((unsigned)k[2]<<16);
    case 2 : a+=((unsigned)k[1]<<8);
    case 1 : a+=k[0];
     /* case 0: nothing left to add */
  }
  mix(a,b,c);
   /*-------------------------------------------- report the result */
  return (c & HASHSIZE);
  }
};


// comparators: all one 1 place
// comparator for structures
bool compf_short (const short *lhs, const short *rhs);
struct comps_short {
  bool operator() (const short *lhs, const short *rhs) const {
    return compf_short(lhs, rhs);
  }
};

// comparator for hash_entries
bool compf_entries (const struct_en *lhs, const struct_en *rhs);
bool compf_entries2 (const struct_en &lhs, const struct_en &rhs);
struct comps_entries {
  bool operator() (const struct_en *lhs, const struct_en *rhs) const {
    if (lhs->energy!=rhs->energy) return lhs->energy<rhs->energy;
    return compf_short(lhs->structure, rhs->structure);
  }
  bool operator() (const struct_en &lhs, const struct_en &rhs) const {
    if (lhs.energy!=rhs.energy) return lhs.energy<rhs.energy;
    return compf_short(lhs.structure, rhs.structure);
  }
};

bool compf_short_rev (const short *lhs, const short *rhs);
struct comps_short_rev {
  bool operator() (const short *lhs, const short *rhs) const {
    return compf_short_rev(lhs, rhs);
  }
};

// comparator for hash_entries
bool compf_entries_rev (const struct_en *lhs, const struct_en *rhs);
struct comps_entries_rev {
  bool operator() (const struct_en *lhs, const struct_en *rhs) const {
    if (lhs->energy!=rhs->energy) return lhs->energy>rhs->energy;
    return compf_short_rev(lhs->structure, rhs->structure);
  }
  bool operator() (const Structure *lhs, const Structure *rhs) const {
    return *lhs<*rhs;
  }
};

// print stats about hash
void print_stats(std::unordered_map<struct_en, gw_struct, hash_fncts, hash_eq> &structs);
// add stats from hash to output map
void add_stats(std::unordered_map<struct_en, gw_struct, hash_fncts, hash_eq> &structs, std::map<struct_en, int, comps_entries> &output);


// free hash
void free_hash(std::unordered_map<struct_en, gw_struct, hash_fncts, hash_eq> &structs);
//void free_hash(unordered_map<Structure, gw_struct, hash_fncts, hash_eq> &structs);
void free_hash(std::unordered_set<struct_en*, hash_fncts, hash_eq> &structs);
void free_hash(std::unordered_set<Structure*, hash_fncts, hash_eq> &structs);

// entry handling
struct_en *copy_entry(const struct_en *he);
void free_entry(struct_en *he);
#endif
