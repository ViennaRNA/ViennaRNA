#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hash_util.h"

using namespace std;

// does not depend on number - just on dot-bracket notation:
bool compf_short (const short *lhs, const short *rhs) {
  int i=1;
  char l,r;
  while (i<=lhs[0]) {
    l = (lhs[i]==0?'.':(lhs[i]<lhs[lhs[i]]?')':'('));
    r = (rhs[i]==0?'.':(rhs[i]<rhs[rhs[i]]?')':'('));
    if (l != r) {
      //fprintf (stderr, "%c %c %d %d\n", l, r, l, r);
      break;
    }
    i++;
  }
  return (i<=lhs[0] && l<r);
}

bool compf_entries (const struct_en *lhs, const struct_en *rhs) {
  if (lhs->energy!=rhs->energy) return lhs->energy<rhs->energy;
  return compf_short(lhs->structure, rhs->structure);
}
bool compf_entries2 (const struct_en &lhs, const struct_en &rhs)
{
  if (lhs.energy!=rhs.energy) return lhs.energy<rhs.energy;
  return compf_short(lhs.structure, rhs.structure);
}

bool compf_short_rev (const short *lhs, const short *rhs) {
  /*int i=1;
  char l,r;
  while (i<=lhs[0]) {
    l = (lhs[i]==0?'.':(lhs[i]<lhs[lhs[i]]?')':'('));
    r = (rhs[i]==0?'.':(rhs[i]<rhs[rhs[i]]?')':'('));
    if (l != r) break;
    i++;
  }
  return (i<=lhs[0] && l>r);  /// REVERSE ORDERING*/
  return compf_short(rhs, lhs);
}

bool compf_entries_rev (const struct_en *lhs, const struct_en *rhs) {
  if (lhs->energy!=rhs->energy) return lhs->energy>rhs->energy;
  return compf_short_rev(lhs->structure, rhs->structure);
}

void print_stats(unordered_map<struct_en, gw_struct, hash_fncts, hash_eq> &structs)
{
  double mean = 0.0;
  int count = 0;
  double entropy = 0.0;
  unordered_map<struct_en, gw_struct, hash_fncts>::iterator it;
  for (it=structs.begin(); it!=structs.end(); it++) {
    count += it->second.count;
    mean += (it->first.energy)*(it->second.count);
    entropy += it->second.count*log(it->second.count);
  }

  mean /= (double)count*100.0;
  entropy = entropy/(double)count - log(count);

  fprintf(stderr, "Mean  : %.3f (Entrpy: %.3f)\n", mean, entropy);
}

void add_stats(unordered_map<struct_en, gw_struct, hash_fncts, hash_eq> &structs, map<struct_en, int, comps_entries> &output)
{
  unordered_map<struct_en, gw_struct, hash_fncts>::iterator it;
  for (it=structs.begin(); it!=structs.end(); it++) {
    // add stats:
    //fprintf(stderr, "struct: %s %6.2f %d\n", pt_to_str(it->second.he.structure).c_str(), it->second.he.energy/100.0, it->second.count);

    if (output.count(it->second.he) == 0) {
      fprintf(stderr, "ERROR: output does not contain structure it should!!!\n");
      //if (!Opt.pknots) exit(EXIT_FAILURE);
    }
    output[it->second.he] += it->second.count-1;
  }
}

// free hash
void free_hash(unordered_map<struct_en, gw_struct, hash_fncts, hash_eq> &structs)
{
  unordered_map<struct_en, gw_struct, hash_fncts>::iterator it;
  for (it=structs.begin(); it!=structs.end(); it++) {
    free(it->first.structure);
  }
  structs.clear();
}

// free hash
void free_hash(unordered_set<struct_en*, hash_fncts, hash_eq> &structs)
{
  unordered_set<struct_en*, hash_fncts, hash_eq>::iterator it;
  for (it=structs.begin(); it!=structs.end(); it++) {
    free((*it)->structure);
    free(*it);
  }
  structs.clear();
}

// free hash
void free_hash(unordered_set<Structure*, hash_fncts, hash_eq> &structs)
{
  unordered_set<Structure*, hash_fncts, hash_eq>::iterator it;
  for (it=structs.begin(); it!=structs.end(); it++) {
    //free((*it)->str);
    delete *it;
  }
  structs.clear();
}
