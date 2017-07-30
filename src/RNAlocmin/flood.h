#ifndef __FLOOD_H
#define __FLOOD_H

#include "hash_util.h"
#include "globals.h"
extern "C" {
  #include "move_set.h"
  }


// hash_entry comparator function
bool compare_vect (const struct_en &lhs, const struct_en &rhs);
bool compare_vect (const Structure &lhs, const Structure &rhs);

// flood the structure - return one below saddle structure (should be freed then) energy of saddle is in "saddle_en"
  // minh_total - if set to true then flood also down and try to find energetically lower LM
  // maxh = height of flood (0 = infinity)
  // if returns NULL - in saddle_en is fail status - 1 for maxh reached, 0 otherwise
struct_en* flood(const struct_en &str, SeqInfo &sqi, int &saddle_en, int maxh = 0, bool pknots = false, bool minh_total = false);

#endif
