#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <queue>

extern "C" {
  //#include "utils.h"
  //#include "move_set_inside.h"
}

#include "flood.h"
#include "RNAlocmin.h"

using namespace std;

// global priority queue for stuff in flooding (does not hold memory - memory is in hash)
priority_queue<struct_en*, vector<struct_en*>, comps_entries_rev> neighs;
priority_queue<Structure*, vector<Structure*>, comps_entries_rev> neighs2;
int energy_lvl;
bool debugg;
int top_lvl;
int min_lvl;
bool minh_total;
bool found_exit;
// hash for the flooding
unordered_set<struct_en*, hash_fncts, hash_eq> hash_flood (HASHSIZE);
unordered_set<struct_en*, hash_fncts, hash_eq>::iterator it_hash;

unordered_set<Structure*, hash_fncts, hash_eq> hash_flood2 (HASHSIZE);
unordered_set<Structure*, hash_fncts, hash_eq>::iterator it_hash2;


void copy_se(struct_en *dest, const struct_en *src) {
  copy_arr(dest->structure, src->structure);
  dest->energy = src->energy;
}

struct_en *allocopy_se(const struct_en *src) {
  struct_en *dest = (struct_en*)malloc(sizeof(struct_en));
  dest->structure = allocopy(src->structure);
  dest->energy = src->energy;
  return dest;
}

// function to do on all the items...
int flood_func(struct_en *input, struct_en *output)
{
  // have we seen him?
  it_hash = hash_flood.find(input);
  if (it_hash != hash_flood.end()) {
    // nothing to do with already processed structure
    if (debugg) fprintf(stderr,     "   already seen: %s %.2f\n", pt_to_str(input->structure).c_str(), input->energy/100.0);
    return 0;
  } else {
    // found escape? (its energy is lower than our energy lvl and we havent seen it)
    if (input->energy < energy_lvl) {
      // if minh_total, then continue:
      if (minh_total) {
        // if we are lower than our min_lvl, we have found the exit:
        if (input->energy < min_lvl) {
          copy_se(output, input);
          found_exit = true;
          if (debugg) fprintf(stderr,   "    escape(min): %s %.2f\n", pt_to_str(input->structure).c_str(), input->energy/100.0);
          return 1;
        } else {
          //add it:
          if (debugg) fprintf(stderr, "    adding(min): %s %.2f\n", pt_to_str(input->structure).c_str(), input->energy/100.0);
          // just add it to the queue... and to hash
          struct_en *he_tmp = (struct_en*)space(sizeof(struct_en));
          he_tmp->structure = allocopy(input->structure);
          he_tmp->energy = input->energy;
          neighs.push(he_tmp);
          hash_flood.insert(he_tmp);
          return 0;
        }
      }

      // ends flood and return it as a structure to walk down
      copy_se(output, input);
      found_exit = true;
      if (debugg) fprintf(stderr,   "       escape  : %s %.2f\n", pt_to_str(input->structure).c_str(), input->energy/100.0);
      return 1;
    } else {
      if (input->energy > top_lvl) {
        if (debugg) fprintf(stderr, "energy too high: %s %.2f\n", pt_to_str(input->structure).c_str(), input->energy/100.0);
        return 0;
      } else {
        if (debugg) fprintf(stderr, "       adding  : %s %.2f\n", pt_to_str(input->structure).c_str(), input->energy/100.0);
        // just add it to the queue... and to hash
        struct_en *he_tmp = (struct_en*)space(sizeof(struct_en));
        he_tmp->structure = allocopy(input->structure);
        he_tmp->energy = input->energy;
        neighs.push(he_tmp);
        hash_flood.insert(he_tmp);
        return 0;
      }
    }
  }
}

// function to do on all the items...
int flood_func2(Structure *input, Structure *output)
{
  // have we seen him?
  it_hash2 = hash_flood2.find(input);
  if (it_hash2 != hash_flood2.end()) {
    // nothing to do with already processed structure
    if (debugg) fprintf(stderr,     "   already seen: %s %.2f\n", pt_to_str(input->str).c_str(), input->energy/100.0);
    return 0;
  } else {
    // found escape? (its energy is lower than our energy lvl and we havent seen it)
    if (input->energy < energy_lvl) {
      // if minh_total, then continue:
      if (minh_total) {
        // if we are lower than our min_lvl, we have found the exit:
        if (input->energy < min_lvl) {
          // ends flood and return it as a structure to walk down
          *output = *input;
          found_exit = true;
          if (debugg) fprintf(stderr,   "    escape(min): %s %.2f\n", pt_to_str(input->str).c_str(), input->energy/100.0);
          return 1;
        } else {
          //add it:
          if (debugg) fprintf(stderr, "    adding(min): %s %.2f\n", pt_to_str(input->str).c_str(), input->energy/100.0);
          // just add it to the queue... and to hash
          Structure *str_tmp = new Structure(*input);
          neighs2.push(str_tmp);
          hash_flood2.insert(str_tmp);
          return 0;
        }
      }


      // ends flood and return it as a structure to walk down
      *output = *input;
      found_exit = true;
      if (debugg) fprintf(stderr,   "       escape  : %s %.2f\n", pt_to_str(input->str).c_str(), input->energy/100.0);
      return 1;
    } else {
      if (input->energy > top_lvl) {
        if (debugg) fprintf(stderr, "energy too high: %s %.2f\n", pt_to_str(input->str).c_str(), input->energy/100.0);
        return 0;
      } else {
        if (debugg) fprintf(stderr, "       adding  : %s %.2f\n", pt_to_str(input->str).c_str(), input->energy/100.0);
        // just add it to the queue... and to hash
        Structure *str_tmp = new Structure(*input);
        neighs2.push(str_tmp);
        hash_flood2.insert(str_tmp);
        return 0;
      }
    }
  }
}

struct_en* flood(const struct_en &he, SeqInfo &sqi, int &saddle_en, int maxh, bool pknots, bool flood_total)
{
  int count = 0;
  debugg = Opt.verbose_lvl>2;

  struct_en *res = NULL;

  // save deg.first + create deg
  bool first = Opt.first;
  Opt.first = true;

  // if minh specified, assign top_lvl and flood_total
  if (maxh>0) {
    top_lvl = he.energy + maxh;
    minh_total = flood_total;
    min_lvl = he.energy;
  } else {
    top_lvl = 1e9;
  }

  ///#### PKNOTS
  if (pknots) {
    // init priority queue
    while (!neighs2.empty()) {
      //fprintf(stderr, "-neighs size: %d\n", (int)neighs.size());
      neighs2.pop();
    }

    // init hash
    free_hash(hash_flood2);
    found_exit = false;

    // add the first structure to hash, get its adress and add it to priority queue
    {
      Structure *he_tmp = new Structure(he.structure, he.energy);
      neighs2.push(he_tmp);
      hash_flood2.insert(he_tmp);
    }

    // FLOOOD!
    while ((int)hash_flood2.size() < Opt.floodMax) {
      // should not be empty (only when maxh specified)
      if (neighs2.empty()) break;

      // get structure
      Structure *he_top = neighs2.top();
      neighs2.pop();
      energy_lvl = he_top->energy;

      if (Opt.verbose_lvl>2) fprintf(stderr, "  neighbours of: %s %.2f (%d)\n", pt_to_str(he_top->str).c_str(), he_top->energy/100.0, (int)neighs2.size());

      int verbose = Opt.verbose_lvl<2?0:Opt.verbose_lvl-2;
      he_top->energy = browse_neighs_pk_pt(sqi.seq, he_top, sqi.s0, sqi.s1, Opt.shift, verbose, flood_func2);

      if (found_exit && Opt.verbose_lvl>2) fprintf(stderr, "sad= %6.2f    : %s %.2f\n", saddle_en/100.0, pt_to_str(he_top->str).c_str(), he_top->energy/100.0);

      // did we find exit from basin?
      if (found_exit) {
        res = (struct_en*)malloc(sizeof(struct_en));
        res->structure = allocopy(he_top->str);
        res->energy = he_top->energy;
        break;
      }

      count++;
    }

    // return status in saddle_en :/
    if (!found_exit) {
      saddle_en = (neighs2.empty() ? 1 : 0);
    }

    // destroy queue
    while (!neighs2.empty()) {
      //fprintf(stderr, "-neighs size: %d\n", (int)neighs.size());
      neighs2.pop();
    }

    // destroy hash
    free_hash(hash_flood2);
  } else {  /// ######## NOT PKNOTS!!!
    // init priority queue
    while (!neighs.empty()) {
      //fprintf(stderr, "-neighs size: %d\n", (int)neighs.size());
      neighs.pop();
    }

    // init hash
    free_hash(hash_flood);
    found_exit = false;


    // add the first structure to hash, get its adress and add it to priority queue
    {
      struct_en *he_tmp = allocopy_se(&he);
      neighs.push(he_tmp);
      hash_flood.insert(he_tmp);
    }

    // FLOOOD!
    while ((int)hash_flood.size() < Opt.floodMax) {
      // should not be empty (only when maxh specified)
      if (neighs.empty()) break;

      // get structure
      struct_en *he_top = neighs.top();
      neighs.pop();
      energy_lvl = he_top->energy;

      if (Opt.verbose_lvl>2) fprintf(stderr, "  neighbours of: %s %.2f\n", pt_to_str(he_top->structure).c_str(), he_top->energy/100.0);

      int verbose = Opt.verbose_lvl<2?0:Opt.verbose_lvl-2;
      he_top->energy = browse_neighs_pt(sqi.seq, he_top->structure, sqi.s0, sqi.s1, verbose, Opt.shift, Opt.noLP, flood_func);

      if (found_exit && Opt.verbose_lvl>2) fprintf(stderr, "sad= %6.2f    : %s %.2f\n", saddle_en/100.0, pt_to_str(he_top->structure).c_str(), he_top->energy/100.0);

      // did we find exit from basin?
      if (found_exit) {
        res = allocopy_se(he_top);
        break;
      }

      count++;
    }

    // return status in saddle_en :/
    if (!found_exit) {
      saddle_en = (neighs.empty() ? 1 : 0);
    }

    // destroy queue
    while (!neighs.empty()) {
      //fprintf(stderr, "-neighs size: %d\n", (int)neighs.size());
      neighs.pop();
    }

    // destroy hash
    free_hash(hash_flood);
  }  /// #### END OF PKNOTS BRANCH

  // restore deg options
  Opt.first = first;

  // return found? structure
  return res;
}
