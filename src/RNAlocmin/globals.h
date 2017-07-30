#ifndef __GLOBALS_H
#define __GLOBALS_H

#include <set>
#include <vector>

extern "C" {
  #include "RNAlocmin_cmdline.h"
}
#include "hash_util.h"
#include "move_set_pk.h"

// some global counters
static int num_moves = 0;
static int seq_len;

// structure for sequence related stuff:

class SeqInfo {
public:
  char *seq;
  short *s0;
  short *s1;

  SeqInfo();
  ~SeqInfo();
  void Init(char *seq);
};

// cute options singleton class
class Options {
  // options
public:
  float minh;   // leave out shallow minima (should be relativelly small)
  bool noLP;    // no lone pairs
  bool EOM;     // use energy_of_move
  bool first;   // use first descent, not deepest
  bool rand;    // use random walk, not deepest
  bool shift;   // use shifts?
  int verbose_lvl; // level of verbosity
  int floodMax; // cap for flooding
  bool neighs;  // use neighborhood routines?

  bool pknots; // flag for pseudoknots.

public:
  Options();

  // return 0 if success
  int Init(gengetopt_args_info &args_info);
};
/*
// structure for degeneracy
class Degen {
  // for degeneracy (structures with equal energies)
public:
  int current;    // all structures here have this energy
  set<short*, comps_short> processed;
  set<short*, comps_short> unprocessed;

public:
  Degen();
  ~Degen();
  void Clear();
};*/

// some good functions


// some singleton objects
//extern Degen Deg;
extern Options Opt;
//extern Encoded Enc;

#endif
