#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include <stack>

extern "C" {
  #include "pair_mat.h"
  #include "fold.h"
}

#include "globals.h"
#include "move_set_pk.h"

// some singleton objects
Options Opt;

SeqInfo::SeqInfo()
{
  seq = NULL;
  s0 = s1 = NULL;

  update_fold_params();
}

SeqInfo::~SeqInfo()
{
  if (s0) free(s0);
  if (s1) free(s1);
  if (seq) free(seq);
}

void SeqInfo::Init(char *seq) {
  this->seq = (char*) malloc(sizeof(char)*strlen(seq)+1);
  strcpy(this->seq, seq);
  make_pair_matrix();
  update_fold_params();
  s0 = encode_sequence(seq, 0);
  s1 = encode_sequence(seq, 1);
}

Options::Options()
{
}

int Options::Init(gengetopt_args_info &args_info)
{
  int ret = 0;

  if (args_info.min_num_arg<0) {
    fprintf(stderr, "Number of local minima should be non-negative integer (min-num)\n");
    ret = -1;
  }

  if (args_info.find_num_given && args_info.find_num_arg<=0) {
    fprintf(stderr, "Number of local minima should be positive integer (find-num)\n");
    ret = -1;
  }

  if (args_info.verbose_lvl_arg<0 || args_info.verbose_lvl_arg>4) {
    if (args_info.verbose_lvl_arg<0) args_info.verbose_lvl_arg = 0;
    else args_info.verbose_lvl_arg = 4;
    fprintf(stderr, "WARNING: level of verbosity is not in range (0-4), setting it to %d\n", args_info.verbose_lvl_arg);
  }

  if (args_info.temp_arg<-273.15) {
    fprintf(stderr, "Temperature cannot be below absolute zero\n");
    ret = -1;
  }

  if (args_info.floodMax_arg<0) {
    fprintf(stderr, "Flood cap must be non-negative\n");
    ret = -1;
  }

  if (args_info.floodPortion_arg<0.0 || args_info.floodPortion_arg>1.0) {
    args_info.floodPortion_arg = (args_info.floodPortion_arg<0.0 ? 0.0 : 1.0);
    fprintf(stderr, "WARNING: floodPortion is not in range (0.0-1.0), setting it to %.1f\n", args_info.floodPortion_arg);
  }

  if (args_info.depth_arg<=0) {
    fprintf(stderr, "Depth of findpath search should be positive integer\n");
    ret = -1;
  }

  if (args_info.minh_arg<0.0) {
    fprintf(stderr, "Depth of findpath search should be non-negative number\n");
    ret = -1;
  }

  if (args_info.numIntervals_arg<0) {
    fprintf(stderr, "Number of intervals should be non-negative number\n");
    ret = -1;
  }

  if (args_info.eRange_given && args_info.eRange_arg<0.0) {
    fprintf(stderr, "Energy range should be non-negative number\n");
    ret = -1;
  }

  if (args_info.dangles_arg<0 || args_info.dangles_arg>3) {
    fprintf(stderr, "Dangle treatment constant should be 0, 1, 2, or 3\n");
    ret = -1;
  }

  if (ret ==-1) return -1;

  // adjust options
  minh = (int)(args_info.minh_arg*100);
  noLP = args_info.noLP_flag;
  EOM = !args_info.useEOS_flag;
  first = args_info.walk_arg[0]=='F';
  rand = args_info.walk_arg[0]=='R';
  shift = args_info.move_arg[0]=='S';
  verbose_lvl = args_info.verbose_lvl_arg;
  floodMax = args_info.floodMax_arg;
  pknots = args_info.pseudoknots_flag;
  neighs = args_info.neighborhood_flag;

  return ret;
}

/*Degen::Degen()
{
  current = 0;
}

Degen::~Degen()
{
  Clear();
}

void Degen::Clear()
{
  // unprocessed
  for (set<short*, comps_short>::iterator it=unprocessed.begin(); it!=unprocessed.end(); it++) {
    if (*it) free(*it);
  }
  unprocessed.clear();

  // processed
  for (set<short*, comps_short>::iterator it=processed.begin(); it!=processed.end(); it++) {
    if (*it) free(*it);
  }
  processed.clear();
}*/
