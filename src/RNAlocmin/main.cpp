#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <string.h>

#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <memory>

extern "C" {
  #include "fold.h"
  #include "findpath.h"
  #include "RNAlocmin_cmdline.h"
  #include "utils.h"
  #include "read_epars.h"
  #include "fold_vars.h"
  #include "params.h"
  #include "move_set_inside.h"
}

#include "findpath_pk.h"

#include "hash_util.h"
#include "globals.h"
#include "RNAlocmin.h"
#include "flood.h"
#include "RNAlocmin.h"
#include "move_set_pk.h"
#include "neighbourhood.h"

#include "barrier_tree.h"

using namespace std;

// dirty hack to "allegiance" stuff
 //to be filled by 1st run
map<struct_en, struct_en, comps_entries> str_to_LM;
vector<struct_en> structures;
 // to be filled after
map<struct_en, int, comps_entries> LM_to_LMnum;
static bool allegiance = false;


inline bool isSeq(char *p)
{
  // check first two chars - should be enough
  char seq[]="ACGTUactgu";
  int ok = 0;
  for (unsigned int i=0; i<strlen(seq); i++) {
    if (p[0]==seq[i]) {
      ok++;
    }
    if (p[1]==seq[i]) {
      ok++;
    }
    if (ok==2) break;
  }

  return (ok==2);
}

inline bool isStruct(char *p)
{
  // check first two chars - should be enough
  if ((p[0]=='.' || p[0]=='(') && (p[1]=='.' || p[1]=='(' || p[1]=='[')) return true;
  else return false;
}

inline bool isEnergy(char *p, float &energy)
{
  if (sscanf(p, "%f", &energy) == 1) return true;
  else return false;
}

inline int en_fltoi(float en)
{
  if (en < 0.0) return (int)(en*100 - 0.5);
  else return (int)(en*100 + 0.5);
}

struct barr_info { // info taken from barriers
  int father;
  int e_diff;
  int bsize;
  int fbsize;
  // unused, unmodified
  float fen;
  int grad;
  float feng;

  barr_info &operator+(barr_info &second) {
    bsize += second.bsize;
    if (second.father!=father) father = -1;
    fbsize += second.fbsize;
    grad += second.grad;
    return *this;
  }
};

// functions that are down in file ;-)
char *read_seq(char *seq_arg, char **name_out);
int move(unordered_map<struct_en, gw_struct, hash_fncts, hash_eq> &structs, map<struct_en, int, comps_entries> &output, set<struct_en, comps_entries> &output_shallow, SeqInfo &sqi, bool pure_output);
char *read_previous(char *previous, map<struct_en, int, comps_entries> &output);
char *read_barr(char *previous, map<struct_en, barr_info, comps_entries> &output);

int main(int argc, char **argv)
{
  clock_t clck1 = clock();

  // testing
  //test();
  //exit(EXIT_SUCCESS);

  // parse arguments
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    fprintf(stderr, "ERROR: argument parsing problem.\n");
    exit(EXIT_FAILURE);
  }

  //check for bad arguments
  if (Opt.Init(args_info) !=0 ) {
    fprintf(stderr, "ERROR: one or more bad arguments, exiting...\n");
    exit(EXIT_FAILURE);
  }

  // adjust temperature
  if (args_info.temp_given) {
    temperature = args_info.temp_arg;
  }

  // read parameter file
  if (args_info.paramFile_given) {
    read_parameter_file(args_info.paramFile_arg);
  }

  // dangle setup
  if (args_info.dangles_given) {
    dangles = args_info.dangles_arg;
    model_detailsT  md;
    set_model_details(&md);
    md.dangles = dangles;
  }

  // adjust dependencies:
  if (args_info.barr_name_given) {args_info.bartree_flag = true;}
  if (args_info.rates_file_given) {args_info.rates_flag = true;}

  // degeneracy setup
  if (args_info.degeneracy_off_flag) {
    degeneracy_handling(0);
    Neighborhood::SwitchOffDegen();
  }

  //try_pk();
  //exit(0);

  // keep track of structures & statistics
  map<struct_en, int, comps_entries> output; // structures plus energies to output (+ how many hits has this minima)
  map<struct_en, barr_info, comps_entries> output_barr; // structures plus energies to output (+ barr_info)
  set<struct_en, comps_entries> output_shallow; // shallow structures (if minh specified)
  char *seq = NULL;
  char *name = NULL;

  // read previous LM or/and sequence
  if (args_info.fix_barriers_given) {
    seq = read_barr(args_info.fix_barriers_arg, output_barr);
  } else {
    if (args_info.previous_given) {
      seq = read_previous(args_info.previous_arg, output);
    } else {
      seq = read_seq(args_info.seq_arg, &name);
    }
  }

  if (args_info.verbose_lvl_arg>1) fprintf(stderr, "%s\n", seq);

  seq_len = strlen(seq);

  // allegiance:
  FILE *alleg = NULL;
  if (args_info.allegiance_given) {
    alleg = fopen(args_info.allegiance_arg, "w");
    if (alleg) {
      allegiance = true;
      fprintf(alleg, "       %s\n", seq);
    }
  }

  // time?
  if (args_info.verbose_lvl_arg>0) {
    fprintf(stderr, "Time to initialize: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
    clck1 = clock();
  }

  // ########################## begin main loop - reads structures from RNAsubopt and process them
  int count = max(output.size(), output_barr.size());  //num of local minima
  SeqInfo sqi;
  sqi.Init(seq);

  int not_canonical = 0;

  if (!args_info.fix_barriers_given) {

    // if direct output:
    if (args_info.just_output_flag) printf("%s\n", seq);

    // hash
    unordered_map<struct_en, gw_struct, hash_fncts, hash_eq> structs (HASHSIZE); // structures to minima map
    while ((!args_info.find_num_given || count != args_info.find_num_arg) && !args_info.just_read_flag) {
      int res = move(structs, output, output_shallow, sqi, args_info.just_output_flag);

      // print out
      //if (Opt.verbose_lvl>0 && num_moves%10000==0) fprintf(stderr, "processed %d, minima %d, time %f secs.\n", num_moves, count, (clock()-clck1)/(double)CLOCKS_PER_SEC);
      if (Opt.verbose_lvl>0 && num_moves%(Opt.pknots?1000:10000)==0 && num_moves!=0) fprintf(stderr, "processed %d, minima %d, time %f secs.\n", num_moves, (int)output.size(), (clock()-clck1)/(double)CLOCKS_PER_SEC);

      // evaluate results
      if (res==0)   continue; // same structure has been processed already
      if (res==-1)  break; // error or end
      if (res==-2)  not_canonical++;
      if (res==1)   count=output.size();
    }

    if (args_info.just_output_flag) {
      // end it

      // release resources
      if (seq!=NULL) free(seq);
      if (name!=NULL) free(name);
      cmdline_parser_free(&args_info);
      freeP();
      free_arrays();

      // time?
      if (args_info.verbose_lvl_arg>0) {
        fprintf(stderr, "Printing results + freeing args: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
        clck1 = clock();
      }
      return 0;
    }


    // ########################## end main loop - reads structures from RNAsubopt and process them

    //int num_of_structures = hash_size();
    if (args_info.verbose_lvl_arg>0) print_stats(structs);
    add_stats(structs, output);

    // time?
    if (args_info.verbose_lvl_arg>0) {
      fprintf(stderr, "Main loop (deepest descent from RNAsubopt): %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
      clck1 = clock();
    }

    // vectors for rates computation
    int num = ((count > args_info.min_num_arg && args_info.min_num_arg!=0) ? args_info.min_num_arg : count);
    vector<string> output_str;
    output_str.resize(num);
    vector<struct_en> output_he;
    struct_en h;
    h.structure = NULL;
    output_he.resize(num, h);
    vector<int> output_en;
    output_en.resize(num);
    vector<int> output_num;
    output_num.resize(num);

    // threshold for flooding
    int threshold;

    int i=0;
    int ii=0;
    for (map<struct_en, int, comps_entries>::iterator it=output.begin(); it!=output.end(); it++) {
      ii++;
      // if not enough minima
      if (i<num) {
        // first check if the output is not shallow
        if (Opt.minh>0) {
          int saddle;
          struct_en *escape = flood(it->first, sqi, saddle, Opt.minh, args_info.pseudoknots_flag, !args_info.minh_lite_flag);

          if (args_info.verbose_lvl_arg>0 && ii%100 == 0) {
            fprintf(stderr, "non-shallow remained: %d / %d; time: %.2f secs.\n", i, ii, (clock()-clck1)/(double)CLOCKS_PER_SEC);
          }

          // shallow
          if (escape) {
            if (args_info.verbose_lvl_arg>1) {
              fprintf(stderr, "shallow: %s %6.2f (saddle: %s %6.2f)\n", pt_to_str_pk(it->first.structure).c_str(), it->first.energy/100.0, pt_to_str_pk(escape->structure).c_str(), escape->energy/100.0);
            }
            free(escape->structure);
            free(escape);
            free(it->first.structure);
            continue;
          }
        }

        if (args_info.verbose_lvl_arg > 2) fprintf(stderr, "%4d %s %6.2f %6d\n", i+1, pt_to_str_pk(it->first.structure).c_str(), it->first.energy/100.0, it->second);
        output_str[i]=pt_to_str_pk(it->first.structure);
        output_num[i]=it->second;
        output_he[i]=it->first;
        output_en[i]=it->first.energy;
        i++;

        // allegiance:
        if (allegiance) LM_to_LMnum[it->first] = i;

      } else { // we have enough minima
        free(it->first.structure);
      }
    }
    output.clear();

    // allegiance:
    if (allegiance) {
      for (int i=0; i<(int)structures.size(); i++) {
        fprintf(alleg, "%6d %s %6.2f %6d\n", i+1, pt_to_str_pk(structures[i].structure).c_str(), structures[i].energy/100.0, LM_to_LMnum[str_to_LM[structures[i]]]);
        //free(structures[i].structure);
      }
      structures.clear();
      LM_to_LMnum.clear();
      str_to_LM.clear();
    }

    // erase possible NULL elements...
    /*if (args_info.noSort_flag) {
      for (int j=output_he.size()-1; j>=0; j--) {
        if (output_he[j].structure==NULL) {
          output_he.erase(output_he.begin()+j);
          output_str.erase(output_str.begin()+j);
          output_en.erase(output_en.begin()+j);
          output_num.erase(output_num.begin()+j);
        }
      }
    } else {*/
      output_he.resize(i);
      output_str.resize(i);
      output_en.resize(i);
      output_num.resize(i);
    //}

    // time?
    if (Opt.minh>0 && args_info.verbose_lvl_arg>0) {
      fprintf(stderr, "Discarding shallow minima: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
      clck1 = clock();
    }

    // array of energy barriers
    float *energy_barr = NULL;
    bool *findpath_barr = NULL;

    // find saddles - fill energy barriers
    if (args_info.rates_flag || args_info.bartree_flag || args_info.barrier_file_given) {
      // threshold for flooding
      vector<int> tmp = output_num;
      sort(tmp.begin(), tmp.end());
      int thr = num*args_info.floodPortion_arg;
      thr--;
      threshold = (thr<0 ? 0 : tmp[thr]);

      // nodes
      std::unique_ptr<nodeT []> nodes(new nodeT[num]);
      energy_barr = (float*) malloc(num*num*sizeof(float));
      for (int i=0; i<num*num; i++) energy_barr[i]=1e10;
      findpath_barr = (bool*) malloc(num*num*sizeof(bool));
      for (int i=0; i<num*num; i++) findpath_barr[i]=false;

      // fill nodes
      for (int i=0; i<num; i++) {
        nodes[i].father = -1;
        nodes[i].height = output_en[i]/100.0;
        nodes[i].label = NULL;
        nodes[i].color = 0.0;
        nodes[i].saddle_height = 1e10;
      }

      int flooded = 0;
      // init union-findset
      init_union(num);
      // first try to flood the highest bins
      for (int i=num-1; i>=0; i--) {
        // flood only if low number of walks ended there
        if (output_num[i]<=threshold && Opt.floodMax>0) {
          //copy_arr(Enc.pt, output_he[i].structure);
          if (args_info.verbose_lvl_arg>2) fprintf(stderr,   "flooding  (%3d): %s %.2f\n", i+1, output_str[i].c_str(), output_he[i].energy/100.0);

          int saddle;
          struct_en *he = flood(output_he[i], sqi, saddle, Opt.minh, args_info.pseudoknots_flag);

          // print info
          if (args_info.verbose_lvl_arg>1) {
            if (he) {
              fprintf(stderr, "below     (%3d): %s %.2f\n"
                              "en: %7.2f  is: %s %.2f\n", i,
                      output_str[i].c_str(), output_he[i].energy/100.0, saddle/100.0,
                      pt_to_str_pk(he->structure).c_str(), he->energy/100.0);
            } else {
              fprintf(stderr, "unsucesful(%3d): %s %.2f\n", i,
                      output_str[i].c_str(), output_he[i].energy/100.0);
            }
          }
          // if flood succesfull - walk down to find father minima
          if (he) {
            // walk down
            move_set(*he, sqi);

            // now check if we have the minimum already (hopefuly yes ;-) )
            vector<struct_en>::iterator it;
            it = lower_bound(output_he.begin(), output_he.end(), *he, compf_entries2);

            if (args_info.verbose_lvl_arg>1) fprintf(stderr, "minimum: %s %.2f\n", pt_to_str_pk(he->structure).c_str(), he->energy/100.0);
            // we dont need it again

            hash_eq heq;
            if (it!=output_he.end() && heq(&*it, he)) {
              int pos = (int)(it-output_he.begin());
              if (args_info.verbose_lvl_arg>1) fprintf(stderr, "found father at pos: %d\n", pos);

              flooded++;
              energy_barr[i*num+pos] = energy_barr[pos*num+i] = saddle/100.0;

              // union set
              //fprintf(stderr, "join: %d %d\n", min(i, pos), max(i, pos));
              union_set(min(i, pos), max(i, pos));
            }
            free(he->structure);
            free(he);
          }
        }
      }

      // time?
      if (args_info.verbose_lvl_arg>0) {
        fprintf(stderr, "Flood(%d(%d)/%d): %.2f secs.\n", flooded, (int)(num*args_info.floodPortion_arg), num, (clock() - clck1)/(double)CLOCKS_PER_SEC);
        clck1 = clock();
      }

      // for others, just do findpath
      int findpath = 0;
      set<int> to_findpath;
      for (int i=0; i<num; i++) to_findpath.insert(find(i));

      if (args_info.verbose_lvl_arg>1) {
        fprintf(stderr, "Minima left to findpath (their father = -1): ");
        for (set<int>::iterator it=to_findpath.begin(); it!=to_findpath.end(); it++) {
          fprintf(stderr, "%d ", *it);
        }
        fprintf(stderr, "\n");
      }

      // findpath:
      for (set<int>::iterator it=to_findpath.begin(); it!=to_findpath.end(); it++) {
        set<int>::iterator it2=it;
        it2++;
        for (; it2!=to_findpath.end(); it2++) {
          if (args_info.pseudoknots_flag) energy_barr[(*it2)*num+(*it)] = energy_barr[(*it)*num+(*it2)] = find_saddle_pk(seq, output_str[*it].c_str(), output_str[*it2].c_str(), args_info.depth_arg)/100.0;
          else energy_barr[(*it2)*num+(*it)] = energy_barr[(*it)*num+(*it2)] = find_saddle(seq, output_str[*it].c_str(), output_str[*it2].c_str(), args_info.depth_arg)/100.0;
          findpath_barr[(*it2)*num+(*it)] = findpath_barr[(*it)*num+(*it2)] = true;
          if (args_info.verbose_lvl_arg>0 && findpath %10000==0){
            fprintf(stderr, "Findpath:%7d/%7d\n", findpath, (int)(to_findpath.size()*(to_findpath.size()-1)/2));
          }
          findpath++;
        }
      }

      // debug output
      if (args_info.verbose_lvl_arg>2) {
        fprintf(stderr, "Energy barriers:\n");
        //bool symmetric = true;
        for (int i=0; i<num; i++) {
          for (int j=0; j<num; j++) {
            fprintf(stderr, "%8.2g%c ", energy_barr[i*num+j], (findpath_barr[i*num+j]?'~':' '));
            //if (energy_barr[i*num+j] != energy_barr[j*num+i]) symmetric = false;
          }
          fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");
        //fprintf(stderr, "%s", (symmetric? "":"non-symmetric energy barriers!!\n"));
      }

      // time?
      if (args_info.verbose_lvl_arg>0) {
        fprintf(stderr, "Findpath(%d/%d): %.2f secs.\n", findpath, num*(num-1)/2, (clock() - clck1)/(double)CLOCKS_PER_SEC);
        clck1 = clock();
      }

      // create rates for treekin
      if (args_info.rates_flag) {
        print_rates(args_info.rates_file_arg, args_info.temp_arg, num, energy_barr, output_en);
      }

      // saddles for evaluation
      if (args_info.barrier_file_given) {
        print_rates(args_info.barrier_file_arg, args_info.temp_arg, num, energy_barr, output_en, true);
      }

      // generate barrier tree?
      if (args_info.bartree_flag) {

        //PS_tree_plot(nodes, num, "tst.ps");

        // make tree (fill missing nodes)
        make_tree(num, energy_barr, findpath_barr, nodes.get());

        // plot it!
        PS_tree_plot(nodes.get(), num, args_info.barr_name_arg);
      }

      // time?
      if (args_info.verbose_lvl_arg>0) {
        fprintf(stderr, "Rates + barrier tree generation: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
        clck1 = clock();
      }

      // ############### actual output
      // printf output with fathers! maybe
      printf("     %s\n", seq);
      for (unsigned int i=0; i<output_str.size(); i++) {
        if (args_info.eRange_given) {
          if ((output_en[i] - output_en[0]) >  args_info.eRange_arg*100 ) {
            break;
          }
        }
        printf("%4d %s %6.2f", i+1, output_str[i].c_str(), output_en[i]/100.0);
        if (args_info.bartree_flag) printf(" %4d %6.2f\n", nodes[i].father+1, nodes[i].saddle_height-nodes[i].height);
        printf(" %6d\n", output_num[i]);
      }

    } else {

      // printf output without fathers!
      printf("     %s\n", seq);
      for (unsigned int i=0; i<output_str.size(); i++) {
        if (args_info.eRange_given) {
          if ((output_en[i] - output_en[0]) >  args_info.eRange_arg*100 ) {
            break;
          }
        }
        printf("%4d %s %6.2f %6d", i+1, output_str[i].c_str(), output_en[i]/100.0, output_num[i]);
        printf("\n");
      }
    }

    // release smth
    for(unsigned int i=0; i<output_he.size(); i++) {
      free(output_he[i].structure);
    }
    // free hash
    free_hash(structs);

    // release res:
    if (energy_barr!=NULL) free(energy_barr);
    if (findpath_barr!=NULL) free(findpath_barr);

  } else { // fix-barrier part (just print output):
    int mfe = output_barr.begin()->first.energy;
    printf("     %s\n", seq);
    int i=1;
    for (map<struct_en, barr_info, comps_entries>::iterator it=output_barr.begin(); it!=output_barr.end(); it++) {
      if (args_info.eRange_given) {
        if ((it->first.energy - mfe) >  args_info.eRange_arg*100 ) {
          break;
        }
      }
      const struct_en &he = it->first;
      barr_info &bi = it->second;
      printf("%4d %s %6.2f", i, pt_to_str_pk(he.structure).c_str(), he.energy/100.0);
      printf(" %4d %6.2f %6d %6d %10.6f %6d %10.6f\n", bi.father+1, bi.e_diff/100.0, bi.bsize, bi.fbsize, bi.fen, bi.grad, bi.feng);
      i++;
    }
  }

  // release resources
  if (seq!=NULL) free(seq);
  if (name!=NULL) free(name);
  cmdline_parser_free(&args_info);
  freeP();
  free_arrays();

  // time?
  if (args_info.verbose_lvl_arg>0) {
    fprintf(stderr, "Printing results + freeing args: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
    clck1 = clock();
  }

  return 0;
}

char *read_previous(char *previous, map<struct_en, int, comps_entries> &output)
{
  char *seq = NULL;
  FILE *fprev;
  fprev = fopen(previous, "r");
  if (fprev == NULL) {
    fprintf(stderr, "Cannot open file \"%s\".\n", previous);
    exit(EXIT_FAILURE);
  }
  char *line = my_getline(fprev);
  char *p = strtok(line, " \t\n");
  while (p) {
    if (isSeq(p)) {
      seq = (char*)malloc((strlen(p)+1)*sizeof(char));
      strcpy(seq, p);
      break;
    }
    p = strtok(NULL, " \t\n");
  }
  free(line);

  if (seq == NULL) {
    fprintf(stderr, "Couldn't find sequence on first line of file \"%s\"\n", previous);
    fclose(fprev);
    exit(EXIT_FAILURE);
  }
  // read previously found minima:
  struct_en he;
  while ((line = my_getline(fprev))) {
    p = strtok(line, " \t\n");

    he.structure = NULL;
    he.energy = INT_MAX;

    // read the stuff
    int num;
    sscanf(p, "%d", &num);
    p = strtok(NULL, " \t\n");
    if (p && isStruct(p)) {
      he.structure = make_pair_table_PK(p);
    }
    p = strtok(NULL, " \t\n");
    if (p && he.structure && he.energy==INT_MAX) {
      float en;
      if (sscanf(p, "%f", &en)==1) {
        he.energy = en_fltoi(en);
      }
    }
    // 3 aternatives: read father, e_diff and count; or read only count; or read nothing
    int father;
    float e_diff;
    int count;
    p = strtok(NULL, " \t\n");
    if (p && he.structure && he.energy!=INT_MAX) {
      sscanf(p, "%d", &father);
      p = strtok(NULL, " \t\n");
      if (p && he.structure && he.energy!=INT_MAX) {
        sscanf(p, "%f", &e_diff);
        p = strtok(NULL, " \t\n");
        if (p && he.structure && he.energy!=INT_MAX) {
          sscanf(p, "%d", &count);
          p = strtok(NULL, " \t\n");
          if (p && he.structure && he.energy!=INT_MAX) {
            count = 0; // because we are reading barriers file without info about count...
          }
        } else count = 0;

      } else {
        count = father;
      }
    } else count = 0;
    // insert
    output[he]=count;
    free(line);
  }

  fclose(fprev);

  return seq;
}

char *read_barr(char *barr_arg, map<struct_en, barr_info, comps_entries> &output)
{
  char *seq = NULL;
  FILE *fbarr;
  fbarr = fopen(barr_arg, "r");
  if (fbarr == NULL) {
    fprintf(stderr, "Cannot open file \"%s\".\n", barr_arg);
    exit(EXIT_FAILURE);
  }
  char *line = my_getline(fbarr);
  char *p = strtok(line, " \t\n");
  while (p) {
    if (isSeq(p)) {
      seq = (char*)malloc((strlen(p)+1)*sizeof(char));
      strcpy(seq, p);
      break;
    }
    p = strtok(NULL, " \t\n");
  }
  free(line);

  if (seq == NULL) {
    fprintf(stderr, "Couldn't find sequence on first line of file \"%s\"\n", barr_arg);
    fclose(fbarr);
    exit(EXIT_FAILURE);
  }
  // read previously found minima:
  while ((line = my_getline(fbarr))) {
    p = strtok(line, " \t\n");

    struct_en he;
    he.structure = NULL;
    he.energy = INT_MAX;

    // read the stuff
    int num;
    sscanf(p, "%d", &num);
    p = strtok(NULL, " \t\n");
    if (p && isStruct(p)) {
      he.structure = make_pair_table_PK(p);
    }
    p = strtok(NULL, " \t\n");
    if (p && he.structure && he.energy==INT_MAX) {
      float en;
      if (sscanf(p, "%f", &en)==1) {
        he.energy = en_fltoi(en);
      }
    }
    barr_info bi;
    int i =0;
    float en;
    while (p && he.structure && he.energy!=INT_MAX) {
      p = strtok(NULL, " \t\n");
      i++;
      switch (i){
        case 1: sscanf(p, "%d", &bi.father); break;
        case 2: sscanf(p, "%f", &en); bi.e_diff = en_fltoi(en); break;
        case 3: sscanf(p, "%d", &bi.bsize); break;
        case 4: sscanf(p, "%d", &bi.fbsize); break;
        case 5: sscanf(p, "%f", &bi.fen); break;
        case 6: sscanf(p, "%d", &bi.grad); break;
        case 7: sscanf(p, "%f", &bi.feng); break;
      }
    }

    // try to move it:
    SeqInfo sqi;
    sqi.Init(seq);
    he.energy = Opt.pknots ? energy_of_struct_pk(seq, he.structure, sqi.s0, sqi.s1, 0): energy_of_structure_pt(seq, he.structure, sqi.s0, sqi.s1, 0);
    int last_en = he.energy;
    move_set(he, sqi);
    /*
    //fprintf(stderr, "%f\n", he.energy);

    while ((Opt.rand? move_rand(he) : move_set(he))!=0) {
      Deg.Clear();
    }
    Deg.Clear();*/

    // print changes:
    if (last_en != he.energy) {
      fprintf(stderr, "%6.2f -> %6.2f %s\n", last_en/100.0, he.energy/100.0, pt_to_str_pk(he.structure).c_str());
    }

    // if we have it already:
    map<struct_en, barr_info, comps_entries>::iterator it;
    if ((it = output.find(he))!=output.end()) {
      it->second = it->second + bi;
      it->second.e_diff += last_en - he.energy;
      free(he.structure);
    } else {
      output[he] = bi;
    }

    free(line);
  }

  fclose(fbarr);
  return seq;
}

char *read_seq(char *seq_arg, char **name_out)
{
  char *name = NULL;
  char *seq = NULL;
  bool have_seq = false;
  // read sequence
  FILE *fseq;
  fseq = fopen(seq_arg, "r");
  if (fseq == NULL) {
    fprintf(stderr, "WARNING: Cannot open file \"%s\".\n", seq_arg);
    // try to read it from 1st line of input:
    seq = my_getline(stdin);
    int j = 0;
    for (unsigned int i=0; i<strlen(seq); i++) {
      if (seq[i]=='A' || seq[i]=='C' || seq[i]=='T' || seq[i]=='G' || seq[i]=='a' || seq[i]=='c' || seq[i]=='t' || seq[i]=='g' || seq[i]=='U' || seq[i]=='u') seq[j++] = seq[i];
    }
    seq[j]='\0';
    seq = (char*)realloc(seq, sizeof(char)*(j+1));
    if (j<5) {
      fprintf(stderr, "ERROR: Sequence not found in input nor in \"%s\" file.\n", seq_arg);
      exit(EXIT_FAILURE);
    } else {
      have_seq = true;
    }
  }
  if (!have_seq) {
    name = my_getline(fseq);
    if (name == NULL) {
      fprintf(stderr, "ERROR: File \"%s\" empty.\n", seq_arg);
      fclose(fseq);
      exit(EXIT_FAILURE);
    }
    seq = my_getline(fseq);
    if (seq == NULL || name[0]!='>') {
      //fprintf(stderr, "WARNING: File \"%s\" not in FASTA format. Using it as a sequence.\n", seq_arg);
      if (seq!=NULL) free(seq);
      seq = name;
      name = NULL;
    }
    // seq on more lines??
    char *seq2;
    while ((seq2=my_getline(fseq))!=NULL && isSeq(seq2)) {
      seq = (char*) realloc(seq, sizeof(char)*(strlen(seq)+strlen(seq2)+1));
      strcpy(seq+strlen(seq), seq2);
      free(seq2);
      //fprintf(stderr, "%s %d\n", seq, (int)strlen(seq));
    }

    fclose(fseq);
    if (name) (*name_out) = name;
  }

  // covert seq's T's to U's
  for (int i=0; i<(int)strlen(seq); i++) {
    if (seq[i]=='T') seq[i]='U';
    if (seq[i]==' ') seq[i]='\0';
  }
  return seq;
}


int move(unordered_map<struct_en, gw_struct, hash_fncts, hash_eq> &structs, map<struct_en, int, comps_entries> &output, set<struct_en, comps_entries> &output_shallow, SeqInfo &sqi, bool pure_output)
{
  // read a line
  char *line = my_getline(stdin);
  if (line == NULL) return -1;
  if (line[0]=='>') {
    free(line);
    return 0;
  }

  float energy=1e10;

  // process lines
  char *p = line;
  char *sep = " \t\n";
  char *temp = NULL;

  bool struct_found = false;
  bool energy_found = false;

  // read the structure
  p = strtok(line, sep);
  while(p!=NULL && !(struct_found && energy_found)) {
    //fprintf(stderr, "%s\n", p);
    if (isStruct(p)) {
      if (struct_found) fprintf(stderr, "WARNING: On line \"%s\" two structure-like strings found!\n", line);
      else {
        temp = p;
      }
      struct_found = true;
    } else {
      if (isEnergy(p, energy)) {
        energy_found = true;
      }
    }

    p = strtok(NULL, sep);
  }
  p = temp;
  if (!struct_found) {
    fprintf(stderr, "WARNING: On line \"%s\" no structure-like string found!\n", line);
    free(line);
    return 0;
  }

  // count moves
  num_moves++;

  // find length of structure
  int len=0;
  while (p[len]!='\0' && p[len]!=' ') len++;

  if (len!=seq_len) {
    fprintf(stderr, "WARNING: Unequal lengths:\n(structure) %s\n (sequence) %s\n", p, sqi.seq);
    free(line);
    return -0;
  }

  // make make_pair
  struct_en str;
  str.structure = Opt.pknots? make_pair_table_PK(p):make_pair_table(p);

  // only H,K,L,M types allowed:
  if (!str.structure) {
    free(line);
    return 0;
  } else {
    str.energy = Opt.pknots? energy_of_struct_pk(sqi.seq, str.structure, sqi.s0, sqi.s1, Opt.verbose_lvl>3):energy_of_structure_pt(sqi.seq, str.structure, sqi.s0, sqi.s1, 0);
    free(line);
  }

  // if pure, just do descend and print it:
  if (pure_output) {
    //is it canonical (noLP)
    if (Opt.noLP && find_lone_pair(str.structure)!=-1) {
      if (Opt.verbose_lvl>0) fprintf(stderr, "WARNING: structure \"%s\" has lone pairs, skipping...\n", pt_to_str_pk(str.structure).c_str());
      free(str.structure);
      return -2;
    }

    //debugging
    if (Opt.verbose_lvl>1) fprintf(stderr, "proc(pure): %d %s\n", num_moves, pt_to_str_pk(str.structure).c_str());

    // descend
    int gw_length = move_set(str, sqi);
    // only some types of PK allowed!!!
    if (Opt.pknots && str.energy == INT_MAX) {
      free(str.structure);
      return 0;
    }

    if (Opt.verbose_lvl>2) fprintf(stderr, "\n  %s %d %d\n", pt_to_str_pk(str.structure).c_str(), str.energy, gw_length);
    printf("%s %6.2f %4d\n", pt_to_str_pk(str.structure).c_str(), str.energy/100.0, gw_length);
    free(str.structure);
    return 1;
  }

  // check if it was before
  unordered_map<struct_en, gw_struct, hash_fncts, hash_eq>::iterator it_s = structs.find(str);

  // if it was - release memory + get another
  if (it_s != structs.end()) {
    it_s->second.count++;
    free(str.structure);
    return 0;
  } else {
    // find energy only if not in input (not working - does energy_of_move require energy_of_struct run first???)
    //str.energy = Enc.Energy(str);
    //printf("\n%s\n%s %7.2f\n\n", Enc.seq, pt_to_str_pk(str.structure).c_str(), str.energy/100.0);
    /*if (1 || !energy_found) str.energy = Enc.Energy(str);
    else str.energy = (int)(energy*100.0+(energy<0.0 ? -0.5 : 0.5));*/

    // allegiance hack:
    struct_en he_str = str;
    if (allegiance) {
      structures.push_back(he_str);
    }

    //is it canonical (noLP)
    if (Opt.noLP && find_lone_pair(str.structure)!=-1) {
      if (Opt.verbose_lvl>0) fprintf(stderr, "WARNING: structure \"%s\" has lone pairs, skipping...\n", pt_to_str_pk(str.structure).c_str());
      free(str.structure);
      return -2;
    }

    // copy it anew
    struct_en old = str;
    str.structure = allocopy(str.structure);

    //debugging
    if (Opt.verbose_lvl>1) fprintf(stderr, "processing: %d %s\n", num_moves, pt_to_str_pk(str.structure).c_str());

    // descend
    move_set(str, sqi);
    // only some types of PK allowed!!!
    if (Opt.pknots && str.energy == INT_MAX) {
      free(str.structure);
      free(old.structure);
      return 0;
    }

    // insert into hash (memory is here only on left side)
    gw_struct &lm = structs[old];
    lm.count = 1;
    /*
    int i;
    while ((i = (Opt.rand? move_rand(str) : move_set(str)))!=0) {
      Deg.Clear();
    }
    Deg.Clear();*/

    if (Opt.verbose_lvl>2) fprintf(stderr, "\n  %s %d\n", pt_to_str_pk(str.structure).c_str(), str.energy);

    // save for output
    map<struct_en, int, comps_entries>::iterator it;
    if ((it = output.find(str)) != output.end()) {
      it->second++;
      lm.he = it->first;
      free(str.structure);
      // allegiance hack:
      if (allegiance) str_to_LM[he_str] = it->first;
    } else {
      //str.num = output.size();
      lm.he = str;
      output.insert(make_pair(str, 1));
      // allegiance hack:
      if (allegiance) str_to_LM[he_str] = str;
    }
  }

  return 1;
}
