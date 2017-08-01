#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <unordered_map>
#include <stack>
#include <algorithm>
#include <vector>
#include <ctime>

extern "C" {
  #include "pair_mat.h"
  #include "fold.h"
  #include "loop_energies.h"

  #include "move_set_inside.h"
}

#include "pknots.h"

static float time_eos = 0.0;
static paramT *P = NULL;

void freeP()
{
  if (P) free(P);
  P = NULL;
}

float get_eos_time()
{
  return time_eos;
}

// ################ issue: TODO GU pair on end adds 0.5 penalty....

inline int Ext_Loop(int p, int q, short *s0, paramT *P)
{
  int energy = E_ExtLoop(pair[s0[p]][s0[q]], -1, -1, P);
  return energy;
}

inline int M_Loop(int p, int q, short *s0, paramT *P)
{
  int energy = E_MLstem(pair[s0[p]][s0[q]], -1, -1, P);
  return energy;
}

using namespace std;
/*
void copy_arr(short *dest, const short *src)
{
  if (!src || !dest) {
    fprintf(stderr, "Empty pointer in copying\n");
    return;
  }
  memcpy(dest, src, sizeof(short)*(src[0]+1));
}

short *allocopy(const short *src)
{
  short *res = (short*) malloc(sizeof(short)*(src[0]+1));
  copy_arr(res, src);
  return res;
}

char *allocopy(const char *src)
{
  char *res = (char*) malloc(sizeof(char)*(strlen(src)+1));
  strcpy(res, src);
  return res;
}*/

short *make_pair_table_PK(const char *str)
{
  const int maxp = 4;
 /* const char ptl[] = {"([{<"};
  const char ptr[] = {")]}>"};*/

  // stack of end points of each parenthesis
  vector<vector<int> > pars(maxp);

  int length = strlen(str);
  short *structure = (short*) malloc(sizeof(short)*(length+1));
  for (int i=1; i<=length; i++) {
    switch (str[i-1]) {
    case '.': structure[i] = 0; break;
    case '(': pars[0].push_back(i); break;
    case ')': structure[i] = pars[0].back(); structure[pars[0].back()] = i; pars[0].pop_back(); break;
    case '[': pars[1].push_back(i); break;
    case ']': structure[i] = pars[1].back(); structure[pars[1].back()] = i; pars[1].pop_back(); break;
    case '{': pars[2].push_back(i); break;
    case '}': structure[i] = pars[2].back(); structure[pars[2].back()] = i; pars[2].pop_back(); break;
    case '<': pars[3].push_back(i); break;
    case '>': structure[i] = pars[3].back(); structure[pars[3].back()] = i; pars[3].pop_back(); break;
    }
  }

  pars.clear();

  structure[0] = length;

  return structure;
}

string pt_to_str_pk(const short *str)
{
  // we can do with 4 types now:
  const int maxp = 4;
  const char ptl[] = {"([{<"};
  const char ptr[] = {")]}>"};

  // result:
  string res;
  res.reserve(str[0]);

  // stack of end points of each parenthesis
  vector<stack<int> > pars(maxp);

  // temprary structure
  vector<int> tmp(str[0]+1, 0);

  // precomputation
  for (int i=1; i<=str[0]; i++) {
    if (str[i]>i) {
      int j=0;
      //fprintf(stderr, "%d %d \n", (int)pars[j].size(), pars[j].empty());
      while (j<maxp && !pars[j].empty() && pars[j].top()<str[i]) {
        j++;
      }
      if (j==maxp) {
        fprintf(stderr, "Cannot print it with %d types of parentheses!!!\n", maxp);
        throw;
        return res;
      } else {
        // jth parenthesis:
        pars[j].push(str[i]);
        tmp[i] = j;
        tmp[str[i]] = j;
      }
    } else if (str[i]>0 && str[i]<i) {
      pars[tmp[i]].pop();
    }
  }

  // filling the result:
  for (int i=1; i<=str[0]; i++) {
    if (str[i]==0) res += '.';
    else if (str[i]>i) res += ptl[tmp[i]];
    else res += ptr[tmp[i]];
  }
  res+='\0';

  return res;
}

char* pt_to_chars_pk(const short *str)
{
  // we can do with 4 types now:
  const int maxp = 4;
  const char ptl[] = {"([{<"};
  const char ptr[] = {")]}>"};

  // result:
  char *res = (char*) malloc((str[0]+1)*sizeof(char));

  // stack of end points of each parenthesis
  vector<stack<int> > pars(maxp);

  // temprary structure
  vector<int> tmp(str[0]+1, 0);

  // precomputation
  for (int i=1; i<=str[0]; i++) {
    if (str[i]>i) {
      int j=0;
      //fprintf(stderr, "%d %d \n", (int)pars[j].size(), pars[j].empty());
      while (j<maxp && !pars[j].empty() && pars[j].top()<str[i]) {
        j++;
      }
      if (j==maxp) {
        fprintf(stderr, "Cannot print it with %d types of parentheses!!!\n", maxp);
        free(res);
        return NULL;
      } else {
        // jth parenthesis:
        pars[j].push(str[i]);
        tmp[i] = j;
        tmp[str[i]] = j;
      }
    } else if (str[i]>0 && str[i]<i) {
      pars[tmp[i]].pop();
    }
  }

  // filling the result:
  for (int i=1; i<=str[0]; i++) {
    if (str[i]==0) res[i-1] = '.';
    else if (str[i]>i) res[i-1] = ptl[tmp[i]];
    else res[i-1] = ptr[tmp[i]];
  }
  res[str[0]] ='\0';

  return res;
}
void pt_to_chars_pk(const short *str, char *dest)
{
  // we can do with 4 types now:
  const int maxp = 4;
  const char ptl[] = {"([{<"};
  const char ptr[] = {")]}>"};

  // stack of end points of each parenthesis
  vector<stack<int> > pars(maxp);

  // temprary structure
  vector<int> tmp(str[0]+1, 0);

  // precomputation
  for (int i=1; i<=str[0]; i++) {
    if (str[i]>i) {
      int j=0;
      //fprintf(stderr, "%d %d \n", (int)pars[j].size(), pars[j].empty());
      while (j<maxp && !pars[j].empty() && pars[j].top()<str[i]) {
        j++;
      }
      if (j==maxp) {
        fprintf(stderr, "Cannot print it with %d types of parentheses!!!\n", maxp);
        return ;
      } else {
        // jth parenthesis:
        pars[j].push(str[i]);
        tmp[i] = j;
        tmp[str[i]] = j;
      }
    } else if (str[i]>0 && str[i]<i) {
      pars[tmp[i]].pop();
    }
  }

  // filling the result:
  for (int i=1; i<=str[0]; i++) {
    if (str[i]==0) dest[i-1] = '.';
    else if (str[i]>i) dest[i-1] = ptl[tmp[i]];
    else dest[i-1] = ptr[tmp[i]];
  }
  dest[str[0]] ='\0';
}

int loop_energy_pk(int begin, short *str, short *s0, short *s1, paramT *P, bool &multiloop) {

  int end = str[begin];
  if (begin==0) end++;
  int alone = end - begin - 1;
  int energy = 0;
  int alone_upto = begin;
  int loops = 0;

  //fprintf(stderr, "le_pk %d %d\n", begin, multiloop);

  for (int i=begin+1; i<end; i++) {
    ///TODO dangles == 0
    if (str[i]>i && str[i]>alone_upto && str[i]<end) { // '(' ,but not nested '('
      if (multiloop) {
        energy += M_Loop(i, str[i], s0, P);
        loops ++;
        // and subtract alone positions
        if (alone_upto<i) {
          alone -= str[i] - i +1;
        } else {
          alone -= str[i] - alone_upto;
        }
        alone_upto = str[i];
      } else {
        energy += Ext_Loop(i, str[i], s0, P);
        //fprintf(stderr, "le_pk EXT %d %d %d\n", i, str[i], energy);
        alone_upto = str[i];
      }

      //if (tree->debug) fprintf(stderr, "eval EN: %3d %3d =>  %7.2f  (%7.2f -> %7.2f)\n", it->second->begin, it->second->end, (energy - old_energy)/100.0, old_energy/100.0, energy/100.0);
    }
  }

  // case that we don't have a multiloop
  if (loops == 1 && multiloop) {
    multiloop = false;
    return loop_energy(str, s0, s1, begin);
  }

  if (multiloop) {
    energy += M_Loop(begin, end, s0, P);
    energy += P->MLclosing;
    /* logarithmic ML loop energy if logML */
    if (logML && (alone>6)) {
      energy += 6*P->MLbase+(int)(P->lxc*log((double)alone/6.));
    } else {
      energy += (alone*P->MLbase);
    }
  }

  return energy;
}

struct LE_ret {
  int energy;
  int ending;  // my ending - only the end of the simple bp
  int pk_ending; // ending of whole pknot
  int pk_loops; // number of pknot loops
};

void try_pk() {
  energy_of_struct_pk("CCCCCCCCCCGGCCCCGGGGGGGGG",
                      "(((..((...)).[[[)))...]]]", 4);
}

// checks if we are inside of pknot or inside something else (return true if inside pknot)
bool OKStack(vector<int> &stack, int pk[4]) {
  if (stack.size() == 0) return true;
  int data = stack.back();
  return data==pk[0] || data==pk[1] || data==pk[2] || data==pk[3];
}

bool OKStacks(vector<int> stacks[3], int pk[4]) {
  return OKStack(stacks[0], pk) && OKStack(stacks[1], pk) && OKStack(stacks[2], pk);
}

LE_ret loop_energy_rec(int i, short *str, short *s0, short *s1, paramT *P, Helpers &hlps, int encircling = 0)
{
  // first find the type and count alone points:
  int begin = i;
  int true_begin = begin;
  int bpairs = 0;
  int last_in_begin = hlps.last_open;
  hlps.last_open = i;

  hlps.str_toleft[i] = i;
  hlps.str_torght[i] = i;

  // for pk counting
  int opening = 1;
  int max_op = 1;

  int last_pk = 0;
  int pk[4] = {i,0,0,0};
  int pk_outer[4] = {i,0,0,0};

  int reg_inloops = 0;
  int next_i;

  //int ending = str[begin];
  i++;

  #define END() str[pk_outer[last_pk]]

  LE_ret ret;
  ret.energy = 0;
  ret.ending = str[begin];
  ret.pk_loops = 1;
  ret.pk_ending = ret.ending;

  // debug:
  //fprintf(stderr, "%4d called loop_energy_rec %s\n", begin, pt_to_str_pk(str).c_str());

  // if we have done it before:
  /*if (0 && str_type[i]!=ROOT) {
    ret.ending = ending;
    ret.pknot = str_type[i]>=P_H?true:false;
    ret.energy = str_energy[i];
    return ret;
  }*/

  // and start crawling to the right:
  while (i < END()) {
    next_i = i;
    if (str[i] == 0) {
    } else
    if (str[i] > i) {
      if (str[i] < END()) { // regular inloop - skip it
        reg_inloops++;
        if (last_pk>0 && (i < str[pk[last_pk-1]] && str[i] > str[pk[last_pk-1]])) { // we had wrong pk
          bool multiloop = true;
          int energy_tmp = loop_energy_pk(pk[last_pk], str, s0, s1, P, multiloop);
          hlps.str_type[pk[last_pk]] = multiloop?N_M:N_S;
          hlps.str_energy[pk[last_pk]] = energy_tmp;

          ret.pk_ending = max(ret.pk_ending, (int)str[pk[last_pk]]);
          //ending = str[pk[last_pk]];
          ret.energy += energy_tmp;
          pk[last_pk] = i;
        } else { // do everything inbetween i, str[i] (+ maybe more)
          LE_ret ret_tmp = loop_energy_rec(i, str, s0, s1, P, hlps, max(encircling, ret.ending));
          next_i = ret_tmp.pk_ending;
          bpairs++;
          // if we are encircling a pseudoknot, we have to add 2+ loops:
          if (ret_tmp.pk_loops > 1 && ret_tmp.pk_ending<ret.ending) {
            reg_inloops += ret_tmp.pk_loops-1;
          } else { // or we are extending a pknot, then we have to keep the pknot info for a multiloop above.
            ret.pk_ending = max(ret.pk_ending, ret_tmp.pk_ending);
            ret.pk_loops = max(ret.pk_loops, ret_tmp.pk_loops);
          }
          ret.energy += ret_tmp.energy;
        }

      } else if (str[i] > END()) { // now we know we are in pknot / we found another part of pknot
        hlps.str_torght[hlps.last_open] = hlps.last_open;
        hlps.last_open = 0;
        // pknot identifying
        last_pk++;
        opening++;
        if (max_op < opening) max_op = opening;
        pk[last_pk] = i;
        pk_outer[last_pk] = i;
        ret.pk_ending = max(ret.pk_ending, (int)str[i]);
      }
      // now we can do the left & right stacking:
      if (hlps.last_open != 0) {
        hlps.str_toleft[i] = hlps.last_open;
        hlps.str_torght[hlps.last_open] = i;
      } else {
        hlps.str_toleft[i] = i;
        hlps.str_torght[i] = i;
      }

      hlps.last_open = i;
    } else
    if (str[i] < i) { // should happen only in pknot:
      hlps.str_torght[hlps.last_open] = hlps.last_open;
      hlps.last_open = 0;
      for (int j=0; j<last_pk; j++) if (str[i]==pk[j]) {
        opening--;
        break;
      }
      if (str[i]<true_begin) true_begin = str[i];
    }
    i = next_i +1;
  }

  hlps.last_open = last_in_begin;
  int energy;

  // return value 3 options:
  BPAIR_TYPE type = P_H;
  if (last_pk == 0) {
    if (reg_inloops <= 1) {
      energy = loop_energy(str, s0, s1, begin);
      ret.energy += energy;
      type = N_S;
    } else {
      bool multiloop = true;
      energy = loop_energy_pk(begin, str, s0, s1, P, multiloop);
      ret.energy += energy;
      type = N_M;
    }

  } else {

    if (last_pk==3) type = P_M;
    if (last_pk==2 && max_op==2) type = P_K;
    if (last_pk==2 && max_op==3) type = P_L;
    if (encircling > ret.pk_ending) {
      energy = beta1mp_pen[type];
    } else energy = beta1_pen[type];

    // compute alone:
    int alone = 0;
    vector<int> stacks[3];
    for (int j=true_begin; j<=ret.pk_ending; j++) {
      if (str[j] == 0) {
        if (OKStacks(stacks, pk)) {
          alone++;
        }
      } else if (str[j] > j) {
        int k=0;
        while (stacks[k].size() > 0 && str[j]>str[stacks[k].back()]) k++;
        stacks[k].push_back(j);
      } else if (str[j] < j) {
        int k=0;
        while (k<3 && (stacks[k].size() == 0 || str[j] != stacks[k].back())) k++;
        if (k<3) stacks[k].pop_back();
      }
    }
    hlps.beta1 = energy;
    energy += beta2_pen[type]*(bpairs+last_pk+1) + beta3_pen[type]*(alone);
    hlps.beta2 = bpairs+last_pk+1;
    hlps.beta3 = alone;

    // weird 0.1:
    if (type == P_L || type == P_K || type == P_M) energy -= 10;

    ret.energy += energy;
    ret.pk_loops = max(last_pk+1, ret.pk_loops);
  }

  // assign a type:
  for (int j=0; j<=last_pk; j++) {
    hlps.str_type[pk[j]] = type;
  }
  hlps.str_energy[pk[0]] = energy;

  //fprintf(stderr, "%4d  ended loop_energy_rec\n", begin);

  return ret;
}

Helpers::Helpers(int length)
{
  Create(length);
}

void Helpers::Create(int length)
{
  str_energy.resize(length + 1, 0);
  str_type.resize(length+1, ROOT);
  str_toleft.resize(length+1, 0);
  str_torght.resize(length+1, 0);
  beta1 = beta2 = beta3 = last_open = 0;
}

int energy_of_struct_pk(const char *seq, char *structure, int verbose)
{
  if (P == NULL) {
    make_pair_matrix();
    update_fold_params();
    P = scale_parameters();
  }
  short *str = make_pair_table_PK(structure);
  int res = energy_of_struct_pk(seq, str, verbose);
  free(str);
  return res;
}

int energy_of_struct_pk(const char *seq, short *structure, int verbose)
{
  if (P == NULL) {
    make_pair_matrix();
    update_fold_params();
    P = scale_parameters();
  }
  short *s0 = encode_sequence(seq, 0);
  short *s1 = encode_sequence(seq, 1);

  int res = energy_of_struct_pk(seq, structure, s0, s1, verbose);

  free(s1);
  free(s0);

  return res;
}

int energy_of_struct_pk(const char *seq, short *structure, short *s0, short *s1, int verbose)
{
  clock_t time = clock();

  if (P == NULL) {
    make_pair_matrix();
    update_fold_params();
    P = scale_parameters();
  }
  short *str = structure;

  // some debug/helper arrays:
  Helpers hlps(str[0]);

  int energy = 0;

  // go through the structure:
  for (int i=1; i<=str[0]; i++) {
    if (str[i]>i && hlps.str_type[i]==ROOT) {  //'('
      // found new loop, evaluate the energy
      LE_ret ret = loop_energy_rec(i, str, s0, s1, P, hlps);
      i = ret.ending;
      energy += ret.energy;
      // just verbose
      if (verbose && hlps.beta1 != 0) {
        fprintf(stderr, "found pknot: %4d %4d %4d\n", hlps.beta1, hlps.beta2, hlps.beta3);
      }
    }
  }

  //add the energy of external loop
  bool multiloop = false;
  int ext = loop_energy_pk(0, str, s0, s1, P, multiloop);
  energy += ext;

  string type;
  for (int i=1; i<=str[0]; i++) {
    type += bpair_type_sname[hlps.str_type[i]];
  }

  if (verbose) fprintf(stderr, "%s\n%s", seq, pt_to_str_pk(structure).c_str());
  if (verbose) fprintf(stderr, " %7.2f\n%s\n", energy/100.0, type.c_str());

  if (verbose) {
    int summ = ext;
    fprintf(stderr, "%7.2f EXT\n", ext/100.0);
    for (int i=1; i<=str[0]; i++) {
      fprintf(stderr, "%7.2f %5d %5d\n", hlps.str_energy[i]/100.0, hlps.str_toleft[i], hlps.str_torght[i]);
      summ += hlps.str_energy[i];
    }

    fprintf(stderr, "%7.2f summed\n", summ/100.0);
  }

  float time_tmp = (clock()-time)/(double)CLOCKS_PER_SEC;
  //fprintf(stderr, "time_tmp = %10g\n", time_tmp);
  time_eos += time_tmp;

  return energy;
}


Pseudoknot::Pseudoknot()
{
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) imat[i][j] = false;

  size = 0;
}
/*
// constructor
Pseudoknot::Pseudoknot(int left, int right, int k)
{
  assert(left>0);
  if (left > right) swap(left, right);
  if (k && left > k) swap(left, k);
  if (k && k > right) swap(k, right);

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) imat[i][j] = false;

  int pos_right = k?2:1;

  parts[0].insert(left);
  parts[pos_right].insert(right);

  imat[0][pos_right] = 1;

  if (k) {
    parts[1].insert(k);
    imat[0][1] = 1;
    imat[1][pos_right] = 1;
  }
}*/

// contains some bpair?
bool Pseudoknot::Contains(int left)
{
  return parts[0].count(left) || parts[1].count(left) || parts[2].count(left)|| parts[3].count(left);
}

/*bool cross(int left, int right, int left2, int right2) {
  return (left<left2 && right<right2 && left2<right) ||
         (left>left2 && right>right2 && right2>left);
}

bool Pseudoknot::Cross(short *str, int left, int right)
{
  int left2 = *parts[0].begin();
  if (cross(left, right, left2, str[left2])) return true;
  left2 = *parts[1].begin();
  if (cross(left, right, left2, str[left2])) return true;
  if (parts[2].empty()) return false;
  left2 = *parts[2].begin();
  if (cross(left, right, left2, str[left2])) return true;
  if (parts[3].empty()) return false;
  left2 = *parts[3].begin();
  if (cross(left, right, left2, str[left2])) return true;
  return false;
}*/


  // can we insert it inside?
int Pseudoknot::Inside(short *str, int left, int right)
{
  for (int i=0; i<4; i++) {
    if (parts[i].empty()) return -1;

    set<int>::iterator it = parts[i].upper_bound(left);

    // check if we are inside :   it--(..left[..it( ... it)..right]..it--)
    if ((it == parts[i].end() || str[*it]<right) &&
        (it == parts[i].begin() || str[*(--it)]>right)) {
          return i;
    }
  }

  return -1;
}

bool IsViable(int size, bool imat[4][4]) {
/*PH    PK    PL    PM
  0100  0100  0110  0110
  0000  0010  0010  0011
  0000  0000  0000  0001
  0000  0000  0000  0000*/

  if (size == 2) return true;
  if (size == 3) return true;
  if (size == 4) {
    return imat[0][1] && imat[0][2] && !imat[0][3] &&
           imat[1][2] && imat[1][3] &&
           imat[2][3];
  }
  return false;
}

Pseudoknot::Pseudoknot(const Pseudoknot &pknot)
{
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) imat[i][j] = pknot.imat[i][j];
    parts[i] = pknot.parts[i];
  }
  size = pknot.size;
}

bool Pseudoknot::CanInsert(int left, vector<int> &numbers, bool insert)
{
  if (size==4) return false;

  bool imat2[4][4];
  int new_pos = 0;
  while (!parts[new_pos].empty() && *parts[new_pos].begin()<left) new_pos++;
  for (int i=0; i<size; i++) {
    for (int j=i+1; j<size; j++) {
      int ii = i >= new_pos ? i+1:i;
      int jj = j >= new_pos ? j+1:j;
      imat2[ii][jj] = imat[i][j];
    }
  }

  for (int i=0; i<size+1; i++) imat2[min(new_pos, i)][max(new_pos, i)] = false;
  for (unsigned int i=0; i<numbers.size(); i++) {
    int number = numbers[i]>=new_pos ? numbers[i]+1:numbers[i];
    imat2[min(new_pos, number)][max(new_pos,number)] = true;
  }

  bool viable = IsViable(size+1, imat2);

  if (viable && insert) {
    for (int i=3; i>new_pos; i--) {
      swap(parts[i], parts[i-1]);
    }
    parts[new_pos].insert(left);
    size++;

    for (int i=0; i<size; i++) {
      for (int j=i+1; j<size; j++) {
        imat[i][j] = imat2[i][j];
      }
    }
  }

  return viable;
}

bool Pseudoknot::Delete(int left)
{
  for (int i=0; i<4; i++) {
  // try to find it
    set<int>::iterator it = parts[i].find(left);
    if (it != parts[i].end()) {
      // delete it
      parts[i].erase(it);
      // change type?
      if (parts[i].empty()) {
        size--;
        if (size==1) return true; // cancelled H type
        if (size==2 && i==1 && !imat[0][2]) { // cancelled K type
          size--;
          return true;
        }
        for (int j=0; j<size; j++) {
          int jj = j>=i?j+1:j;
          for (int k=j+1; k<size; k++) {
            int kk = k>=i?k+1:k;
            imat[j][k] = imat[jj][kk];
          }
          parts[j] = parts[jj];
        }
        parts[size].clear();
        imat[0][size] = imat[1][size] = imat[2][size] = false;
        return true;
      }
      return true;
    }
  }

  return false;
}

Pseudoknot &Pseudoknot::operator=(const Pseudoknot &pknot)
{
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) this->imat[i][j] = pknot.imat[i][j];
    this->parts[i] = pknot.parts[i];
  }
  this->size = pknot.size;

  return *this;
}

Structure::Structure(int length)
{
  this->str = (short*) malloc(sizeof(short)*(length+1));
  for (int i=1; i<=length; i++) this->str[i] = 0;
  this->str[0] = length;

  this->energy = 0;
}

Structure::Structure(short *structure, int energy)
{
  int length = structure[0];
  this->str = (short*) malloc(sizeof(short)*(length+1));
  for (int i=1; i<=length; i++) this->str[i] = 0;
  this->str[0] = length;

  // assign all bpairs:
  for (int i=1; i<=structure[0]; i++) {
    if (structure[i]>i) ViableInsert(i, structure[i], true);
  }

  this->energy = energy;
}

Structure::Structure(char *structure, int energy)
{
  int length = strlen(structure);
  this->str = (short*) malloc(sizeof(short)*(length+1));
  for (int i=1; i<=length; i++) this->str[i] = 0;
  this->str[0] = length;

  // assign all bpairs:
  short *str_tmp = make_pair_table_PK(structure);
  for (int i=1; i<=str_tmp[0]; i++) {
    if (str_tmp[i]>i) ViableInsert(i, str_tmp[i], true);
  }
  free(str_tmp);

  this->energy = energy;
}

Structure::Structure(const char *seq, short *structure, short *s0, short *s1)
{
  int length = structure[0];
  this->str = (short*) malloc(sizeof(short)*(length+1));
  for (int i=1; i<=length; i++) this->str[i] = 0;
  this->str[0] = length;

  // assign all bpairs:
  for (int i=1; i<=structure[0]; i++) {
    if (structure[i]>i) ViableInsert(i, structure[i], true);
  }

  energy = energy_of_struct_pk(seq, str, s0, s1, 0);
}

Structure::Structure(const char *seq, char *structure, short *s0, short *s1)
{
  int length = strlen(structure);
  this->str = (short*) malloc(sizeof(short)*(length+1));
  for (int i=1; i<=length; i++) this->str[i] = 0;
  this->str[0] = length;

  // assign all bpairs:
  short *str_tmp = make_pair_table_PK(structure);
  for (int i=1; i<=str_tmp[0]; i++) {
    if (str_tmp[i]>i) ViableInsert(i, str_tmp[i], true);
  }
  free(str_tmp);

  energy = energy_of_struct_pk(seq, str, s0, s1, 0);
}

Structure::Structure(const Structure &second)
{
  energy = second.energy;
  str = allocopy(second.str);

  pknots = second.pknots;

  bpair_pknot = second.bpair_pknot;
}

Structure::~Structure()
{
  free(str);
}

bool const Structure::operator<(const Structure &second) const
{
  if (energy != second.energy) return energy<second.energy;
  int i=1;
  char l=0,r=0;
  while (i<=str[0]) {
    l = (str[i]==0?'.':(str[i]<str[str[i]]?')':'('));
    r = (second.str[i]==0?'.':(second.str[i]<second.str[second.str[i]]?')':'('));
    if (l != r) return l<r;
    //if (str[i] != second.str[i]) return str[i] > second.str[i];
    i++;
  }

  return false;
}

bool const Structure::operator==(const Structure &second) const
{
  if (energy != second.energy) return false;
  int i=1;
  while (i<=str[0]) {
    if (str[i] != second.str[i]) return false;
    i++;
  }
  return true;
}

Structure &Structure::operator=(const Structure &second)
{
  energy = second.energy;
  copy_arr(str, second.str);

  pknots = second.pknots;

  bpair_pknot = second.bpair_pknot;

  return *this;
}

// checks if we can add it to the 'cross' vector (if it is compatible with )
bool CrossOuter(vector<int> &cross, int point, short *str) {
  if (cross.empty()) return true;
  int left = cross.back();
  int right = str[left];
  if (point<left && str[point]>right) return true;
  return false;
}

bool CrossInner(vector<int> &cross, int point, short *str) {
  if (cross.empty()) return true;
  int left = cross.back();
  int right = str[left];
  if (point>left && str[point]<right) return true;
  return false;
}

// is the bpair viable?
INS_FLAG Structure::ViableInsert(int left, int right, bool insert)
{
  if (str[left] || str[right]) return NO_INS;

  INS_FLAG res = NO_INS;

  // get crossings:
  vector<set<int> > crossings;
  vector<vector<int> > cross_tmp;
  for (int i=left+1; i<right; i++) { // can be faster:
    if (str[i] == 0) continue;
    if (str[i] < left) { // ')'
      int k=0;
      if ((int)cross_tmp.size()<=k) cross_tmp.resize(k+1);
      if ((int)crossings.size()<=k) crossings.resize(k+1);
      while (!CrossOuter(cross_tmp[k], str[i], str)) {
        k++;
        if ((int)cross_tmp.size()<=k) cross_tmp.resize(k+1);
        if ((int)crossings.size()<=k) crossings.resize(k+1);
      }
      cross_tmp[k].push_back(str[i]);
      crossings[k].insert(str[i]);
    } else if (str[i] > right) { // '('
      int k=0;
      if ((int)cross_tmp.size()<=k) cross_tmp.resize(k+1);
      if ((int)crossings.size()<=k) crossings.resize(k+1);
      while (!CrossInner(cross_tmp[k], i, str)) {
        k++;
        if ((int)cross_tmp.size()<=k) cross_tmp.resize(k+1);
        if ((int)crossings.size()<=k) crossings.resize(k+1);
      }
      cross_tmp[k].push_back(i);
      crossings[k].insert(i);
    }
  }
  // sort the sets:
  sort(crossings.begin(), crossings.end());

  // find if it crosses a pknot:
  Pseudoknot *cross_pk = NULL;
  int cross_index = -1;
  if (pknots.size()>0) {
    for (unsigned int i=0; i<crossings.size(); i++) {
      set<int>::iterator it = crossings[i].begin();



      int index = bpair_pknot[*it];
      Pseudoknot *pk_tmp = Pknot_index(index);
      for (it++; it!=crossings[i].end(); it++) {
        Pseudoknot *pk_now = Pknot_bpair(*it);
        if (pk_tmp != pk_now) return NO_INS; // split the stack!
      }
      if (pk_tmp != NULL) {
        if (cross_pk == NULL) {
          cross_pk = pk_tmp;
          cross_index = index;
        } else if (cross_pk != pk_tmp) return NO_INS; // we have found 2 colliding pknots
      }
    }
  }

  /*for (unsigned int i=0; i<pknots.size(); i++) {
    if (pknots[i].Cross(str, left, right)) {
      if (cross_pk == NULL) cross_pk = &pknots[i];
      else return false;
    }
  }*/

  if (cross_pk) {
    // check if we don't split the parts
    vector<int> numbers;
    for (unsigned int i=0; i<crossings.size(); i++) {
      for (int j=0; j<cross_pk->size; j++) {
        if (*cross_pk->parts[j].begin() == *crossings[i].begin()) {
          if (cross_pk->parts[j] == crossings[i]) {
            numbers.push_back(j);
            // stack is correct
            break;
          }
          // stack is split
          return NO_INS;
        }
      }
      // each crossing we must find otherwise wrong.
      if (numbers.size()==i) return NO_INS;
    }

    // now we have the numbers which we border (0-3), the rest should be checked if we can insert it in it
    int inside = cross_pk->Inside(str, left, right);
    if (inside>=0) {
      // check if neighbourhood fits  // essentially numbers[x] == Imat(x, inside)
      int k=0;
      for (int i=0; i<cross_pk->size; i++) {
        if (cross_pk->Imat(inside, i)) {
          if (k>(int)numbers.size()-1) return NO_INS;
          if (numbers[k]==i) k++;
          else return NO_INS;
        }
      }
      // lastly check if we do not cross more than we should
      if (k != (int)numbers.size()) return NO_INS;

      // neighbourhood fits -> we can insert it in the parts[inside]
      if (insert) cross_pk->parts[inside].insert(left);
      res = INSIDE_PK;
    } else {
      // we cannot insert it, but maybe we can change the type:
      bool ci = cross_pk->CanInsert(left, numbers, insert);
      if (ci) res = CHNG_PK;
      else res = NO_INS;
    }

  } else {
    // now we can insert it since we do not cross any pk.
    if (insert) {
      switch (crossings.size()) {
      case 0:
        break;
      case 1: {
        Pseudoknot pknot;
        int cross_first = *crossings[0].begin();
        int bpair = (left > cross_first)? 1:0;
        pknot.parts[bpair].insert(left);
        pknot.parts[bpair==0?1:0] = crossings[0];
        pknot.imat[0][1] = true;
        pknot.size = 2;
        pknots.push_back(pknot);
        cross_index = pknots.size()-1;
        for (set<int>::iterator it=crossings[0].begin(); it!=crossings[0].end(); it++) bpair_pknot[*it] = pknots.size()-1;
        break;
        }
      case 2: {
        Pseudoknot pknot;
        int bpair = (left > *crossings[0].begin())? ((left > *crossings[1].begin())? 2:1):0;
        pknot.parts[bpair].insert(left);
        pknot.parts[bpair<=0?1:0] = crossings[0];
        pknot.parts[bpair<=1?2:1] = crossings[1];
        pknot.imat[0][1] = true;
        pknot.imat[1][2] = true;
        pknot.size = 3;
        pknots.push_back(pknot);
        cross_index = pknots.size()-1;
        for (set<int>::iterator it=crossings[0].begin(); it!=crossings[0].end(); it++) bpair_pknot[*it] = pknots.size()-1;
        for (set<int>::iterator it=crossings[1].begin(); it!=crossings[1].end(); it++) bpair_pknot[*it] = pknots.size()-1;
        break;
        }
      default: assert(false);
      }
    }
    switch (crossings.size()) {
      case 0: res = REG_INS; break;
      case 1:
      case 2: res = CREATE_PK; break;
    }
  }

  if (res != NO_INS && insert) {
    bpair_pknot[left] = cross_index;
    str[left] = right;
    str[right] = left;
  }
  return res;
}

bool Structure::Delete(int left)
{
  if (left<0) left = -left;
  if (!str[left]) return false;

  int right = str[left];

  for (unsigned int i=0; i<pknots.size(); i++) {
    if (pknots[i].Delete(left)) {
      if (pknots[i].size == 1) {
        for (int j=0; j<4; j++) {
          for (set<int>::iterator it = pknots[i].parts[j].begin(); it!=pknots[i].parts[j].end(); it++) {
            bpair_pknot[*it] = -1;
          }
        }
        // erase the pknot and fix bpair_pknot
        pknots.erase(pknots.begin()+i);
        for (map<int, int>::iterator it=bpair_pknot.begin(); it!=bpair_pknot.end(); it++) {
          if (it->second >= (int)i) it->second--;
        }
      }
      break;
    }
  }

  bpair_pknot.erase(left);
  str[left] = 0;
  str[right] = 0;
  return true;
}

int Structure::UndoMove()
{
  if (undo_l == 0) {
    fprintf(stderr, "ERROR: Undoing non-existent move!\n");
    return 0;
  }

  int left = undo_l;
  int right = undo_r;

  int dif_en = - energy + undo_en;
  energy = undo_en;

  // switch?
  if (left>0 && right>0 && (str[left]>0 || str[right]>0)) {
    Shift(left, right);
  } else {
    if (left<0) {
      Delete(left);
    } else {
      Insert(left, right);
    }
  }

  undo_l = 0;
  undo_r = 0;
  undo_en = 0;

  return dif_en;
}

int Structure::MakeMove(const char *seq, short *s0, short *s1, int left, int right)
{
  undo_en = energy;
  bool change = false;

  // switch?
  if (left>0 && right>0 && (str[left]>0 || str[right]>0)) {
    if (str[left]>0) {
      undo_l = left;
      undo_r = str[left];
    } else {
      undo_r = right;
      undo_l = str[right];
    }
    change = (Shift(left, right) != NO_INS);
  } else {
    if (left<0) {
      undo_l = -left;
      undo_r = -right;
      change = Delete(left);
    } else {
      undo_l = -left;
      undo_r = -right;
      change = (Insert(left, right) != NO_INS);
    }
  }

  if (change) energy = energy_of_struct_pk(seq, str, s0, s1, 0);
  return energy - undo_en;
}

INS_FLAG Structure::CanShift(int left, int right)
{
  if ((str[left]>0 && str[right]>0)  ||
      (str[left]==0 && str[right]==0)) return NO_INS;

  int left_orig = str[left]>0? left : str[right];
  int right_orig = str[right]>0? right : str[left];

  Delete(left_orig);

  INS_FLAG flag = CanInsert(left, right);

  Insert(left_orig, right_orig);

  return flag;
}

INS_FLAG Structure::Shift(int left, int right)
{
  if ((str[left]>0 && str[right]>0)  ||
      (str[left]==0 && str[right]==0)) return NO_INS;

  int left_orig = str[left]>0? left : str[right];
  int right_orig = str[right]>0? right : str[left];

  Delete(left_orig);

  INS_FLAG flag = Insert(left, right);

  if (flag==NO_INS) Insert(left_orig, right_orig);

  return flag;
}

INS_FLAG Structure::CanInsert(int left, int right)
{
  return ViableInsert(left, right, false);
}

INS_FLAG Structure::Insert(int left, int right)
{
  return ViableInsert(left, right, true);
}

int Contains_PK(short *str)
{
  stack<int> pairing;
  pairing.push(str[0]+1);
  for (int i=1; i<=str[0]; i++) {
    if (str[i]==0) continue;
    if (str[i]>i) { //'('
      if (str[i]>pairing.top()) return i;
      pairing.push(str[i]);
    }
    else { //')'
      pairing.pop();
    }
  }
  return 0;
}
