#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <queue>


#include "findpath_pk.h"
#include "move_set_pk.h"

extern "C" {
  #include "pair_mat.h"
}

using namespace std;

/**
 *  \brief
 */
struct move_fp {
  int left;  /* i,j>0 insert; i,j<0 delete */
  int right;
  int when;  /* 0 if still available, else resulting distance from start */
  int dE;     // change in energy

  move_fp(int i, int j) {
    if (abs(i)>abs(j)) swap(i, j);
    left = i;
    right = j;
    dE = 0;
    when = 0;
  }

  bool const operator==(const move_fp &second) const {
    return (second.left == left) && (second.right == right);
  }

  bool const operator<(const move_fp &second) const {
    if (second.left == left) return right < second.right;
    else return left < second.left;
  }
};

class intermediate_pk;

// compare intermeds by a structure.
struct compare_intermeds {
  bool operator()(const intermediate_pk &left, const intermediate_pk &right) const;
};

struct compare_struct {
  bool operator()(const short *lhs, const short *rhs) const {
    int i=1;
    short l=0,r=0;
    while (i<=lhs[0]) {
      l = lhs[i];
      r = rhs[i];
      if (l != r) {
        //fprintf (stderr, "%s\n%s\n[%d] %c %c %d %d\n", pt_to_str_pk(lhs).c_str(), pt_to_str_pk(rhs).c_str(), i, l, r, l, r);
        break;
      }
      i++;
    }
    return (i<=lhs[0] && l<r);
  }
};

class intermediate_pk {
public:
  short *structure;      /**<  \brief  pair table */
  int Sen;              /**<  \brief  saddle energy so far */
  int energy;          /**<  \brief  current energy */
  int dist;           /** distance to the str1 */

  vector<move_fp> moves_done;
  set<move_fp> moves_todo;

  Structure pknot; // maybe for faster computation

  // comparator for priority queue
  bool const operator<(const intermediate_pk &second) const {
    if (dist == second.dist) {
      if (Sen == second.Sen) {
        if (energy == second.energy) {
          compare_struct cs;
          return cs(structure, second.structure);
        } else return energy>second.energy;
      } else return Sen>second.Sen;
    } else return dist>second.dist;
  }

  intermediate_pk &operator=(const intermediate_pk &second)
  {
    structure = second.structure;
    Sen = second.Sen;
    energy = second.energy;
    dist = second.dist;

    pknot = second.pknot;

    moves_done = second.moves_done;
    moves_todo = second.moves_todo;

    return *this;
  }

  //intermediate_pk ()

  intermediate_pk(int length):
    pknot(length)
  {
    structure = NULL;
    Sen = energy = 0;
    dist = INT_MAX;
  }

  intermediate_pk(const intermediate_pk &second):
    pknot(second.pknot)
  {
    structure = second.structure;
    Sen = second.Sen;
    energy = second.energy;
    dist = second.dist;

    moves_done = second.moves_done;
    moves_todo = second.moves_todo;
  }

  intermediate_pk(short *str, int energy, vector<move_fp> todo):
    pknot(str, energy)
  {
    structure = allocopy(str);

    Sen = this->energy = energy;
    for (unsigned int i=0; i<todo.size(); i++) moves_todo.insert(todo[i]);
    dist = 0;
  }

  intermediate_pk(intermediate_pk &prev, set<move_fp>::iterator it, const char *seq, short *s0, short *s1, short *new_struct = NULL):
    pknot(prev.pknot) {
    if (new_struct) structure = new_struct;
    else {
      structure = allocopy(prev.structure);
      if (it->left>0) {
        structure[it->left] = it->right;
        structure[it->right] = it->left;
      } else if (it->left<0) {
        structure[-it->left] = 0;
        structure[-it->right] = 0;
      }
    }

    int energy_chng = pknot.MakeMove(seq, s0, s1, it->left, it->right);

    energy = prev.energy + energy_chng;
    Sen = max(prev.Sen, energy);
    dist = prev.dist+1;

    moves_todo = prev.moves_todo;
    moves_todo.erase(*it);

    moves_done = prev.moves_done;
    moves_done.push_back(*it);
    moves_done[moves_done.size()-1].dE = energy_chng;
  }
};

// compare intermeds by a structure.

bool compare_intermeds::operator()(const intermediate_pk &left, const intermediate_pk &right) const {
  int i=1;
  short *lhs = left.structure;
  short *rhs = right.structure;
  char l=0,r=0;
  while (i<=lhs[0]) {
    l = lhs[i];
    r = rhs[i];
    if (l != r) {
      //fprintf (stderr, "%c %c %d %d\n", l, r, l, r);
      break;
    }
    i++;
  }
  return (i<=lhs[0] && l<r);
}

class Findpath
{
public:
  int maxkeep;

  const char *seq;
  short *s0;
  short *s1;
  int verbose_lvl;

  // hash for already seen structures -> best Saddle energy
  map<short *, int, compare_struct> structs_visited;

  // current state of computation
  priority_queue<intermediate_pk> pqueue;

  // result after computing the whole stuff
  intermediate_pk result;

  Findpath(const char *seq, int verbose_lvl = 0):
    result(strlen(seq))
  {
    this->seq = seq;
    maxkeep = 10;

    make_pair_matrix();
    s0 = encode_sequence(seq, 0);
    s1 = encode_sequence(seq, 1);
    this->verbose_lvl = verbose_lvl;
  }

  ~Findpath() {
    if (result.structure) free(result.structure);
    if (s0) free(s0);
    if (s1) free(s1);
    freeP();
  }

  int ComputeSaddle(short *str1, short *str2);
  path_pk *GetPath(short *str1, short *str2, bool chars = true);

  vector<move_fp> GetMoves(const short *str1, const short *str2);

  bool Insert(intermediate_pk &prev, set<move_fp>::iterator it, int verbose_lvl);

  void SetMaxKeep(int maxkeep) {this->maxkeep = maxkeep;}
};

vector<move_fp> Findpath::GetMoves(const short *str1, const short *str2)
{
  vector<move_fp> result;

  // checks
  if (!str1 || !str2) {
    fprintf(stderr, "ERROR: NULL pointer!\n");
    return result;
  }
  if (str1[0] != str2[0]) {
    fprintf(stderr, "ERROR: unequal lengths!\n");
    return result;
  }

  // find unequal bpairs!
  for (int i=1; i<=str1[0]; i++) {
    // str1 && str2
    if (str1[i]==str2[i]) continue;
    if (str2[i]!=0 && str2[i]>i) {
      result.push_back(move_fp(i, str2[i]));
    }
    if (str1[i]!=0 && str1[i]>i) {
      result.push_back(move_fp(-i, -str1[i]));
    }
  }

  return result;
}


bool Findpath::Insert(intermediate_pk &prev, set<move_fp>::iterator it, int verbose_lvl)
{
  bool inserted = false;
  // temporarily move the structure:
  short *structure = allocopy(prev.structure);
  if (it->left>0) {
    structure[it->left] = it->right;
    structure[it->right] = it->left;
  } else if (it->left<0) {
    structure[-it->left] = 0;
    structure[-it->right] = 0;
  }


  // get the energy
  intermediate_pk next(prev, it, seq, s0, s1, structure);
  int energy_chng = -prev.energy + next.pknot.energy;

  //if (verbose_lvl > 1) fprintf(stderr, "TIN: %s %6.2f %6.2f %d ", pt_to_str_pk(next.structure).c_str(), next.Sen/100.0, next.energy/100.0, next.dist);

  // check if we have encoutered it:
  map<short*, int, compare_struct>::iterator sit;

  /*for (sit=structs_visited.begin(); sit!=structs_visited.end(); sit++) {
    fprintf(stderr, "---: %s %7.2f\n", pt_to_str_pk(sit->first).c_str(), sit->second/100.0);
  }*/

  if ((sit = structs_visited.find(structure))!=structs_visited.end()) {
    // better:
    // update
    if (sit->second > prev.energy + energy_chng) {
      sit->second = prev.energy + energy_chng;
      inserted = true;
    }
    // worse or the same:
    // do nothing.
  } else {
    inserted = true;
    structs_visited[allocopy(structure)] = prev.energy + energy_chng;
  }

    // insert it
  if (inserted) {
    pqueue.push(next);
    if (verbose_lvl > 1) fprintf(stderr, "INS: %s %6.2f %6.2f %d\n", pt_to_str_pk(next.structure).c_str(), next.Sen/100.0, next.energy/100.0, next.dist);
    //if (verbose_lvl > 1) fprintf(stderr, "CHc: %s %6.2f\n", pt_to_str_pk(next.pknot.str).c_str(), next.pknot.energy/100.0);

    //if (verbose_lvl > 1) fprintf(stderr, "PRV: %s %6.2f %6.2f %d\n", pt_to_str_pk(prev.structure).c_str(), prev.Sen/100.0, prev.energy/100.0, prev.dist);
    //if (verbose_lvl > 1) fprintf(stderr, "CHc: %s %6.2f\n", pt_to_str_pk(prev.pknot.str).c_str(), prev.pknot.energy/100.0);

  } else {
    free(structure);
    //if (verbose_lvl > 1) fprintf(stderr, "\n");
  }

  // return if inserted
  return inserted;
}

bool MoveStr(short *structure, int left, int right)
{
  if (left<0 && structure[-left] == -right && structure[-right] == -left) {
    structure[-left] = 0;
    structure[-right] = 0;
    return true;
  }
  if (left >0 && structure[left]==0 && structure[right]==0) {
    structure[left] = right;
    structure[right] = left;
    return true;
  }
  return false;
}

int Findpath::ComputeSaddle(short *str1, short *str2)
{
  // compute the distance and all the moves.
  vector<move_fp> moves = GetMoves(str1, str2);

  int dist = (int)moves.size();

  intermediate_pk inter_tmp(str1, energy_of_struct_pk(seq, str1, s0, s1, verbose_lvl), moves);
  pqueue.push(inter_tmp);
  structs_visited[allocopy(inter_tmp.structure)] = inter_tmp.energy;

  // for each distance do
  for (int i=0; i<dist; i++) {

    if (verbose_lvl > 1) fprintf(stderr, "STR: distance %4d\n", i);


    // for each in maxkeep do:
    int cnt = 0;
    while (cnt < maxkeep) {

      int size = pqueue.size();

      //if (size == 0 && i == dist-1) break;

      if (size == 0) fprintf(stderr, "%s\n%s\n", pt_to_str_pk(str1).c_str(), pt_to_str_pk(str2).c_str());

      intermediate_pk inter = pqueue.top();
      //if (verbose_lvl > 1) fprintf(stderr, "ATT: %s %6.2f %6.2f %d\n", pt_to_str_pk(inter.structure).c_str(), inter.Sen/100.0, inter.energy/100.0, inter.dist);
      //if (verbose_lvl > 1) fprintf(stderr, "CHc: %s %6.2f\n", pt_to_str_pk(inter.pknot.str).c_str(), inter.pknot.energy/100.0);


      if (inter.dist ==i) pqueue.pop();
      else break;

      // if we are going to proceed the structure with lower energy than optimal:
      int stored_en = structs_visited[inter.structure];
      if (stored_en != inter.energy) {
        free(inter.structure);
        continue;
      }

      if (verbose_lvl > 1) fprintf(stderr, "GET: %s %6.2f %6.2f %d%c\n", pt_to_str_pk(inter.structure).c_str(), inter.Sen/100.0, inter.energy/100.0, inter.dist, inter.pknot.pknots.size()>0?'P':'-');
      //if (verbose_lvl > 1) fprintf(stderr, "CHc: %s %6.2f\n", pt_to_str_pk(inter.pknot.str).c_str(), inter.pknot.energy/100.0);


      // apply movement
      for (set<move_fp>::iterator it=inter.moves_todo.begin(); it!=inter.moves_todo.end(); it++) {



        if (verbose_lvl > 1) {
          short *st = allocopy(inter.structure);
          if (MoveStr(st, it->left, it->right)) {
            fprintf(stderr, "TRY: %s %6.2f %6.2f %d\n", pt_to_str_pk(st).c_str(), inter.Sen/100.0, inter.energy/100.0, inter.dist);
          }
          free(st);
        }
        // check if we can move it
        if (it->left<0 || inter.pknot.CanInsert(it->left, it->right)) { // can be faster

          //fprintf(stderr, "%s %d %d\n", pt_to_str_pk(inter.structure).c_str(), it->left, it->right);

          //int energy_chng = energy_of_move_PK(inter.pk, inter.structure, Enc.s0, Enc.s1, it->left, it->right, Enc.verbose_lvl);
          //if (structs_visited.find()) // faster

          // insert new one into the queue ;-)
          Insert(inter, it, verbose_lvl);
          /*intermediate next (inter, it, energy_chng);
          pqueue.push(next);*/

        }
      }
      // if we want to have it like the findpath...
      ///if (!Contains_PK(inter.structure))
      cnt++;

      free(inter.structure);

    }
    //clean:
    while (pqueue.top().dist == i) {
      free(pqueue.top().structure);
      pqueue.pop();
    }
  }

  result = pqueue.top();
  pqueue.pop();

  //clean again:
  while (!pqueue.empty()) {
    free(pqueue.top().structure);
    pqueue.pop();
  }

  for (auto it=structs_visited.begin(); it!=structs_visited.end(); it++) {
    free(it->first);
  }

  return result.Sen;
}

char *allocopy(const char *src)
{
  char *res = (char*) malloc(sizeof(char)*((int)strlen(src)+1));
  strcpy(res, src);
  return res;
}

path_pk *Findpath::GetPath(short *str1, short *str2, bool chars)
{
  ComputeSaddle(str1, str2);

  path_pk *res = (path_pk*) malloc(sizeof(path_pk)*(result.moves_done.size()+2));

  res[result.moves_done.size()+1].structure = NULL;
  res[result.moves_done.size()+1].s = NULL;
  res[result.moves_done.size()].structure = allocopy(result.structure);
  res[result.moves_done.size()].s = chars?pt_to_chars_pk(result.structure):NULL;
  res[result.moves_done.size()].en = result.energy;

  for (int i=result.moves_done.size()-1; i>=0; i--) {
    res[i].structure = allocopy(res[i+1].structure);
    int left = result.moves_done[i].left;
    int right = result.moves_done[i].right;
    if (left > 0) {
      res[i].structure[left] = 0;
      res[i].structure[right] = 0;
    } else {
      res[i].structure[-left] = -right;
      res[i].structure[-right] = -left;
    }
    res[i].s = chars?pt_to_chars_pk(res[i].structure):NULL;
    res[i].en = res[i+1].en - result.moves_done[i].dE; // minus the difference
  }

  return res;
}

int verbosity = 0;

int find_saddle_pk(const char *seq,
                    const char *struc1,
                    const char *struc2,
                    int max)
{
  Findpath fp(seq, verbosity);
  fp.SetMaxKeep(max);
  short *str1 = make_pair_table_PK(struc1);
  short *str2 = make_pair_table_PK(struc2);
  int res = fp.ComputeSaddle(str1, str2);
  free(str1);
  free(str2);
  return res;
}

path_pk* get_path_pk( const char *seq,
                  const char *s1,
                  const char* s2,
                  int maxkeep)
{
  Findpath fp(seq, verbosity);
  fp.SetMaxKeep(maxkeep);
  short *str1 = make_pair_table_PK(s1);
  short *str2 = make_pair_table_PK(s2);
  path_pk *res = fp.GetPath(str1, str2);
  free(str1);
  free(str2);
  return res;
}

path_pk* get_path_light_pk( const char *seq,
                  short *s1,
                  short* s2,
                  int maxkeep)
{
  Findpath fp(seq, verbosity);
  fp.SetMaxKeep(maxkeep);
  path_pk *res = fp.GetPath(s1, s2, false);
  return res;
}

void free_path_pk(path_pk *path)
{
  path_pk *old = path;
  while (path->structure) {
    if (path->s) free(path->s);
    free(path->structure);
    path++;
  }
  free(old);
}



