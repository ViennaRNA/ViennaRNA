#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#include "neighbourhood.h"
#include "RNAlocmin.h"
#include "hash_util.h"

extern "C" {
  #include "move_set_inside.h"
  #include "fold.h"
  #include "pair_mat.h"
}

#define MINGAP 3

// declare static members:
char *Neighborhood::seq = NULL;
short *Neighborhood::s0 = NULL;
short *Neighborhood::s1 = NULL;
int Neighborhood::debug = 0;

// for degeneracy:
int Neighborhood::energy_deg = 0;
bool Neighborhood::deal_degen = 1;
std::vector<Neighborhood*> Neighborhood::degen_todo;
std::vector<Neighborhood*> Neighborhood::degen_done;

void error_message(char *str, int i = -1, int j = -1, int k = -1, int l = -1)
{
  fprintf(stderr, str, i, j, k, l);
  exit(EXIT_FAILURE);
}

Neigh::Neigh(int i, int j, int energy)
{
  this->i = i;
  this->j = j;

  energy_change = energy;
}

Neigh::Neigh():Neigh(0, 0)
{
}

bool const Neigh::operator<(const Neigh &second) const
{
  // typ = insert, nothing, delete
  int typ_us = i>0?0:(i==0?1:2);
  int typ_them = second.i>0?0:(second.i==0?1:2);

  // if types of neighbor are same, then just compare the numbers in them
  if (typ_us == typ_them) return i<second.i;
  else return typ_us < typ_them;
}

Loop::Loop(int i, int j)
{
  this->left = i;
  this->right = j;

  energy = INT_MAX;
}

/*Loop::Loop(Loop &second)
{
  *this = second;
}*/

inline bool
compat(char a, char b){

  if (a=='A' && b=='U') return true;
  if (a=='C' && b=='G') return true;
  if (a=='G' && b=='U') return true;
  if (a=='U' && b=='A') return true;
  if (a=='G' && b=='C') return true;
  if (a=='U' && b=='G') return true;
  /* and with T's*/
  if (a=='A' && b=='T') return true;
  if (a=='T' && b=='A') return true;
  if (a=='G' && b=='T') return true;
  if (a=='T' && b=='G') return true;
  return false;
}

int Loop::GenNeighs(char *seq, short *pt)
{
  neighs.clear();
  int res = -1;

  for (int i=left+1; i<right; i++) {
    if (pt[i]>i) {   // '('
      if (res == -1) res = i;
      i = pt[i];
      continue;
    }
    for (int j=i+1; j<right; j++) {
      if (pt[j]>j) {
        j = pt[j];
        continue;
      }
      if (j-i>MINGAP && pt[j] == 0 && compat(seq[i-1], seq[j-1])) {
        neighs.push_back(Neigh(i,j));
      }
    }
  }

  return res;
}

int Loop::EvalLoop(short *pt, short *s0, short *s1, bool inside)
{
  energy = loop_energy(pt, s0, s1, left);

  if (inside) {
    for (int i=0; i<(int)neighs.size(); i++) {
      pt[neighs[i].i] = neighs[i].j;
      pt[neighs[i].j] = neighs[i].i;
      neighs[i].energy_change = -energy + loop_energy(pt, s0, s1, neighs[i].i) + loop_energy(pt, s0, s1, left);
      pt[neighs[i].i] = 0;
      pt[neighs[i].j] = 0;
    }
  }

  if (Neighborhood::debug) fprintf(stderr, "EvalLoop %s (%3d, %3d) = %4d\n", pt_to_str(pt).c_str(), left, right, energy);

  return energy;
}

void debug_loops(std::vector<Loop*> &loops)
{
  for (int i=0; i<(int)loops.size(); i++) {
    if (loops[i] && loops[i]->left != i) fprintf(stderr, "WARNING: loops[i]->left=%d i=%d", loops[i]->left, i);
  }
}

Neighborhood::Neighborhood(char *seq_in, short *s0, short *s1, short *pt, bool eval)
{
  this->pt = allocopy(pt);
  this->seq = seq_in;
  /*
  seq = (char*)malloc((strlen(seq_in)+1)*sizeof(char));
  strcpy(seq, seq_in);*/

  this->s0 = s0;
  this->s1 = s1;

  energy = INT_MAX;

   // create array of loops:
  loops.resize(pt[0]+1);
  for (int i=0; i<(int)loops.size(); i++) loops[i] = NULL;

  // generate the external loop
  Loop *newone = new Loop(0, pt[0]+1);
  loops[0] = newone;
  int i = newone->GenNeighs(seq, pt);

  // generate the neighbourhood (inserts)
  if (i != -1) {
    for (; i<pt[0]; i++) {
      if (pt[i] != 0 && pt[i]>i) {
        Loop *newone = new Loop(i, pt[i]);
        loops[i] = newone;
        int k = newone->GenNeighs(seq, pt);

        // jump to next -- either inside or outside
        if (k!=-1) i = k-1;
        else i = pt[i];
      }
    }
  }

  if (eval) EvalNeighs(true);
  StartEnumerating();
  //debug_loops(loops);
}

Neighborhood::Neighborhood(const Neighborhood &second)
{
  pt = NULL;
  HardCopy(second);
}

Neighborhood::~Neighborhood()
{
  Free();
}

void Neighborhood::Free()
{
  if (debug && pt) fprintf(stderr, "Free     %s %6.2f\n", pt_to_str(pt).c_str(), energy/100.0);
  top_loop.clear();

  if (pt) free(pt);
  for (int i=0; i<(int)loops.size(); i++) {
    if (loops[i]) {
      delete loops[i];
      loops[i] = NULL;
    }
  }
}

void Neighborhood::HardCopy(const Neighborhood &second)
{
  Free();
  this->pt = allocopy(second.pt);
  this->energy = second.energy;
  this->loopnum = second.loopnum;
  this->neighnum = second.neighnum;
  this->top_loop = second.top_loop;

  if (debug) fprintf(stderr, "HardCopy %s %6.2f\n", pt_to_str(pt).c_str(), energy/100.0);

  loops.resize(second.loops.size(), NULL);
  for (int i=0; i<(int)second.loops.size(); i++) {
    if (second.loops[i]) loops[i] = new Loop(*second.loops[i]);
  }
  debug_loops(loops);
}

void Neighborhood::ClearStatic()
{
  if (seq) { seq = NULL; }
  if (s0) { s0 = NULL; }
  if (s1) { s1 = NULL; }
  ClearDegen();
}

bool const Neighborhood::operator<(const Neighborhood &second) const
{
  if (second.energy != energy) return energy<second.energy;
  else return compf_short(pt, second.pt);
}

inline int find_enclosing(short *pt, int i)
{
  int beg = i-1;
  for (; beg>0; beg--) {
    if (pt[beg]==0) continue;
    if (pt[beg]>beg) break;
    if (pt[beg]<beg) beg = pt[beg];
  }
  return beg;
}

int Neighborhood::AddBase(int i, int j, bool reeval)
{
  // find enclosing loop (can be better)
  int beg = find_enclosing(pt, i);

  int size = -loops[beg]->neighs.size();

  // insert it + generate new neighbors
  if (loops[i]) error_message("Loop %3d already set!!!", i);
  Loop* newloop = new Loop(i,j);
  loops[i] = newloop;
  newloop->GenNeighs(seq, pt);
  pt[i] = j;
  pt[j] = i;
  int energy_chng = 0;
  if (reeval) energy_chng += newloop->EvalLoop(pt, s0, s1, true);

  // delete the neighbors that are wrong now (can be better)
  if (reeval) energy_chng -= loops[beg]->energy;
  loops[beg]->GenNeighs(seq, pt);
  if (reeval) energy_chng += loops[beg]->EvalLoop(pt, s0, s1, true);

  // update energy
  energy += energy_chng;

  // calculate size
  size += loops[i]->neighs.size() + loops[beg]->neighs.size();

  return size;
}

int Neighborhood::RemBase(int i, int j, bool reeval)
{
  // loop exists?
  if (loops[i] == NULL) error_message("There is no loop at point %d!!!\n", i);
  if (loops[i]->left!=i || loops[i]->right!=j) error_message("Different end: removing (%d, %d); exists (%d, %d)\n", i, j, loops[i]->left, loops[j]->right);

  // find the upper loop:
  int upper = find_enclosing(pt, i);
  int size = -loops[upper]->neighs.size() - loops[i]->neighs.size();

  // delete this one
  int energy_chng = 0;
  if (reeval) energy_chng -= loops[i]->energy;
  delete loops[i];
  loops[i] = NULL;
  pt[i] = 0;
  pt[j] = 0;

  // recompute the upper one:
  if (reeval) energy_chng -= loops[upper]->energy;
  loops[upper]->GenNeighs(seq, pt);
  if (reeval) energy_chng += loops[upper]->EvalLoop(pt, s0, s1, true);
  size += loops[upper]->neighs.size();

  // energy assign
  energy += energy_chng;

  return size;
}

int Neighborhood::ApplyNeigh(Neigh &neigh, bool reeval)
{
  if (neigh.i>0) return AddBase(neigh.i, neigh.j, reeval);
  else return RemBase(-neigh.i, -neigh.j, reeval);
}

int Neighborhood::PrintNeighs()
{
  int res = 0;
  for (int i=0; i<(int)loops.size(); i++) {
    if (loops[i]) {
      if (i!=0) res++; // one delete move per loop
      fprintf(stdout, "Loop %3d %3d - %5d (%d neighbors):\n", loops[i]->left, loops[i]->right, loops[i]->energy, (int)loops[i]->neighs.size());
      for (int j=0; j<(int)loops[i]->neighs.size(); j++) {
        fprintf(stdout, "  %3d %3d %5d\n", loops[i]->neighs[j].i, loops[i]->neighs[j].j, loops[i]->neighs[j].energy_change);
        res++;
      }
    }
  }
  return res;
}

int Neighborhood::PrintEnum(bool inserts_first)
{
  int res = 0;
  Neigh tmp(0,0,0);
  StartEnumerating(inserts_first);
  while (NextNeighbor(tmp, inserts_first)) {
    fprintf(stderr, "  %3d %3d %5d\n", tmp.i, tmp.j, tmp.energy_change);
    res++;
  }
  return res;
}

void Neighborhood::PrintStr()
{
  fprintf(stdout, "%s %6.2f\n", pt_to_str(pt).c_str(), energy/100.0);
}

int Neighborhood::EvalNeighs(bool full)
{
  energy = 0;
  for (int i=0; i<(int)loops.size(); i++) {
    if (loops[i]) energy += loops[i]->EvalLoop(pt, s0, s1, full);
  }

  return energy;
}

int Neighborhood::RemEnergy(short *pt, int loop, int last_loop)
{
  // find last loop if not provided
  if (last_loop == -1) {
    for (last_loop = loop-1; last_loop>0; last_loop--) {
      if (pt[last_loop] == 0) continue;
      if (pt[last_loop] > last_loop) break;
    }
  }

  // resolve energy:
  pt[loops[loop]->left] = 0;
  pt[loops[loop]->right] = 0;
  int change = -loops[loop]->energy - loops[last_loop]->energy + loop_energy(pt, s0, s1, loops[last_loop]->left);
  pt[loops[loop]->left] = loops[loop]->right;
  pt[loops[loop]->right] = loops[loop]->left;

  return change;
}

std::string Neighborhood::GetPT(Neigh &next)
{
  std::string ret;
  // delete:
  if (next.i<0) {
    pt[-next.i] = 0;
    pt[-next.j] = 0;
    ret = pt_to_str(pt);
    pt[-next.i] = -next.j;
    pt[-next.j] = -next.i;
  } else { // insert
    pt[next.i] = next.j;
    pt[next.j] = next.i;
    ret = pt_to_str(pt);
    pt[next.i] = 0;
    pt[next.j] = 0;
  }
  return ret;
}

int Neighborhood::MoveLowest(bool first, bool reeval)
{
  int lowest = 0;

  // debug:
  if (debug) fprintf(stderr, "MoveLows %s %6.2f\n", pt_to_str(pt).c_str(), energy/100.0);
  if (debug>1) PrintEnum();

  StartEnumerating();
  Neigh next;
  bool lowest_found = false;
  Neigh lowest_n;
  while (NextNeighbor(next)) { // linear -- with a better data class, we would be able to pick the lowest energy first and have it in O(1) (or O(log n)) -- the majority of time comes from energy evaluation, so maybe not a reasonto do that...
    // degeneracy!
    if (lowest == 0 && next.energy_change == 0 && deal_degen) {
      if (debug) fprintf(stderr, "FndEqual %s %6.2f (%3d, %3d)\n", GetPT(next).c_str(), (next.energy_change+energy)/100.0, next.i, next.j);
      AddDegen(next);
    }

    // two options: either we have ound the lower one, or we have found the same energetically, but lower lexikografically
    if (next.energy_change < lowest ||
        (next.energy_change == lowest && (lowest > 0 || !deal_degen) && next < lowest_n)) {
      if (debug) fprintf(stderr, "FndLower %s %6.2f (%3d, %3d)\n", GetPT(next).c_str(), (next.energy_change+energy)/100.0, next.i, next.j);
      ClearDegen();
      lowest = next.energy_change;
      lowest_n = next;
      lowest_found = true;
      if (first) break;
    }
  }

  // solve degen:
  if (deal_degen && (degen_done.size() + degen_todo.size() > 0)) return SolveDegen(false, reeval, lowest, first);


  // apply it: (in case of no degeneracy)
  if (lowest_found) {
    ApplyNeigh(lowest_n);
  }

  if (lowest_found && lowest == 0) return 1;
  else return lowest;
}

int Neighborhood::SolveDegen(bool random, bool reeval, int lowest, bool first)
{
  // resolve degeneracy
  if (degen_todo.size() > 0) {
    if (energy == energy_deg) {
      degen_done.push_back(new Neighborhood(*this));
      if (debug) fprintf(stderr, "AddDoneD %s %6.2f\n", pt_to_str(pt).c_str(), (energy)/100.0);
    }
    Neighborhood *todo = degen_todo[0];
    degen_todo.erase(degen_todo.begin());
    int degen_en = random?todo->MoveRandom(reeval):todo->MoveLowest(first, reeval);
    HardCopy(*todo);
    delete todo;  // maybe can be better
    //ClearDegen();
    return degen_en+lowest;
  }

  // now chose the lowest one:
  if (degen_done.size() > 0) {
    // chose the lowest one lexicographically:
    Neighborhood *res = this;
    if (debug) fprintf(stderr, "LwstLexT %s %6.2f\n", pt_to_str(pt).c_str(), (energy)/100.0);
    for (int i=0; i<(int)degen_done.size(); i++) {
      if (debug) fprintf(stderr, "LwstLex  %s %6.2f\n", pt_to_str(degen_done[i]->pt).c_str(), (degen_done[i]->energy)/100.0);
      if (*degen_done[i] < *res) res = degen_done[i];
    }

    if (debug) fprintf(stderr, "LwstLexW %s %6.2f\n", pt_to_str(res->pt).c_str(), (res->energy)/100.0);

    int diff_en = energy - res->energy;
    if (this != res) HardCopy(*res);
    ClearDegen();
    return diff_en;
  }

  fprintf(stderr, "###### WARNING! we should not be here at all! (SolveDegen)");
  return 0; // will never happen hopefully
}

int Neighborhood::MoveRandom(bool reeval)
{
  srand(time(NULL));

  // debug:
  if (debug) fprintf(stderr, "MoveRND  %s %6.2f\n", pt_to_str(pt).c_str(), energy/100.0);
  if (debug>1) PrintEnum();

  int lowers = 0;
  int equals = 0;
  StartEnumerating();
  Neigh next;
  while (NextNeighbor(next)) {
    if (next.energy_change < 0) lowers++;
    if (next.energy_change == 0) equals++;
  }

  // if found any lowers, then draw random and go there
  if (lowers>0) {
    ClearDegen();
    int rnd = rand() % lowers;
    StartEnumerating();
    while (NextNeighbor(next)) {
      if (next.energy_change < 0) {
        if (rnd > 0) rnd--;
        else break;
      }
    }

    ApplyNeigh(next, reeval);
    return next.energy_change;
  }

  // else just add all equals to degen and do degen stuff:
  if (equals == 0) return 0;
  else {
    if (deal_degen) {
      StartEnumerating();
      while (NextNeighbor(next)) {
        if (next.energy_change == 0) AddDegen(next);
      }
    }
  }

  // solve degen:
  if (deal_degen && (degen_done.size() + degen_todo.size() > 0)) return SolveDegen(true, reeval);

  return 0;
}

void Neighborhood::StartEnumerating(bool inserts_first)
{
  top_loop.clear();
  loopnum = 0;
  neighnum = -1; // if neighnum == -1 we do deletes
  deletes = false;
  IncreaseCount(inserts_first);
}

inline void Neighborhood::MoveLoop()
{
  neighnum = -1;
  top_loop.push_back(loopnum);
  loopnum++;
  while (loopnum<(int)loops.size() && loops[loopnum]==NULL) loopnum++;
  while (top_loop.size()>0 && loops[top_loop.back()]->right < loopnum) top_loop.pop_back();
}

void Neighborhood::IncreaseCount(bool inserts_first)
{
  if (!deletes) {
    // increase the count
    neighnum++;
    if (neighnum >= (int)loops[loopnum]->neighs.size()) {
      MoveLoop();

    }
  } else {
    MoveLoop();
  }

  // now if the inserts are first:
  if (inserts_first) {
    if ((int)loops.size() <= loopnum && !deletes) {
      // switch to deletes:
      deletes = true;
      loopnum = 0;
      MoveLoop();
    }
    if (!deletes && neighnum==-1) IncreaseCount(inserts_first);
  }
}

bool Neighborhood::NextNeighbor(Neigh &res, bool inserts_first, bool with_energy)
{
  // first get end:
  if ((int)loops.size() <= loopnum) return false;

  // deletes and inserts:
  if (neighnum == -1) {
    // deletes
    res = Neigh(-loops[loopnum]->left, -loops[loopnum]->right, with_energy?RemEnergy(pt, loopnum, top_loop.back()):INT_MAX);
  } else {
    // inserts
    res = loops[loopnum]->neighs[neighnum];
  }

  // increae the enumeration count
  IncreaseCount(inserts_first);

  //if (loop.size() <= loopnum) return false;
  return true;
}

void Neighborhood::ClearDegen()
{
  // debug:
  if (debug && (degen_done.size() + degen_todo.size() > 0)) fprintf(stderr, "ClrDegen (%d, %d)\n", (int)degen_todo.size(), (int)degen_done.size());

  for (int i=0; i<(int)degen_done.size(); i++) {
    delete degen_done[i];
  }

  for (int i=0; i<(int)degen_todo.size(); i++) {
    delete degen_todo[i];
  }

  degen_done.clear();
  degen_todo.clear();
}

void Neighborhood::SwitchOffDegen()
{
  deal_degen = false;
}

bool Neighborhood::AddDegen(Neigh &neigh)
{
  int res = false;

  // check if already there:
  ApplyNeigh(neigh);

  // debug:
  if (debug) fprintf(stderr, "AddDegen %s %6.2f (%3d, %3d)\n", pt_to_str(pt).c_str(), energy/100.0, neigh.i, neigh.j);

  // assign energy if first:
  if (degen_done.size() == 0 && degen_todo.size() == 0) energy_deg = energy;

  // check:
  if (energy_deg != energy) {
    fprintf(stderr, "WARNING: energies do not match in AddDegen (%d != %d)\n", energy_deg, energy);
  }

  // search him in degen_*
  for (int i=0; i<(int)degen_todo.size(); i++) {
    if (*degen_todo[i] == *this) {
      res = true;
      break;
    }
  }
  if (!res) {
    for (int i=0; i<(int)degen_done.size(); i++) {
      if (*degen_done[i] == *this) {
        res = true;
        break;
      }
    }
  }

  // add if not found
  if (!res)  {
    degen_todo.push_back(new Neighborhood(*this));
    // debug:
    if (debug) fprintf(stderr, "AddTodoD %s %6.2f (%3d, %3d)\n", pt_to_str(pt).c_str(), energy/100.0, neigh.i, neigh.j);
  }

  // return state
  Neigh to_apply(-neigh.i, -neigh.j, -neigh.energy_change);
  ApplyNeigh(to_apply);

  return !res;
}

extern "C" {
  #include "utils.h"
}

void test()
{
//char num[] = "123456789012345678901234567890123456";
  char seq[] = "UAGAGGUUCGGUGUUGAACGUGUAAAGAAGUAAUGAUUCGUCUUGUUACCAGCAUGGUAUUCCCGCUUCCUGAUGUAUCUGGGAAUAAAGCAAAUAGUAC";
  char str0[] = "((.(.(((((....))))).).))..................((((...........(((((((((........))....)))))))..)))).......";
  char str1[] = "((.(.(((((....))))).).))..................(((((..........(((((((((........))....))))))).))))).......";
  short *pt0 = make_pair_table(str0);
  short *pt1 = make_pair_table(str1);

  make_pair_matrix(); // dunno if needed
  update_fold_params();

  short *s0 = encode_sequence(seq, 0);
  short *s1 = encode_sequence(seq, 1);

  Neighborhood nh0(seq, s0, s1, pt0);
  Neighborhood nh1(seq, s0, s1, pt1);

  fprintf(stderr, "%d\n", nh0 < nh1);
  //nh.EvalNeighs(true);
  nh0.PrintEnum(false);
  nh0.PrintEnum();

  //int size = nh.AddBase(15,25,true);
  //fprintf(stderr, "adding %d %d added %d neighbours\n", 15,25,size);
  //nh.PrintNeighs();
  //nh.PrintEnum();
  //nh.PrintStr();
  //while(nh.MoveLowest());
  //nh.PrintStr();

  free(pt0);
  free(pt1);
  free(s0);
  free(s1);
  free_arrays();
  Neighborhood::ClearStatic();
}
