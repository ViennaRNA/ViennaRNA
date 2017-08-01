#ifndef __PKNOTS_H
#define __PKNOTS_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <vector>
#include <set>
#include <map>
#include <string>

enum BPAIR_TYPE {N_S, N_M, P_H, P_K, P_L, P_M, ROOT};
const char bpair_type_name[][5] = {"S", "M", "P_H", "P_K", "P_L", "P_M", "ROOT"};
const char bpair_type_sname[] = "smHKLM_";
const int beta1_pen[] =   {0,    0,  960, 1260, 1460, 1760,    0};  // simple pk penalty
const int beta1mp_pen[] = {0,    0, 1500, 1800, 2000, 2300,    0};  // pk in a multiloop
const int beta2_pen[] =   {0,    0,   10,   10,   10,   10,    0};  // penalty for each bpair forming a pknot
const int beta3_pen[] =   {0,    0,   10,   10,   10,   10,    0};  // penalty for a nucleotide not paired

short *make_pair_table_PK(const char *str);
std::string pt_to_str_pk(const short *str);
char* pt_to_chars_pk(const short *str);
void pt_to_chars_pk(const short *str, char *dest);
//short *allocopy(const short *src);
//char *allocopy(const char *src);
//void copy_arr(short *dest, const short *src);

int Contains_PK(short *str);

//debug
float get_eos_time();

class Pseudoknot
{
public:
  bool imat[4][4];
  std::set<int> parts[4];
  int size;
public:
  // contains some bpair?
  bool Contains(int left);

  // crosses this bpair?
  //bool Cross(short *str, int left, int right);

  // can we insert it inside?
  int Inside(short *str, int left, int right);

  // can we insert it normally?
  bool CanInsert(int left, std::vector<int> &numbers, bool insert = false);

  // delete a bpair:
  bool Delete(int left);

  inline bool Imat(int a, int b) {
    return a>b?imat[b][a]:imat[a][b];
  }

  // constructor
  Pseudoknot();
  Pseudoknot(const Pseudoknot &pknot);

  Pseudoknot &operator=(const Pseudoknot &pknot);
};

enum INS_FLAG {NO_INS, REG_INS, INSIDE_PK, CREATE_PK, CHNG_PK};

class Structure {

public:
  // pknots:
  std::vector<Pseudoknot> pknots;

  // map bpairs to pknots or nothing if the bpair is standard
  std::map<int, int> bpair_pknot;

public:
  // old-school structure
  short *str;
  //char *seq;

  // energy:
  int energy;

private:
  int undo_l;
  int undo_r;
  int undo_en;

public:
  // constructor
  Structure(const char *seq, short *structure, short *s0, short *s1);
  Structure(const char *seq, char *structure, short *s0, short *s1);
  Structure(short *structure, int energy);
  Structure(char *structure, int energy);
  Structure(const Structure &second);
  Structure(int length);
  ~Structure();

  inline Pseudoknot *Pknot_index(int index) {
    if (index == -1) return NULL;
    return &pknots[index];
  }

  inline Pseudoknot *Pknot_bpair(int index) {
    return Pknot_index(bpair_pknot[index]);
  }

  bool const operator<(const Structure &second) const;
  bool const operator==(const Structure &second) const;
  Structure &operator=(const Structure &second);
private:
  // is the bpair viable?
  INS_FLAG ViableInsert(int left, int right, bool insert = false);
  INS_FLAG Insert(int left, int right);
  INS_FLAG Shift(int left, int right);
  bool Delete(int left);
public:
  // do a move
  int MakeMove(const char *seq, short *s0, short *s1, int left, int right);
  int UndoMove();

  INS_FLAG CanInsert(int left, int right);
  INS_FLAG CanShift(int left, int right);
};

class Helpers {
public:
  std::vector<int> str_energy;
  std::vector<BPAIR_TYPE> str_type;
  std::vector<int> str_torght;
  std::vector<int> str_toleft;
  int last_open;

  int beta2;
  int beta3;
  int beta1;

  void Create(int length);
  Helpers(int length);
};

int energy_of_struct_pk(const char *seq, char *structure, int verbose = 0);
int energy_of_struct_pk(const char *seq, short *structure, int verbose = 0);
int energy_of_struct_pk(const char *seq, short *structure, short *s0, short *s1, int verbose = 0);
void freeP(); // cleanup - run after last call of energy_of_struct - not necessary


void try_pk();

#endif
