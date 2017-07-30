#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "move_set_pk.h"

extern "C" {
  #include "pair_mat.h"
}

/* maximum degeneracy value - if degeneracy is greater than this, program segfaults*/
#define MAX_DEGEN 100
#define MINGAP 3

//#define space(a) malloc(a)

#define bool int
#define true 1
#define false 0

int compare(Structure *lhs, Structure *rhs)
{
  return (*lhs)<(*rhs);
}

int find_min(Structure *arr[MAX_DEGEN], int begin, int end) {
  Structure *min = arr[begin];
  int min_num = begin;
  int i;

  for (i=begin+1; i<end; i++) {
    if (compare(arr[i], min)) {
      min = arr[i];
      min_num = i;
    }
  }
  return min_num;
}

int equals(const Structure *first, const Structure *second)
{
  return (*first)==(*second);
}

/* ############################## DECLARATION #####################################*/
/* private functions & declarations*/

static int cnt_move = 0;
int count_move() {return cnt_move;}

void print_str_pk(FILE *out, short *str);

/*check if base is lone*/
int lone_base(short *pt, int i);

/* internal struct with moves, sequence, degeneracy and options*/
typedef struct _Encoded {
  // sequence
  const char  *seq;
  short *s0;
  short *s1;

  /* moves*/
  int   bp_left;
  int   bp_right;

  /* options*/
  //int noLP;
  int verbose_lvl;
  int first;
  int shift;
  int all_neighs; // sould be on for shifts!

  /* degeneracy*/
  int begin_unpr;
  int begin_pr;
  int end_unpr;
  int end_pr;
  Structure *processed[MAX_DEGEN];
  Structure *unprocessed[MAX_DEGEN];
  int current_en;

  /* moves in random (needs to be freed afterwards)*/
  int *moves_from;
  int *moves_to;
  int num_moves;

  /* function for flooding */
  int (*funct) (Structure*, Structure*);
} Encoded;

/* frees all things allocated by degeneracy...*/
void free_degen(Encoded *Enc)
{
  int i;
  for (i=Enc->begin_unpr; i<Enc->end_unpr; i++) {
    if (Enc->unprocessed[i]) {
      delete Enc->unprocessed[i];
      Enc->unprocessed[i]=NULL;
    }
  }
  for (i=Enc->begin_pr; i<Enc->end_pr; i++) {
    if (Enc->processed[i]) {
      delete Enc->processed[i];
      Enc->processed[i]=NULL;
    }
  }
  Enc->begin_pr=0;
  Enc->begin_unpr=0;
  Enc->end_pr=0;
  Enc->end_unpr=0;
}

/* ############################## IMPLEMENTATION #####################################*/


void print_str_pk(FILE *out, short *str)
{
  fprintf(out, "%s", pt_to_str_pk(str).c_str());
}

inline void do_move(short *pt, int bp_left, int bp_right)
{
  /* delete*/
  if (bp_left<0) {
    pt[-bp_left]=0;
    pt[-bp_right]=0;
  } else { /* insert*/
    pt[bp_left]=bp_right;
    pt[bp_right]=bp_left;
  }
}

/* done with all structures along the way to deepest*/
int update_deepest(Encoded *Enc, Structure *str, Structure *min)
{
  /* apply move + get its energy*/
  int tmp_en;
  tmp_en = str->energy + str->MakeMove(Enc->seq, Enc->s0, Enc->s1, Enc->bp_left, Enc->bp_right);

  /* use f_point if we have it */
  if (Enc->funct) {
    int end = Enc->funct(str, min);

    // undo moves
    str->UndoMove();
    Enc->bp_left=0;
    Enc->bp_right=0;

    return (end?1:0);
  }

  if (Enc->verbose_lvl>1) { fprintf(stderr, "  "); print_str_pk(stderr, str->str); fprintf(stderr, " %d\n", tmp_en); }

  /* better deepest*/
  if (tmp_en < min->energy) {

    *min = *str;

    /* delete degeneracy*/
    free_degen(Enc);

    /* undo moves*/
    str->UndoMove();
    Enc->bp_left=0;
    Enc->bp_right=0;
    return 1;
  }

  /* degeneracy*/
  if ((str->energy == min->energy) && (Enc->current_en == min->energy)) {
    int found = 0;
    int i;
    for (i=Enc->begin_pr; i<Enc->end_pr; i++) {
      if (equals(Enc->processed[i], str)) {
        found = 1;
        break;
      }
    }
    for (i=Enc->begin_unpr; !found && i<Enc->end_unpr; i++) {
      if (equals(Enc->unprocessed[i], str)) {
        found = 1;
        break;
      }
    }

    if (!found) {
      //fprintf(stderr, "%s %6.2f\n", pt_to_str_pk(str->str).c_str(), str->energy);
      Enc->unprocessed[Enc->end_unpr] = new Structure(*str);
      if (Enc->end_unpr == MAX_DEGEN) {
        fprintf(stderr, "ERROR: Degeneracy too high!!! %s\n", Enc->seq);
        for (int j=Enc->begin_unpr; j<Enc->end_unpr; j++) {
          fprintf(stderr, "D%3d %s %6.2f\n", j, pt_to_str_pk(Enc->unprocessed[j]->str).c_str(), Enc->unprocessed[j]->energy/100.0);
        }
        fprintf(stderr, "\n");
        for (int j=Enc->begin_pr; j<Enc->end_pr; j++) {
          fprintf(stderr, "D%3d %s %6.2f\n", j, pt_to_str_pk(Enc->processed[j]->str).c_str(), Enc->processed[j]->energy/100.0);
        }
				exit(EXIT_FAILURE);
      }
      Enc->end_unpr++;
    }
  }

  /* undo moves*/
  str->UndoMove();
  Enc->bp_left=0;
  Enc->bp_right=0;
  return 0;
}


/* deletions move set*/
int deletions_pk(Encoded *Enc, Structure *str, Structure *minim)
{
  int cnt = 0;
  short *pt = str->str;
  int len = pt[0];
  int i;

  for (i=1; i<=len; i++) {
    if (pt[i] && pt[i]>i) {  /* '('*/
      Enc->bp_left=-i;
      Enc->bp_right=-pt[i];


      cnt += update_deepest(Enc, str, minim);
      /* in case useFirst is on and structure is found, end*/
      if (Enc->first && cnt > 0) return cnt;
    }
  }
  return cnt;
}

  /* compatible base pair?*/
inline bool compat(char a, char b) {
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

/* try insert base pair (i,j)*/
inline bool try_insert(const short *pt, const char *seq, int i, int j)
{
  if (i<=0 || j<=0 || i>pt[0] || j>pt[0]) return false;
  return (j-i>MINGAP && pt[j]==0 && pt[i]==0 && compat(seq[i-1], seq[j-1]));
}

/* try switch base pair (i,j)*/
inline bool try_shift(const short *pt, const char *seq, int i, int j)
{
  if (i<=0 || j<=0 || i>pt[0] || j>pt[0]) return false;
  return (j-i>MINGAP && compat(seq[i-1], seq[j-1]) &&
          ((pt[i]==0 && pt[j]!=0) || (pt[i]!=0 && pt[j]==0)) );
}

// try insert base pair (i,j)
inline bool try_insert_seq(const char *seq, int i, int j)
{
  if (i<=0 || j<=0) return false;
  return (j-i>MINGAP && compat(seq[i-1], seq[j-1]));
}

/* insertions move set*/
int insertions_pk(Encoded *Enc, Structure *str, Structure *minim)
{
  int cnt = 0;
  short *pt = str->str;
  int len = pt[0];
  int i,j;

  for (i=1; i<=len; i++) {
    if (pt[i]==0) {
      for (j=i+1; j<=len; j++) {
        /* end if found closing bracket*/
        //if (pt[j]!=0 && pt[j]<j) break;  //')'
        /*if (pt[j]!=0 && pt[j]>j) {       //'('
          j = pt[j];
          continue;
        }*/
        /* if conditions are met, do insert*/
        if (try_insert(pt, Enc->seq, i, j)) {
          INS_FLAG iflag = str->CanInsert(i, j);
          if (iflag != NO_INS && (Enc->all_neighs || iflag == REG_INS || iflag == INSIDE_PK)) {
            Enc->bp_left=i;
            Enc->bp_right=j;

            //fprintf(stderr, "CI: %4d %4d\n", i, j);


            cnt += update_deepest(Enc, str, minim);
            /* in case useFirst is on and structure is found, end*/
            if (Enc->first && cnt > 0) return cnt;
          } /*else  {
            str->str[i] = j;
            str->str[j] = i;
            fprintf(stdout, "%s    BAD\n", pt_to_str_pk(str->str).c_str());
            str->str[i] = 0;
            str->str[j] = 0;
          }*/
        }
      }
    }
  }
  return cnt;
}

/* shifts move set*/
int shifts_pk(Encoded *Enc, Structure *str, Structure *minim)
{
  int cnt = 0;
  short *pt = str->str;
  int len = pt[0];
  int i,j;

  for (i=1; i<=len; i++) {
    if (pt[i]==0) continue;

    // first '('
    if (pt[i]>i) {

      for (j=i+1; j<=len; j++) {

        if (try_shift(pt, Enc->seq, i, j)) {
          INS_FLAG iflag = str->CanShift(i, j);
          if (iflag != NO_INS && (Enc->all_neighs || iflag == REG_INS || iflag == INSIDE_PK)) {
            Enc->bp_left=i;
            Enc->bp_right=j;

            //fprintf(stderr, "CS: %4d %4d\n", i, j);


            cnt += update_deepest(Enc, str, minim);
            /* in case useFirst is on and structure is found, end*/
            if (Enc->first && cnt > 0) return cnt;
          }
        }
      }
    }

    // then ')'
    if (pt[i]<i) {
      for (j=1; j<i; j++) {

        if (try_shift(pt, Enc->seq, i, j)) {
          INS_FLAG iflag = str->CanShift(i, j);
          if (iflag != NO_INS && (Enc->all_neighs || iflag == REG_INS || iflag == INSIDE_PK)) {
            Enc->bp_left=i;
            Enc->bp_right=j;

            //fprintf(stderr, "CS: %4d %4d\n", i, j);


            cnt += update_deepest(Enc, str, minim);
            /* in case useFirst is on and structure is found, end*/
            if (Enc->first && cnt > 0) return cnt;
          }
        }
      }
    }

  }
  return cnt;
}

/* move to deepest (or first) neighbour*/
int move_set(Encoded *Enc, Structure *str_in)
{
  /* count how many times called*/
  cnt_move++;

  /* count better neighbours*/
  int cnt = 0;

  /* deepest descent*/
  Structure *str = new Structure(*str_in);
  Structure *min = new Structure(*str);
  Enc->current_en = str->energy;

  if (Enc->verbose_lvl>1) { fprintf(stderr, "  start of MS:\n  "); print_str_pk(stderr, str->str); fprintf(stderr, " %d\n\n", str->energy); }

  /* if using first dont do all of them*/
  bool end = false;
  /* insertions*/
  if (!end) cnt += insertions_pk(Enc, str, min);
  if (Enc->first && cnt>0) end = true;
  if (Enc->verbose_lvl>1) fprintf(stderr, "\n");

  /* deletions*/
  if (!end) cnt += deletions_pk(Enc, str, min);
  if (Enc->first && cnt>0) end = true;

  /* shifts*/
  if (Enc->shift) {
    if (!end) cnt += shifts_pk(Enc, str, min);
    if (Enc->first && cnt>0) end = true;
  }

  /* if degeneracy occurs, solve it!*/
  if (!end && (Enc->end_unpr - Enc->begin_unpr)>0) {
    Enc->processed[Enc->end_pr] = str;
    Enc->end_pr++;
    str = Enc->unprocessed[Enc->begin_unpr];
    Enc->unprocessed[Enc->begin_unpr]=NULL;
    Enc->begin_unpr++;
    delete min;
    cnt += move_set(Enc, str);
  } else {
    /* write output to str*/
    *str = *min;
    delete min;
  }

  /* resolve degeneracy in local minima*/
  if ((Enc->end_pr - Enc->begin_pr)>0) {
    Enc->processed[Enc->end_pr]=str;
    Enc->end_pr++;

    int min = find_min(Enc->processed, Enc->begin_pr, Enc->end_pr);
    Structure *tmp = Enc->processed[min];
    Enc->processed[min] = Enc->processed[Enc->begin_pr];
    Enc->processed[Enc->begin_pr] = tmp;
    str = Enc->processed[Enc->begin_pr];
    Enc->begin_pr++;
    free_degen(Enc);
  }

  if (Enc->verbose_lvl>1 && !(Enc->first)) { fprintf(stderr, "\n  end of MS:\n  "); print_str_pk(stderr, str->str); fprintf(stderr, " %d\n\n", str->energy); }

  *str_in = *str;
  delete str;

  return cnt;
}

void construct_moves(Encoded *Enc, short *structure)
{
  /* generate all possible moves (less than n^2)*/
  Enc->num_moves = 0;
  int i;
  for (i=1; i<=structure[0]; i++) {
    if (structure[i]!=0) {
      if (structure[i]<i) continue;
      Enc->moves_from[Enc->num_moves]=-i;
      Enc->moves_to[Enc->num_moves]=-structure[i];
      Enc->num_moves++;
      //fprintf(stderr, "add  d(%d, %d)\n", i, str.structure[i]);
    } else {
      int j;
      for (j=i+1; j<=structure[0]; j++) {
        //fprintf(stderr, "check (%d, %d)\n", i, j);
        if (structure[j]==0) {
          if (try_insert_seq(Enc->seq,i,j)) {
            Enc->moves_from[Enc->num_moves]=i;
            Enc->moves_to[Enc->num_moves]=j;
            Enc->num_moves++;
            //fprintf(stderr, "add  i(%d, %d)\n", i, j);
            continue;
          }
        } /*else if (structure[j]>j) { // '('
          j = structure[j];
        } else break;*/
      }
    }
  }

  /* permute them */
  for (i=0; i<Enc->num_moves-1; i++) {
    int rnd = rand();
    rnd = rnd % (Enc->num_moves-i) + i;
    int swp;
    swp = Enc->moves_from[i];
    Enc->moves_from[i]=Enc->moves_from[rnd];
    Enc->moves_from[rnd]=swp;
    swp = Enc->moves_to[i];
    Enc->moves_to[i]=Enc->moves_to[rnd];
    Enc->moves_to[rnd]=swp;
  }
}

int move_rset(Encoded *Enc, Structure *str_in)
{
  /* count how many times called*/
  cnt_move++;

  /* count better neighbours*/
  int cnt = 0;

  /* deepest descent*/
  Structure *str = new Structure(*str_in);
  Structure *min = new Structure(*str);
  Enc->current_en = str->energy;

  if (Enc->verbose_lvl>1) { fprintf(stderr, "  start of MR:\n  "); print_str_pk(stderr, str->str); fprintf(stderr, " %d\n\n", str->energy); }

  // construct and permute possible moves
  construct_moves(Enc, str->str);

  /* find first the lower one*/
  int i;
  for (i=0; i<Enc->num_moves; i++) {
    Enc->bp_left = Enc->moves_from[i];
    Enc->bp_right = Enc->moves_to[i];
    cnt = update_deepest(Enc, str, min);
    if (cnt) break;
  }

  /* if degeneracy occurs, solve it!*/
  if (!cnt && (Enc->end_unpr - Enc->begin_unpr)>0) {
    Enc->processed[Enc->end_pr] = str;
    Enc->end_pr++;
    str = Enc->unprocessed[Enc->begin_unpr];
    Enc->unprocessed[Enc->begin_unpr]=NULL;
    Enc->begin_unpr++;
    delete min;
    cnt += move_rset(Enc, str);
  } else {
    /* write output to str*/

    *str = *min;
    delete min;
  }

  /* resolve degeneracy in local minima*/
  if ((Enc->end_pr - Enc->begin_pr)>0) {
    Enc->processed[Enc->end_pr]=str;
    Enc->end_pr++;

    int min = find_min(Enc->processed, Enc->begin_pr, Enc->end_pr);
    Structure *tmp = Enc->processed[min];
    Enc->processed[min] = Enc->processed[Enc->begin_pr];
    Enc->processed[Enc->begin_pr] = tmp;
    str = Enc->processed[Enc->begin_pr];
    Enc->begin_pr++;
    free_degen(Enc);
  }

  *str_in = *str;
  delete str;

  return cnt;
}

/*check if base is lone*/
/*int lone_base(short *pt, int i)
{
  if (i<=0 || i>pt[0]) return 0;
  // is not a base pair
  if (pt[i]==0) return 0;

  // base is lone:
  if (i-1>0) {
    // is base pair and is the same bracket
    if (pt[i-1]!=0 && ((pt[i-1]<pt[pt[i-1]]) == (pt[i]<pt[pt[i]]))) return 0;
  }

  if (i+1<=pt[0]) {
    if (pt[i+1]!=0 && ((pt[i-1]<pt[pt[i-1]]) == (pt[i]<pt[pt[i]]))) return 0;
  }

  return 1;
}*/

int move_standard_pk_pt(const char *seq,
                  Structure *str,
                  short *s0,
                  short *s1,
                  enum MOVE_TYPE type,
                  int shifts,
                  int verbosity_level)
{
  switch (type){
  case GRADIENT: move_gradient_pk(seq, str, s0, s1, shifts, verbosity_level); break;
  case FIRST: move_first_pk(seq, str, s0, s1, shifts, verbosity_level); break;
  case ADAPTIVE: move_adaptive_pk(seq, str, s0, s1, shifts, verbosity_level); break;
  }

  return str->energy;
}

int move_standard_pk(const char *seq,
                  char *struc,
                  enum MOVE_TYPE type,
                  int shifts,
                  int verbosity_level)
{
  make_pair_matrix();
  short *s0 = encode_sequence(seq, 0);
  short *s1 = encode_sequence(seq, 1);

  Structure *str = new Structure(seq, struc, s0, s1);

  int energy = move_standard_pk_pt(seq, str, s0, s1, type, shifts, verbosity_level);

  free(s0);
  free(s1);

  pt_to_chars_pk(str->str, struc);
  delete str;

  return energy;
}

int move_gradient_pk(const char *seq,
                  Structure *str,
                  short *s0,
                  short *s1,
                  int shifts,
                  int verbosity_level)
{
  cnt_move = 0;

  Encoded enc;
  enc.seq = seq;
  enc.s0 = s0;
  enc.s1 = s1;

  /* moves*/
  enc.bp_left=0;
  enc.bp_right=0;

  /* options*/
  enc.verbose_lvl=verbosity_level;
  enc.first=0;
  enc.all_neighs=shifts;
  enc.shift = shifts;

  /* degeneracy*/
  enc.begin_unpr=0;
  enc.begin_pr=0;
  enc.end_unpr=0;
  enc.end_pr=0;
  enc.current_en=0;

  // function
  enc.funct=NULL;

  int i;
  for (i=0; i<MAX_DEGEN; i++) enc.processed[i]=enc.unprocessed[i]=NULL;

  /*struct_en str;
  str.structure = allocopy(ptable);
  str.energy = energy_of_structure_pt(enc.seq, str.structure, enc.s0, enc.s1, 0);
*/

  while (move_set(&enc, str)!=0) {
    free_degen(&enc);
  }
  free_degen(&enc);

  return str->energy;
}

int move_first_pk(const char *seq,
                  Structure *str,
                  short *s0,
                  short *s1,
                  int shifts,
                  int verbosity_level)
{
  cnt_move = 0;

  Encoded enc;
  enc.seq = seq;
  enc.s0 = s0;
  enc.s1 = s1;

  /* moves*/
  enc.bp_left=0;
  enc.bp_right=0;

  /* options*/
  enc.verbose_lvl=verbosity_level;
  enc.first=1;
  enc.all_neighs=shifts;
  enc.shift = shifts;

  /* degeneracy*/
  enc.begin_unpr=0;
  enc.begin_pr=0;
  enc.end_unpr=0;
  enc.end_pr=0;
  enc.current_en=0;

  // function
  enc.funct=NULL;

  int i;
  for (i=0; i<MAX_DEGEN; i++) enc.processed[i]=enc.unprocessed[i]=NULL;

  while (move_set(&enc, str)!=0) {
    free_degen(&enc);
  }
  free_degen(&enc);

  return str->energy;
}

int move_adaptive_pk(const char *seq,
                  Structure *str,
                  short *s0,
                  short *s1,
                  int shifts,
                  int verbosity_level)
{
  srand(time(NULL));

  cnt_move = 0;

  Encoded enc;
  enc.seq = seq;
  enc.s0 = s0;
  enc.s1 = s1;

  /* moves*/
  enc.bp_left=0;
  enc.bp_right=0;

  /* options*/
  enc.verbose_lvl=verbosity_level;
  enc.first=1;
  enc.all_neighs=shifts;
  enc.shift = shifts;

  /* degeneracy*/
  enc.begin_unpr=0;
  enc.begin_pr=0;
  enc.end_unpr=0;
  enc.end_pr=0;
  enc.current_en=0;

  // function
  enc.funct=NULL;

  // allocate memory for moves
  int length = str->str[0];
  enc.moves_from = (int*) space(length*length*sizeof(int));
  enc.moves_to = (int*) space(length*length*sizeof(int));

  int i;
  for (i=0; i<MAX_DEGEN; i++) enc.processed[i]=enc.unprocessed[i]=NULL;

  while (move_rset(&enc, str)!=0) {
    free_degen(&enc);
  }
  free_degen(&enc);

  free(enc.moves_from);
  free(enc.moves_to);

  return str->energy;
}

int browse_neighs_pk( const char *seq,
                   char *struc,
                   int shifts,
                   int verbosity_level,
                   int (*funct) (Structure*, Structure*))

{
  make_pair_matrix();
  short *s0 = encode_sequence(seq, 0);
  short *s1 = encode_sequence(seq, 1);

  Structure *str = new Structure(seq, struc, s0, s1);

  int res = browse_neighs_pk_pt(seq, str, s0, s1, shifts, verbosity_level, funct);

  free(s0);
  free(s1);

  delete str;

  return res;
}

int browse_neighs_pk_pt(const char *seq,
                  Structure *str,
                  short *s0,
                  short *s1,
                  int shifts,
                  int verbosity_level,
                  int (*funct) (Structure*, Structure*))
{
  cnt_move = 0;

  Encoded enc;
  enc.seq = seq;
  enc.s0 = s0;
  enc.s1 = s1;

  /* moves*/
  enc.bp_left=0;
  enc.bp_right=0;

  /* options*/
  enc.verbose_lvl=verbosity_level;
  enc.first=1;
  enc.all_neighs=1;
  enc.shift = shifts;

  /* degeneracy*/
  enc.begin_unpr=0;
  enc.begin_pr=0;
  enc.end_unpr=0;
  enc.end_pr=0;
  enc.current_en=0;

  // function
  enc.funct=funct;

  int i;
  for (i=0; i<MAX_DEGEN; i++) enc.processed[i]=enc.unprocessed[i]=NULL;

  move_set(&enc, str);
  free_degen(&enc);

  return str->energy;
}

/*
// sample usage:
int main() {
  char seq[20] = "ACCCCCCTCTGTAGGGGGA";
  char str[20] = ".((.(.........).)).";

  // move to the local minimum and display it
  int energy = move_standard(seq, str, GRADIENT, 0, 0, 0);
  fprintf(stdout, "%s %6.2f\n\n", str, energy/100.0);

  //now create an array of every structure in neighbourhood of str structure
  struct_en *list = NULL;
  int list_length = 0;

  int get_list(struct_en *new_one, struct_en *old_one)
  {
    // enlarge the list
    list_length++;
    list = (struct_en*) realloc(list, list_length*sizeof(struct_en));

    // copy the structure
    list[list_length-1].energy = new_one->energy;
    list[list_length-1].structure = allocopy(new_one->structure);

    // we want to continue -> return 0
    return 0;
  }
  browse_neighs(seq, str, 0, 0, 0, get_list);

  // print them and free the memory:
  int i;
  for (i=0; i<list_length; i++) {
    print_stren(stdout, &list[i]);
    free(list[i].structure);
  }
  free(list);

  return 0;
}*/


