#ifndef __MOVE_SET_PK_H
#define __MOVE_SET_PK_H

#include "pknots.h"

extern "C" {
  #include "move_set_inside.h"
}

//enum MOVE_TYPE {GRADIENT, FIRST, ADAPTIVE};

/* walking methods (verbose_lvl 0-2, shifts = use shift moves? noLP = no lone pairs? (not compatible with shifts))
    input:    seq - sequence
              ptable - structure encoded with make_pair_table() from pair_mat.h
              s, s1 - sequence encoded with encode_sequence from pair_mat.h
    methods:  deepest - lowest energy structure is used
              first - first found lower energy structure is used
              rand - random lower energy structure is used
    returns local minima structure in ptable and its energy in 10kcal/mol as output */

int move_gradient_pk(const char *seq,
                  Structure *str,
                  short *s0,
                  short *s1,
                  int shifts,
                  int verbosity_level);
int move_first_pk(const char *seq,
                  Structure *str,
                  short *s0,
                  short *s1,
                  int shifts,
                  int verbosity_level);
int move_adaptive_pk(const char *seq,
                  Structure *str,
                  short *s0,
                  short *s1,
                  int shifts,
                  int verbosity_level);

/* standardized method that encapsulates above "_pt" methods
  input:  seq - sequence
          struc - structure in dot-bracket notation
          type - type of move selection according to MOVE_TYPE enum
  return: energy of LM
          structure of LM in struc in bracket-dot notation
*/
int move_standard_pk(const char *seq,
                  char *struc,
                  enum MOVE_TYPE type,
                  int shifts,
                  int verbosity_level);

int move_standard_pk_pt(const char *seq,
                  Structure *str,
                  short *s0,
                  short *s1,
                  enum MOVE_TYPE type,
                  int shifts,
                  int verbosity_level);

/* browse_neighbours and perform "funct" function on each of them (used mainly for user specified flooding)
    input:    seq - sequence
              ptable - structure encoded with make_pair_table() from pair_mat.h
              s, s1 - sequence encoded with encode_sequence from pair_mat.h
              funct - function (structure from neighbourhood, structure from input) toperform on every structure in neigbourhood (if the function returns non-zero, the iteration through neighbourhood stops.)
    returns energy of the structure funct sets as second argument*/
int browse_neighs_pk_pt(const char *seq,
                   Structure  *str,
                   short *s0,
                   short *s1,
                   int shifts,
                   int verbosity_level,
                   int (*funct) (Structure*, Structure*));

int browse_neighs_pk(const char *seq,
                   char *struc,
                   int shifts,
                   int verbosity_level,
                   int (*funct) (Structure*, Structure*));

#endif



