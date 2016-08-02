#ifndef __MOVE_SET_H
#define __MOVE_SET_H

/**
 *  @brief  Data structure for energy_of_move()
 */
typedef struct _struct_en{
  int energy;        /* energy in 10kcal/mol*/
  short *structure;  /* structure in energy_of_move format*/
} struct_en;

/* prints structure*/
void print_stren(FILE *out, struct_en *str);
void print_str(FILE *out, short *str);

/* copying functions*/
void copy_arr(short *dest, short *src); /*just copy*/
short *allocopy(short *src);            /*copy and make space*/

enum MOVE_TYPE {GRADIENT, FIRST, ADAPTIVE};

/* walking methods (verbose_lvl 0-2, shifts = use shift moves? noLP = no lone pairs? (not compatible with shifts))
    input:    seq - sequence
              ptable - structure encoded with make_pair_table() from pair_mat.h
              s, s1 - sequence encoded with encode_sequence from pair_mat.h
    methods:  deepest - lowest energy structure is used
              first - first found lower energy structure is used
              rand - random lower energy structure is used
    returns local minima structure in ptable and its energy in 10kcal/mol as output */

int move_gradient( char *seq,
                  short *ptable,
                  short *s,
                  short *s1,
                  int verbosity_level,
                  int shifts,
                  int noLP);
int move_first( char *seq,
                short *ptable,
                short *s,
                short *s1,
                int verbosity_level,
                int shifts,
                int noLP);
int move_adaptive(  char *seq,
                short *ptable,
                short *s,
                short *s1,
                int verbosity_level);

/* standardized method that encapsulates above "_pt" methods
  input:  seq - sequence
          struc - structure in dot-bracket notation
          type - type of move selection according to MOVE_TYPE enum
  return: energy of LM
          structure of LM in struc in bracket-dot notation
*/
int move_standard(char *seq,
                  char *struc,
                  enum MOVE_TYPE type,
                  int verbosity_level,
                  int shifts,
                  int noLP);


/* browse_neighbours and perform funct function on each of them (used mainly for user specified flooding)
    input:    seq - sequence
              ptable - structure encoded with make_pair_table() from pair_mat.h
              s, s1 - sequence encoded with encode_sequence from pair_mat.h
              funct - function (structure from neighbourhood, structure from input) toperform on every structure in neigbourhood (if the function returns non-zero, the iteration through neighbourhood stops.)
    returns energy of the structure funct sets as second argument*/
int browse_neighs_pt( char *seq,
                   short *ptable,
                   short *s,
                   short *s1,
                   int verbosity_level,
                   int shifts,
                   int noLP,
                   int (*funct) (struct_en*, struct_en*));

int browse_neighs( char *seq,
                   char *struc,
                   int verbosity_level,
                   int shifts,
                   int noLP,
                   int (*funct) (struct_en*, struct_en*));

#endif



