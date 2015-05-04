#ifndef __MOVE_SET_H
#define __MOVE_SET_H

/* used data structure*/
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

/* walking methods (verbose_lvl 0-2, shifts = use shift moves? noLP = no lone pairs? (not compatible with shifts))
    input:    seq - sequence
              ptable - structure encoded with make_pair_table() from pair_mat.h
              s, s1 - sequence encoded with encode_sequence from pair_mat.h
    methods:  deepest - lowest energy structure is used
              first - first found lower energy structure is used
              rand - random lower energy structure is used
    returns local minima structure in ptable and its energy in 10kcal/mol as output */

int move_deepest( char *seq,
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
int move_rand(  char *seq,
                short *ptable,
                short *s,
                short *s1,
                int verbosity_level);


/* browse_neighbours and do funct function on each of them (used mainly for user specified flooding)
    input:    seq - sequence
              ptable - structure encoded with make_pair_table() from pair_mat.h
              s, s1 - sequence encoded with encode_sequence from pair_mat.h
              funct - function (moved structure, current structure (or altered by funct)) to do with every structure in neigbourhood
    returns energy of the structure funct sets as second argument*/
int browse_neighs( char *seq,
                   short *ptable,
                   short *s,
                   short *s1,
                   int verbosity_level,
                   int shifts,
                   int noLP,
                   int (*funct) (struct_en*, struct_en*));

#endif
