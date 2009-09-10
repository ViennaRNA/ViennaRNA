#ifndef __VIENNA_RNA_PACKAGE_FOLD_H__
#define __VIENNA_RNA_PACKAGE_FOLD_H__

extern  int logML;      /* if nonzero use logarithmic ML energy in energy_of_struct */
extern  int uniq_ML;    /* do ML decomposition uniquely (for subopt) */
extern  int cut_point;  /* set to first pos of second seq for cofolding */
extern  int eos_debug;  /* verbose info from energy_of_struct */
extern  int circ;

/* function from fold.c */
float fold(const char *sequence, char *structure);
/* calculate mfe-structure of sequence */
float energy_of_struct(const char *string, const char *structure);
int   energy_of_struct_pt(const char *string, short *ptable, short *s, short *s1);
/* calculate energy of string on structure */
void  free_arrays(void);           /* free arrays for mfe folding */
void  initialize_fold(int length); /* allocate arrays for folding */
void  update_fold_params(void);    /* recalculate parameters */
char  *backtrack_fold_from_pair(char *sequence, int i, int j);
int   loop_energy(short * ptable, short *s, short *s1, int i);
void  export_fold_arrays(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);

/* some circfold related functions...	*/
float circfold(const char *string, char *structure);
float energy_of_circ_struct(const char *string, const char *structure);
void  export_circfold_arrays(int *Fc_p, int *FcH_p, int *FcI_p, int *FcM_p, int **fM2_p, int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);


/* finally moved the loop energy function declarations to this header...  */
/* BUT: The functions only exist for backward compatibility reasons!      */
/* You better include "loop_energies.h" and call the functions:           */
/* E_Hairpin() and E_IntLoop() which are threadsafe as they get a pointer */
/* to the energy parameter datastructure as additional argument           */

int   LoopEnergy(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1);
int   HairpinE(int size, int type, int si1, int sj1, const char *string);

#endif
