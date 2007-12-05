/* function from fold.c */
extern float  fold(const char *sequence, char *structure);
/* calculate mfe-structure of sequence */
extern float  energy_of_struct(const char *string, const char *structure);
/* calculate energy of string on structure */
extern void   free_arrays(void);           /* free arrays for mfe folding */
extern void   initialize_fold(int length); /* allocate arrays for folding */
extern void   update_fold_params(void);    /* recalculate parameters */
extern char  *backtrack_fold_from_pair(char *sequence, int i, int j);
extern int loop_energy(short * ptable, short *s, short *s1, int i);
extern void		export_fold_arrays(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);

/* some circfold related functions...	*/
extern	float	circfold(const char *string, char *structure);
extern	float	energy_of_circ_struct(const char *string, const char *structure);
extern	void	export_circfold_arrays(int *Fc_p, int *FcH_p, int *FcI_p, int *FcM_p, int **fM2_p, int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);
