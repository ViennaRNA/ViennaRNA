/* to use floats instead of doubles in pf_fold() comment next line */
#define LARGE_PF 
#ifdef  LARGE_PF
#define FLT_OR_DBL double
#else
#define FLT_OR_DBL float
#endif

extern int  noGU;           /* GU not allowed at all */
extern int  no_closingGU;   /* GU allowed only inside stacks */
extern int  tetra_loop;     /* Fold with specially stable 4-loops */
extern int  energy_set;     /* 0 = BP; 1=any mit GC; 2=any mit AU-parameter */
extern int  dangles;	    /* use dangling end energies (not in part_func!) */
/*@null@*/
extern char *nonstandards;  /* contains allowed non standard bases */
extern double temperature;   /* rescale parameters to this temperature */
extern int  james_rule;     /* interior loops of size 2 get energy 0.8Kcal and
			       no mismatches, default 1 */
extern int  logML;          /* use logarithmic multiloop energy function */
extern int  cut_point;      /* first position of 2nd strand for co-folding */

struct bond {               /* base pair */
   int i;
   int j;
};
typedef struct bond bondT;            
extern bondT  *base_pair; /* list of base pairs */

extern FLT_OR_DBL *pr;          /* base pairing prob. matrix */
extern int   *iindx;            /* pr[i,j] -> pr[iindx[i]-j] */
extern double pf_scale;         /* scaling factor to avoid float overflows*/
extern int    fold_constrained; /* fold with constraints */
extern int    do_backtrack;     /* calculate pair prob matrix in part_func() */
extern int    noLonelyPairs;    /* avoid helices of length 1 */
extern char backtrack_type;     /* usually 'F'; 'C' require (1,N) to be bonded;
				   'M' seq is part of a multi loop */
char * option_string(void);
