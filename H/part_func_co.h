/* functions from part_func.c */
/* calculate partition function and base pair probabilities */
#define LARGE_PF
#ifdef  LARGE_PF
#define FLT_OR_DBL double
#else
#define FLT_OR_DBL float
#endif
extern int mirnatog; /*toggles no intrabp in 2nd mol*/

extern float  co_pf_fold(char *sequence, char *structure); /* calculate partition function and base pair probabilities */
extern void   init_co_pf_fold(int length);
extern void   free_co_pf_arrays(void);
extern void   update_co_pf_params(int length); /*recalculate energy parameters */
extern char   co_bppm_symbol(float *x);    /* string representation of structure */
extern void   compute_probabilities(double *FEAB,double *FEAA,
				     double *FEBB,double *FEA,
				     double *FEB, struct plist  *prAB,
				     struct plist  *prAA, struct plist  *prBB,
				     struct plist  *prA, struct plist  *prB,
				     int Alength,int Blength);

extern float *get_monomerefreeenergies();

typedef struct ConcEnt {
  double A0;    /*start concentration A*/
  double B0;    /*start concentration B*/
  double ABc;   /*End concentration AB*/
  double AAc;
  double BBc;
  double Ac;
  double Bc;
} ConcEnt;



typedef struct pairpro{
  struct plist *AB;
  struct plist *AA;
  struct plist *A;
  struct plist *B;
  struct plist *BB;
}pairpro;


extern struct ConcEnt  *get_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double * startconc);

extern struct plist *get_plist(struct plist *pl, int length, double cut_off);

extern struct plist *get_mfe_plist(struct plist *pl);
extern int make_probsum(int length, char *name);
