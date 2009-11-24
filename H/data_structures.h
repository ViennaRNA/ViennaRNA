#ifndef __VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H__
#define __VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H__

#include "energy_const.h"
#include "list.h"

/* to use floats instead of doubles in pf_fold() comment next line */
#define LARGE_PF
#ifdef  LARGE_PF
#define FLT_OR_DBL double
#else
#define FLT_OR_DBL float
#endif

#ifndef NBASES
#define NBASES 8
#endif

#ifndef MAXALPHA
#define MAXALPHA 20         /* maximal length of alphabet */
#endif

#define MAXDOS            1000 /* maximum density of states discretization for subopt */

/*
* ############################################################
* Here are the type definitions of various datastructures
* shared among the Vienna RNA Package
* ############################################################ 
*/

/* this datastructure is used as input parameter in functions of PS_dot.c */
typedef struct plist {
  unsigned int i;
  unsigned int j;
  float p;
} plist;

/* this datastructure is used as input parameter in functions of PS_dot.c */
typedef struct cpair {
  int i,j,mfe;
  float p, hue, sat;
} cpair;

/*
* this is a workarround for the SWIG Perl Wrapper RNA plot function
* that returns an array of type COORDINATE
*/
typedef struct {
  float X; /* X coords */
  float Y; /* Y coords */
} COORDINATE;

typedef struct sect {
  unsigned int  i;
  unsigned int  j;
  int ml;
} sect; /* stack of partial structures for backtracking */

typedef struct bond {               /* base pair */
   unsigned int i;
   unsigned int j;
} bondT;

typedef struct bondTEn {               /* base pair with associated energy*/
   int i;
   int j;
   int energy;
} bondTEn;

typedef struct paramT{
  int id;
  int stack[NBPAIRS+1][NBPAIRS+1];
  int hairpin[31];
  int bulge[MAXLOOP+1];
  int internal_loop[MAXLOOP+1];
  int mismatchI[NBPAIRS+1][5][5];
  int mismatchH[NBPAIRS+1][5][5];
  int mismatchM[NBPAIRS+1][5][5];
  int dangle5[NBPAIRS+1][5];
  int dangle3[NBPAIRS+1][5];
  int int11[NBPAIRS+1][NBPAIRS+1][5][5];
  int int21[NBPAIRS+1][NBPAIRS+1][5][5][5];
  int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
  int F_ninio[5];
  double lxc;
  int MLbase;
  int MLintern[NBPAIRS+1];
  int MLclosing;
  int TerminalAU;
  int DuplexInit;
  int TETRA_ENERGY[200];
  char Tetraloops[1401];
  int Triloop_E[40];
  char Triloops[241];
  double temperature;
}  paramT;

typedef struct pf_paramT{
  int id;
  double expstack[NBPAIRS+1][NBPAIRS+1];
  double exphairpin[31]; 
  double expbulge[MAXLOOP+1];
  double expinternal[MAXLOOP+1];
  double expmismatchI[NBPAIRS+1][5][5];
  double expmismatchH[NBPAIRS+1][5][5];
  double expmismatchM[NBPAIRS+1][5][5];
  double expdangle5[NBPAIRS+1][5];
  double expdangle3[NBPAIRS+1][5];
  double expint11[NBPAIRS+1][NBPAIRS+1][5][5];
  double expint21[NBPAIRS+1][NBPAIRS+1][5][5][5];
  double expint22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
  double expninio[5][MAXLOOP+1];
  double lxc;
  double expMLbase; 
  double expMLintern[NBPAIRS+1];
  double expMLclosing;
  double expTermAU;
  double expDuplexInit;
  double exptetra[40];
  char Tetraloops[1401];
  double expTriloop[40];
  char Triloops[241];
  double temperature;
}  pf_paramT;

typedef struct fold_variables{
  int     allocated;
  int     *f5;              /* energy of 5' end */
  int     *f3;              /* energy of 3' end */
  int     *fc;              /* energy of 3' end */
  int     *c;               /* energy array, given that i-j pair */
  int     *fML;             /* multi-loop auxiliary energy array */
  int     *fM1;             /* second ML array, for subopt and circfolding */
  int     *Fmi;             /* holds row i of fML (avoids jumps in memory) */
  int     *DMLi;            /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
  int     *DMLi1;           /*             MIN(fML[i+1,k]+fML[k+1,j])  */
  int     *DMLi2;           /*             MIN(fML[i+2,k]+fML[k+1,j])  */
  int     *cc;              /* linear array for calculating canonical structures */
  int     *cc1;             /*   "     "        */
  int     Fc;
  int     FcH;         
  int     FcI;
  int     FcM;
  int     *fM2;

} fold_variables;

typedef struct pfold_variables{
  int         allocated;
  FLT_OR_DBL  *q;
  FLT_OR_DBL  *qb;
  FLT_OR_DBL  *qm;
  FLT_OR_DBL  *qm1;
  FLT_OR_DBL  *qqm;
  FLT_OR_DBL  *qqm1;
  FLT_OR_DBL  *qq;
  FLT_OR_DBL  *qq1;
  FLT_OR_DBL  *pr;
  FLT_OR_DBL  *prml;
  FLT_OR_DBL  *prm_l;
  FLT_OR_DBL  *prm_l1;
  FLT_OR_DBL  *q1k;
  FLT_OR_DBL  *qln;
  FLT_OR_DBL  qo;
  FLT_OR_DBL  qho;
  FLT_OR_DBL  qio;
  FLT_OR_DBL  qmo;
  FLT_OR_DBL  *qm2;
} pfold_variables;


typedef struct void_matrices{
  int allocated;
  int type;
  void **f;
  void **f3;
  void **f5;
  void **c;
  void **ml;
  void **m1;
  void **m2;
  void *Fc;
  void *FcH;
  void *FcI;
  void *FcM;
} void_matrices;



typedef struct common_parameters{
  unsigned int options;                 /* this is the options vector that may replace several integer switches in the future */
  int     noGU;
  int     no_closingGU;
  
  int     tetra_loop;                   /* Fold with specially stable 4-loops */
  int     energy_set;                   /* 0 = BP; 1=any mit GC; 2=any mit AU-parameter */
  int     dangles;                      /* use dangling end energies (not in part_func!) */
  char    *nonstandards;                /* contains allowed non standard bases */
  double  temperature;                  /* rescale parameters to this temperature */
  int     james_rule;                   /* interior loops of size 2 get energy 0.8Kcal and
                                           no mismatches, default 1 */
  int     logML;                        /* use logarithmic multiloop energy function */
  int     do_backtrack;                 /* calculate pair prob matrix in part_func() */
  int     noLonelyPairs;                /* avoid helices of length 1 */
  int     uniq_ML;                      /* do ML decomposition uniquely (for subopt or circ, mfe and pf algorithms depend on it) */

  int     min_hairpin;
  int     eos_debug;                    /* verbose info from energy_of_struct */

  short   alias[MAXALPHA+1];
  int     pair[MAXALPHA+1][MAXALPHA+1];
  int     rtype[8];                     /* rtype[pair[i][j]]:=pair[j][i] */
  char    Law_and_Order[10];
  int     BP_pair[NBASES][NBASES];
  int     circ;                         /* assume molecule to be circular */
  int     mirnatog;  
  int     particular_type;              /* the type of datastructure that contains the particalur parameters */
  void    *particular_parameters;       /* a pointer to the datastructure with all particular parameters */

  int     particular_array_type;        /* the type of the array datastructure */
  void    *particular_arrays;           /* a pointer to the datastructure that contains all the matrices for rcalculations */
  
} common_parameters;

typedef struct mfe_parameters{
  char  *sequence;        /* the sequence */
  char  *structure;       /* the structure constraint and mfe structure output */
  char  *ptype;           /* precomputed array of pair types */
  int   seq_length;       /* the length of thesequence */
  short *S, *S1;
  int   init_length;
  int   cut_point;        /* first position of 2nd strand for co-folding */
  int   zuker;            /* flag for cofold whether to compute zuker subopt double matrices or regular cofold */
  int   do_backtrack;     
  int   fold_constrained; /* fold with constraints */
  char  backtrack_type;   /* usually 'F'; 'C' require (1,N) to be bonded;
                              'M' seq is part of a multi loop */
  int   *BP;              /* contains the structure constrainsts: BP[i]
                             -1: | = base must be paired
                             -2: < = base must be paired with j<i
                             -3: > = base must be paired with j>i
                             -4: x = base must not pair
                             positive int: base is paired with int      */

  paramT  *P;
  int     bonus;
  bondT   *base_pair;     /* list of base pairs */

  int   *index;           /* an index for moving in the matrices */
} mfe_parameters;

typedef struct ali_mfe_parameters{
  int   n_seq;
  char  **sequences;      /* the sequences */
  char  **names;          /* the names of the sequences */
  char  *structure;       /* the structure constraint and mfe structure output */
  int   seq_length;       /* the length of thesequence */
  int   *pscore;           /* precomputed array of pair types */
  short **S, **S3,**S5;   /* S5[s][i] holds next base 5' of i in sequence s,
                             S3[s][i] holds next base 3' of i in sequence s */
  int   oldAliEn;         /* use old alifold energies (with gaps) */
  int   ribo;             /* use ribosum instead of classic covariance term */
  char   *RibosumFile;
  char  **Ss;
  unsigned short **a2s;
  double cv_fact;
  double nc_fact;
  int   init_length;
  int   cut_point;        /* first position of 2nd strand for co-folding */
  int   do_backtrack;     
  int   fold_constrained; /* fold with constraints */
  char  backtrack_type;   /* usually 'F'; 'C' require (1,N) to be bonded;
                              'M' seq is part of a multi loop */
  int   *BP;              /* contains the structure constrainsts: BP[i]
                             -1: | = base must be paired
                             -2: < = base must be paired with j<i
                             -3: > = base must be paired with j>i
                             -4: x = base must not pair
                             positive int: base is paired with int      */

  paramT  *P;
  int     bonus;
  bondT   *base_pair;     /* list of base pairs */

  int   *index;           /* an index for moving in the matrices */
} ali_mfe_parameters;

typedef struct pf_parameters{
  char  *sequence;        /* the sequence */
  char  *structure;       /* the structure constraint and mfe structure output */
  char  *ptype;           /* precomputed array of pair types */
  int   seq_length;       /* the length of thesequence */
  short *S, *S1;
  int   init_length;
  double  init_temp;
  int   cut_point;        /* first position of 2nd strand for co-folding */
  int   do_backtrack;     /* calculate pair prob matrix in part_func() */
  int   fold_constrained; /* fold with constraints */
  char  backtrack_type;   /* usually 'F'; 'C' require (1,N) to be bonded;
                              'M' seq is part of a multi loop */
  int   *index;           /* an index for moving in the matrices */
  FLT_OR_DBL  *scale;
  FLT_OR_DBL  *exphairpin;
  FLT_OR_DBL  *expMLbase;
  FLT_OR_DBL  pf_scale;
  pf_paramT   *pf_params;   /* holds all [unscaled] pf energy parameters */
  int         *jindx;       /* index for moving in qm1*/
  int   *BP;              /* contains the structure constrainsts: BP[i]
                             -1: | = base must be paired
                             -2: < = base must be paired with j<i
                             -3: > = base must be paired with j>i
                             -4: x = base must not pair
                             positive int: base is paired with int      */

} pf_parameters;

typedef struct subopt_parameters{
  char  *sequence;        /* the sequence */
  char  *structure;       /* the structure constraint and mfe structure output */
  char  *ptype;           /* precomputed array of pair types */
  int   seq_length;       /* the length of thesequence */
  short *S, *S1;
  int   init_length;
  double  init_temp;
  int   cut_point;        /* first position of 2nd strand for co-folding */
  int   do_backtrack;     /* calculate pair prob matrix in part_func() */
  int   fold_constrained; /* fold with constraints */
  char  backtrack_type;   /* usually 'F'; 'C' require (1,N) to be bonded;
                              'M' seq is part of a multi loop */
  int   *index;           /* an index for moving in the matrices */

  int   *BP;              /* contains the structure constrainsts: BP[i]
                             -1: | = base must be paired
                             -2: < = base must be paired with j<i
                             -3: > = base must be paired with j>i
                             -4: x = base must not pair
                             positive int: base is paired with int      */
  paramT  *P;
  int     subopt_sorted;                /* sort output by energy */
  LIST    *Stack;
  int     nopush;
  int     delta;
  int     best_energy;                  /* best_energy = remaining energy */
  int     minimal_energy;               /* minimum free energy */
  int     element_energy;               /* internal energy of a structural element */
  int     threshold;                    /* minimal_energy + delta */
  double  print_energy;                 /* printing threshold for use with logML */
  int     density_of_states[MAXDOS+1];
} subopt_parameters;

/*
* ############################################################
* SUBOPT data structures
* ############################################################
*/
typedef struct {
    char *structure;
    LIST *Intervals;
    int partial_energy;
    /* int best_energy;   */ /* best attainable energy */
} STATE;

typedef struct {
  unsigned int i;
  unsigned int j;
} PAIR;

typedef struct {
    unsigned int i;
    unsigned int j;
    int array_flag;
} INTERVAL;

typedef struct {
  float energy;                            /* energy of structure */
  char *structure;
} SOLUTION;

/*
* ############################################################
* COFOLD data structures
* ############################################################
*/
typedef struct cofoldF {
  /* free energies for: */
  double F0AB; /* null model without DuplexInit */
  double FAB;  /* all states with DuplexInit corretion */
  double FcAB; /* true hybrid states only */
  double FA;   /* monomer A */
  double FB;   /* monomer B */
} cofoldF;

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

typedef struct {
   unsigned i;        /* i,j in [0, n-1] */
   unsigned j;
   float p;      /* probability */
   float ent;    /* pseudo entropy for p(i,j) = S_i + S_j - p_ij*ln(p_ij) */
   short bp[8];  /* frequencies of pair_types */
   char comp;    /* 1 iff pair is in mfe structure */
}  pair_info;


/* some data structures for findpath.c */

typedef struct move {
  int i;  /* i,j>0 insert; i,j<0 delete */
  int j;
  int when;  /* 0 if still available, else resulting distance from start */
  int E;
} move_t;

typedef struct intermediate {
  short *pt;     /* pair table */
  int Sen;       /* saddle energy so far */
  int curr_en;   /* current energy */
  move_t *moves; /* remaining moves to target */
} intermediate_t;

typedef struct path {
  double en;
  char *s;
} path_t;

#endif
