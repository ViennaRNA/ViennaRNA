#ifndef __VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H__
#define __VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H__

#include "energy_const.h"

/**
 *  \file data_structures.h
 *  \brief All datastructures and typedefs shared among the Vienna RNA Package can be found here
 */

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

/**
 *  \brief this datastructure is used as input parameter in functions of PS_dot.h and others
 */
typedef struct plist {
  int i;
  int j;
  float p;
} plist;

/* this datastructure is used as input parameter in functions of PS_dot.c */
typedef struct cpair {
  int i,j,mfe;
  float p, hue, sat;
} cpair;

/**
 *  \brief this is a workarround for the SWIG Perl Wrapper RNA plot function
 *  that returns an array of type COORDINATE
 */
typedef struct {
  float X; /* X coords */
  float Y; /* Y coords */
} COORDINATE;

typedef struct sect {
  int  i;
  int  j;
  int ml;
} sect; /* stack of partial structures for backtracking */

typedef struct bondT {               /* base pair */
   unsigned int i;
   unsigned int j;
} bondT;

typedef struct bondTEn {               /* base pair with associated energy*/
   int i;
   int j;
   int energy;
} bondTEn;

/**
 *  The datastructure that contains temperature scaled energy parameters.
 */
typedef struct{
  int id;
  int stack[NBPAIRS+1][NBPAIRS+1];
  int hairpin[31];
  int bulge[MAXLOOP+1];
  int internal_loop[MAXLOOP+1];
  int mismatchExt[NBPAIRS+1][5][5];
  int mismatchI[NBPAIRS+1][5][5];
  int mismatch1nI[NBPAIRS+1][5][5];
  int mismatch23I[NBPAIRS+1][5][5];
  int mismatchH[NBPAIRS+1][5][5];
  int mismatchM[NBPAIRS+1][5][5];
  int dangle5[NBPAIRS+1][5];
  int dangle3[NBPAIRS+1][5];
  int int11[NBPAIRS+1][NBPAIRS+1][5][5];
  int int21[NBPAIRS+1][NBPAIRS+1][5][5][5];
  int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
  int ninio[5];
  double lxc;
  int MLbase;
  int MLintern[NBPAIRS+1];
  int MLclosing;
  int TerminalAU;
  int DuplexInit;
  int Tetraloop_E[200];
  char Tetraloops[1401];
  int Triloop_E[40];
  char Triloops[241];
  int Hexaloop_E[40];
  char Hexaloops[1801];
  int TripleC;
  int MultipleCA;
  int MultipleCB;
  double temperature;
}  paramT;

/**
 *  The datastructure that contains temperature scaled Boltzmann weights of the energy parameters.
 */
typedef struct{
  int     id;
  double  expstack[NBPAIRS+1][NBPAIRS+1];
  double  exphairpin[31];
  double  expbulge[MAXLOOP+1];
  double  expinternal[MAXLOOP+1];
  double  expmismatchExt[NBPAIRS+1][5][5];
  double  expmismatchI[NBPAIRS+1][5][5];
  double  expmismatch23I[NBPAIRS+1][5][5];
  double  expmismatch1nI[NBPAIRS+1][5][5];
  double  expmismatchH[NBPAIRS+1][5][5];
  double  expmismatchM[NBPAIRS+1][5][5];
  double  expdangle5[NBPAIRS+1][5];
  double  expdangle3[NBPAIRS+1][5];
  double  expint11[NBPAIRS+1][NBPAIRS+1][5][5];
  double  expint21[NBPAIRS+1][NBPAIRS+1][5][5][5];
  double  expint22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
  double  expninio[5][MAXLOOP+1];
  double  lxc;
  double  expMLbase;
  double  expMLintern[NBPAIRS+1];
  double  expMLclosing;
  double  expTermAU;
  double  expDuplexInit;
  double  exptetra[40];
  double  exptri[40];
  double  exphex[40];
  char    Tetraloops[1401];
  double  expTriloop[40];
  char    Triloops[241];
  char    Hexaloops[1801];
  double  expTripleC;
  double  expMultipleCA;
  double  expMultipleCB;
  double  temperature;
  double  kT;
}  pf_paramT;



/*
* ############################################################
* SUBOPT data structures
* ############################################################
*/

typedef struct {
  int i;
  int j;
} PAIR;

typedef struct {
    int i;
    int j;
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

/**
 *  \brief A base pair info structure
 * 
 *  for each base pair (i,j) the structure lists: its probability
 *  'p', an entropy-like measure for its well-definedness 'ent',
 *  and in 'bp[]' the frequency of each type of pair. 'bp[0]'
 *  contains the number of non-compatible sequences, 'bp[1]' the
 *  number of CG pairs, etc.
 */
typedef struct {
   unsigned i;        /* i,j in [0, n-1] */
   unsigned j;
   float p;      /* probability */
   float ent;    /* pseudo entropy for p(i,j) = S_i + S_j - p_ij*ln(p_ij) */
   short bp[8];  /* frequencies of pair_types */
   char comp;    /* 1 iff pair is in mfe structure */
} pair_info;


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

/*
* ############################################################
* RNAup data structures
* ############################################################
*/
typedef struct pu_contrib { /* contributions to prob_unpaired in */
  double **H; /* hairpin loops */
  double **I; /* interior loops */
  double **M; /* multi loops */
  double **E; /* exterior loop */
  int length; /* length of the input sequence */
  int w;      /* longest unpaired region */
} pu_contrib;

typedef struct interact { /* contributions to prob_unpaired in */
  double *Pi; /* probabilities of interaction */
  double *Gi; /* free energies of interaction */
  double Gikjl; /* full free energy for interaction between [k,i] k<i
                   in longer seq and [j,l] j<l in shorter seq */
  double Gikjl_wo; /* Gikjl without contributions for prob_unpaired */
  int i; /* k<i in longer seq */
  int k; /* k<i in longer seq */
  int j; /*j<l in shorter seq */
  int l; /*j<l in shorter seq */
  int length; /* length of longer sequence */
} interact;

typedef struct pu_out { /* collect all free_energy of beeing unpaired
                           values for output */
  int len;        /* sequence length */
  int u_vals;     /* number of different -u values */
  int contribs;   /* [-c "SHIME"] */
  char **header;  /* header line */
  double **u_values; /* (differnet -u values * [-c "SHIME"]) * seq len */
} pu_out;

typedef struct constrain { /* constrains for cofolding */
  int *indx;
  char *ptype;
} constrain;

/*
* ############################################################
* RNAduplex data structures
* ############################################################
*/

typedef struct {
  int i;
  int j;
  char *structure;
  float energy;
} duplexT;

/*
* ############################################################
* PKplex data structures
* ############################################################
*/

typedef struct dupVar{
  int i;
  int j;
  int end;
  char *structure;
  double energy;
  int offset;
  double dG1;
  double dG2;
  double ddG;
  int tb;
  int te;
  int qb;
  int qe;
} dupVar;



/*
* ############################################################
* 2Dfold data structures
* ############################################################
*/
typedef struct{
  int k;
  int l;
  float en;
  char *s;
} TwoDfold_solution;

typedef struct{
  paramT          *P;
  int             do_backtrack;
  char            *ptype;   /* precomputed array of pair types */
  char            *sequence;
  short           *S, *S1;
  unsigned int    maxD1;
  unsigned int    maxD2;


  unsigned int    *mm1;         /* maximum matching matrix, reference struct 1 disallowed */
  unsigned int    *mm2;         /* maximum matching matrix, reference struct 2 disallowed */

  int             *my_iindx;    /* index for moving in quadratic distancy dimsensions */

  double          temperature;

  unsigned int    *referenceBPs1; /* matrix containing number of basepairs of reference structure1 in interval [i,j] */
  unsigned int    *referenceBPs2; /* matrix containing number of basepairs of reference structure2 in interval [i,j] */
  unsigned int    *bpdist;        /* matrix containing base pair distance of reference structure 1 and 2 on interval [i,j] */

  short           *reference_pt1;
  short           *reference_pt2;
  int             circ;
  int             dangles;
  unsigned int    seq_length;

  int             ***E_F5;
  int             ***E_F3;
  int             ***E_C;
  int             ***E_M;
  int             ***E_M1;
  int             ***E_M2;

  int             **E_Fc;
  int             **E_FcH;
  int             **E_FcI;
  int             **E_FcM;

  int             **l_min_values;
  int             **l_max_values;
  int             *k_min_values;
  int             *k_max_values;

  int             **l_min_values_m;
  int             **l_max_values_m;
  int             *k_min_values_m;
  int             *k_max_values_m;

  int             **l_min_values_m1;
  int             **l_max_values_m1;
  int             *k_min_values_m1;
  int             *k_max_values_m1;

  int             **l_min_values_f;
  int             **l_max_values_f;
  int             *k_min_values_f;
  int             *k_max_values_f;

  int             **l_min_values_f3;
  int             **l_max_values_f3;
  int             *k_min_values_f3;
  int             *k_max_values_f3;

  int             **l_min_values_m2;
  int             **l_max_values_m2;
  int             *k_min_values_m2;
  int             *k_max_values_m2;

  int             *l_min_values_fc;
  int             *l_max_values_fc;
  int             k_min_values_fc;
  int             k_max_values_fc;

  int             *l_min_values_fcH;
  int             *l_max_values_fcH;
  int             k_min_values_fcH;
  int             k_max_values_fcH;

  int             *l_min_values_fcI;
  int             *l_max_values_fcI;
  int             k_min_values_fcI;
  int             k_max_values_fcI;

  int             *l_min_values_fcM;
  int             *l_max_values_fcM;
  int             k_min_values_fcM;
  int             k_max_values_fcM;

  /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
  int             *E_F5_rem;
  int             *E_F3_rem;
  int             *E_C_rem;
  int             *E_M_rem;
  int             *E_M1_rem;
  int             *E_M2_rem;

  int             E_Fc_rem;
  int             E_FcH_rem;
  int             E_FcI_rem;
  int             E_FcM_rem;

#ifdef COUNT_STATES
  unsigned long             ***N_F5;
  unsigned long             ***N_C;
  unsigned long             ***N_M;
  unsigned long             ***N_M1;
#endif
} TwoDfold_vars;

typedef struct{
  int k;
  int l;
  FLT_OR_DBL  q;
} TwoDpfold_solution;

typedef struct{

  unsigned int    alloc;
  char            *ptype;   /* precomputed array of pair types */
  char            *sequence;
  short           *S, *S1;
  double          temperature;      /* temperature in last call to scale_pf_params */
  double          init_temp;      /* temperature in last call to scale_pf_params */
  unsigned int    maxD1;
  unsigned int    maxD2;

  FLT_OR_DBL  *scale;
  FLT_OR_DBL  pf_scale;
  pf_paramT   *pf_params;     /* holds all [unscaled] pf parameters */

  int             *my_iindx;         /* index for moving in quadratic distancy dimsensions */
  int             *jindx;         /* index for moving in the triangle matrix qm1 */

  unsigned int    *referenceBPs1;    /* matrix containing number of basepairs of reference structure1 in interval [i,j] */
  unsigned int    *referenceBPs2;    /* matrix containing number of basepairs of reference structure2 in interval [i,j] */
  short           *reference_pt1;
  short           *reference_pt2;
  unsigned int    *mm1;         /* maximum matching matrix, reference struct 1 disallowed */
  unsigned int    *mm2;         /* maximum matching matrix, reference struct 2 disallowed */
  unsigned int    *bpdist;      /* matrix containing base pair distance of reference structure 1 and 2 on interval [i,j] */
  int             circ;
  int             dangles;
  unsigned int    seq_length;

  FLT_OR_DBL      ***Q;
  FLT_OR_DBL      ***Q_B;
  FLT_OR_DBL      ***Q_M;
  FLT_OR_DBL      ***Q_M1;
  FLT_OR_DBL      ***Q_M2;

  FLT_OR_DBL      **Q_c;
  FLT_OR_DBL      **Q_cH;
  FLT_OR_DBL      **Q_cI;
  FLT_OR_DBL      **Q_cM;

  int             **l_min_values;
  int             **l_max_values;
  int             *k_min_values;
  int             *k_max_values;

  int             **l_min_values_b;
  int             **l_max_values_b;
  int             *k_min_values_b;
  int             *k_max_values_b;

  int             **l_min_values_m;
  int             **l_max_values_m;
  int             *k_min_values_m;
  int             *k_max_values_m;

  int             **l_min_values_m1;
  int             **l_max_values_m1;
  int             *k_min_values_m1;
  int             *k_max_values_m1;

  int             **l_min_values_m2;
  int             **l_max_values_m2;
  int             *k_min_values_m2;
  int             *k_max_values_m2;

  int             *l_min_values_qc;
  int             *l_max_values_qc;
  int             k_min_values_qc;
  int             k_max_values_qc;

  int             *l_min_values_qcH;
  int             *l_max_values_qcH;
  int             k_min_values_qcH;
  int             k_max_values_qcH;

  int             *l_min_values_qcI;
  int             *l_max_values_qcI;
  int             k_min_values_qcI;
  int             k_max_values_qcI;

  int             *l_min_values_qcM;
  int             *l_max_values_qcM;
  int             k_min_values_qcM;
  int             k_max_values_qcM;

  /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
  FLT_OR_DBL      *Q_rem;
  FLT_OR_DBL      *Q_B_rem;
  FLT_OR_DBL      *Q_M_rem;
  FLT_OR_DBL      *Q_M1_rem;
  FLT_OR_DBL      *Q_M2_rem;

  FLT_OR_DBL      Q_c_rem;
  FLT_OR_DBL      Q_cH_rem;
  FLT_OR_DBL      Q_cI_rem;
  FLT_OR_DBL      Q_cM_rem;

} TwoDpfold_vars;

#endif
