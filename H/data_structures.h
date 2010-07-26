#ifndef __VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H__
#define __VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H__

#include "energy_const.h"

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
  int i;
  int j;
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
*** The datastructure that contains temperature scaled energy parameters.
**/
typedef struct paramT{
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
  double temperature;
}  paramT;

/**
*** The datastructure that contains temperature scaled Boltzmann weights of the energy parameters.
**/
typedef struct pf_paramT{
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


#endif
