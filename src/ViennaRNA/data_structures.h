#ifndef __VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H__
#define __VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H__

#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/params.h>


/**
 *  \file data_structures.h
 *
 *  \addtogroup   data_structures   Common Data Structures and Preprocessor Macros
 *  @{
 *
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

/**
 *  \brief Maximum density of states discretization for subopt
 */
#define MAXDOS                1000

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
  int type;
} plist;

/**
 *  \brief this datastructure is used as input parameter in functions of PS_dot.c
 */
typedef struct cpair {
  int i,j,mfe;
  float p, hue, sat;
} cpair;


/**
 *  \brief  Stack of partial structures for backtracking
 */
typedef struct sect {
  int  i;
  int  j;
  int ml;
} sect;

/**
 *  \brief  Base pair
 */
typedef struct bondT {
   unsigned int i;
   unsigned int j;
} bondT;

/**
 *  \brief  Base pair with associated energy
 */
typedef struct bondTEn {
   int i;
   int j;
   int energy;
} bondTEn;


/*
* ############################################################
* SUBOPT data structures
* ############################################################
*/

/**
 *  \brief  Base pair data structure used in subopt.c
 */
typedef struct {
  int i;
  int j;
} PAIR;

/**
 *  \brief  Sequence interval stack element used in subopt.c
 */
typedef struct {
    int i;
    int j;
    int array_flag;
} INTERVAL;

/**
 *  \brief  Solution element from subopt.c
 */
typedef struct {
  float energy;       /**< \brief Free Energy of structure in kcal/mol */
  char *structure;    /**< \brief Structure in dot-bracket notation */
} SOLUTION;

/*
* ############################################################
* COFOLD data structures
* ############################################################
*/

/**
 *  \brief
 */
typedef struct cofoldF {
  /* free energies for: */
  double F0AB;  /**< \brief Null model without DuplexInit */
  double FAB;   /**< \brief all states with DuplexInit correction */
  double FcAB;  /**< \brief true hybrid states only */
  double FA;    /**< \brief monomer A */
  double FB;    /**< \brief monomer B */
} cofoldF;

/**
 *  \brief
 */
typedef struct ConcEnt {
  double A0;    /**< \brief start concentration A */
  double B0;    /**< \brief start concentration B */
  double ABc;   /**< \brief End concentration AB */
  double AAc;
  double BBc;
  double Ac;
  double Bc;
} ConcEnt;

/**
 *  \brief
 */
typedef struct pairpro{
  struct plist *AB;
  struct plist *AA;
  struct plist *A;
  struct plist *B;
  struct plist *BB;
} pairpro;

/**
 *  \brief A base pair info structure
 *
 *  For each base pair (i,j) with i,j in [0, n-1] the structure lists:
 *  - its probability 'p'
 *  - an entropy-like measure for its well-definedness 'ent'
 *  - the frequency of each type of pair in 'bp[]'
 *    + 'bp[0]' contains the number of non-compatible sequences
 *    + 'bp[1]' the number of CG pairs, etc.
 */
typedef struct {
   unsigned i;    /**<  \brief  nucleotide position i */
   unsigned j;    /**<  \brief  nucleotide position j */
   float p;       /**< \brief  Probability */
   float ent;     /**< \brief  Pseudo entropy for \f$ p(i,j) = S_i + S_j - p_ij*ln(p_ij) \f$ */
   short bp[8];   /**< \brief  Frequencies of pair_types */
   char comp;     /**< \brief  1 iff pair is in mfe structure */
} pair_info;


/*
* ############################################################
* FINDPATH data structures
* ############################################################
*/

/**
 *  \brief
 */
typedef struct move {
  int i;  /* i,j>0 insert; i,j<0 delete */
  int j;
  int when;  /* 0 if still available, else resulting distance from start */
  int E;
} move_t;

/**
 *  \brief
 */
typedef struct intermediate {
  short *pt;      /**<  \brief  pair table */
  int Sen;        /**<  \brief  saddle energy so far */
  int curr_en;    /**<  \brief  current energy */
  move_t *moves;  /**<  \brief  remaining moves to target */
} intermediate_t;

/**
 *  \brief
 */
typedef struct path {
  double en;
  char *s;
} path_t;

/*
* ############################################################
* RNAup data structures
* ############################################################
*/

/**
 *  \brief contributions to p_u
 */
typedef struct pu_contrib {
  double **H; /**<  \brief  hairpin loops */
  double **I; /**<  \brief  interior loops */
  double **M; /**<  \brief  multi loops */
  double **E; /**<  \brief  exterior loop */
  int length; /**<  \brief  length of the input sequence */
  int w;      /**<  \brief  longest unpaired region */
} pu_contrib;

/**
 *  \brief
 */
typedef struct interact {
  double *Pi;       /**<  \brief  probabilities of interaction */
  double *Gi;       /**<  \brief  free energies of interaction */
  double Gikjl;     /**<  \brief  full free energy for interaction between [k,i] k<i
                                  in longer seq and [j,l] j<l in shorter seq */
  double Gikjl_wo;  /**<  \brief  Gikjl without contributions for prob_unpaired */
  int i;            /**<  \brief  k<i in longer seq */
  int k;            /**<  \brief  k<i in longer seq */
  int j;            /**<  \brief  j<l in shorter seq */
  int l;            /**<  \brief  j<l in shorter seq */
  int length;       /**<  \brief  length of longer sequence */
} interact;

/**
 *  \brief  Collection of all free_energy of beeing unpaired values for output
 */
typedef struct pu_out {
  int len;            /**<  \brief  sequence length */
  int u_vals;         /**<  \brief  number of different -u values */
  int contribs;       /**<  \brief  [-c "SHIME"] */
  char **header;      /**<  \brief  header line */
  double **u_values;  /**<  \brief  (the -u values * [-c "SHIME"]) * seq len */
} pu_out;

/**
 *  \brief  constraints for cofolding
 */
typedef struct constrain{
  int *indx;
  char *ptype;
} constrain;

/*
* ############################################################
* RNAduplex data structures
* ############################################################
*/

/**
 *  \brief
 */
typedef struct {
  int i;
  int j;
  int end;
  char *structure;
  double energy;
  double energy_backtrack;
  double opening_backtrack_x;
  double opening_backtrack_y;
  int offset;
  double dG1;
  double dG2;
  double ddG;
  int tb;
  int te;
  int qb;
  int qe;
} duplexT;

/*
* ############################################################
* RNAsnoop data structures
* ############################################################
*/

/**
 *  \brief
 */
typedef struct node {
  int k;
  int energy;
  struct node *next;
} folden;

/**
 *  \brief
 */
typedef struct {
  int i;
  int j;
  int u;
  char *structure;
  float energy;
  float Duplex_El;
  float Duplex_Er;
  float Loop_E;
  float Loop_D;
  float pscd;
  float psct;
  float pscg;
  float Duplex_Ol;
  float Duplex_Or;
  float Duplex_Ot;
  float fullStemEnergy;
} snoopT;


/*
* ############################################################
* PKplex data structures
* ############################################################
*/

/**
 *  \brief
 */
typedef struct dupVar{
  int i;
  int j;
  int end;
  char *pk_helix;
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
  int inactive;
  int processed;
} dupVar;

/**
 * @}
 */


/*
* ############################################################
* VRNA fold compound related functions
* ############################################################
*/

/**
 *  \addtogroup   basic_data_structures  Basic Data Structures for Structure Prediction and Evaluation
 *  @{
 *
 *  \brief  This module provides interfaces that deal with the most basic data structures used
 *          in structure predicting and energy evaluating function of the RNAlib.
 *
 *          Throughout the RNAlib, a data structure, the #vrna_fold_compound, is used to group
 *          information and data that is required for structure prediction and energy evaluation.
 *          Here, you'll find interface functions to create, modify, and delete #vrna_fold_compound
 *          data structures.
 */

typedef struct{
  unsigned int allocated; /* flag keeper for fast evaluation which matrices have been allocated */
  unsigned int length;
  int     *c;   /* energy array, given that i-j pair */
  int     *f5;  /* energy of 5' end */
  int     *f3;  /* energy of 3' end */
  int     *fc;  /* energy from i to cutpoint (and vice versa if i>cut) */
  int     *fML; /* multi-loop auxiliary energy array */
  int     *fM1; /* second ML array, only for subopt */
  int     *fM2; /* fM2 = multiloop region with exactly two stems, extending to 3' end */
  int     *ggg; /* energies of g-quadruplexes */
  int     Fc;   /* parts of the exterior loop energies for circfolding */
  int     FcH;
  int     FcI;
  int     FcM;
} vrna_mx_mfe_t;

typedef struct{
  unsigned int allocated;
  unsigned int length;
  FLT_OR_DBL  *q;
  FLT_OR_DBL  *qb;
  FLT_OR_DBL  *qm;
  FLT_OR_DBL  *qm1;
  FLT_OR_DBL  *probs;
  FLT_OR_DBL  *q1k;
  FLT_OR_DBL  *qln;
  FLT_OR_DBL  *G;

  FLT_OR_DBL  qo;
  FLT_OR_DBL  *qm2;
  FLT_OR_DBL  qho;
  FLT_OR_DBL  qio;
  FLT_OR_DBL  qmo;

  FLT_OR_DBL  *scale;
  FLT_OR_DBL  *expMLbase;
} vrna_mx_pf_t;

typedef struct{

  unsigned int  type;            /**<  \brief The type of the fold_compound */
  unsigned int  length;
  int           cutpoint;               /**<  \brief  The position of the (cofold) cutpoint within the provided sequence.
                                      If there is no cutpoint, this field will be set to -1
                                */
  union {
    struct {
      char  *sequence;
      short *sequence_encoding;
      short *sequence_encoding2;
      char  *ptype;                 /**<  \brief Pair type array
                                     
                                          Contains the numerical encoding of the pair type for each pair (i,j) used
                                          in MFE, Partition function and Evaluation computations.
                                          \note This array is always indexed via jindx, in contrast to previously
                                          different indexing between mfe and pf variants!
                                          \see  get_indx(), vrna_get_ptypes()
                                    */
      char  *ptype_pf_compat;       /**<  \brief  ptype array indexed via iindx
                                          \deprecated  This attribute will vanish in the future!
                                          It's meant for backward compatibility only!
                                    */
      struct vrna_scT *sc;
    };
    struct {
      char  **sequences;
      unsigned int    n_seq;
      char            *cons_seq;
      short           *S_cons;
      short           **S;
      short           **S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
      short           **S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
      char            **Ss;
      unsigned short  **a2s;
      int             *pscore;     /* precomputed array of pair types */
      struct vrna_scT **scs;
      int             oldAliEn;
    };
  };

  struct vrna_hcT   *hc;

  vrna_mx_mfe_t     *matrices;
  vrna_mx_pf_t      *exp_matrices;

  struct paramT     *params;
  struct pf_paramT  *exp_params;

  int               *iindx;
  int               *jindx;

} vrna_fold_compound;

/**
 *  \brief fold_compound_type Single Sequence
 */
#define   VRNA_VC_TYPE_SINGLE     1

/**
 *  \brief fold_compound_type Sequence Alignment
 */
#define   VRNA_VC_TYPE_ALIGNMENT  2

/* the definitions below should be used for functions that return/receive/destroy fold compound data structures */

/**
 *  \brief  Option flag to specify requirement of Minimum Free Energy (MFE) DP matrices
 *          and corresponding set of energy parameters
 *
 *  \see vrna_get_fold_compound(), vrna_get_fold_compound_ali(), #VRNA_OPTION_EVAL_ONLY
 */
#define VRNA_OPTION_MFE             1

/**
 *  \brief  Option flag to specify requirement of Partition Function (PF) DP matrices
 *          and corresponding set of Boltzmann factors
 *
 *  \see vrna_get_fold_compound(), vrna_get_fold_compound_ali(), #VRNA_OPTION_EVAL_ONLY
 */
#define VRNA_OPTION_PF              2

#define VRNA_OPTION_HYBRID          4

/**
 *  \brief  Option flag to specify that neither MFE, nor PF DP matrices are required
 *
 *  Use this flag in conjuntion with #VRNA_OPTION_MFE, and #VRNA_OPTION_PF to save
 *  memory for a #vrna_fold_compound obtained from vrna_get_fold_compound(), or vrna_get_fold_compound_ali()
 *  in cases where only energy evaluation but no structure prediction is required.
 *
 *  \see vrna_get_fold_compound(), vrna_get_fold_compound_ali(), vrna_eval_structure()
 */
#define VRNA_OPTION_EVAL_ONLY       8



/**
 *  \brief  Retrieve a #vrna_fold_compound data structure for single sequences and hybridizing sequences
 *
 *  This function provides an easy interface to obtain a prefilled #vrna_fold_compound by passing a single
 *  sequence, or two contatenated sequences as input. For the latter, sequences need to be seperated by
 *  an '&' character like this: \verbatim char *sequence = "GGGG&CCCC"; \endverbatim
 *
 *  The optional parameter 'md_p' can be used to specify the model details for computations on the #vrna_fold_compounds
 *  content. The third parameter 'options' is used to specify the DP matrix requirements and the corresponding set
 *  of energy parameters. Use the macros:
 *
 *  - #VRNA_OPTION_MFE
 *  - #VRNA_OPTION_PF
 *  - #VRNA_OPTION_EVAL_ONLY
 *
 *  to specify the required type of computations that will be performed with the #vrna_fold_compound.
 *
 *  \note The sequence string must be uppercase, and should contain only RNA (resp. DNA) alphabet depending
 *        on what energy parameter set is used
 *
 *  \see  vrna_get_fold_compound_ali(), #vrna_md_t, #VRNA_OPTION_MFE, #VRNA_OPTION_PF, #VRNA_OPTION_EVAL_ONLY
 *
 *  \param    sequence    A single sequence, or two concatenated sequences seperated by an '&' character
 *  \param    md_p        An optional set of model details
 *  \param    options     The options for DP matrices memory allocation
 *  \return               A prefilled vrna_fold_compound that can be readily used for computations
 */
vrna_fold_compound *vrna_get_fold_compound( const char *sequence,
                                            vrna_md_t *md_p,
                                            unsigned int options);

/**
 *  \brief  Retrieve a #vrna_fold_compound data structure for sequence alignments
 *
 *  This function provides an easy interface to obtain a prefilled #vrna_fold_compound by passing an
 *  alignment of sequences.
 *
 *  The optional parameter 'md_p' can be used to specify the model details for computations on the #vrna_fold_compounds
 *  content. The third parameter 'options' is used to specify the DP matrix requirements and the corresponding set
 *  of energy parameters. Use the macros:
 *
 *  - #VRNA_OPTION_MFE
 *  - #VRNA_OPTION_PF
 *  - #VRNA_OPTION_EVAL_ONLY
 *
 *  to specify the required type of computations that will be performed with the #vrna_fold_compound.
 *
 *  \note The sequence strings must be uppercase, and should contain only RNA (resp. DNA) alphabet including
 *        gap characters depending on what energy parameter set is used.
 *
 *  \see  vrna_get_fold_compound(), #vrna_md_t, #VRNA_OPTION_MFE, #VRNA_OPTION_PF, #VRNA_OPTION_EVAL_ONLY,
 *        read_clustal()
 *
 *  \param    sequences   A sequence alignment including 'gap' characters
 *  \param    md_p        An optional set of model details
 *  \param    options     The options for DP matrices memory allocation
 *  \return               A prefilled vrna_fold_compound that can be readily used for computations
 */
vrna_fold_compound *vrna_get_fold_compound_ali( const char **sequences,
                                                vrna_md_t *md_p,
                                                unsigned int options);

/**
 *  \brief  Free memory occupied by a #vrna_fold_compound
 *
 *  \see vrna_get_fold_compound(), vrna_get_fold_compound_ali(), vrna_free_mfe_matrices(), vrna_free_pf_matrices()
 *
 *  \param  vc  The #vrna_fold_compound that is to be erased from memory
 */
void vrna_free_fold_compound(vrna_fold_compound *vc);

/**
 *  \brief  Free memory occupied by the Minimum Free Energy (MFE) Dynamic Programming (DP) matrices
 *
 *  \see vrna_get_fold_compound(), vrna_get_fold_compound_ali(), vrna_free_fold_compound(), vrna_free_pf_matrices()
 *
 *  \param  vc  The #vrna_fold_compound storing the MFE DP matrices that are to be erased from memory
 */
void vrna_free_mfe_matrices(vrna_fold_compound *vc);

/**
 *  \brief  Free memory occupied by the Partition Function (PF) Dynamic Programming (DP) matrices
 *
 *  \see vrna_get_fold_compound(), vrna_get_fold_compound_ali(), vrna_free_fold_compound(), vrna_free_mfe_matrices()
 *
 *  \param  vc  The #vrna_fold_compound storing the PF DP matrices that are to be erased from memory
 */
void vrna_free_pf_matrices(vrna_fold_compound *vc);

/**
 *  @}
 */

#endif
