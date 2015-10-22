#ifndef VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H
#define VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H

/**
 *  @file data_structures.h
 *
 *  @addtogroup   data_structures   Common Data Structures and Preprocessor Macros
 *  @{
 *
 *  @brief All datastructures and typedefs shared among the Vienna RNA Package can be found here
 */

/* below are several typedef's we use throughout the ViennaRNA library */

/** @brief Typename for the fold_compound data structure #vrna_fc
 *  @ingroup basic_data_structures
 */
typedef struct vrna_fc        vrna_fold_compound;

/** @brief Typename for the base pair repesenting data structure #vrna_basepair */
typedef struct vrna_basepair  vrna_basepair_t;

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

#ifdef VRNA_BACKWARD_COMPAT

typedef struct vrna_basepair  PAIR;
typedef struct vrna_subopt_solution   SOLUTION;
typedef struct vrna_subopt_interval   INTERVAL;

#endif

#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/dp_matrices.h>
#include <ViennaRNA/constraints.h>

/*
* ############################################################
* Here are the type definitions of various datastructures
* shared among the Vienna RNA Package
* ############################################################
*/

/**
 *  @brief this datastructure is used as input parameter in functions of PS_dot.h and others
 */
typedef struct plist {
  int i;
  int j;
  float p;
  int type;
} plist;

/**
 *  @brief this datastructure is used as input parameter in functions of PS_dot.c
 */
typedef struct cpair {
  int i,j,mfe;
  float p, hue, sat;
} cpair;


/**
 *  @brief  Stack of partial structures for backtracking
 */
typedef struct sect {
  int  i;
  int  j;
  int ml;
} sect;

/**
 *  @brief  Base pair
 */
typedef struct bondT {
   unsigned int i;
   unsigned int j;
} bondT;

/**
 *  @brief  Base pair with associated energy
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
 *  @brief  Base pair data structure used in subopt.c
 */
struct vrna_basepair {
  int i;
  int j;
};

/**
 *  @brief  Sequence interval stack element used in subopt.c
 */
struct vrna_subopt_interval {
    int i;
    int j;
    int array_flag;
};

/**
 *  @brief  Solution element from subopt.c
 */
struct vrna_subopt_solution {
  float energy;       /**< @brief Free Energy of structure in kcal/mol */
  char *structure;    /**< @brief Structure in dot-bracket notation */
};

/*
* ############################################################
* COFOLD data structures
* ############################################################
*/

/**
 *  @brief
 */
typedef struct cofoldF {
  /* free energies for: */
  double F0AB;  /**< @brief Null model without DuplexInit */
  double FAB;   /**< @brief all states with DuplexInit correction */
  double FcAB;  /**< @brief true hybrid states only */
  double FA;    /**< @brief monomer A */
  double FB;    /**< @brief monomer B */
} cofoldF;

/**
 *  @brief
 */
typedef struct ConcEnt {
  double A0;    /**< @brief start concentration A */
  double B0;    /**< @brief start concentration B */
  double ABc;   /**< @brief End concentration AB */
  double AAc;
  double BBc;
  double Ac;
  double Bc;
} ConcEnt;

/**
 *  @brief
 */
typedef struct pairpro{
  struct plist *AB;
  struct plist *AA;
  struct plist *A;
  struct plist *B;
  struct plist *BB;
} pairpro;

/**
 *  @brief A base pair info structure
 *
 *  For each base pair (i,j) with i,j in [0, n-1] the structure lists:
 *  - its probability 'p'
 *  - an entropy-like measure for its well-definedness 'ent'
 *  - the frequency of each type of pair in 'bp[]'
 *    + 'bp[0]' contains the number of non-compatible sequences
 *    + 'bp[1]' the number of CG pairs, etc.
 */
typedef struct {
   unsigned i;    /**<  @brief  nucleotide position i */
   unsigned j;    /**<  @brief  nucleotide position j */
   float p;       /**< @brief  Probability */
   float ent;     /**< @brief  Pseudo entropy for @f$ p(i,j) = S_i + S_j - p_ij*ln(p_ij) @f$ */
   short bp[8];   /**< @brief  Frequencies of pair_types */
   char comp;     /**< @brief  1 iff pair is in mfe structure */
} pair_info;


/*
* ############################################################
* FINDPATH data structures
* ############################################################
*/

/**
 *  @brief
 */
typedef struct move {
  int i;  /* i,j>0 insert; i,j<0 delete */
  int j;
  int when;  /* 0 if still available, else resulting distance from start */
  int E;
} move_t;

/**
 *  @brief
 */
typedef struct intermediate {
  short *pt;      /**<  @brief  pair table */
  int Sen;        /**<  @brief  saddle energy so far */
  int curr_en;    /**<  @brief  current energy */
  move_t *moves;  /**<  @brief  remaining moves to target */
} intermediate_t;

/**
 *  @brief
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
 *  @brief contributions to p_u
 */
typedef struct pu_contrib {
  double **H; /**<  @brief  hairpin loops */
  double **I; /**<  @brief  interior loops */
  double **M; /**<  @brief  multi loops */
  double **E; /**<  @brief  exterior loop */
  int length; /**<  @brief  length of the input sequence */
  int w;      /**<  @brief  longest unpaired region */
} pu_contrib;

/**
 *  @brief
 */
typedef struct interact {
  double *Pi;       /**<  @brief  probabilities of interaction */
  double *Gi;       /**<  @brief  free energies of interaction */
  double Gikjl;     /**<  @brief  full free energy for interaction between [k,i] k<i
                                  in longer seq and [j,l] j<l in shorter seq */
  double Gikjl_wo;  /**<  @brief  Gikjl without contributions for prob_unpaired */
  int i;            /**<  @brief  k<i in longer seq */
  int k;            /**<  @brief  k<i in longer seq */
  int j;            /**<  @brief  j<l in shorter seq */
  int l;            /**<  @brief  j<l in shorter seq */
  int length;       /**<  @brief  length of longer sequence */
} interact;

/**
 *  @brief  Collection of all free_energy of beeing unpaired values for output
 */
typedef struct pu_out {
  int len;            /**<  @brief  sequence length */
  int u_vals;         /**<  @brief  number of different -u values */
  int contribs;       /**<  @brief  [-c "SHIME"] */
  char **header;      /**<  @brief  header line */
  double **u_values;  /**<  @brief  (the -u values * [-c "SHIME"]) * seq len */
} pu_out;

/**
 *  @brief  constraints for cofolding
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
 *  @brief
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
 *  @brief
 */
typedef struct node {
  int k;
  int energy;
  struct node *next;
} folden;

/**
 *  @brief
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
 *  @brief
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
 *  @addtogroup   basic_data_structures  Basic Data Structures for Structure Prediction and Evaluation
 *  @{
 *
 *  @brief  This module provides interfaces that deal with the most basic data structures used
 *          in structure predicting and energy evaluating function of the RNAlib.
 *
 *          Throughout the RNAlib, a data structure, the #vrna_fold_compound, is used to group
 *          information and data that is required for structure prediction and energy evaluation.
 *          Here, you'll find interface functions to create, modify, and delete #vrna_fold_compound
 *          data structures.
 */

/**
 *  @brief  An enumerator that is used to specify the type of a #vrna_fold_compound
 */
typedef enum {
  VRNA_VC_TYPE_SINGLE,    /**< Type is suitable for single, and hybridizing sequences */
  VRNA_VC_TYPE_ALIGNMENT  /**< Type is suitable for sequence alignments (consensus structure prediction) */
} vrna_vc_t;


/**
 *  @brief  The most basic data structure required by many functions throughout the RNAlib
 *
 *  @note   Please read the documentation of this data structure carefully! Some attributes are only available for
 *  specific types this data structure can adopt.
 *
 *  @warning  Reading/Writing from/to attributes that are not within the scope of the current type usually result
 *  in undefined behavior!
 *
 *  @see  #vrna_fold_compound.type, vrna_get_fold_compound(), vrna_get_fold_compound_ali(), vrna_free_fold_compound(),
 *        #VRNA_VC_TYPE_SINGLE, #VRNA_VC_TYPE_ALIGNMENT
 */
struct vrna_fc{

  /**
      @name Common data fields
      @{
   */
  vrna_vc_t         type;           /**<  @brief  The type of the #vrna_fold_compound.
                                      @details Currently possible values are #VRNA_VC_TYPE_SINGLE, and #VRNA_VC_TYPE_ALIGNMENT
                                      @warning Do not edit this attribute, it will be automagically set by
                                            the corresponding get() methods for the #vrna_fold_compound.
                                            The value specified in this attribute dictates the set of other
                                            attributes to use within this data structure.
                                    */
  unsigned int      length;         /**<  @brief  The length of the sequence (or sequence alignment) */
  int               cutpoint;       /**<  @brief  The position of the (cofold) cutpoint within the provided sequence.
                                      If there is no cutpoint, this field will be set to -1
                                    */

  vrna_hc_t         *hc;            /**<  @brief  The hard constraints data structure used for structure prediction */

  vrna_mx_mfe_t     *matrices;      /**<  @brief  The MFE DP matrices */
  vrna_mx_pf_t      *exp_matrices;  /**<  @brief  The PF DP matrices  */

  vrna_param_t      *params;        /**<  @brief  The precomputed free energy contributions for each type of loop */
  vrna_exp_param_t  *exp_params;    /**<  @brief  The precomputed free energy contributions as Boltzmann factors  */

  int               *iindx;         /**<  @brief  DP matrix accessor  */
  int               *jindx;         /**<  @brief  DP matrix accessor  */

  /**
      @}
   */

#if __STDC_VERSION__ >= 201112L
    /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif

  /**
      @name Data fields available for single/hybrid structure prediction
      @{
   */
      char  *sequence;              /**<  @brief  The input sequence string
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_SINGLE @endverbatim
                                    */
      short *sequence_encoding;     /**<  @brief  Numerical encoding of the input sequence
                                          @see    vrna_sequence_encode()
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_SINGLE @endverbatim
                                    */
      short *sequence_encoding2;
      char  *ptype;                 /**<  @brief  Pair type array
                                     
                                          Contains the numerical encoding of the pair type for each pair (i,j) used
                                          in MFE, Partition function and Evaluation computations.
                                          @note This array is always indexed via jindx, in contrast to previously
                                          different indexing between mfe and pf variants!
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_SINGLE @endverbatim
                                          @see    vrna_get_indx(), vrna_get_ptypes()
                                    */
      char  *ptype_pf_compat;       /**<  @brief  ptype array indexed via iindx
                                          @deprecated  This attribute will vanish in the future!
                                          It's meant for backward compatibility only!
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_SINGLE @endverbatim
                                    */
      vrna_sc_t *sc;                /**<  @brief  The soft constraints for usage in structure prediction and evaluation
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_SINGLE @endverbatim
                                    */

  /**
      @}
   */

#if __STDC_VERSION__ >= 201112L
    /* C11 support for unnamed unions/structs */
    };
    struct {
#endif

  /**
      @name Data fields for consensus structure prediction
      @{
   */
      char  **sequences;            /**<  @brief  The aligned sequences
                                          @note   The end of the alignment is indicated by a NULL pointer in the second dimension
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_ALIGNMENT @endverbatim
                                    */
      unsigned int    n_seq;        /**<  @brief  The number of sequences in the alignment
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_ALIGNMENT @endverbatim
                                    */
      char            *cons_seq;    /**<  @brief  The consensus sequence of the aligned sequences
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_ALIGNMENT @endverbatim
                                    */
      short           *S_cons;      /**<  @brief  Numerical encoding of the consensus sequence
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_ALIGNMENT @endverbatim
                                    */
      short           **S;          /**<  @brief  Numerical encoding of the sequences in the alignment
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_ALIGNMENT @endverbatim
                                    */
      short           **S5;         /**<  @brief    S5[s][i] holds next base 5' of i in sequence s
                                          @warning  Only available if @verbatim type==VRNA_VC_TYPE_ALIGNMENT @endverbatim
                                    */
      short           **S3;         /**<  @brief    Sl[s][i] holds next base 3' of i in sequence s
                                          @warning  Only available if @verbatim type==VRNA_VC_TYPE_ALIGNMENT @endverbatim
                                    */
      char            **Ss;
      unsigned short  **a2s;
      int             *pscore;      /**<  @brief  Precomputed array of pair types expressed as pairing scores
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_ALIGNMENT @endverbatim
                                    */
      vrna_sc_t       **scs;        /**<  @brief  A set of soft constraints (for each sequence in the alignment)
                                          @warning   Only available if @verbatim type==VRNA_VC_TYPE_ALIGNMENT @endverbatim
                                    */
      int             oldAliEn;

  /**
      @}
   */
#if __STDC_VERSION__ >= 201112L
    };
  };
#endif

  /**
   *  @name Additional data fields for Distance Class Partitioning
   *
   *  These data fields are typically populated with meaningful data only if used in the context of Distance Class Partitioning
   *  @{
   */
  unsigned int    maxD1;          /**<  @brief  Maximum allowed base pair distance to first reference */
  unsigned int    maxD2;          /**<  @brief  Maximum allowed base pair distance to second reference */
  short           *reference_pt1; /**<  @brief  A pairtable of the first reference structure */
  short           *reference_pt2; /**<  @brief  A pairtable of the second reference structure */

  unsigned int    *referenceBPs1; /**<  @brief  Matrix containing number of basepairs of reference structure1 in interval [i,j] */
  unsigned int    *referenceBPs2; /**<  @brief  Matrix containing number of basepairs of reference structure2 in interval [i,j] */
  unsigned int    *bpdist;        /**<  @brief  Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j] */

  unsigned int    *mm1;           /**<  @brief  Maximum matching matrix, reference struct 1 disallowed */
  unsigned int    *mm2;           /**<  @brief  Maximum matching matrix, reference struct 2 disallowed */
  
  /**
      @}
   */

  /**
   *  @name Additional data fields for local folding
   *
   *  These data fields are typically populated with meaningful data only if used in the context of local folding
   *  @{
   */
  int             window_size;    /**<  @brief  window size for local folding sliding window approach */
  char            **ptype_local;  /**<  @brief  Pair type array (for local folding) */
  /**
      @}
   */

};


/* the definitions below should be used for functions that return/receive/destroy fold compound data structures */

/**
 *  @brief  Option flag to specify requirement of Minimum Free Energy (MFE) DP matrices
 *          and corresponding set of energy parameters
 *
 *  @see vrna_get_fold_compound(), vrna_get_fold_compound_ali(), #VRNA_OPTION_EVAL_ONLY
 */
#define VRNA_OPTION_MFE             1

/**
 *  @brief  Option flag to specify requirement of Partition Function (PF) DP matrices
 *          and corresponding set of Boltzmann factors
 *
 *  @see vrna_get_fold_compound(), vrna_get_fold_compound_ali(), #VRNA_OPTION_EVAL_ONLY
 */
#define VRNA_OPTION_PF              2

#define VRNA_OPTION_HYBRID          4

/**
 *  @brief  Option flag to specify that neither MFE, nor PF DP matrices are required
 *
 *  Use this flag in conjuntion with #VRNA_OPTION_MFE, and #VRNA_OPTION_PF to save
 *  memory for a #vrna_fold_compound obtained from vrna_get_fold_compound(), or vrna_get_fold_compound_ali()
 *  in cases where only energy evaluation but no structure prediction is required.
 *
 *  @see vrna_get_fold_compound(), vrna_get_fold_compound_ali(), vrna_eval_structure()
 */
#define VRNA_OPTION_EVAL_ONLY       8

#define VRNA_OPTION_WINDOW          16


/**
 *  @brief  Retrieve a #vrna_fold_compound data structure for single sequences and hybridizing sequences
 *
 *  This function provides an easy interface to obtain a prefilled #vrna_fold_compound by passing a single
 *  sequence, or two contatenated sequences as input. For the latter, sequences need to be seperated by
 *  an '&' character like this: @verbatim char *sequence = "GGGG&CCCC"; @endverbatim
 *
 *  The optional parameter 'md_p' can be used to specify the model details for successive computations
 *  based on the content of the generated #vrna_fold_compound.
 *  The third parameter 'options' is used to specify the DP matrix requirements and the corresponding set
 *  of energy parameters. Use the macros:
 *
 *  - #VRNA_OPTION_MFE
 *  - #VRNA_OPTION_PF
 *  - #VRNA_OPTION_WINDOW
 *  - #VRNA_OPTION_EVAL_ONLY
 *
 *  to specify the required type of computations that will be performed with the #vrna_fold_compound.
 *
 *  @note The sequence string must be uppercase, and should contain only RNA (resp. DNA) alphabet depending
 *        on what energy parameter set is used
 *
 *  @see  vrna_get_fold_compound_ali(), #vrna_md_t, #VRNA_OPTION_MFE, #VRNA_OPTION_PF,
 *        #VRNA_OPTION_EVAL_ONLY, #VRNA_OPTION_WINDOW
 *
 *  @param    sequence    A single sequence, or two concatenated sequences seperated by an '&' character
 *  @param    md_p        An optional set of model details
 *  @param    options     The options for DP matrices memory allocation
 *  @return               A prefilled vrna_fold_compound that can be readily used for computations
 */
vrna_fold_compound *vrna_get_fold_compound( const char *sequence,
                                            vrna_md_t *md_p,
                                            unsigned int options);

/**
 *  @brief  Retrieve a #vrna_fold_compound data structure for sequence alignments
 *
 *  This function provides an easy interface to obtain a prefilled #vrna_fold_compound by passing an
 *  alignment of sequences.
 *
 *  The optional parameter 'md_p' can be used to specify the model details for successive computations
 *  based on the content of the generated #vrna_fold_compound.
 *  The third parameter 'options' is used to specify the DP matrix requirements and the corresponding set
 *  of energy parameters. Use the macros:
 *
 *  - #VRNA_OPTION_MFE
 *  - #VRNA_OPTION_PF
 *  - #VRNA_OPTION_EVAL_ONLY
 *
 *  to specify the required type of computations that will be performed with the #vrna_fold_compound.
 *
 *  @note The sequence strings must be uppercase, and should contain only RNA (resp. DNA) alphabet including
 *        gap characters depending on what energy parameter set is used.
 *
 *  @see  vrna_get_fold_compound(), #vrna_md_t, #VRNA_OPTION_MFE, #VRNA_OPTION_PF, #VRNA_OPTION_EVAL_ONLY,
 *        read_clustal()
 *
 *  @param    sequences   A sequence alignment including 'gap' characters
 *  @param    md_p        An optional set of model details
 *  @param    options     The options for DP matrices memory allocation
 *  @return               A prefilled vrna_fold_compound that can be readily used for computations
 */
vrna_fold_compound *vrna_get_fold_compound_ali( const char **sequences,
                                                vrna_md_t *md_p,
                                                unsigned int options);

vrna_fold_compound *vrna_get_fold_compound_2D(const char *sequence,
                                              const char *s1,
                                              const char *s2,
                                              vrna_md_t *md_p,
                                              unsigned int options);

/**
 *  @brief  Free memory occupied by a #vrna_fold_compound
 *
 *  @see vrna_get_fold_compound(), vrna_get_fold_compound_ali(), vrna_free_mfe_matrices(), vrna_free_pf_matrices()
 *
 *  @param  vc  The #vrna_fold_compound that is to be erased from memory
 */
void vrna_free_fold_compound(vrna_fold_compound *vc);

/**
 *  @}
 */

#endif
