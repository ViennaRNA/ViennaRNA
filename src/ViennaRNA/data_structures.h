#ifndef VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H
#define VIENNA_RNA_PACKAGE_DATA_STRUCTURES_H

/**
 *  @file     data_structures.h
 *  @ingroup  data_structures
 *  @brief    Various data structures and pre-processor macros
 */

/**
 *  @addtogroup   data_structures
 *  @brief All datastructures and typedefs shared among the Vienna RNA Package can be found here
 *
 *  @{
 *  @ingroup data_structures
 */

/* below are several convenience typedef's we use throughout the ViennaRNA library */

/** @brief Typename for the fold_compound data structure #vrna_fc_s
 *  @ingroup fold_compound
 */
typedef struct vrna_fc_s vrna_fold_compound_t;

/** @brief Typename for the base pair repesenting data structure #vrna_basepair_s */
typedef struct vrna_basepair_s vrna_basepair_t;

/** @brief Typename for the base pair list repesenting data structure #vrna_plist_s */
typedef struct vrna_plist_s vrna_plist_t;

/** @brief Typename for the base pair stack repesenting data structure #vrna_bp_stack_s */
typedef struct vrna_bp_stack_s vrna_bp_stack_t;

/** @brief Typename for data structure #vrna_cpair_s */
typedef struct vrna_cpair_s vrna_cpair_t;

/** @brief Typename for stack of partial structures #vrna_sect_s */
typedef struct vrna_sect_s vrna_sect_t;

typedef struct vrna_data_linear_s vrna_data_lin_t;

typedef struct vrna_color_s vrna_color_t;

/** @brief Typename for floating point number in partition function computations */
#ifdef  USE_FLOAT_PF
typedef float FLT_OR_DBL;
#else
typedef double FLT_OR_DBL;
#endif

/**
 *  @brief Callback to free memory allocated for auxiliary user-provided data
 *
 *  @ingroup fold_compound
 *  This type of user-implemented function usually deletes auxiliary data structures.
 *  The user must take care to free all the memory occupied by the data structure passed.
 *
 *  @param data    The data that needs to be free'd
 */
typedef void (vrna_callback_free_auxdata)(void *data);

/**
 *  @brief Callback to perform specific user-defined actions before, or after recursive computations
 *
 *  @ingroup fold_compound
 *  @see #VRNA_STATUS_MFE_PRE, #VRNA_STATUS_MFE_POST, #VRNA_STATUS_PF_PRE, #VRNA_STATUS_PF_POST
 *  @param status   The status indicator
 *  @param data     The data structure that was assigned with vrna_fold_compound_add_auxdata()
 *  @param status   The status indicator
 */
typedef void (vrna_callback_recursion_status)(unsigned char status,
                                              void          *data);

/**
 *  @brief  Status message indicating that MFE computations are about to begin
 *
 *  @ingroup fold_compound
 *  @see  #vrna_fold_compound_t.stat_cb, vrna_callback_recursion_status(), vrna_mfe(), vrna_fold(), vrna_circfold(),
 *        vrna_alifold(), vrna_circalifold(), vrna_cofold()
 */
#define VRNA_STATUS_MFE_PRE     (unsigned char)1

/**
 *  @brief  Status message indicating that MFE computations are finished
 *
 *  @ingroup fold_compound
 *  @see  #vrna_fold_compound_t.stat_cb, vrna_callback_recursion_status(), vrna_mfe(), vrna_fold(), vrna_circfold(),
 *        vrna_alifold(), vrna_circalifold(), vrna_cofold()
 */
#define VRNA_STATUS_MFE_POST    (unsigned char)2

/**
 *  @brief  Status message indicating that Partition function computations are about to begin
 *
 *  @ingroup fold_compound
 *  @see  #vrna_fold_compound_t.stat_cb, vrna_callback_recursion_status(), vrna_pf()
 */
#define VRNA_STATUS_PF_PRE      (unsigned char)3

/**
 *  @brief  Status message indicating that Partition function computations are finished
 *
 *  @ingroup fold_compound
 *  @see  #vrna_fold_compound_t.stat_cb, vrna_callback_recursion_status(), vrna_pf()
 */
#define VRNA_STATUS_PF_POST     (unsigned char)4


#define VRNA_PLIST_TYPE_BASEPAIR      0
#define VRNA_PLIST_TYPE_GQUAD         1
#define VRNA_PLIST_TYPE_H_MOTIF       2
#define VRNA_PLIST_TYPE_I_MOTIF       3
#define VRNA_PLIST_TYPE_UD_MOTIF      4


/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT


#ifdef VRNA_BACKWARD_COMPAT

/* the following typedefs are for backward compatibility only */

/**
 *  @brief Old typename of #vrna_basepair_s
 *  @deprecated Use #vrna_basepair_t instead!
 */
typedef struct vrna_basepair_s PAIR;

/**
 *  @brief Old typename of #vrna_plist_s
 *  @deprecated Use #vrna_plist_t instead!
 */
typedef struct vrna_plist_s plist;

/**
 *  @brief Old typename of #vrna_cpair_s
 *  @deprecated Use #vrna_cpair_t instead!
 */
typedef struct vrna_cpair_s cpair;

/**
 *  @brief Old typename of #vrna_sect_s
 *  @deprecated Use #vrna_sect_t instead!
 */
typedef struct vrna_sect_s sect;

/**
 *  @brief Old typename of #vrna_bp_stack_s
 *  @deprecated Use #vrna_bp_stack_t instead!
 */
typedef struct vrna_bp_stack_s bondT;

#endif

#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/dp_matrices.h>
#include <ViennaRNA/constraints.h>
#include <ViennaRNA/grammar.h>
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"

/*
 * ############################################################
 * Here are the type definitions of various datastructures
 * shared among the Vienna RNA Package
 * ############################################################
 */

/**
 *  @brief  Base pair data structure used in subopt.c
 */
struct vrna_basepair_s {
  int i;
  int j;
};

/**
 *  @brief this datastructure is used as input parameter in functions of PS_dot.h and others
 */
struct vrna_plist_s {
  int   i;
  int   j;
  float p;
  int   type;
};

/**
 *  @brief this datastructure is used as input parameter in functions of PS_dot.c
 */
struct vrna_cpair_s {
  int   i, j, mfe;
  float p, hue, sat;
};

struct vrna_color_s {
  float hue;
  float sat;
  float bri;
};

struct vrna_data_linear_s {
  unsigned int  position;
  float         value;
  vrna_color_t  color;
};


/**
 *  @brief  Stack of partial structures for backtracking
 */
struct vrna_sect_s {
  int i;
  int j;
  int ml;
};

/**
 *  @brief  Base pair stack element
 */
struct vrna_bp_stack_s {
  unsigned int  i;
  unsigned int  j;
};


/*
 * ############################################################
 * RNAup data structures
 * ############################################################
 */

/**
 *  @brief contributions to p_u
 */
typedef struct pu_contrib {
  double  **H;    /**<  @brief  hairpin loops */
  double  **I;    /**<  @brief  interior loops */
  double  **M;    /**<  @brief  multi loops */
  double  **E;    /**<  @brief  exterior loop */
  int     length; /**<  @brief  length of the input sequence */
  int     w;      /**<  @brief  longest unpaired region */
} pu_contrib;

/**
 *  @brief  interaction data structure for RNAup
 */
typedef struct interact {
  double  *Pi;      /**<  @brief  probabilities of interaction */
  double  *Gi;      /**<  @brief  free energies of interaction */
  double  Gikjl;    /**<  @brief  full free energy for interaction between [k,i] k<i
                     *            in longer seq and [j,l] j<l in shorter seq */
  double  Gikjl_wo; /**<  @brief  Gikjl without contributions for prob_unpaired */
  int     i;        /**<  @brief  k<i in longer seq */
  int     k;        /**<  @brief  k<i in longer seq */
  int     j;        /**<  @brief  j<l in shorter seq */
  int     l;        /**<  @brief  j<l in shorter seq */
  int     length;   /**<  @brief  length of longer sequence */
} interact;

/**
 *  @brief  Collection of all free_energy of beeing unpaired values for output
 */
typedef struct pu_out {
  int     len;        /**<  @brief  sequence length */
  int     u_vals;     /**<  @brief  number of different -u values */
  int     contribs;   /**<  @brief  [-c "SHIME"] */
  char    **header;   /**<  @brief  header line */
  double  **u_values; /**<  @brief  (the -u values * [-c "SHIME"]) * seq len */
} pu_out;

/**
 *  @brief  constraints for cofolding
 */
typedef struct constrain {
  int   *indx;
  char  *ptype;
} constrain;

/*
 * ############################################################
 * RNAduplex data structures
 * ############################################################
 */

/**
 *  @brief  Data structure for RNAduplex
 */
typedef struct {
  int     i;
  int     j;
  int     end;
  char    *structure;
  double  energy;
  double  energy_backtrack;
  double  opening_backtrack_x;
  double  opening_backtrack_y;
  int     offset;
  double  dG1;
  double  dG2;
  double  ddG;
  int     tb;
  int     te;
  int     qb;
  int     qe;
} duplexT;

/*
 * ############################################################
 * RNAsnoop data structures
 * ############################################################
 */

/**
 *  @brief  Data structure for RNAsnoop (fold energy list)
 */
typedef struct node {
  int         k;
  int         energy;
  struct node *next;
} folden;

/**
 *  @brief  Data structure for RNAsnoop
 */
typedef struct {
  int   i;
  int   j;
  int   u;
  char  *structure;
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
 *  @brief  Data structure used in RNApkplex
 */
typedef struct dupVar {
  int     i;
  int     j;
  int     end;
  char    *pk_helix;
  char    *structure;
  double  energy;
  int     offset;
  double  dG1;
  double  dG2;
  double  ddG;
  int     tb;
  int     te;
  int     qb;
  int     qe;
  int     inactive;
  int     processed;
} dupVar;

/**
 *  @brief  Dummy symbol to check whether the library was build using C11/C++11 features
 *
 *  By default, several data structures of our new v3.0 API use C11/C++11 features, such
 *  as unnamed unions, unnamed structs. However, these features can be deactivated at
 *  compile time to allow building the library and executables with compilers that do not
 *  support these features.
 *
 *  Now, the problem arises that once our static library is compiled and a third-party
 *  application is supposed to link against it, it needs to know, at compile time, how to
 *  correctly address particular data structures. This is usually implicitely taken care of
 *  through the API exposed in our header files. Unfortunately, we had some preprocessor directives
 *  in our header files that changed the API depending on the capabilities of the compiler
 *  the third-party application is build with. This in turn prohibited the use of an RNAlib
 *  compiled without C11/C++11 support in a program that compiles/links with enabled C11/C++11
 *  support and vice-versa.
 *
 *  Therefore, we introduce this dummy symbol which can be used to check, whether the
 *  static library was build with C11/C++11 features.
 *
 *  @note If the symbol is present, the library was build with enabled C11/C++11 features support
 *  and no action is required. However, if the symbol is missing in RNAlib >= 2.2.9, programs
 *  that link to RNAlib must define a pre-processor identifier @em VRNA_DISABLE_C11_FEATURES before
 *  including any ViennaRNA Package header file, for instance by adding a @em CPPFLAG
 *  @code
 * CPPFLAGS+=-DVRNA_DISABLE_C11_FEATURES
 *  @endcode
 *
 *  @since v2.2.9
 */
#ifndef VRNA_DISABLE_C11_FEATURES
void vrna_C11_features(void);


#endif

/**
 * @}
 */


/*
 * ############################################################
 * VRNA fold compound related functions
 * ############################################################
 */

/**
 *  @addtogroup   fold_compound   The Fold Compound
 *  @{
 *
 *  @brief  This module provides interfaces that deal with the most basic data structure used
 *          in structure predicting and energy evaluating function of the RNAlib.
 *
 *          Throughout the entire RNAlib, the #vrna_fold_compound_t, is used to group
 *          information and data that is required for structure prediction and energy evaluation.
 *          Here, you'll find interface functions to create, modify, and delete #vrna_fold_compound_t
 *          data structures.
 */

/**
 *  @brief  An enumerator that is used to specify the type of a #vrna_fold_compound_t
 */
typedef enum {
  VRNA_FC_TYPE_SINGLE,      /**< Type is suitable for single, and hybridizing sequences */
  VRNA_FC_TYPE_COMPARATIVE  /**< Type is suitable for sequence alignments (consensus structure prediction) */
} vrna_fc_type_e;


/**
 *  @brief  The most basic data structure required by many functions throughout the RNAlib
 *
 *  @note   Please read the documentation of this data structure carefully! Some attributes are only available for
 *  specific types this data structure can adopt.
 *
 *  @warning  Reading/Writing from/to attributes that are not within the scope of the current type usually result
 *  in undefined behavior!
 *
 *  @see  #vrna_fold_compound_t.type, vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_fold_compound_free(),
 *        #VRNA_FC_TYPE_SINGLE, #VRNA_FC_TYPE_COMPARATIVE
 */
struct vrna_fc_s {
  /**
   *  @name Common data fields
   *  @{
   */
  vrna_fc_type_e type;              /**<  @brief  The type of the #vrna_fold_compound_t.
                                     * @details Currently possible values are #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE
                                     * @warning Do not edit this attribute, it will be automagically set by
                                     *      the corresponding get() methods for the #vrna_fold_compound_t.
                                     *      The value specified in this attribute dictates the set of other
                                     *      attributes to use within this data structure.
                                     */
  unsigned int      length;         /**<  @brief  The length of the sequence (or sequence alignment) */
  int               cutpoint;       /**<  @brief  The position of the (cofold) cutpoint within the provided sequence.
                                     * If there is no cutpoint, this field will be set to -1
                                     */

  unsigned int      *strand_number; /**<  @brief  The strand number a particular nucleotide is associated with */

  vrna_hc_t         *hc;            /**<  @brief  The hard constraints data structure used for structure prediction */

  vrna_mx_mfe_t     *matrices;      /**<  @brief  The MFE DP matrices */
  vrna_mx_pf_t      *exp_matrices;  /**<  @brief  The PF DP matrices  */

  vrna_param_t      *params;        /**<  @brief  The precomputed free energy contributions for each type of loop */
  vrna_exp_param_t  *exp_params;    /**<  @brief  The precomputed free energy contributions as Boltzmann factors  */

  int               *iindx;         /**<  @brief  DP matrix accessor  */
  int               *jindx;         /**<  @brief  DP matrix accessor  */

  /**
   *  @}
   *
   *  @name User-defined data fields
   *  @{
   */
  vrna_callback_recursion_status  *stat_cb;       /**<  @brief  Recursion status callback (usually called just before, and
                                                   *            after recursive computations in the library
                                                   *    @see    vrna_callback_recursion_status(), vrna_fold_compound_add_callback()
                                                   */

  void                            *auxdata;       /**<  @brief  A pointer to auxiliary, user-defined data
                                                   *    @see vrna_fold_compound_add_auxdata(), #vrna_fold_compound_t.free_auxdata
                                                   */

  vrna_callback_free_auxdata      *free_auxdata;  /**<  @brief A callback to free auxiliary user data whenever the fold_compound itself is free'd
                                                   *    @see  #vrna_fold_compound_t.auxdata, vrna_callback_free_auxdata()
                                                   */

  /**
   *  @}
   *
   *  @name Secondary Structure Decomposition (grammar) related data fields
   *  @{
   */

  /* data structure to adjust additional structural domains, such as G-quadruplexes */
  vrna_sd_t     *domains_struc;             /**<  @brief  Additional structured domains */

  /* data structure to adjust additional contributions to unpaired stretches, e.g. due to protein binding */
  vrna_ud_t     *domains_up;                /**<  @brief  Additional unstructured domains */

  /* auxiliary (user-defined) extension to the folding grammar */
  vrna_gr_aux_t *aux_grammar;

  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif

  /**
   *  @name Data fields available for single/hybrid structure prediction
   *  @{
   */
  char *sequence;                   /**<  @brief  The input sequence string
                                     *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                     */
  short *sequence_encoding;         /**<  @brief  Numerical encoding of the input sequence
                                     *    @see    vrna_sequence_encode()
                                     *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                     */
  short *sequence_encoding2;
  char *ptype;                      /**<  @brief  Pair type array
                                     *
                                     *    Contains the numerical encoding of the pair type for each pair (i,j) used
                                     *    in MFE, Partition function and Evaluation computations.
                                     *    @note This array is always indexed via jindx, in contrast to previously
                                     *    different indexing between mfe and pf variants!
                                     *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                     *    @see    vrna_idx_col_wise(), vrna_ptypes()
                                     */
  char *ptype_pf_compat;            /**<  @brief  ptype array indexed via iindx
                                     *    @deprecated  This attribute will vanish in the future!
                                     *    It's meant for backward compatibility only!
                                     *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                     */
  vrna_sc_t *sc;                    /**<  @brief  The soft constraints for usage in structure prediction and evaluation
                                     *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                     */

  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /**
   *  @name Data fields for consensus structure prediction
   *  @{
   */
  char            **sequences;      /**<  @brief  The aligned sequences
                                     *    @note   The end of the alignment is indicated by a NULL pointer in the second dimension
                                     *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                     */
  unsigned int    n_seq;            /**<  @brief  The number of sequences in the alignment
                                     *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                     */
  char            *cons_seq;        /**<  @brief  The consensus sequence of the aligned sequences
                                     *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                     */
  short           *S_cons;          /**<  @brief  Numerical encoding of the consensus sequence
                                     *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                     */
  short           **S;              /**<  @brief  Numerical encoding of the sequences in the alignment
                                     *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                     */
  short           **S5;             /**<  @brief    S5[s][i] holds next base 5' of i in sequence s
                                     *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                     */
  short           **S3;             /**<  @brief    Sl[s][i] holds next base 3' of i in sequence s
                                     *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                     */
  char            **Ss;
  unsigned short  **a2s;
  int             *pscore;            /**<  @brief  Precomputed array of pair types expressed as pairing scores
                                       *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                       */
  int             **pscore_local;     /**<  @brief  Precomputed array of pair types expressed as pairing scores
                                       *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                       */
  short           *pscore_pf_compat;  /**<  @brief  Precomputed array of pair types expressed as pairing scores indexed via iindx
                                       *    @deprecated  This attribute will vanish in the future!
                                       *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                       */
  vrna_sc_t       **scs;              /**<  @brief  A set of soft constraints (for each sequence in the alignment)
                                       *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                       */
  int             oldAliEn;

  /**
   *  @}
   */
#ifndef VRNA_DISABLE_C11_FEATURES
};
};
#endif

  /**
   *  @name Additional data fields for Distance Class Partitioning
   *
   *  These data fields are typically populated with meaningful data only if used in the context of Distance Class Partitioning
   *  @{
   */
  unsigned int  maxD1;            /**<  @brief  Maximum allowed base pair distance to first reference */
  unsigned int  maxD2;            /**<  @brief  Maximum allowed base pair distance to second reference */
  short         *reference_pt1;   /**<  @brief  A pairtable of the first reference structure */
  short         *reference_pt2;   /**<  @brief  A pairtable of the second reference structure */

  unsigned int  *referenceBPs1;   /**<  @brief  Matrix containing number of basepairs of reference structure1 in interval [i,j] */
  unsigned int  *referenceBPs2;   /**<  @brief  Matrix containing number of basepairs of reference structure2 in interval [i,j] */
  unsigned int  *bpdist;          /**<  @brief  Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j] */

  unsigned int  *mm1;             /**<  @brief  Maximum matching matrix, reference struct 1 disallowed */
  unsigned int  *mm2;             /**<  @brief  Maximum matching matrix, reference struct 2 disallowed */

  /**
   *  @}
   */

  /**
   *  @name Additional data fields for local folding
   *
   *  These data fields are typically populated with meaningful data only if used in the context of local folding
   *  @{
   */
  int   window_size;              /**<  @brief  window size for local folding sliding window approach */
  char  **ptype_local;            /**<  @brief  Pair type array (for local folding) */
  /**
   *  @}
   */
};


/* the definitions below should be used for functions that return/receive/destroy fold compound data structures */

/**
 *  @brief  Option flag to specify default settings/requirements
 */
#define VRNA_OPTION_DEFAULT         0U

/**
 *  @brief  Option flag to specify requirement of Minimum Free Energy (MFE) DP matrices
 *          and corresponding set of energy parameters
 *
 *  @see vrna_fold_compound(), vrna_fold_compound_comparative(), #VRNA_OPTION_EVAL_ONLY
 */
#define VRNA_OPTION_MFE             1U

/**
 *  @brief  Option flag to specify requirement of Partition Function (PF) DP matrices
 *          and corresponding set of Boltzmann factors
 *
 *  @see vrna_fold_compound(), vrna_fold_compound_comparative(), #VRNA_OPTION_EVAL_ONLY
 */
#define VRNA_OPTION_PF              2U

/**
 *  @brief  Option flag to specify requirement of dimer DP matrices
 */
#define VRNA_OPTION_HYBRID          4U

/**
 *  @brief  Option flag to specify that neither MFE, nor PF DP matrices are required
 *
 *  Use this flag in conjuntion with #VRNA_OPTION_MFE, and #VRNA_OPTION_PF to save
 *  memory for a #vrna_fold_compound_t obtained from vrna_fold_compound(), or vrna_fold_compound_comparative()
 *  in cases where only energy evaluation but no structure prediction is required.
 *
 *  @see vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_eval_structure()
 */
#define VRNA_OPTION_EVAL_ONLY       8U

/**
 *  @brief  Option flag to specify requirement of DP matrices for local folding approaches
 */
#define VRNA_OPTION_WINDOW          16U

/**
 *  @brief  Retrieve a #vrna_fold_compound_t data structure for single sequences and hybridizing sequences
 *
 *  This function provides an easy interface to obtain a prefilled #vrna_fold_compound_t by passing a single
 *  sequence, or two contatenated sequences as input. For the latter, sequences need to be seperated by
 *  an '&' character like this: @verbatim char *sequence = "GGGG&CCCC"; @endverbatim
 *
 *  The optional parameter @p md_p can be used to specify the model details for successive computations
 *  based on the content of the generated #vrna_fold_compound_t. Passing NULL will instruct the function
 *  to use default model details.
 *  The third parameter @p options may be used to specify dynamic programming (DP) matrix requirements.
 *  Use the macros:
 *
 *  - #VRNA_OPTION_MFE
 *  - #VRNA_OPTION_PF
 *  - #VRNA_OPTION_WINDOW
 *  - #VRNA_OPTION_EVAL_ONLY
 *  - #VRNA_OPTION_DEFAULT
 *
 *  to specify the required type of computations that will be performed with the #vrna_fold_compound_t.
 *
 *  If you just need the folding compound serving as a container for your data, you can simply pass
 *  #VRNA_OPTION_DEFAULT to the @p option parameter. This creates a #vrna_fold_compound_t without DP
 *  matrices, thus saving memory. Subsequent calls of any structure prediction function will then take
 *  care of allocating the memory required for the DP matrices.
 *  If you only intend to evaluate structures instead of actually predicting them, you may use the
 *  #VRNA_OPTION_EVAL_ONLY macro. This will seriously speedup the creation of the #vrna_fold_compound_t.
 *
 *  @note The sequence string must be uppercase, and should contain only RNA (resp. DNA) alphabet depending
 *        on what energy parameter set is used
 *
 *  @see  vrna_fold_compound_free(), vrna_fold_compound_comparative(), #vrna_md_t, #VRNA_OPTION_MFE,
 *        #VRNA_OPTION_PF, #VRNA_OPTION_EVAL_ONLY, #VRNA_OPTION_WINDOW
 *
 *  @param    sequence    A single sequence, or two concatenated sequences seperated by an '&' character
 *  @param    md_p        An optional set of model details
 *  @param    options     The options for DP matrices memory allocation
 *  @return               A prefilled vrna_fold_compound_t that can be readily used for computations
 */
vrna_fold_compound_t *
vrna_fold_compound(const char   *sequence,
                   vrna_md_t    *md_p,
                   unsigned int options);


/**
 *  @brief  Retrieve a #vrna_fold_compound_t data structure for sequence alignments
 *
 *  This function provides an easy interface to obtain a prefilled #vrna_fold_compound_t by passing an
 *  alignment of sequences.
 *
 *  The optional parameter @p md_p can be used to specify the model details for successive computations
 *  based on the content of the generated #vrna_fold_compound_t. Passing NULL will instruct the function
 *  to use default model details.
 *  The third parameter @p options may be used to specify dynamic programming (DP) matrix requirements.
 *  Use the macros:
 *
 *  - #VRNA_OPTION_MFE
 *  - #VRNA_OPTION_PF
 *  - #VRNA_OPTION_EVAL_ONLY
 *  - #VRNA_OPTION_DEFAULT
 *
 *  to specify the required type of computations that will be performed with the #vrna_fold_compound_t.
 *
 *  If you just need the folding compound serving as a container for your data, you can simply pass
 *  #VRNA_OPTION_DEFAULT to the @p option parameter. This creates a #vrna_fold_compound_t without DP
 *  matrices, thus saving memory. Subsequent calls of any structure prediction function will then take
 *  care of allocating the memory required for the DP matrices.
 *  If you only intend to evaluate structures instead of actually predicting them, you may use the
 *  #VRNA_OPTION_EVAL_ONLY macro. This will seriously speedup the creation of the #vrna_fold_compound_t.
 *
 *  @note The sequence strings must be uppercase, and should contain only RNA (resp. DNA) alphabet including
 *        gap characters depending on what energy parameter set is used.
 *
 *  @see  vrna_fold_compound_free(), vrna_fold_compound(), #vrna_md_t, #VRNA_OPTION_MFE, #VRNA_OPTION_PF,
 *        #VRNA_OPTION_EVAL_ONLY, read_clustal()
 *
 *  @param    sequences   A sequence alignment including 'gap' characters
 *  @param    md_p        An optional set of model details
 *  @param    options     The options for DP matrices memory allocation
 *  @return               A prefilled vrna_fold_compound_t that can be readily used for computations
 */
vrna_fold_compound_t *
vrna_fold_compound_comparative(const char   **sequences,
                               vrna_md_t    *md_p,
                               unsigned int options);


vrna_fold_compound_t *
vrna_fold_compound_TwoD(const char    *sequence,
                        const char    *s1,
                        const char    *s2,
                        vrna_md_t     *md_p,
                        unsigned int  options);


int
vrna_fold_compound_prepare(vrna_fold_compound_t *vc,
                           unsigned int         options);


/**
 *  @brief  Free memory occupied by a #vrna_fold_compound_t
 *
 *  @see vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_mx_mfe_free(), vrna_mx_pf_free()
 *
 *  @param  vc  The #vrna_fold_compound_t that is to be erased from memory
 */
void
vrna_fold_compound_free(vrna_fold_compound_t *vc);


/**
 *  @brief  Add auxiliary data to the #vrna_fold_compound_t
 *
 *  This function allows one to bind arbitrary data to a #vrna_fold_compound_t which may later on be used
 *  by one of the callback functions, e.g. vrna_callback_recursion_status(). To allow for proper cleanup
 *  of the memory occupied by this auxiliary data, the user may also provide a pointer to a cleanup function
 *  that free's the corresponding memory. This function will be called automatically when the #vrna_fold_compound_t
 *  is free'd with vrna_fold_compound_free().
 *
 *  @note Before attaching the arbitrary data pointer, this function will call the vrna_callback_free_auxdata()
 *        on any pre-existing data that is already attached.
 *
 *  @see vrna_callback_free_auxdata()
 *  @param  vc    The fold_compound the arbitrary data pointer should be associated with
 *  @param  data  A pointer to an arbitrary data structure
 *  @param  f     A pointer to function that free's memory occupied by the arbitrary data (May be NULL)
 */
void vrna_fold_compound_add_auxdata(vrna_fold_compound_t        *vc,
                                    void                        *data,
                                    vrna_callback_free_auxdata  *f);


/**
 *  @brief  Add a recursion status callback to the #vrna_fold_compound_t
 *
 *  Binding a recursion status callback function to a #vrna_fold_compound_t allows one to perform
 *  arbitrary operations just before, or after an actual recursive computations, e.g. MFE prediction,
 *  is performed by the RNAlib. The callback function will be provided with a pointer to its
 *  #vrna_fold_compound_t, and a status message. Hence, it has complete access to all variables that
 *  incluence the recursive computations.
 *
 *  @see  vrna_callback_recursion_status(), #vrna_fold_compound_t,
 *        #VRNA_STATUS_MFE_PRE, #VRNA_STATUS_MFE_POST, #VRNA_STATUS_PF_PRE, #VRNA_STATUS_PF_POST
 *
 *  @param  vc    The fold_compound the callback function should be attached to
 *  @param  f     The pointer to the recursion status callback function
 */
void vrna_fold_compound_add_callback(vrna_fold_compound_t           *vc,
                                     vrna_callback_recursion_status *f);


/**
 *  @}
 */

#endif
