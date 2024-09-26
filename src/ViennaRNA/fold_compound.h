#ifndef VIENNA_RNA_PACKAGE_FOLD_COMPOUND_H
#define VIENNA_RNA_PACKAGE_FOLD_COMPOUND_H

#ifdef VRNA_WARN_DEPRECATED
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

/**
 *  @file     fold_compound.h
 *  @ingroup  fold_compound
 *  @brief    The Basic Fold Compound API
 */

/**
 *  @addtogroup   fold_compound   The Fold Compound
 *  @{
 */

/**
 *  @brief Typename for the fold_compound data structure #vrna_fc_s
 */
typedef struct vrna_fc_s vrna_fold_compound_t;

/**
 *  @brief Callback to free memory allocated for auxiliary user-provided data
 *
 *  This type of user-implemented function usually deletes auxiliary data structures.
 *  The user must take care to free all the memory occupied by the data structure passed.
 *
 *  @callback
 *  @parblock
 *  This callback is supposed to free memory occupied by an auxiliary data structure.
 *  It will be called when the #vrna_fold_compound_t is erased from memory through a
 *  call to vrna_fold_compound_free() and will be passed the address of memory previously
 *  bound to the #vrna_fold_compound_t via vrna_fold_compound_add_auxdata().
 *  @endparblock
 *
 *  @see vrna_fold_compound_add_auxdata(), vrna_fold_compound_free(), vrna_fold_compound_add_callback()
 *
 *  @param data    The data that needs to be free'd
 */
typedef void (*vrna_auxdata_free_f)(void *data);


/**
 *  @brief  Callback to prepare user-provided data based on some event
 *
 *  @param  fc          The fold compound the event was emitted from
 *  @param  data        The user-defined data that may be prepared according to the event
 *  @param  event       The event
 *  @param  event_data  Some additional data corresponding to the event
 *  @return             non-zero if the preparation was successful, 0 otherwise
 */
typedef int (*vrna_auxdata_prepare_f)(vrna_fold_compound_t  *fc,
                                      void                  *data,
                                      unsigned int          event,
                                      void                  *event_data);


/**
 *  @brief Callback to free memory allocated for auxiliary user-provided data
 *
 *  @deprecated   Use vrna_auxdata_free_f(void *data) instead!
 */
DEPRECATED(typedef void (vrna_callback_free_auxdata)(void *data),
           "Use vrna_auxdata_free_f instead!");

/**
 *  @brief Callback to perform specific user-defined actions before, or after recursive computations
 *
 *  @callback
 *  @parblock
 *  This function will be called to notify a third-party implementation about the status of
 *  a currently ongoing recursion. The purpose of this callback mechanism is to provide users
 *  with a simple way to ensure pre- and post conditions for auxiliary mechanisms attached to
 *  our implementations.
 *  @endparblock
 *
 *  @see vrna_fold_compound_add_auxdata(), vrna_fold_compound_add_callback(), vrna_mfe(),
 *       vrna_pf(),
 *       #VRNA_STATUS_MFE_PRE, #VRNA_STATUS_MFE_POST, #VRNA_STATUS_PF_PRE, #VRNA_STATUS_PF_POST
 *
 *  @param fc       The fold compound emitting the status
 *  @param status   The status indicator
 *  @param data     The data structure that was assigned with vrna_fold_compound_add_auxdata()
 */
typedef void (*vrna_recursion_status_f)(vrna_fold_compound_t *fc,
                                        unsigned char status,
                                        void          *data);

DEPRECATED(typedef void (vrna_callback_recursion_status)(unsigned char status,
                                              void          *data),
           "Use vrna_recursion_status_f instead!");


/**
 *  @brief  Status message indicating that MFE computations are about to begin
 *
 *  @see  #vrna_fold_compound_t.stat_cb, vrna_recursion_status_f(), vrna_mfe(), vrna_fold(), vrna_circfold(),
 *        vrna_alifold(), vrna_circalifold(), vrna_cofold()
 */
#define VRNA_STATUS_MFE_PRE     (unsigned char)1

/**
 *  @brief  Status message indicating that MFE computations are finished
 *
 *  @see  #vrna_fold_compound_t.stat_cb, vrna_recursion_status_f(), vrna_mfe(), vrna_fold(), vrna_circfold(),
 *        vrna_alifold(), vrna_circalifold(), vrna_cofold()
 */
#define VRNA_STATUS_MFE_POST    (unsigned char)2

/**
 *  @brief  Status message indicating that Partition function computations are about to begin
 *
 *  @see  #vrna_fold_compound_t.stat_cb, vrna_recursion_status_f(), vrna_pf()
 */
#define VRNA_STATUS_PF_PRE      (unsigned char)3

/**
 *  @brief  Status message indicating that Partition function computations are finished
 *
 *  @see  #vrna_fold_compound_t.stat_cb, vrna_recursion_status_f(), vrna_pf()
 */
#define VRNA_STATUS_PF_POST     (unsigned char)4


#include <ViennaRNA/model.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/sequences/sequence.h>
#include <ViennaRNA/datastructures/dp_matrices.h>
#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/constraints/soft.h>
#include <ViennaRNA/grammar/basic.h>
#include <ViennaRNA/structured_domains.h>
#include <ViennaRNA/unstructured_domains.h>

#ifdef VRNA_WITH_SVM
#include <ViennaRNA/zscore/basic.h>
#endif


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
 *          specific types this data structure can adopt.
 *
 *  @warning  Reading/Writing from/to attributes that are not within the scope of the current type usually result
 *            in undefined behavior!
 *
 *  @see  #vrna_fold_compound_t.type, vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_fold_compound_free(),
 *        #VRNA_FC_TYPE_SINGLE, #VRNA_FC_TYPE_COMPARATIVE
 */
struct vrna_fc_s {
  /**
   *  @name Common data fields
   *  @{
   */
  const vrna_fc_type_e type;        /**<  @brief  The type of the #vrna_fold_compound_t.
                                     * @details Currently possible values are #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE
                                     * @warning Do not edit this attribute, it will be automagically set by
                                     *      the corresponding get() methods for the #vrna_fold_compound_t.
                                     *      The value specified in this attribute dictates the set of other
                                     *      attributes to use within this data structure.
                                     */
  unsigned int      length;         /**<  @brief  The length of the sequence (or sequence alignment) */
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
  DEPRECATED(int cutpoint,
             "Use strand_* members instead");
                                    /**<  @brief  The position of the (cofold) cutpoint within the provided sequence.
                                     * If there is no cutpoint, this field will be set to -1
                                     */
#endif
  unsigned int      *strand_number; /**<  @brief  The strand number a particular nucleotide is associated with */
  unsigned int      *strand_order;  /**<  @brief  The strand order, i.e. permutation of current concatenated sequence */
  unsigned int      *strand_order_uniq; /**<  @brief  The strand order array where identical sequences have the same ID */
  unsigned int      *strand_start;  /**<  @brief  The start position of a particular strand within the current concatenated sequence */
  unsigned int      *strand_end;    /**<  @brief  The end (last) position of a particular strand within the current concatenated sequence */

  unsigned int      strands;        /**<  @brief  Number of interacting strands */
  vrna_seq_t        *nucleotides;   /**<  @brief  Set of nucleotide sequences */
  vrna_msa_t        *alignment;     /**<  @brief  Set of alignments */

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
  vrna_recursion_status_f   stat_cb;       /**<  @brief  Recursion status callback (usually called just before, and
                                            *            after recursive computations in the library
                                            *    @see    vrna_recursion_status_f(), vrna_fold_compound_add_callback()
                                            */

  void                      *auxdata;      /**<  @brief  A pointer to auxiliary, user-defined data
                                            *    @see vrna_fold_compound_add_auxdata(), #vrna_fold_compound_t.free_auxdata
                                            */

  vrna_auxdata_free_f       free_auxdata;  /**<  @brief A callback to free auxiliary user data whenever the fold_compound itself is free'd
                                            *    @see  #vrna_fold_compound_t.auxdata, vrna_auxdata_free_f()
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
  vrna_gr_aux_t aux_grammar;               /**<  @brief  Additional decomposition grammar rules */

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
      short *encoding5;
      short *encoding3;

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
      char          **sequences;        /**<  @brief  The aligned sequences
                                         *    @note   The end of the alignment is indicated by a NULL pointer in the second dimension
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      unsigned int  n_seq;              /**<  @brief  The number of sequences in the alignment
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      char          *cons_seq;          /**<  @brief  The consensus sequence of the aligned sequences
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         *S_cons;            /**<  @brief  Numerical encoding of the consensus sequence
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         **S;                /**<  @brief  Numerical encoding of the sequences in the alignment
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         **S5;               /**<  @brief    S5[s][i] holds next base 5' of i in sequence s
                                         *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         **S3;               /**<  @brief    Sl[s][i] holds next base 3' of i in sequence s
                                         *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
  char          **Ss;
  unsigned int  **a2s;
      int           *pscore;              /**<  @brief  Precomputed array of pair types expressed as pairing scores
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
      int           **pscore_local;       /**<  @brief  Precomputed array of pair types expressed as pairing scores
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
      short         *pscore_pf_compat;    /**<  @brief  Precomputed array of pair types expressed as pairing scores indexed via iindx
                                           *    @deprecated  This attribute will vanish in the future!
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
      vrna_sc_t     **scs;                /**<  @brief  A set of soft constraints (for each sequence in the alignment)
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
  int           oldAliEn;

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
#ifdef VRNA_WITH_SVM
  vrna_zsc_dat_t  zscore_data;    /**<  @brief  Data structure with settings for z-score computations */
#endif

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
#define VRNA_OPTION_MFE             (1 << 0)

/**
 *  @brief  Option flag to specify requirement of Partition Function (PF) DP matrices
 *          and corresponding set of Boltzmann factors
 *
 *  @see vrna_fold_compound(), vrna_fold_compound_comparative(), #VRNA_OPTION_EVAL_ONLY
 */
#define VRNA_OPTION_PF              (1 << 1)

/**
 *  @brief  Option flag to specify requirement of dimer DP matrices
 */
#define VRNA_OPTION_HYBRID          (1 << 2)

/**
 *  @brief  Option flag to specify that neither MFE, nor PF DP matrices are required
 *
 *  Use this flag in conjuntion with #VRNA_OPTION_MFE, and #VRNA_OPTION_PF to save
 *  memory for a #vrna_fold_compound_t obtained from vrna_fold_compound(), or vrna_fold_compound_comparative()
 *  in cases where only energy evaluation but no structure prediction is required.
 *
 *  @see vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_eval_structure()
 */
#define VRNA_OPTION_EVAL_ONLY       (1 << 3)

/**
 *  @brief  Option flag to specify requirement of DP matrices for local folding approaches
 */
#define VRNA_OPTION_WINDOW          (1 << 4)


#define VRNA_OPTION_F5              (1 << 5)
#define VRNA_OPTION_F3              (1 << 6)
#define VRNA_OPTION_WINDOW_F5       (VRNA_OPTION_WINDOW | VRNA_OPTION_F5)
#define VRNA_OPTION_WINDOW_F3       (VRNA_OPTION_WINDOW | VRNA_OPTION_F3)

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
 *
 *  #### Options ####
 *  * #VRNA_OPTION_DEFAULT  - @copybrief #VRNA_OPTION_DEFAULT
 *  * #VRNA_OPTION_MFE      - @copybrief #VRNA_OPTION_MFE
 *  * #VRNA_OPTION_PF       - @copybrief #VRNA_OPTION_PF
 *  * #VRNA_OPTION_WINDOW   - @copybrief #VRNA_OPTION_WINDOW
 *
 *  The above options may be OR-ed together.
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
 *  @see  vrna_fold_compound_free(), vrna_fold_compound_comparative(), #vrna_md_t
 *
 *  @param    sequence    A single sequence, or two concatenated sequences seperated by an '&' character
 *  @param    md_p        An optional set of model details
 *  @param    options     The options for DP matrices memory allocation
 *  @return               A prefilled vrna_fold_compound_t ready to be used for computations (may be @p NULL on error)
 */
vrna_fold_compound_t *
vrna_fold_compound(const char       *sequence,
                   const vrna_md_t  *md_p,
                   unsigned int     options);


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
 *
 *  #### Options ####
 *  * #VRNA_OPTION_DEFAULT  - @copybrief #VRNA_OPTION_DEFAULT
 *  * #VRNA_OPTION_MFE      - @copybrief #VRNA_OPTION_MFE
 *  * #VRNA_OPTION_PF       - @copybrief #VRNA_OPTION_PF
 *  * #VRNA_OPTION_WINDOW   - @copybrief #VRNA_OPTION_WINDOW
 *
 *  The above options may be OR-ed together.
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
 *  @return               A prefilled vrna_fold_compound_t ready to be used for computations (may be @p NULL on error)
 */
vrna_fold_compound_t *
vrna_fold_compound_comparative(const char   **sequences,
                               vrna_md_t    *md_p,
                               unsigned int options);


vrna_fold_compound_t *
vrna_fold_compound_comparative2(const char                **sequences,
                                const char                **names,
                                const unsigned char       *orientation,
                                const unsigned long long  *start,
                                const unsigned long long  *genome_size,
                                vrna_md_t                 *md_p,
                                unsigned int              options);


vrna_fold_compound_t *
vrna_fold_compound_TwoD(const char    *sequence,
                        const char    *s1,
                        const char    *s2,
                        vrna_md_t     *md_p,
                        unsigned int  options);


int
vrna_fold_compound_prepare(vrna_fold_compound_t *fc,
                           unsigned int         options);


/**
 *  @brief  Free memory occupied by a #vrna_fold_compound_t
 *
 *  @see vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_mx_mfe_free(), vrna_mx_pf_free()
 *
 *  @param  fc  The #vrna_fold_compound_t that is to be erased from memory
 */
void
vrna_fold_compound_free(vrna_fold_compound_t *fc);


/**
 *  @brief  Add auxiliary data to the #vrna_fold_compound_t
 *
 *  This function allows one to bind arbitrary data to a #vrna_fold_compound_t which may later on be used
 *  by one of the callback functions, e.g. vrna_recursion_status_f(). To allow for proper cleanup
 *  of the memory occupied by this auxiliary data, the user may also provide a pointer to a cleanup function
 *  that free's the corresponding memory. This function will be called automatically when the #vrna_fold_compound_t
 *  is free'd with vrna_fold_compound_free().
 *
 *  @note Before attaching the arbitrary data pointer, this function will call the vrna_auxdata_free_f()
 *        on any pre-existing data that is already attached.
 *
 *  @see vrna_auxdata_free_f()
 *  @param  fc    The fold_compound the arbitrary data pointer should be associated with
 *  @param  data  A pointer to an arbitrary data structure
 *  @param  f     A pointer to function that free's memory occupied by the arbitrary data (May be NULL)
 */
void
vrna_fold_compound_add_auxdata(vrna_fold_compound_t *fc,
                               void                 *data,
                               vrna_auxdata_free_f  f);


/**
 *  @brief  Add a recursion status callback to the #vrna_fold_compound_t
 *
 *  Binding a recursion status callback function to a #vrna_fold_compound_t allows one to perform
 *  arbitrary operations just before, or after an actual recursive computations, e.g. MFE prediction,
 *  is performed by the RNAlib. The callback function will be provided with a pointer to its
 *  #vrna_fold_compound_t, and a status message. Hence, it has complete access to all variables that
 *  incluence the recursive computations.
 *
 *  @see  vrna_recursion_status_f(), #vrna_fold_compound_t,
 *        #VRNA_STATUS_MFE_PRE, #VRNA_STATUS_MFE_POST, #VRNA_STATUS_PF_PRE, #VRNA_STATUS_PF_POST
 *
 *  @param  fc    The fold_compound the callback function should be attached to
 *  @param  f     The pointer to the recursion status callback function
 */
void
vrna_fold_compound_add_callback(vrna_fold_compound_t    *fc,
                                vrna_recursion_status_f f);


/**
 *  @}
 */

#endif
