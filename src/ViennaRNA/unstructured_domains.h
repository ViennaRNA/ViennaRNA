#ifndef VIENNA_RNA_PACKAGE_UNSTRUCTURED_DOMAIN_H
#define VIENNA_RNA_PACKAGE_UNSTRUCTURED_DOMAIN_H

#ifdef VRNA_WARN_DEPRECATED
# if defined(DEPRECATED)
#   undef DEPRECATED
# endif
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
 *  @file unstructured_domains.h
 *  @ingroup domains_up
 *  @brief    Functions to modify unstructured domains, e.g. to incorporate ligands binding to unpaired stretches
 */

/** @brief Typename for the ligand binding extension data structure #vrna_unstructured_domain_s
 *  @ingroup domains_up
 */
typedef struct vrna_unstructured_domain_s vrna_ud_t;

typedef struct vrna_unstructured_domain_motif_s vrna_ud_motif_t;

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/structures/problist.h>

/**
 *  @brief Callback to retrieve binding free energy of a ligand bound to an unpaired sequence segment
 *
 *  @ingroup domains_up
 *
 *  @callback
 *  @parblock
 *  This function will be called to determine the additional energy contribution of a specific unstructured
 *  domain, e.g. the binding free energy of some ligand.
 *  @endparblock
 *
 *  @param  fc        The current #vrna_fold_compound_t
 *  @param  i         The start of the unstructured domain (5' end)
 *  @param  j         The end of the unstructured domain (3' end)
 *  @param  loop_type The loop context of the unstructured domain
 *  @param  data      Auxiliary data
 *  @return           The auxiliary energy contribution in deka-cal/mol
 */
typedef int (*vrna_ud_f)(vrna_fold_compound_t  *fc,
                                      int                   i,
                                      int                   j,
                                      unsigned int          loop_type,
                                      void                  *data);

DEPRECATED(typedef int (vrna_callback_ud_energy)(vrna_fold_compound_t  *fc,
                                      int                   i,
                                      int                   j,
                                      unsigned int          loop_type,
                                      void                  *data),
           "Use vrna_ud_f instead!");

/**
 *  @brief Callback to retrieve Boltzmann factor of the binding free energy of a ligand bound to an unpaired sequence segment
 *  @ingroup domains_up
 *
 *  @callback
 *  @parblock
 *  This function will be called to determine the additional energy contribution of a specific unstructured
 *  domain, e.g. the binding free energy of some ligand (Partition function variant, i.e. the Boltzmann factors
 *  instead of actual free energies).
 *  @endparblock
 *
 *  @param  fc        The current #vrna_fold_compound_t
 *  @param  i         The start of the unstructured domain (5' end)
 *  @param  j         The end of the unstructured domain (3' end)
 *  @param  loop_type The loop context of the unstructured domain
 *  @param  data      Auxiliary data
 *  @return           The auxiliary energy contribution as Boltzmann factor
 */
typedef FLT_OR_DBL (*vrna_ud_exp_f)(vrna_fold_compound_t *fc,
                                                 int                  i,
                                                 int                  j,
                                                 unsigned int         loop_type,
                                                 void                 *data);

DEPRECATED(typedef FLT_OR_DBL (vrna_callback_ud_exp_energy)(vrna_fold_compound_t *fc,
                                                 int                  i,
                                                 int                  j,
                                                 unsigned int         loop_type,
                                                 void                 *data),
          "Use vrna_ud_exp_f instead!");

/**
 *  @brief Callback for pre-processing the production rule of the ligand binding to unpaired stretches feature
 *
 *  @ingroup domains_up
 *
 *  @callback
 *  @parblock
 *  The production rule for the unstructured domain grammar extension
 *  @endparblock
 */
typedef void (*vrna_ud_production_f)(vrna_fold_compound_t *fc,
                                           void                 *data);

DEPRECATED(typedef void (vrna_callback_ud_production)(vrna_fold_compound_t *fc,
                                           void                 *data),
          "Use vrna_ud_production_f instead!");

/**
 *  @brief Callback for pre-processing the production rule of the ligand binding to unpaired stretches feature (partition function variant)
 *
 *  @ingroup domains_up
 *
 *  @callback
 *  @parblock
 *  The production rule for the unstructured domain grammar extension (Partition function variant)
 *  @endparblock
 */
typedef void (*vrna_ud_exp_production_f)(vrna_fold_compound_t *fc,
                                               void                 *data);

DEPRECATED(typedef void (vrna_callback_ud_exp_production)(vrna_fold_compound_t *fc,
                                               void                 *data),
           "Use vrna_ud_exp_production_f instead!");


/**
 *  @brief Callback to store/add equilibrium probability for a ligand bound to an unpaired sequence segment
 *  @ingroup domains_up
 *
 *  @callback
 *  @parblock
 *  A callback function to store equilibrium probabilities for the unstructured domain feature
 *  @endparblock
 */
typedef void (*vrna_ud_add_probs_f)(vrna_fold_compound_t  *fc,
                                          int                   i,
                                          int                   j,
                                          unsigned int          loop_type,
                                          FLT_OR_DBL            exp_energy,
                                          void                  *data);

DEPRECATED(typedef void (vrna_callback_ud_probs_add)(vrna_fold_compound_t  *fc,
                                          int                   i,
                                          int                   j,
                                          unsigned int          loop_type,
                                          FLT_OR_DBL            exp_energy,
                                          void                  *data),
          "Use vrna_ud_add_probs_f instead!");

/**
 *  @brief Callback to retrieve equilibrium probability for a ligand bound to an unpaired sequence segment
 *  @ingroup domains_up
 *
 *  @callback
 *  @parblock
 *  A callback function to retrieve equilibrium probabilities for the unstructured domain feature
 *  @endparblock
 */
typedef FLT_OR_DBL (*vrna_ud_get_probs_f)(vrna_fold_compound_t  *fc,
                                                int                   i,
                                                int                   j,
                                                unsigned int          loop_type,
                                                int                   motif,
                                                void                  *data);

DEPRECATED(typedef FLT_OR_DBL (vrna_callback_ud_probs_get)(vrna_fold_compound_t  *fc,
                                                int                   i,
                                                int                   j,
                                                unsigned int          loop_type,
                                                int                   motif,
                                                void                  *data),
            "Use vrna_ud_get_probs_f instead!");


/**
 *  @brief Flag to indicate ligand bound to unpiared stretch in the exterior loop
 *  @ingroup domains_up
 */
#define VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP    1U

/**
 *  @brief Flag to indicate ligand bound to unpaired stretch in a hairpin loop
 *  @ingroup domains_up
 */
#define VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP     2U

/**
 *  @brief Flag to indicate ligand bound to unpiared stretch in an internal loop
 *  @ingroup domains_up
 */
#define VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP    4U

/**
 *  @brief Flag to indicate ligand bound to unpiared stretch in a multibranch loop
 *  @ingroup domains_up
 */
#define VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP     8U

/**
 *  @brief Flag to indicate ligand binding without additional unbound nucleotides (motif-only)
 *  @ingroup domains_up
 */
#define VRNA_UNSTRUCTURED_DOMAIN_MOTIF       16U

/**
 *  @brief Flag to indicate ligand bound to unpiared stretch in any loop (convenience macro)
 *  @ingroup domains_up
 */
#define VRNA_UNSTRUCTURED_DOMAIN_ALL_LOOPS   (VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | \
                                              VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP | \
                                              VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP | \
                                              VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP)

/**
 *  @brief  Data structure to store all functionality for ligand binding
 *  @ingroup domains_up
 */
struct vrna_unstructured_domain_s {
  /*
   **********************************
   * Keep track of all motifs added
   **********************************
   */
  int           uniq_motif_count;                   /**<  @brief The unique number of motifs of different lengths */
  unsigned int  *uniq_motif_size;                   /**<  @brief An array storing a unique list of motif lengths */

  int           motif_count;                        /**<  @brief Total number of distinguished motifs */
  char          **motif;                            /**<  @brief Motif sequences */
  char          **motif_name;                       /**<  @brief Motif identifier/name */
  unsigned int  *motif_size;                        /**<  @brief Motif lengths */
  double        *motif_en;                          /**<  @brief Ligand binding free energy contribution */
  unsigned int  *motif_type;                        /**<  @brief Type of motif, i.e. loop type the ligand binds to */

  /*
   **********************************
   * Grammar extension for ligand
   * binding
   **********************************
   */
  vrna_ud_production_f     prod_cb;        /**<  @brief Callback to ligand binding production rule, i.e. create/fill DP free energy matrices
                                                   *    @details This callback will be executed right before the actual secondary structure decompositions,
                                                   *    and, therefore, any implementation must not interleave with the regular DP matrices.
                                                   */
  vrna_ud_exp_production_f exp_prod_cb;    /**<  @brief Callback to ligand binding production rule, i.e. create/fill DP partition function matrices */
  vrna_ud_f         energy_cb;      /**<  @brief Callback to evaluate free energy of ligand binding to a particular unpaired stretch */
  vrna_ud_exp_f     exp_energy_cb;  /**<  @brief Callback to evaluate Boltzmann factor of ligand binding to a particular unpaired stretch */
  void                            *data;          /**<  @brief Auxiliary data structure passed to energy evaluation callbacks */
  vrna_auxdata_free_f      free_data;      /**<  @brief Callback to free auxiliary data structure */
  vrna_ud_add_probs_f      probs_add;      /**<  @brief Callback to store/add outside partition function */
  vrna_ud_get_probs_f      probs_get;      /**<  @brief Callback to retrieve outside partition function */
};


struct vrna_unstructured_domain_motif_s {
  int start;
  int number;
};


/**
 *  @brief Detect unstructured domains in centroid structure
 *
 *  Given a centroid structure and a set of unstructured domains compute
 *  the list of unstructured domain motifs present in the centroid.
 *  Since we do not explicitly annotate unstructured domain motifs in
 *  dot-bracket strings, this function can be used to check for the
 *  presence and location of unstructured domain motifs under the
 *  assumption that the dot-bracket string is the centroid structure
 *  of the equiibrium ensemble.
 *
 *  @see vrna_centroid()
 *
 *  @ingroup domains_up
 *
 *  @param  fc        The fold_compound data structure with pre-computed equilibrium probabilities and model settings
 *  @param  structure The centroid structure in dot-bracket notation
 *  @return           A list of unstructured domain motifs (possibly NULL). The last element terminates the list with
 *                    @p start=0, @p number=-1
 */
vrna_ud_motif_t *
vrna_ud_motifs_centroid(vrna_fold_compound_t  *fc,
                        const char            *structure);


/**
 *  @brief Detect unstructured domains in MEA structure
 *
 *  Given an MEA structure and a set of unstructured domains compute
 *  the list of unstructured domain motifs present in the MEA structure.
 *  Since we do not explicitly annotate unstructured domain motifs in
 *  dot-bracket strings, this function can be used to check for the
 *  presence and location of unstructured domain motifs under the
 *  assumption that the dot-bracket string is the MEA structure
 *  of the equiibrium ensemble.
 *
 *  @see MEA()
 *
 *  @ingroup domains_up
 *
 *  @param  fc                The fold_compound data structure with pre-computed equilibrium probabilities and model settings
 *  @param  structure         The MEA structure in dot-bracket notation
 *  @param  probability_list  The list of probabilities to extract the MEA structure from
 *  @return                   A list of unstructured domain motifs (possibly NULL). The last element terminates the list
 *                            with @p start=0, @p number=-1
 */
vrna_ud_motif_t *
vrna_ud_motifs_MEA(vrna_fold_compound_t *fc,
                   const char           *structure,
                   vrna_ep_t            *probability_list);


/**
 *  @brief Detect unstructured domains in MFE structure
 *
 *  Given an MFE structure and a set of unstructured domains compute
 *  the list of unstructured domain motifs present in the MFE structure.
 *  Since we do not explicitly annotate unstructured domain motifs in
 *  dot-bracket strings, this function can be used to check for the
 *  presence and location of unstructured domain motifs under the
 *  assumption that the dot-bracket string is the MFE structure
 *  of the equiibrium ensemble.
 *
 *  @see vrna_mfe()
 *
 *  @ingroup domains_up
 *
 *  @param  fc        The fold_compound data structure with model settings
 *  @param  structure The MFE structure in dot-bracket notation
 *  @return           A list of unstructured domain motifs (possibly NULL). The last element terminates the list with @p start=0, @p number=-1
 */
vrna_ud_motif_t *
vrna_ud_motifs_MFE(vrna_fold_compound_t *fc,
                   const char           *structure);


/**
 *  @brief  Add an unstructured domain motif, e.g. for ligand binding
 *
 *  This function adds a ligand binding motif and the associated binding free energy
 *  to the #vrna_ud_t attribute of a #vrna_fold_compound_t. The motif data
 *  will then be used in subsequent secondary structure predictions. Multiple calls
 *  to this function with different motifs append all additional data to a list of
 *  ligands, which all will be evaluated. Ligand motif data can be removed from the
 *  #vrna_fold_compound_t again using the vrna_ud_remove() function. The loop
 *  type parameter allows one to limit the ligand binding to particular loop type,
 *  such as the exterior loop, hairpin loops, internal loops, or multibranch loops.
 *
 *  @see  #VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP, #VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
 *        #VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP, #VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
 *        #VRNA_UNSTRUCTURED_DOMAIN_ALL_LOOPS, vrna_ud_remove()
 *
 *  @ingroup domains_up
 *
 *  @param  fc          The #vrna_fold_compound_t data structure the ligand motif should be bound to
 *  @param  motif       The sequence motif the ligand binds to
 *  @param  motif_en    The binding free energy of the ligand in kcal/mol
 *  @param  motif_name  The name/id of the motif (may be @p NULL)
 *  @param  loop_type   The loop type the ligand binds to
 *
 */
void  vrna_ud_add_motif(vrna_fold_compound_t  *fc,
                        const char            *motif,
                        double                motif_en,
                        const char            *motif_name,
                        unsigned int          loop_type);


/**
 *  @brief  Get a list of unique motif sizes that start at a certain position within the sequence
 *
 */
int *vrna_ud_get_motif_size_at(vrna_fold_compound_t *fc,
                               int                  i,
                               unsigned int         loop_type);


int *
vrna_ud_get_motifs_at(vrna_fold_compound_t  *fc,
                      int                   i,
                      unsigned int          loop_type);


vrna_ud_motif_t *
vrna_ud_detect_motifs(vrna_fold_compound_t  *fc,
                      const char            *structure);


/**
 *  @brief Remove ligand binding to unpaired stretches
 *
 *  This function removes all ligand motifs that were bound to a #vrna_fold_compound_t using
 *  the vrna_ud_add_motif() function.
 *
 *  @ingroup domains_up
 *
 *  @param fc The #vrna_fold_compound_t data structure the ligand motif data should be removed from
 */
void  vrna_ud_remove(vrna_fold_compound_t *fc);


/**
 *  @brief  Attach an auxiliary data structure
 *
 *  This function binds an arbitrary, auxiliary data structure for user-implemented ligand binding.
 *  The optional callback @p free_cb will be passed the bound data structure whenever the #vrna_fold_compound_t
 *  is removed from memory to avoid memory leaks.
 *
 *  @see  vrna_ud_set_prod_rule_cb(), vrna_ud_set_exp_prod_rule_cb(),
 *        vrna_ud_remove()
 *
 *  @ingroup domains_up
 *
 *  @param  fc      The #vrna_fold_compound_t data structure the auxiliary data structure should be bound to
 *  @param  data    A pointer to the auxiliary data structure
 *  @param  free_cb A pointer to a callback function that free's memory occupied by @p data
 */
void  vrna_ud_set_data(vrna_fold_compound_t       *fc,
                       void                       *data,
                       vrna_auxdata_free_f free_cb);


/**
 *  @brief Attach production rule callbacks for free energies computations
 *
 *  Use this function to bind a user-implemented grammar extension for unstructured
 *  domains.
 *
 *  The callback @p e_cb needs to evaluate the free energy contribution @f$f(i,j)@f$ of
 *  the unpaired segment @f$[i,j]@f$. It will be executed in each of the regular secondary
 *  structure production rules. Whenever the callback is passed the #VRNA_UNSTRUCTURED_DOMAIN_MOTIF
 *  flag via its @p loop_type parameter the contribution of any ligand that consecutively
 *  binds from position @f$i@f$ to @f$j@f$ (the white box) is requested. Otherwise, the callback
 *  usually performs a lookup in the precomputed @p B matrices. Which @p B matrix is
 *  addressed will be indicated by the flags #VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP, #VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP
 *  #VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP, and #VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP. As their names already imply,
 *  they specify exterior loops (@p F production rule), hairpin loops and internal loops
 *  (@p C production rule), and multibranch loops (@p M and @p M1 production rule).
 *
 *  @image xml ligands_up_callback.svg
 *
 *  The @p pre_cb callback will be executed as a pre-processing step right before the
 *  regular secondary structure rules. Usually one would use this callback to fill the
 *  dynamic programming matrices @p U and preparations of the auxiliary data structure
 *  #vrna_unstructured_domain_s.data
 *
 *  @image xml B_prod_rule.svg
 *
 *  @ingroup domains_up
 *
 *  @param  fc      The #vrna_fold_compound_t data structure the callback will be bound to
 *  @param  pre_cb  A pointer to a callback function for the @p B production rule
 *  @param  e_cb    A pointer to a callback function for free energy evaluation
 */
void vrna_ud_set_prod_rule_cb(vrna_fold_compound_t        *fc,
                              vrna_ud_production_f pre_cb,
                              vrna_ud_f     e_cb);


/**
 *  @brief Attach production rule for partition function
 *
 *  This function is the partition function companion of vrna_ud_set_prod_rule_cb().
 *
 *  Use it to bind callbacks to (i) fill the @p U production rule dynamic programming
 *  matrices and/or prepare the #vrna_unstructured_domain_s.data, and (ii) provide a callback
 *  to retrieve partition functions for subsegments @f$ [i,j] @f$.
 *
 *  @image xml   B_prod_rule.svg
 *
 *  @image xml   ligands_up_callback.svg
 *
 *  @ingroup domains_up
 *
 *  @see vrna_ud_set_prod_rule_cb()
 *
 *  @param  fc        The #vrna_fold_compound_t data structure the callback will be bound to
 *  @param  pre_cb    A pointer to a callback function for the @p B production rule
 *  @param  exp_e_cb  A pointer to a callback function that retrieves the partition function
 *                    for a segment @f$[i,j]@f$ that may be bound by one or more ligands.
 */
void  vrna_ud_set_exp_prod_rule_cb(vrna_fold_compound_t             *fc,
                                   vrna_ud_exp_production_f  pre_cb,
                                   vrna_ud_exp_f      exp_e_cb);


void  vrna_ud_set_prob_cb(vrna_fold_compound_t        *fc,
                          vrna_ud_add_probs_f  setter,
                          vrna_ud_get_probs_f  getter);


#endif
