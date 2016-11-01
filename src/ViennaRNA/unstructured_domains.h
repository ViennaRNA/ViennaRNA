#ifndef VIENNA_RNA_PACKAGE_UNSTRUCTURED_DOMAIN_H
#define VIENNA_RNA_PACKAGE_UNSTRUCTURED_DOMAIN_H

/**
 *  @file unstructured_domains.h
 *  @ingroup domains_up
 *  @brief    Functions to modify unstructured domains, e.g. to incorporate ligands binding to unpaired stretches
 */

/**
 *  @addtogroup domains_up
 *
 *  @brief  Add and modify unstructured domains to the RNA folding grammar
 *
 *  This module provides the tools to add and modify unstructured domains to the production rules of the RNA folding grammar.
 *  Usually this functionality is utilized for incorporating ligand binding to unpaired stretches of an RNA.
 *
 *  Unstructured domains appear in the production rules of the RNA folding grammar
 *  whereever new unpaired nucleotides are attached to a growing substructure (see also @cite lorenz:2016b):
 *  @image html   Crecursion.svg
 *  @image latex  Crecursion.eps
 *
 *  The white boxes represent the stretch of RNA bound to the ligand and represented by a more or less specific
 *  sequence motif. The motif itself is considered unable to form basepairs. The additional production rule @em U
 *  is used to precompute the contribution of unpaired stretches possibly bound by one or more ligands. The
 *  auxiliary DP matrix for this production rule is filled right before processing the other (regular) production
 *  rules of the RNA folding grammar.
 *
 *  In a context with @ref domains_struc the grammar is extended as follows:
 *
 *  @image html   GCrecursion.svg
 *  @image latex  GCrecursion.eps
 *
 *  @bug  Although the additional production rule(s) for unstructured domains in the descriptions of this feature
 *        are always treated as 'segments possibly bound to one or more ligands', the current implementation requires
 *        that at least one ligand is bound. The default implementation already takes care of the required changes,
 *        however, upon using callback functions other than the default ones, one has to take care of this fact.
 *        Please also note, that this behavior might change in one of the next releases, such that the decomposition
 *        schemes as shown above comply with the actual implementation.
 *
 *  A default implementation allows one to readily use this feature by simply adding sequence motifs and corresponding
 *  binding free energies with the function vrna_ud_add_motif() (see also @ref ligands_up).
 *
 *  The grammar extension is realized using a callback function that
 *  - evaluates the binding free energy of a ligand to its target sequence segment (white boxes in the figures above), or
 *  - returns the free energy of an unpaired stretch possibly bound by a ligand, stored in the additional @em U DP matrix.
 *
 *  The callback is passed the segment positions, the loop context, and which of the two above mentioned
 *  evaluations are required. A second callback implements the pre-processing step that
 *  prepares the @em U DP matrix by evaluating all possible cases of the additional production rule.
 *  Both callbacks have a default implementation in @em RNAlib, but may be over-written by a
 *  user-implementation, making it fully user-customizable.
 *
 *  For equilibrium probability computations, two additional callbacks exist. One to store/add and one to retrieve the
 *  probability of unstructured domains at particular positions. Our implementation already takes care of computing
 *  the probabilities, but users of the unstructured domain feature are required to provide a mechanism to efficiently
 *  store/add the corresponding values into some external data structure.
 */

/** @brief Typename for the ligand binding extension data structure #vrna_unstructured_domain_s
 *  @ingroup domains_up
 */
typedef struct vrna_unstructured_domain_s  vrna_ud_t;

typedef struct vrna_unstructured_domain_motif_s  vrna_ud_motif_t;

#include <ViennaRNA/data_structures.h>

/**
 *  @brief Callback to retrieve binding free energy of a ligand bound to an unpaired sequence segment
 *  @ingroup domains_up
 */
typedef int (vrna_callback_ud_energy)(vrna_fold_compound_t *vc, int i, int j, unsigned int loop_type, void *data);

/**
 *  @brief Callback to retrieve Boltzmann factor of the binding free energy of a ligand bound to an unpaired sequence segment
 *  @ingroup domains_up
 */
typedef FLT_OR_DBL (vrna_callback_ud_exp_energy)(vrna_fold_compound_t *vc, int i, int j, unsigned int loop_type, void *data);

/**
 *  @brief Callback for pre-processing the production rule of the ligand binding to unpaired stretches feature
 *  @ingroup domains_up
 */
typedef void (vrna_callback_ud_production)(vrna_fold_compound_t *vc, void *data);

/**
 *  @brief Callback for pre-processing the production rule of the ligand binding to unpaired stretches feature (partition function variant)
 *  @ingroup domains_up
 */
typedef void (vrna_callback_ud_exp_production)(vrna_fold_compound_t *vc, void *data);


/**
 *  @brief Callback to store/add equilibrium probability for a ligand bound to an unpaired sequence segment
 *  @ingroup domains_up
 */
typedef void (vrna_callback_ud_probs_add)(vrna_fold_compound_t *vc, int i, int j, unsigned int loop_type, FLT_OR_DBL exp_energy, void *data);

/**
 *  @brief Callback to retrieve equilibrium probability for a ligand bound to an unpaired sequence segment
 *  @ingroup domains_up
 */
typedef FLT_OR_DBL (vrna_callback_ud_probs_get)(vrna_fold_compound_t *vc, int i, int j, unsigned int loop_type, int motif, void *data);


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
 *  @brief Flag to indicate ligand bound to unpiared stretch in an interior loop
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
#define VRNA_UNSTRUCTURED_DOMAIN_ALL_LOOPS   (VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP | VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP)

/**
 *  @brief  Data structure to store all functionality for ligand binding
 *  @ingroup domains_up
 */
struct vrna_unstructured_domain_s {

  /*
  **********************************
    Keep track of all added motifs
  **********************************
  */
  int                             uniq_motif_count; /**<  @brief The unique number of motifs of different lengths */ 
  unsigned int                    *uniq_motif_size; /**<  @brief An array storing a unique list of motif lengths */

  int                             motif_count;      /**<  @brief Total number of distinguished motifs */
  char                            **motif;          /**<  @brief Motif sequences */
  unsigned int                    *motif_size;      /**<  @brief Motif lengths */
  double                          *motif_en;        /**<  @brief Ligand binding free energy contribution */
  unsigned int                    *motif_type;      /**<  @brief Type of motif, i.e. loop type the ligand binds to */

  /*
  **********************************
    Grammar extension for ligand
    binding
  **********************************
  */
  vrna_callback_ud_production     *prod_cb;       /**<  @brief Callback to ligand binding production rule, i.e. create/fill DP free energy matrices
                                                   *    @details This callback will be executed right before the actual secondary structure decompositions,
                                                   *    and, therefore, any implementation must not interleave with the regular DP matrices.
                                                   */
  vrna_callback_ud_exp_production *exp_prod_cb;   /**<  @brief Callback to ligand binding production rule, i.e. create/fill DP partition function matrices */
  vrna_callback_ud_energy         *energy_cb;     /**<  @brief Callback to evaluate free energy of ligand binding to a particular unpaired stretch */
  vrna_callback_ud_exp_energy     *exp_energy_cb; /**<  @brief Callback to evaluate Boltzmann factor of ligand binding to a particular unpaired stretch */
  void                            *data;          /**<  @brief Auxiliary data structure passed to energy evaluation callbacks */
  vrna_callback_free_auxdata      *free_data;     /**<  @brief Callback to free auxiliary data structure */
  vrna_callback_ud_probs_add      *probs_add;     /**<  @brief Callback to store/add outside partition function */
  vrna_callback_ud_probs_get      *probs_get;     /**<  @brief Callback to retrieve outside partition function */
};


struct vrna_unstructured_domain_motif_s {
  int start;
  int number;
};


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
 *  such as the exterior loop, hairpin loops, interior loops, or multibranch loops.
 *
 *  @see  #VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP, #VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
 *  #VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP, #VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP, #VRNA_UNSTRUCTURED_DOMAIN_ALL_LOOPS,
 *  vrna_ud_remove()
 *
 *  @ingroup domains_up
 *
 *  @param  vc        The #vrna_fold_compound_t data structure the ligand motif should be bound to
 *  @param  motif     The sequence motif the ligand binds to
 *  @param  motif_en  The binding free energy of the ligand in kcal/mol
 *  @param  loop_type The loop type the ligand binds to
 *
 */
void  vrna_ud_add_motif(vrna_fold_compound_t *vc,
                        const char *motif,
                        double motif_en,
                        unsigned int loop_type);


/**
 *  @brief  Get a list of unique motif sizes that start at a certain position within the sequence
 *
 */
int *vrna_ud_get_motif_size_at( vrna_fold_compound_t *vc,
                                int i,
                                unsigned int loop_type);


int *
vrna_ud_get_motifs_at(vrna_fold_compound_t *vc,
                      int i,
                      unsigned int loop_type);


vrna_ud_motif_t *
vrna_ud_detect_motifs(vrna_fold_compound_t *vc,
                      const char *structure);


/**
 *  @brief Remove ligand binding to unpaired stretches
 *
 *  This function removes all ligand motifs that were bound to a #vrna_fold_compound_t using
 *  the vrna_ud_add_motif() function.
 *
 *  @ingroup domains_up
 *
 *  @param vc The #vrna_fold_compound_t data structure the ligand motif data should be removed from
 */
void  vrna_ud_remove(vrna_fold_compound_t *vc);

/**
 *  @brief  Attach an auxiliary data structure
 *
 *  This function binds an arbitrary, auxiliary data structure for user-implemented ligand binding.
 *  The optional callback @p free will be passed the bound data structure whenever the #vrna_fold_compound_t
 *  is removed from memory to avoid memory leaks.
 *
 *  @see vrna_ud_set_prod_rule_cb(), vrna_ud_set_exp_prod_rule_cb(),
 *  vrna_ud_remove()
 *
 *  @ingroup domains_up
 *
 *  @param  vc      The #vrna_fold_compound_t data structure the auxiliary data structure should be bound to
 *  @param  data    A pointer to the auxiliary data structure
 *  @param  free_cb A pointer to a callback function that free's memory occupied by @p data
 */
void  vrna_ud_set_data( vrna_fold_compound_t        *vc,
                        void                        *data,
                        vrna_callback_free_auxdata  *free_cb);

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
 *  they specify exterior loops (@p F production rule), hairpin loops and interior loops
 *  (@p C production rule), and multibranch loops (@p M and @p M1 production rule).
 *  
 *  @image html   ligands_up_callback.svg
 *  @image latex  ligands_up_callback.eps
 *
 *  The @p pre_cb callback will be executed as a pre-processing step right before the
 *  regular secondary structure rules. Usually one would use this callback to fill the
 *  dynamic programming matrices @p U and preparations of the auxiliary data structure
 *  #vrna_unstructured_domain_s.data
 *
 *  @image html   B_prod_rule.svg
 *  @image latex  B_prod_rule.eps
 *
 *  @ingroup domains_up
 *
 *  @param  vc      The #vrna_fold_compound_t data structure the callback will be bound to
 *  @param  pre_cb  A pointer to a callback function for the @p B production rule
 *  @param  e_cb    A pointer to a callback function for free energy evaluation
 */
void vrna_ud_set_prod_rule_cb(vrna_fold_compound_t        *vc,
                              vrna_callback_ud_production *pre_cb,
                              vrna_callback_ud_energy     *e_cb);


/**
 *  @brief Attach production rule for partition function
 *
 *  This function is the partition function companion of vrna_ud_set_prod_rule_cb().
 *
 *  Use it to bind callbacks to (i) fill the @p U production rule dynamic programming
 *  matrices and/or prepare the #vrna_unstructured_domain_s.data, and (ii) provide a callback
 *  to retrieve partition functions for subsegments @f$ [i,j] @f$.
 *
 *  @image html   B_prod_rule.svg
 *  @image latex  B_prod_rule.eps
 *
 *  @image html   ligands_up_callback.svg
 *  @image latex  ligands_up_callback.eps
 *
 *  @ingroup domains_up
 *
 *  @see vrna_ud_set_prod_rule_cb()
 *
 *  @param  vc        The #vrna_fold_compound_t data structure the callback will be bound to
 *  @param  pre_cb    A pointer to a callback function for the @p B production rule
 *  @param  exp_e_cb  A pointer to a callback function that retrieves the partition function
 *                    for a segment @f$[i,j]@f$ that may be bound by one or more ligands.
 */
void  vrna_ud_set_exp_prod_rule_cb( vrna_fold_compound_t            *vc,
                                    vrna_callback_ud_exp_production *pre_cb,
                                    vrna_callback_ud_exp_energy     *exp_e_cb);


void  vrna_ud_set_prob_cb(vrna_fold_compound_t        *vc,
                          vrna_callback_ud_probs_add  *setter,
                          vrna_callback_ud_probs_get  *getter);

#endif
