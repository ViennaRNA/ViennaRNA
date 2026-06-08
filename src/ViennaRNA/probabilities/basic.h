#ifndef VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_BASIC_H
#define VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_BASIC_H

#include <ViennaRNA/fold_compound.h>

/**
 *  @file     ViennaRNA/probabilities/basic.h
 *  @ingroup  thermodynamics
 *  @brief    Basic Equilibrium Probability implementations
 *
 *  This file includes the API for the basic equilibrium
 *  probability computations based on the partition function
 *  of an RNA sequence, two concatenated sequences, or a sequence
 *  alignment.
 */


/**
 *  @addtogroup  thermodynamics
 *  @{
 */


/**
 *  @brief  Exterior loop
 */
#define VRNA_LOOP_TYPE_EXT    (unsigned int)(1<<0)

/**
 *  @brief  Hairpin loop
 */
#define VRNA_LOOP_TYPE_HP     (unsigned int)(1<<1)

/**
 *  @brief  Internal loop
 */
#define VRNA_LOOP_TYPE_INT    (unsigned int)(1<<2)

/**
 *  @brief  Multibranch loop
 */
#define VRNA_LOOP_TYPE_MB     (unsigned int)(1<<3)

/**
 *  @brief  Any loop
 */
#define VRNA_LOOP_TYPE_ANY   (VRNA_LOOP_TYPE_EXT | VRNA_LOOP_TYPE_HP | VRNA_LOOP_TYPE_INT | VRNA_LOOP_TYPE_MB)


/**
 *  @brief  Base pair probability flag
 *
 *  Passing this flag to vrna_eqm_probabilities() or vrna_probs_window() activates
 *  callback execution for base pairing probabilities. In turn, the corresponding
 *  callback receives this flag through the @p type argument whenever base pairing
 *  probabilities are provided.
 *
 *  @see  vrna_eqm_probabilities(), vrna_probs_window()
 */
#define VRNA_EQM_PROBS_BPP    (unsigned int)(1 << 12)

/**
 *  @brief  Unpaired probability (accessibility) flag
 *
 *  Passing this flag to vrna_eqm_probabilities() or vrna_probs_window() activates
 *  callback execution for unpaired probabilities, a.k.a. accessibilities. In turn,
 *  the corresponding callback receives this flag through the @p type argument
 *  whenever unpaired probabilities are provided.
 *
 *  @see  vrna_eqm_probabilities(), vrna_probs_window()
 */
#define VRNA_EQM_PROBS_UP     (unsigned int)(1 << 13)

/**
 *  @brief Helix probability flag
 *
 *  Passing this flag to vrna_eqm_probabilities() or vrna_probs_window() activates
 *  callback execution for helix probabilities. In turn, the corresponding callback
 *  receives this flag through the @p type argument whenever stack probabilities are
 *  provided.
 *
 *  @bug  Currently, this flag is a placeholder doing nothing as the corresponding
 *        implementation for helix/stack probability computation is missing.
 *
 *  @see  vrna_eqm_probabilities(), vrna_probs_window()
 */
#define VRNA_EQM_PROBS_HX     (unsigned int)(1 << 14)

/**
 *  @brief  Split probabilities into loop types
 *
 *  Passing this flag to vrna_eqm_probabilities() or vrna_probs_window() indicates that
 *  probabilities should be reported on a loop-type base. 
 *  In turn, the corresponding callback receives the #VRNA_PROBS_WINDOW_UP flag OR-ed
 *  together with the corresponding loop type, i.e.:
 *
 *  * #VRNA_EXT_LOOP  - @copybrief #VRNA_EXT_LOOP
 *  * #VRNA_HP_LOOP   - @copybrief #VRNA_HP_LOOP
 *  * #VRNA_INT_LOOP  - @copybrief #VRNA_INT_LOOP
 *  * #VRNA_MB_LOOP   - @copybrief #VRNA_MB_LOOP
 *  * #VRNA_ANY_LOOP  - @copybrief #VRNA_ANY_LOOP
 *
 *  @see  vrna_eqm_probabilities(), vrna_probs_window(),
 *        #VRNA_EQM_PROBS_UP, #VRNA_EQM_PROBS_BPP
 */
#define VRNA_EQM_PROBS_SPLIT  (unsigned int)(1 << 15)

#define VRNA_EQM_PROBS_MASK   (VRNA_EQM_PROBS_BPP | VRNA_EQM_PROBS_UP | VRNA_EQM_PROBS_HX)

typedef struct {
  unsigned int  mode;
  unsigned int  up_length;
} vrna_eqm_probs_opt_t;

typedef void *vrna_eqm_probs_data_t;

typedef void (*vrna_eqm_probs_f)(vrna_fold_compound_t   *fc,
                                 FLT_OR_DBL             *pr,
                                 unsigned int           j,
                                 unsigned int           a,
                                 unsigned int           b,
                                 unsigned int           type,
                                 vrna_eqm_probs_data_t  data);

int
vrna_eqm_probabilities(vrna_fold_compound_t   *fc,
                       vrna_eqm_probs_f       cb,
                       vrna_eqm_probs_data_t  cb_data,
                       vrna_eqm_probs_opt_t   options);

/* End thermodynamics */
/**@}*/

#endif
