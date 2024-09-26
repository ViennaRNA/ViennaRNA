#ifndef VIENNA_RNA_PACKAGE_GRAMMAR_BASIC_H
#define VIENNA_RNA_PACKAGE_GRAMMAR_BASIC_H

/**
 *  @file     ViennaRNA/grammar/basic.h
 *  @ingroup  grammar
 *  @brief    Basic RNA folding grammar API
 */

/**
 * @addtogroup grammar
 * @{
 */

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/datastructures/array.h>

/**
 *  @brief  A pointer to the auxiliary grammar data structure
 */
typedef struct vrna_gr_aux_s *vrna_gr_aux_t;


/**
 *  @brief  Function prototype for serializing backtracked base pairs and structure elements into a dot-bracket string
 *
 *  This function will be called after backtracking in the MFE predictions to
 *  convert collected base pairs and other information into a dot-bracket-like
 *  structure string.
 *
 *  @callback
 *  @parblock
 *  This callback allows for changing the way how base pairs (and other types of data)
 *  obtained from the default and extended grammar are converted back into a dot-bracket string.
 *  @endparblock
 *
 *  @see  vrna_db_from_bps(), vrna_gr_add_aux()
 *
 *  @param  fc        The fold compound to work on
 *  @param  bp_stack  The base pair stack
 *  @param  data      An arbitrary user-provided data pointer
 *  @return           A '\0'-terminated dot-bracket-like string representing the structure from @p bp_stack (and @p data)
 */
typedef char * (*vrna_gr_serialize_bp_f)(vrna_fold_compound_t *fc,
                                         vrna_bps_t           bp_stack,
                                         void                 *data);


/**
 *  @brief  Prepare the auxiliary grammar rule data
 *
 *  @note   This function is mainly for internal use. Users of the auxiliary grammar API usually
 *          do not need to call this function except for debugging purposes.
 *
 *  @param  fc      The fold compound storing the auxiliary grammar rules
 *  @param  options Options flag(s) that denote what needs to be prepared
 *  @return         non-zero on success, 0 otherwise
 */
unsigned int
vrna_gr_prepare(vrna_fold_compound_t  *fc,
                unsigned int          options);


/**
 *  @brief  Add status event based data preparation callbacks
 *
 *  This function binds additional data structures and corresponding callback
 *  functions for the auxiliary grammar extension API. This might be helpful whenever
 *  certain preparation steps need to be done prior and/or after the actual run of the
 *  prediction algorithms.
 *
 *  @param  fc      The fold compound
 *  @param  cb      The recursion status callback that performs the preparation
 *  @param  data    The data pointer the @p cb callback is working on
 *  @param  prepare_cb  A preparation callback function for parameter @p data
 *  @param  free_cb     A callback to release memory for @p data
 *  @return             The number of status function callbacks bound to the fold compound or 0 on error
 */
unsigned int
vrna_gr_add_status(vrna_fold_compound_t     *fc,
                   vrna_recursion_status_f  cb,
                   void                     *data,
                   vrna_auxdata_prepare_f   prepare_cb,
                   vrna_auxdata_free_f      free_cb);


/**
 *  @brief  Set base pair stack to dot-bracket string serialization function
 *
 *  After backtracking secondary structures, e.g. in MFE predictions, the outside
 *  algorithm usually collects a set of base pairs that then need to be converted
 *  into a dot-bracket string. By default, this conversion is done using the
 *  vrna_db_from_bps() function. However, this function only considers nested
 *  base pairs and no other type of secondary structure elements.
 *
 *  When extending the recursions by additional rules, the default conversion
 *  may not suffice, e.g. because the grammar extension adds 2.5D modules or
 *  pseudoknots. In such cases the user should implement its own dot-bracket
 *  string conversion strategy that may use additional symbols.
 *
 *  This function binds a user-implemented conversion function that must
 *  return a '\0' terminated dot-bracket-like string the same length as the
 *  input sequence. The conversion function will then be used instead of the
 *  default one. In addition to the base pair stack #vrna_bps_t, the user-defined
 *  conversion function may keep track of whatever information is neccessary to
 *  properly convert the backtracked structure into a dot-bracket string. For
 *  that purpose, the @p data pointer can be used, e.g. it can point to the same
 *  data as used in any of the grammar extension rules. The @p prepare_cb and
 *  @p free_cb callbacks can again be used to control preparation and release
 *  of the memory @p data points to. The @p prepare_cb will be called after
 *  all the preparations for the grammar extensions and prior the actual
 *  inside-recursions. The conversion function callback @p cb will be called
 *  after backtracking.
 *
 *  @see  vrna_db_from_bps(), vrna_gr_add_aux()
 *
 *  @param  fc          The fold compound
 *  @param  cb          A pointer to the conversion callback function
 *  @param  data        A pointer to arbitrary data that will be passed through to @p cb (may be @b NULL)
 *  @param  prepare_cb  A function pointer to prepare @p data (may be @b NULL)
 *  @param  free_cb     A function pointer to release memory occupied by @p data (may be @b NULL)
 */
unsigned int
vrna_gr_set_serialize_bp(vrna_fold_compound_t   *fc,
                         vrna_gr_serialize_bp_f cb,
                         void                   *data,
                         vrna_auxdata_prepare_f prepare_cb,
                         vrna_auxdata_free_f    free_cb);


/**
 *  @brief  Remove all auxiliary grammar rules
 *
 *  This function re-sets the fold compound to the default rules
 *  by removing all auxiliary grammar rules
 *
 *  @param  fc  The fold compound
 */
void
vrna_gr_reset(vrna_fold_compound_t *fc);


/**
 * @}
 */

#endif
