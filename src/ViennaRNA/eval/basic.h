#ifndef VIENNA_RNA_PACKAGE_EVAL_BASIC_H
#define VIENNA_RNA_PACKAGE_EVAL_BASIC_H


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
 *  @file     ViennaRNA/eval/basic.h
 *  @ingroup  eval
 *  @brief    Declarations that are essential for many energy evaluation functions
 */


/**
 *  @addtogroup eval
 *  @{
 */


/**
 *  @brief  Quiet level verbosity setting
 */
#define VRNA_VERBOSITY_QUIET     -1


/**
 *  @brief  Default level verbosity setting
 */
#define VRNA_VERBOSITY_DEFAULT    1


/**
 *  @brief  Default options for loop evaluation
 */
#define VRNA_EVAL_LOOP_DEFAULT          0U


/**
 *  @brief  Do not acknowledge hard constraints in loop evaluation
 */
#define VRNA_EVAL_LOOP_NO_HC            1U


/**
 *  @brief  Do not acknowledge soft constraints in loop evaluation
 */
#define VRNA_EVAL_LOOP_NO_SC            2U

/**
 *  @brief  Do not acknowledge any constraints in loop evaluation
 */
#define VRNA_EVAL_LOOP_NO_CONSTRAINTS   (VRNA_EVAL_LOOP_NO_HC | VRNA_EVAL_LOOP_NO_SC)


/**
 * @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup eval_deprecated
 *  @{
 */

/**
 *  @brief first pos of second seq for cofolding
 */
DEPRECATED(extern int  cut_point,
           "Use strand information from vrna_fold_compound_t instead!");

/**
 *  @brief verbose info from energy_of_struct
 */
DEPRECATED(extern int  eos_debug,
           "This variable is obsolete!");

#endif

/**
 * @}
 */

#endif
