#ifndef VIENNA_RNA_PACKAGE_EVAL_BASIC_H
#define VIENNA_RNA_PACKAGE_EVAL_BASIC_H


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
extern int  cut_point;

/**
 *  @brief verbose info from energy_of_struct
 */
extern int  eos_debug;

#endif

/**
 * @}
 */

#endif
