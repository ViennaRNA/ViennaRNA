#ifndef VIENNA_RNA_PACKAGE_EVAL_H
#define VIENNA_RNA_PACKAGE_EVAL_H

#include <stdio.h>
#include <ViennaRNA/data_structures.h>
#include "ViennaRNA/neighbor.h"
#include <ViennaRNA/params.h>   /* for deprecated functions */

#ifdef VRNA_WARN_DEPRECATED
# ifdef __GNUC__
#  define DEPRECATED(func) func __attribute__ ((deprecated))
# else
#  define DEPRECATED(func) func
# endif
#else
# define DEPRECATED(func) func
#endif

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

/**
 *  @file     eval.h
 *  @ingroup  eval
 *  @brief    Functions and variables related to energy evaluation of sequence/structure pairs.
 */


/**
 *  @addtogroup eval
 *  @brief Functions and variables related to free energy evaluation of sequence/structure pairs.
 *
 *  @{
 *  @ingroup  eval
 */

/** @brief set to first pos of second seq for cofolding  */
extern int  cut_point;

/**
 *  @brief verbose info from energy_of_struct
 */
extern int  eos_debug;

/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function allows for energy evaluation of a given pair of structure
 *  and sequence (alignment).
 *  Model details, energy parameters, and possibly soft constraints are used as provided
 *  via the parameter 'vc'. The #vrna_fold_compound_t does not need to contain any DP matrices,
 *  but requires all most basic init values as one would get from a call like this:
 *  @code{.c}
 * vc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_EVAL_ONLY);
 *  @endcode
 *
 *  @note Accepts vrna_fold_compound_t of type #VRNA_FC_TYPE_SINGLE and #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @see  vrna_eval_structure_pt(), vrna_eval_structure_verbose(), vrna_eval_structure_pt_verbose(),
 *        vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_eval_covar_structure()
 *
 *  @param vc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param structure        Secondary structure in dot-bracket notation
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float vrna_eval_structure(vrna_fold_compound_t  *vc,
                          const char            *structure);


/**
 *  @brief Calculate the pseudo energy derived by the covariance scores of a set of aligned sequences
 *
 *  Consensus structure prediction is driven by covariance scores of base pairs in rows of the
 *  provided alignment. This function allows one to retrieve the total amount of this covariance pseudo
 *  energy scores.
 *  The #vrna_fold_compound_t does not need to contain any DP matrices, but requires all most basic
 *  init values as one would get from a call like this:
 *  @code{.c}
 * vc = vrna_fold_compound_comparative(alignment, NULL, VRNA_OPTION_EVAL_ONLY);
 *  @endcode
 *
 *  @note Accepts vrna_fold_compound_t of type #VRNA_FC_TYPE_COMPARATIVE only!
 *
 *  @see  vrna_fold_compound_comparative(), vrna_eval_structure()
 *
 *  @param vc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param structure        Secondary (consensus) structure in dot-bracket notation
 *  @return                 The covariance pseudo energy score of the input structure given the input sequence alignment in kcal/mol
 */
float vrna_eval_covar_structure(vrna_fold_compound_t  *vc,
                                const char            *structure);


/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function allows for energy evaluation of a given sequence/structure pair.
 *  In contrast to vrna_eval_structure() this function assumes default model details
 *  and default energy parameters in order to evaluate the free energy of the secondary
 *  structure. Therefore, it serves as a simple interface function for energy evaluation
 *  for situations where no changes on the energy model are required.
 *
 *  @see vrna_eval_structure(), vrna_eval_structure_pt(), vrna_eval_structure_verbose(), vrna_eval_structure_pt_verbose(),
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param structure        Secondary structure in dot-bracket notation
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float vrna_eval_structure_simple(const char *string,
                                 const char *structure);


/**
 *  @brief Calculate the free energy of an already folded RNA and print contributions on a per-loop base.
 *
 *  This function is a simplyfied version of vrna_eval_structure_v() that uses the @em default
 *  verbosity level.
 * (
 *  @see vrna_eval_structure_pt(), vrna_eval_structure_verbose(), vrna_eval_structure_pt_verbose(),
 *
 *  @param vc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float vrna_eval_structure_verbose(vrna_fold_compound_t  *vc,
                                  const char            *structure,
                                  FILE                  *file);


/**
 *  @brief Calculate the free energy of an already folded RNA and print contributions on a per-loop base.
 *
 *  This function allows for detailed energy evaluation of a given sequence/structure pair.
 *  In contrast to vrna_eval_structure() this function prints detailed energy contributions
 *  based on individual loops to a file handle. If NULL is passed as file handle, this function
 *  defaults to print to stdout. Any positive @p verbosity_level activates potential warning message
 *  of the energy evaluting functions, while values @f$ \ge 1 @f$ allow for detailed control of what
 *  data is printed. A negative parameter @p verbosity_level turns off printing all together.
 *
 *  Model details, energy parameters, and possibly soft constraints are used as provided
 *  via the parameter 'vc'. The fold_compound does not need to contain any DP matrices,
 *  but all the most basic init values as one would get from a call like this:
 *  @code{.c}
 * vc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_EVAL_ONLY);
 *  @endcode
 *
 *  @see vrna_eval_structure_pt(), vrna_eval_structure_verbose(), vrna_eval_structure_pt_verbose(),
 *
 *  @param vc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float vrna_eval_structure_v(vrna_fold_compound_t  *vc,
                            const char            *structure,
                            int                   verbosity_level,
                            FILE                  *file);


/**
 *  @brief Calculate the free energy of an already folded RNA and print contributions per loop.
 *
 *  This function is a simplyfied version of vrna_eval_structure_simple_v() that uses the @em default
 *  verbosity level.
 *
 *  @see  vrna_eval_structure_simple_v(), vrna_eval_structure_verbose(), vrna_eval_structure_pt(),
 *        vrna_eval_structure_verbose(), vrna_eval_structure_pt_verbose()
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float vrna_eval_structure_simple_verbose(const char *string,
                                         const char *structure,
                                         FILE       *file);


/**
 *  @brief Calculate the free energy of an already folded RNA and print contributions per loop.
 *
 *  This function allows for detailed energy evaluation of a given sequence/structure pair.
 *  In contrast to vrna_eval_structure() this function prints detailed energy contributions
 *  based on individual loops to a file handle. If NULL is passed as file handle, this function
 *  defaults to print to stdout. Any positive @p verbosity_level activates potential warning message
 *  of the energy evaluting functions, while values @f$ \ge 1 @f$ allow for detailed control of what
 *  data is printed. A negative parameter @p verbosity_level turns off printing all together.
 *
 *  In contrast to vrna_eval_structure_verbose() this function assumes default model details
 *  and default energy parameters in order to evaluate the free energy of the secondary
 *  structure. Threefore, it serves as a simple interface function for energy evaluation
 *  for situations where no changes on the energy model are required.
 *
 *  @see vrna_eval_structure_verbose(), vrna_eval_structure_pt(), vrna_eval_structure_verbose(), vrna_eval_structure_pt_verbose(),
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float vrna_eval_structure_simple_v(const char *string,
                                   const char *structure,
                                   int        verbosity_level,
                                   FILE       *file);


/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function allows for energy evaluation of a given sequence/structure pair where
 *  the structure is provided in pair_table format as obtained from vrna_ptable().
 *  Model details, energy parameters, and possibly soft constraints are used as provided
 *  via the parameter 'vc'. The fold_compound does not need to contain any DP matrices,
 *  but all the most basic init values as one would get from a call like this:
 *  @code{.c}
 * vc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_EVAL_ONLY);
 *  @endcode
 *
 *  @see vrna_ptable(), vrna_eval_structure(), vrna_eval_structure_pt_verbose()
 *
 *  @param vc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param pt               Secondary structure as pair_table
 *  @return                 The free energy of the input structure given the input sequence in 10cal/mol
 */
int vrna_eval_structure_pt(vrna_fold_compound_t *vc,
                           const short          *pt);


/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  In contrast to vrna_eval_structure_pt() this function assumes default model details
 *  and default energy parameters in order to evaluate the free energy of the secondary
 *  structure. Threefore, it serves as a simple interface function for energy evaluation
 *  for situations where no changes on the energy model are required.
 *
 *  @see vrna_ptable(), vrna_eval_structure_simple(), vrna_eval_structure_pt()
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param pt               Secondary structure as pair_table
 *  @return                 The free energy of the input structure given the input sequence in 10cal/mol
 */
int vrna_eval_structure_pt_simple(const char  *string,
                                  const short *pt);


/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function is a simplyfied version of vrna_eval_structure_simple_v() that uses the @em default
 *  verbosity level.
 *
 *  @see vrna_eval_structure_pt_v(), vrna_ptable(), vrna_eval_structure_pt(), vrna_eval_structure_verbose()
 *
 *  @param vc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param pt               Secondary structure as pair_table
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in 10cal/mol
 */
int vrna_eval_structure_pt_verbose(vrna_fold_compound_t *vc,
                                   const short          *pt,
                                   FILE                 *file);


/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function allows for energy evaluation of a given sequence/structure pair where
 *  the structure is provided in pair_table format as obtained from vrna_ptable().
 *  Model details, energy parameters, and possibly soft constraints are used as provided
 *  via the parameter 'vc'. The fold_compound does not need to contain any DP matrices,
 *  but all the most basic init values as one would get from a call like this:
 *  @code{.c}
 * vc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_EVAL_ONLY);
 *  @endcode
 *  In contrast to vrna_eval_structure_pt() this function prints detailed energy contributions
 *  based on individual loops to a file handle. If NULL is passed as file handle, this function
 *  defaults to print to stdout. Any positive @p verbosity_level activates potential warning message
 *  of the energy evaluting functions, while values @f$ \ge 1 @f$ allow for detailed control of what
 *  data is printed. A negative parameter @p verbosity_level turns off printing all together.
 *
 *  @see vrna_ptable(), vrna_eval_structure_pt(), vrna_eval_structure_verbose()
 *
 *  @param vc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param pt               Secondary structure as pair_table
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in 10cal/mol
 */
int vrna_eval_structure_pt_v(vrna_fold_compound_t *vc,
                             const short          *pt,
                             int                  verbosity_level,
                             FILE                 *file);


/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function is a simplyfied version of vrna_eval_structure_pt_simple_v() that uses the @em default
 *  verbosity level.
 *
 *  @see vrna_eval_structure_pt_simple_v(), vrna_ptable(), vrna_eval_structure_pt_verbose(), vrna_eval_structure_simple()
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param pt               Secondary structure as pair_table
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in 10cal/mol
 */
int vrna_eval_structure_pt_simple_verbose(const char  *string,
                                          const short *pt,
                                          FILE        *file);


/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function allows for energy evaluation of a given sequence/structure pair where
 *  the structure is provided in pair_table format as obtained from vrna_ptable().
 *  Model details, energy parameters, and possibly soft constraints are used as provided
 *  via the parameter 'vc'. The fold_compound does not need to contain any DP matrices,
 *  but all the most basic init values as one would get from a call like this:
 *  @code{.c}
 * vc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_EVAL_ONLY);
 *  @endcode
 *  In contrast to vrna_eval_structure_pt_verbose() this function assumes default model details
 *  and default energy parameters in order to evaluate the free energy of the secondary
 *  structure. Threefore, it serves as a simple interface function for energy evaluation
 *  for situations where no changes on the energy model are required.
 *
 *  @see vrna_ptable(), vrna_eval_structure_pt_v(), vrna_eval_structure_simple()
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param pt               Secondary structure as pair_table
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in 10cal/mol
 */
int vrna_eval_structure_pt_simple_v(const char  *string,
                                    const short *pt,
                                    int         verbosity_level,
                                    FILE        *file);


/**
 * @brief Calculate energy of a loop
 *
 *  @param vc         A vrna_fold_compound_t containing the energy parameters and model details
 *  @param i          position of covering base pair
 *  @param pt         the pair table of the secondary structure
 *  @returns          free energy of the loop in 10cal/mol
 */
int vrna_eval_loop_pt(vrna_fold_compound_t  *vc,
                      int                   i,
                      const short           *pt);


/**
 * @brief Calculate energy of a move (closing or opening of a base pair)
 *
 *  If the parameters m1 and m2 are negative, it is deletion (opening)
 *  of a base pair, otherwise it is insertion (opening).
 *
 *  @see              vrna_eval_move_pt()
 *  @param vc         A vrna_fold_compound_t containing the energy parameters and model details
 *  @param structure  secondary structure in dot-bracket notation
 *  @param m1         first coordinate of base pair
 *  @param m2         second coordinate of base pair
 *  @returns          energy change of the move in kcal/mol
 */
float vrna_eval_move(vrna_fold_compound_t *vc,
                     const char           *structure,
                     int                  m1,
                     int                  m2);


/**
 *
 * @brief Calculate energy of a move (closing or opening of a base pair)
 *
 *  If the parameters m1 and m2 are negative, it is deletion (opening)
 *  of a base pair, otherwise it is insertion (opening).
 *
 *  @see              vrna_eval_move()
 *  @param vc         A vrna_fold_compound_t containing the energy parameters and model details
 *  @param pt         the pair table of the secondary structure
 *  @param m1         first coordinate of base pair
 *  @param m2         second coordinate of base pair
 *  @returns          energy change of the move in 10cal/mol
 */
int vrna_eval_move_pt(vrna_fold_compound_t  *vc,
                      short                 *pt,
                      int                   m1,
                      int                   m2);


int vrna_eval_move_pt_simple(const char *string,
                             short      *pt,
                             int        m1,
                             int        m2);


int
vrna_eval_move_shift_pt(vrna_fold_compound_t  *vc,
                        vrna_move_t           *m,
                        short                 *structure);


#ifdef VRNA_BACKWARD_COMPAT

/**
 *  @brief Calculate the free energy of an already folded RNA using global model detail settings
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  @note OpenMP: This function relies on several global model settings variables and thus is
 *        not to be considered threadsafe. See energy_of_struct_par() for a completely threadsafe
 *        implementation.
 *
 *  @deprecated Use vrna_eval_structure() or vrna_eval_structure_verbose() instead!
 *
 *  @see vrna_eval_structure()
 *
 *  @param string     RNA sequence
 *  @param structure  secondary structure in dot-bracket notation
 *  @param verbosity_level a flag to turn verbose output on/off
 *  @return          the free energy of the input structure given the input sequence in kcal/mol
 */
DEPRECATED(float energy_of_structure(const char *string,
                                     const char *structure,
                                     int        verbosity_level));

/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  @deprecated Use vrna_eval_structure() or vrna_eval_structure_verbose() instead!
 *
 *  @see vrna_eval_structure()
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param parameters       A data structure containing the prescaled energy contributions and the model details.
 *  @param verbosity_level  A flag to turn verbose output on/off
 *  @return                The free energy of the input structure given the input sequence in kcal/mol
 */
DEPRECATED(float energy_of_struct_par(const char    *string,
                                      const char    *structure,
                                      vrna_param_t  *parameters,
                                      int           verbosity_level));

/**
 *  @brief Calculate the free energy of an already folded  circular RNA
 *
 *  @note OpenMP: This function relies on several global model settings variables and thus is
 *        not to be considered threadsafe. See energy_of_circ_struct_par() for a completely threadsafe
 *        implementation.
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  @deprecated Use vrna_eval_structure() or vrna_eval_structure_verbose() instead!
 *
 *  @see vrna_eval_structure()
 *
 *  @param string           RNA sequence
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param verbosity_level  A flag to turn verbose output on/off
 *  @return                The free energy of the input structure given the input sequence in kcal/mol
 */
DEPRECATED(float energy_of_circ_structure(const char  *string,
                                          const char  *structure,
                                          int         verbosity_level));

/**
 *  @brief Calculate the free energy of an already folded circular RNA
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  @deprecated Use vrna_eval_structure() or vrna_eval_structure_verbose() instead!
 *
 *  @see vrna_eval_structure()
 *
 *  @param string           RNA sequence
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param parameters       A data structure containing the prescaled energy contributions and the model details.
 *  @param verbosity_level  A flag to turn verbose output on/off
 *  @return                The free energy of the input structure given the input sequence in kcal/mol
 */
DEPRECATED(float energy_of_circ_struct_par(const char   *string,
                                           const char   *structure,
                                           vrna_param_t *parameters,
                                           int          verbosity_level));


DEPRECATED(float energy_of_gquad_structure(const char *string,
                                           const char *structure,
                                           int        verbosity_level));

DEPRECATED(float energy_of_gquad_struct_par(const char    *string,
                                            const char    *structure,
                                            vrna_param_t  *parameters,
                                            int           verbosity_level));


/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  @note OpenMP: This function relies on several global model settings variables and thus is
 *        not to be considered threadsafe. See energy_of_struct_pt_par() for a completely threadsafe
 *        implementation.
 *
 *  @deprecated Use vrna_eval_structure_pt() or vrna_eval_structure_pt_verbose() instead!
 *
 *  @see vrna_eval_structure_pt()
 *
 *  @param string     RNA sequence
 *  @param ptable     the pair table of the secondary structure
 *  @param s          encoded RNA sequence
 *  @param s1         encoded RNA sequence
 *  @param verbosity_level a flag to turn verbose output on/off
 *  @return          the free energy of the input structure given the input sequence in 10kcal/mol
 */
DEPRECATED(int energy_of_structure_pt(const char  *string,
                                      short       *ptable,
                                      short       *s,
                                      short       *s1,
                                      int         verbosity_level));

/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  @deprecated Use vrna_eval_structure_pt() or vrna_eval_structure_pt_verbose() instead!
 *
 *  @see vrna_eval_structure_pt()
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param ptable           The pair table of the secondary structure
 *  @param s                Encoded RNA sequence
 *  @param s1               Encoded RNA sequence
 *  @param parameters       A data structure containing the prescaled energy contributions and the model details.
 *  @param verbosity_level  A flag to turn verbose output on/off
 *  @return                The free energy of the input structure given the input sequence in 10kcal/mol
 */
DEPRECATED(int energy_of_struct_pt_par(const char   *string,
                                       short        *ptable,
                                       short        *s,
                                       short        *s1,
                                       vrna_param_t *parameters,
                                       int          verbosity_level));


/**
 * @brief Calculate energy of a move (closing or opening of a base pair)
 *
 *  If the parameters m1 and m2 are negative, it is deletion (opening)
 *  of a base pair, otherwise it is insertion (opening).
 *
 *  @deprecated Use vrna_eval_move() instead!
 *
 *  @see vrna_eval_move()
 *
 *  @param string     RNA sequence
 *  @param structure  secondary structure in dot-bracket notation
 *  @param m1         first coordinate of base pair
 *  @param m2         second coordinate of base pair
 *  @returns          energy change of the move in kcal/mol
 */
DEPRECATED(float energy_of_move(const char  *string,
                                const char  *structure,
                                int         m1,
                                int         m2));


/**
 *
 * @brief Calculate energy of a move (closing or opening of a base pair)
 *
 *  If the parameters m1 and m2 are negative, it is deletion (opening)
 *  of a base pair, otherwise it is insertion (opening).
 *
 *  @deprecated Use vrna_eval_move_pt() instead!
 *
 *  @see vrna_eval_move_pt()
 *
 *  @param pt         the pair table of the secondary structure
 *  @param s          encoded RNA sequence
 *  @param s1         encoded RNA sequence
 *  @param m1         first coordinate of base pair
 *  @param m2         second coordinate of base pair
 *  @returns          energy change of the move in 10cal/mol
 */
DEPRECATED(int energy_of_move_pt(short  *pt,
                                 short  *s,
                                 short  *s1,
                                 int    m1,
                                 int    m2));

/**
 * @brief Calculate energy of a loop
 *
 *  @deprecated Use vrna_eval_loop_pt() instead!
 *
 *  @see vrna_eval_loop_pt()
 *
 *  @param ptable     the pair table of the secondary structure
 *  @param s          encoded RNA sequence
 *  @param s1         encoded RNA sequence
 *  @param i          position of covering base pair
 *  @returns          free energy of the loop in 10cal/mol
 */
DEPRECATED(int   loop_energy(short  *ptable,
                             short  *s,
                             short  *s1,
                             int    i));

/**
 *  Calculate the free energy of an already folded RNA
 *
 *  @note This function is not entirely threadsafe! Depending on the state of the global
 *  variable @ref eos_debug it prints energy information to stdout or not...\n
 *
 *  @deprecated This function is deprecated and should not be used in future programs!
 *  Use @ref energy_of_structure() instead!
 *
 *  @see              energy_of_structure, energy_of_circ_struct(), energy_of_struct_pt()
 *  @param string     RNA sequence
 *  @param structure  secondary structure in dot-bracket notation
 *  @return          the free energy of the input structure given the input sequence in kcal/mol
 */
DEPRECATED(float energy_of_struct(const char  *string,
                                  const char  *structure));

/**
 *  Calculate the free energy of an already folded RNA
 *
 *  @note This function is not entirely threadsafe! Depending on the state of the global
 *  variable @ref eos_debug it prints energy information to stdout or not...\n
 *
 *  @deprecated This function is deprecated and should not be used in future programs!
 *  Use @ref energy_of_structure_pt() instead!
 *
 *  @see              make_pair_table(), energy_of_structure()
 *  @param string     RNA sequence
 *  @param ptable     the pair table of the secondary structure
 *  @param s          encoded RNA sequence
 *  @param s1         encoded RNA sequence
 *  @return          the free energy of the input structure given the input sequence in 10kcal/mol
 */
DEPRECATED(int energy_of_struct_pt(const char *string,
                                   short      *ptable,
                                   short      *s,
                                   short      *s1));

/**
 *  Calculate the free energy of an already folded  circular RNA
 *
 *  @note This function is not entirely threadsafe! Depending on the state of the global
 *  variable @ref eos_debug it prints energy information to stdout or not...\n
 *
 *  @deprecated This function is deprecated and should not be used in future programs
 *  Use @ref energy_of_circ_structure() instead!
 *
 *  @see              energy_of_circ_structure(), energy_of_struct(), energy_of_struct_pt()
 *  @param string     RNA sequence
 *  @param structure  secondary structure in dot-bracket notation
 *  @return          the free energy of the input structure given the input sequence in kcal/mol
 */
DEPRECATED(float energy_of_circ_struct(const char *string,
                                       const char *structure));

#endif

/**
 * @}
 */

#endif
