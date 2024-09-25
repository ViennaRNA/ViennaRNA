#ifndef VIENNA_RNA_PACKAGE_EVAL_STRUCTURES_H
#define VIENNA_RNA_PACKAGE_EVAL_STRUCTURES_H

#include <stdio.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/datastructures/char_stream.h>
#include <ViennaRNA/landscape/move.h>
#include <ViennaRNA/params/basic.h>   /* for deprecated functions */
#include <ViennaRNA/eval/basic.h>

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
 *  @file     ViennaRNA/eval/structures.h
 *  @ingroup  eval
 *  @brief    Functions and variables related to energy evaluation of sequence/structure pairs.
 */


/**
 *  @addtogroup eval
 *  @{
 */


/**
 *  @name Basic Energy Evaluation Interface with Dot-Bracket Structure String
 *  @{
 */

/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function allows for energy evaluation of a given pair of structure
 *  and sequence (alignment).
 *  Model details, energy parameters, and possibly soft constraints are used as provided
 *  via the parameter 'fc'. The #vrna_fold_compound_t does not need to contain any DP matrices,
 *  but requires all most basic init values as one would get from a call like this:
 *  @code{.c}
 * fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_EVAL_ONLY);
 *  @endcode
 *
 *  @note Accepts vrna_fold_compound_t of type #VRNA_FC_TYPE_SINGLE and #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @see  vrna_eval_structure_pt(), vrna_eval_structure_verbose(), vrna_eval_structure_pt_verbose(),
 *        vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_eval_covar_structure()
 *
 *  @param fc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param structure        Secondary structure in dot-bracket notation
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float
vrna_eval_structure(vrna_fold_compound_t  *fc,
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
 * fc = vrna_fold_compound_comparative(alignment, NULL, VRNA_OPTION_EVAL_ONLY);
 *  @endcode
 *
 *  @note Accepts vrna_fold_compound_t of type #VRNA_FC_TYPE_COMPARATIVE only!
 *
 *  @see  vrna_fold_compound_comparative(), vrna_eval_structure()
 *
 *  @param fc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param structure        Secondary (consensus) structure in dot-bracket notation
 *  @return                 The covariance pseudo energy score of the input structure given the input sequence alignment in kcal/mol
 */
float
vrna_eval_covar_structure(vrna_fold_compound_t  *fc,
                          const char            *structure);


/**
 *  @brief Calculate the free energy of an already folded RNA and print contributions on a per-loop base.
 *
 *  This function is a simplyfied version of vrna_eval_structure_v() that uses the @em default
 *  verbosity level.
 *
 *  @see vrna_eval_structure_pt(), vrna_eval_structure_verbose(), vrna_eval_structure_pt_verbose(),
 *
 *  @param fc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float
vrna_eval_structure_verbose(vrna_fold_compound_t  *fc,
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
 *  via the parameter 'fc'. The fold_compound does not need to contain any DP matrices,
 *  but all the most basic init values as one would get from a call like this:
 *  @code{.c}
 * fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_EVAL_ONLY);
 *  @endcode
 *
 *  @see vrna_eval_structure_pt(), vrna_eval_structure_verbose(), vrna_eval_structure_pt_verbose(),
 *
 *  @param fc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float
vrna_eval_structure_v(vrna_fold_compound_t  *fc,
                      const char            *structure,
                      int                   verbosity_level,
                      FILE                  *file);


float
vrna_eval_structure_cstr(vrna_fold_compound_t *fc,
                         const char           *structure,
                         int                  verbosity_level,
                         vrna_cstr_t          output_stream);


/* End basic eval interface */
/**@}*/


/**
 *  @name Basic Energy Evaluation Interface with Structure Pair Table
 *  @{
 */

/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function allows for energy evaluation of a given sequence/structure pair where
 *  the structure is provided in pair_table format as obtained from vrna_ptable().
 *  Model details, energy parameters, and possibly soft constraints are used as provided
 *  via the parameter 'fc'. The fold_compound does not need to contain any DP matrices,
 *  but all the most basic init values as one would get from a call like this:
 *  @code{.c}
 * fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_EVAL_ONLY);
 *  @endcode
 *
 *  @see vrna_ptable(), vrna_eval_structure(), vrna_eval_structure_pt_verbose()
 *
 *  @param fc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param pt               Secondary structure as pair_table
 *  @return                 The free energy of the input structure given the input sequence in 10cal/mol
 */
int
vrna_eval_structure_pt(vrna_fold_compound_t *fc,
                       const short          *pt);


/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function is a simplyfied version of vrna_eval_structure_simple_v() that uses the @em default
 *  verbosity level.
 *
 *  @see vrna_eval_structure_pt_v(), vrna_ptable(), vrna_eval_structure_pt(), vrna_eval_structure_verbose()
 *
 *  @param fc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param pt               Secondary structure as pair_table
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in 10cal/mol
 */
int
vrna_eval_structure_pt_verbose(vrna_fold_compound_t *fc,
                               const short          *pt,
                               FILE                 *file);


/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function allows for energy evaluation of a given sequence/structure pair where
 *  the structure is provided in pair_table format as obtained from vrna_ptable().
 *  Model details, energy parameters, and possibly soft constraints are used as provided
 *  via the parameter 'fc'. The fold_compound does not need to contain any DP matrices,
 *  but all the most basic init values as one would get from a call like this:
 *  @code{.c}
 * fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_EVAL_ONLY);
 *  @endcode
 *  In contrast to vrna_eval_structure_pt() this function prints detailed energy contributions
 *  based on individual loops to a file handle. If NULL is passed as file handle, this function
 *  defaults to print to stdout. Any positive @p verbosity_level activates potential warning message
 *  of the energy evaluting functions, while values @f$ \ge 1 @f$ allow for detailed control of what
 *  data is printed. A negative parameter @p verbosity_level turns off printing all together.
 *
 *  @see vrna_ptable(), vrna_eval_structure_pt(), vrna_eval_structure_verbose()
 *
 *  @param fc               A vrna_fold_compound_t containing the energy parameters and model details
 *  @param pt               Secondary structure as pair_table
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in 10cal/mol
 */
int
vrna_eval_structure_pt_v(vrna_fold_compound_t *fc,
                         const short          *pt,
                         int                  verbosity_level,
                         FILE                 *file);


/* End basic eval interface with pair table */
/**@}*/


/**
 *  @name Simplified Energy Evaluation with Sequence and Dot-Bracket Strings
 *  @{
 */

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
float
vrna_eval_structure_simple(const char *string,
                           const char *structure);


/**
 *  @brief  Evaluate the free energy of a sequence/structure pair where the sequence is circular
 *
 *  @see  vrna_eval_structure_simple(), vrna_eval_gquad_structure(), vrna_eval_circ_consensus_structure(),
 *        vrna_eval_circ_structure_v(), vrna_eval_structure()
 *
 *  @param  string    RNA sequence in uppercase letters
 *  @param  structure Secondary structure in dot-bracket notation
 *  @return           The free energy of the structure given the circular input sequence in kcal/mol
 */
float
vrna_eval_circ_structure(const char *string,
                         const char *structure);


/**
 *  @brief  Evaluate the free energy of a sequence/structure pair where the structure may contain G-Quadruplexes
 *
 *  G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences must
 *  be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer G-quadruplex:
 *  @code{.unparsed}
 *  GGAAGGAAAGGAGG
 *  ++..++...++.++
 *  @endcode
 *
 *  @see  vrna_eval_structure_simple(), vrna_eval_circ_structure(), vrna_eval_gquad_consensus_structure(),
 *        vrna_eval_gquad_structure_v(), vrna_eval_structure()
 *
 *  @param  string    RNA sequence in uppercase letters
 *  @param  structure Secondary structure in dot-bracket notation
 *  @return           The free energy of the structure including contributions of G-quadruplexes in kcal/mol
 */
float
vrna_eval_gquad_structure(const char  *string,
                          const char  *structure);


/**
 *  @brief  Evaluate the free energy of a sequence/structure pair where the sequence is circular and
 *          the structure may contain G-Quadruplexes
 *
 *  G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences must
 *  be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer G-quadruplex:
 *  @code{.unparsed}
 *  GGAAGGAAAGGAGG
 *  ++..++...++.++
 *  @endcode
 *
 *  @see  vrna_eval_structure_simple(), vrna_eval_circ_gquad_consensus_structure(),
 *        vrna_eval_circ_gquad_structure_v(), vrna_eval_structure()
 *
 *  @param  string    RNA sequence in uppercase letters
 *  @param  structure Secondary structure in dot-bracket notation
 *  @return           The free energy of the structure including contributions of G-quadruplexes in kcal/mol
 */
float
vrna_eval_circ_gquad_structure(const char *string,
                               const char *structure);


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
float
vrna_eval_structure_simple_verbose(const char *string,
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
 *  @see vrna_eval_structure_verbose(), vrna_eval_structure_pt(), vrna_eval_structure_pt_verbose(),
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float
vrna_eval_structure_simple_v(const char *string,
                             const char *structure,
                             int        verbosity_level,
                             FILE       *file);


/**
 *  @brief  Evaluate free energy of a sequence/structure pair, assume sequence to be circular and
 *          print contributions per loop
 *
 *  This function is the same as vrna_eval_structure_simple_v() but assumes the input sequence
 *  to be circularized.
 *
 *  @see  vrna_eval_structure_simple_v(), vrna_eval_circ_structure(), vrna_eval_structure_verbose()
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float
vrna_eval_circ_structure_v(const char *string,
                           const char *structure,
                           int        verbosity_level,
                           FILE       *file);


/**
 *  @brief  Evaluate free energy of a sequence/structure pair, allow for G-Quadruplexes in the structure
 *          and print contributions per loop
 *
 *  This function is the same as vrna_eval_structure_simple_v() but allows for annotated G-Quadruplexes
 *  in the dot-bracket structure input.
 *
 *  G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences must
 *  be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer G-quadruplex:
 *  @code{.unparsed}
 *  GGAAGGAAAGGAGG
 *  ++..++...++.++
 *  @endcode
 *
 *  @see  vrna_eval_structure_simple_v(), vrna_eval_gquad_structure(), vrna_eval_structure_verbose()
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float
vrna_eval_gquad_structure_v(const char  *string,
                            const char  *structure,
                            int         verbosity_level,
                            FILE        *file);


/**
 *  @brief  Evaluate free energy of a sequence/structure pair, assume sequence to be circular, allow
 *          for G-Quadruplexes in the structure, and print contributions per loop
 *
 *  This function is the same as vrna_eval_structure_simple_v() but assumes the input sequence to
 *  be circular and allows for annotated G-Quadruplexes in the dot-bracket structure input.
 *
 *  G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences must
 *  be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer G-quadruplex:
 *  @code{.unparsed}
 *  GGAAGGAAAGGAGG
 *  ++..++...++.++
 *  @endcode
 *
 *  @param string           RNA sequence in uppercase letters
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
float
vrna_eval_circ_gquad_structure_v(const char *string,
                                 const char *structure,
                                 int        verbosity_level,
                                 FILE       *file);


/* End simplified eval interface */
/**@}*/


/**
 *  @name Simplified Energy Evaluation with Sequence Alignments and Consensus Structure Dot-Bracket String
 *  @{
 */

/**
 *  @brief Calculate the free energy of an already folded RNA sequence alignment
 *
 *  This function allows for energy evaluation for a given multiple sequence alignment
 *  and consensus structure pair.
 *  In contrast to vrna_eval_structure() this function assumes default model details
 *  and default energy parameters in order to evaluate the free energy of the secondary
 *  structure. Therefore, it serves as a simple interface function for energy evaluation
 *  for situations where no changes on the energy model are required.
 *
 *  @note The free energy returned from this function already includes the covariation
 *        pseudo energies that is used fir comparative structure prediction within this
 *        library.
 *
 *  @see  vrna_eval_covar_structure(), vrna_eval_structure(), vrna_eval_structure_pt(),
 *        vrna_eval_structure_verbose(), vrna_eval_structure_pt_verbose()
 *
 *  @param alignment        RNA sequence alignment in uppercase letters and hyphen ('-') to denote gaps
 *  @param structure        Consensus Secondary structure in dot-bracket notation
 *  @return                 The free energy of the consensus structure given the input alignment in kcal/mol
 */
float
vrna_eval_consensus_structure_simple(const char **alignment,
                                     const char *structure);


/**
 *  @brief  Evaluate the free energy of a multiple sequence alignment/consensus structure pair
 *          where the sequences are circular
 *
 *  @note The free energy returned from this function already includes the covariation
 *        pseudo energies that is used fir comparative structure prediction within this
 *        library.
 *
 *  @see  vrna_eval_covar_structure(), vrna_eval_consensus_structure_simple(), vrna_eval_gquad_consensus_structure(),
 *        vrna_eval_circ_structure(), vrna_eval_circ_consensus_structure_v(), vrna_eval_structure()
 *
 *  @param  alignment RNA sequence alignment in uppercase letters
 *  @param  structure Consensus secondary structure in dot-bracket notation
 *  @return           The free energy of the consensus structure given the circular input sequence in kcal/mol
 */
float
vrna_eval_circ_consensus_structure(const char **alignment,
                                   const char *structure);


/**
 *  @brief  Evaluate the free energy of a multiple sequence alignment/consensus structure pair
 *          where the structure may contain G-Quadruplexes
 *
 *  G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences must
 *  be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer G-quadruplex:
 *  @code{.unparsed}
 *  GGAAGGAAAGGAGG
 *  ++..++...++.++
 *  @endcode
 *
 *  @note The free energy returned from this function already includes the covariation
 *        pseudo energies that is used fir comparative structure prediction within this
 *        library.
 *
 *  @see  vrna_eval_covar_structure(), vrna_eval_consensus_structure_simple(), vrna_eval_circ_consensus_structure(),
 *        vrna_eval_gquad_structure(), vrna_eval_gquad_consensus_structure_v(), vrna_eval_structure()
 *
 *  @param  alignment RNA sequence alignment in uppercase letters
 *  @param  structure Consensus secondary structure in dot-bracket notation
 *  @return           The free energy of the consensus structure including contributions of G-quadruplexes in kcal/mol
 */
float
vrna_eval_gquad_consensus_structure(const char  **alignment,
                                    const char  *structure);


/**
 *  @brief  Evaluate the free energy of a multiple sequence alignment/consensus structure pair
 *          where the sequence is circular and the structure may contain G-Quadruplexes
 *
 *  G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences must
 *  be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer G-quadruplex:
 *  @code{.unparsed}
 *  GGAAGGAAAGGAGG
 *  ++..++...++.++
 *  @endcode
 *
 *  @note The free energy returned from this function already includes the covariation
 *        pseudo energies that is used fir comparative structure prediction within this
 *        library.
 *
 *  @see  vrna_eval_covar_structure(), vrna_eval_consensus_structure_simple(), vrna_eval_circ_consensus_structure(),
 *        vrna_eval_gquad_structure(), vrna_eval_circ_gquad_consensus_structure_v(), vrna_eval_structure()
 *
 *  @param  alignment RNA sequence alignment in uppercase letters
 *  @param  structure Consensus secondary structure in dot-bracket notation
 *  @return           The free energy of the consensus structure including contributions of G-quadruplexes in kcal/mol
 */
float
vrna_eval_circ_gquad_consensus_structure(const char **alignment,
                                         const char *structure);


/**
 *  @brief  Evaluate the free energy of a consensus structure for an RNA sequence alignment and print
 *          contributions per loop.
 *
 *  This function is a simplyfied version of vrna_eval_consensus_structure_simple_v() that uses the
 *  @em default verbosity level.
 *
 *  @note The free energy returned from this function already includes the covariation
 *        pseudo energies that is used fir comparative structure prediction within this
 *        library.
 *
 *  @see  vrna_eval_consensus_structure_simple_v(), vrna_eval_structure_verbose(), vrna_eval_structure_pt(),
 *        vrna_eval_structure_pt_verbose()
 *
 *  @param alignment        RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')
 *  @param structure        Consensus secondary structure in dot-bracket notation
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the conensus structure given the aligned input sequences in kcal/mol
 */
float
vrna_eval_consensus_structure_simple_verbose(const char **alignment,
                                             const char *structure,
                                             FILE       *file);


/**
 *  @brief  Evaluate the free energy of a consensus structure for an RNA sequence alignment and print
 *          contributions per loop.
 *
 *  This function allows for detailed energy evaluation of a given sequence alignment/consensus
 *  structure pair. In contrast to vrna_eval_consensus_structure_simple() this function prints
 *  detailed energy contributions based on individual loops to a file handle. If NULL is passed
 *  as file handle, this function defaults to print to stdout. Any positive @p verbosity_level
 *  activates potential warning message of the energy evaluting functions, while values @f$ \ge 1 @f$
 *  allow for detailed control of what data is printed. A negative parameter @p verbosity_level
 *  turns off printing all together.
 *
 *  @note The free energy returned from this function already includes the covariation
 *        pseudo energies that is used fir comparative structure prediction within this
 *        library.
 *
 *  @see vrna_eval_consensus_structure(), vrna_eval_structure()
 *
 *  @param alignment        RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')
 *  @param structure        Consensus secondary structure in dot-bracket notation
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the consensus structure given the sequence alignment in kcal/mol
 */
float
vrna_eval_consensus_structure_simple_v(const char **alignment,
                                       const char *structure,
                                       int        verbosity_level,
                                       FILE       *file);


/**
 *  @brief  Evaluate the free energy of a consensus structure for an alignment of circular RNA sequences
 *          and print contributions per loop.
 *
 *  This function is identical with vrna_eval_consensus_structure_simple_v() but assumed the
 *  aligned sequences to be circular.
 *
 *  @note The free energy returned from this function already includes the covariation
 *        pseudo energies that is used fir comparative structure prediction within this
 *        library.
 *
 *  @see vrna_eval_consensus_structure_simple_v(), vrna_eval_circ_consensus_structure(), vrna_eval_structure()
 *
 *  @param alignment        RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')
 *  @param structure        Consensus secondary structure in dot-bracket notation
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the consensus structure given the sequence alignment in kcal/mol
 */
float
vrna_eval_circ_consensus_structure_v(const char **alignment,
                                     const char *structure,
                                     int        verbosity_level,
                                     FILE       *file);


/**
 *  @brief  Evaluate the free energy of a consensus structure for an RNA sequence alignment, allow for
 *          annotated G-Quadruplexes in the structure and print contributions per loop.
 *
 *  This function is identical with vrna_eval_consensus_structure_simple_v() but allows for annotated
 *  G-Quadruplexes in the consensus structure.
 *
 *  G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences must
 *  be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer G-quadruplex:
 *  @code{.unparsed}
 *  GGAAGGAAAGGAGG
 *  ++..++...++.++
 *  @endcode
 *
 *  @note The free energy returned from this function already includes the covariation
 *        pseudo energies that is used fir comparative structure prediction within this
 *        library.
 *
 *  @see vrna_eval_consensus_structure_simple_v(), vrna_eval_gquad_consensus_structure(), vrna_eval_structure()
 *
 *  @param alignment        RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')
 *  @param structure        Consensus secondary structure in dot-bracket notation
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the consensus structure given the sequence alignment in kcal/mol
 */
float
vrna_eval_gquad_consensus_structure_v(const char  **alignment,
                                      const char  *structure,
                                      int         verbosity_level,
                                      FILE        *file);


/**
 *  @brief  Evaluate the free energy of a consensus structure for an alignment of circular RNA sequences,
 *          allow for annotated G-Quadruplexes in the structure and print contributions per loop.
 *
 *  This function is identical with vrna_eval_consensus_structure_simple_v() but assumes the sequences in
 *  the alignment to be circular and allows for annotated G-Quadruplexes in the consensus structure.
 *
 *  G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences must
 *  be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer G-quadruplex:
 *  @code{.unparsed}
 *  GGAAGGAAAGGAGG
 *  ++..++...++.++
 *  @endcode
 *
 *  @note The free energy returned from this function already includes the covariation
 *        pseudo energies that is used fir comparative structure prediction within this
 *        library.
 *
 *  @see vrna_eval_consensus_structure_simple_v(), vrna_eval_circ_gquad_consensus_structure(), vrna_eval_structure()
 *
 *  @param alignment        RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')
 *  @param structure        Consensus secondary structure in dot-bracket notation
 *  @param verbosity_level  The level of verbosity of this function
 *  @param file             A file handle where this function should print to (may be NULL).
 *  @return                 The free energy of the consensus structure given the sequence alignment in kcal/mol
 */
float
vrna_eval_circ_gquad_consensus_structure_v(const char **alignment,
                                           const char *structure,
                                           int        verbosity_level,
                                           FILE       *file);


/* End simplified comparative eval interface */
/**@}*/


/**
 *  @name Simplified Energy Evaluation with Sequence String and Structure Pair Table
 *  @{
 */

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
int
vrna_eval_structure_pt_simple(const char  *string,
                              const short *pt);


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
int
vrna_eval_structure_pt_simple_verbose(const char  *string,
                                      const short *pt,
                                      FILE        *file);


/**
 *  @brief Calculate the free energy of an already folded RNA
 *
 *  This function allows for energy evaluation of a given sequence/structure pair where
 *  the structure is provided in pair_table format as obtained from vrna_ptable().
 *  Model details, energy parameters, and possibly soft constraints are used as provided
 *  via the parameter 'fc'. The fold_compound does not need to contain any DP matrices,
 *  but all the most basic init values as one would get from a call like this:
 *  @code{.c}
 * fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_EVAL_ONLY);
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
int
vrna_eval_structure_pt_simple_v(const char  *string,
                                const short *pt,
                                int         verbosity_level,
                                FILE        *file);


/* End simplified eval interface with pair table */
/**@}*/

/**
 *  @name Simplified Energy Evaluation with Sequence Alignment and Consensus Structure Pair Table
 *  @{
 */

/**
 *  @brief  Evaluate the Free Energy of a Consensus Secondary Structure given a Sequence Alignment
 *
 *  @note The free energy returned from this function already includes the covariation
 *        pseudo energies that is used fir comparative structure prediction within this
 *        library.
 *
 *  @see  vrna_eval_consensus_structure_simple(), vrna_eval_structure_pt(), vrna_eval_structure(),
 *        vrna_eval_covar_structure()
 *
 *  @param  alignment   RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')
 *  @param  pt          Secondary structure in pair table format
 *  @return             Free energy of the consensus structure in 10cal/mol
 */
int
vrna_eval_consensus_structure_pt_simple(const char  **alignment,
                                        const short *pt);


int
vrna_eval_consensus_structure_pt_simple_verbose(const char  **alignment,
                                                const short *pt,
                                                FILE        *file);


int
vrna_eval_consensus_structure_pt_simple_v(const char  **alignment,
                                          const short *pt,
                                          int         verbosity_level,
                                          FILE        *file);


/* End simplified eval interface with pair table */
/**@}*/


/**
 * @}
 */


/**
 *  @addtogroup eval_loops
 *  @{
 *  @brief  Functions to evaluate the free energy of particular types of loops
 *
 */

/**
 * @brief Calculate energy of a loop
 *
 *  @param fc         A vrna_fold_compound_t containing the energy parameters and model details
 *  @param i          position of covering base pair
 *  @param pt         the pair table of the secondary structure
 *  @returns          free energy of the loop in 10cal/mol
 */
int
vrna_eval_loop_pt(vrna_fold_compound_t  *fc,
                  int                   i,
                  const short           *pt);


/**
 * @brief Calculate energy of a loop
 *
 *  @param fc         A vrna_fold_compound_t containing the energy parameters and model details
 *  @param i          position of covering base pair
 *  @param pt         the pair table of the secondary structure
 *  @param verbosity_level  The level of verbosity of this function
 *  @returns          free energy of the loop in 10cal/mol
 */
int
vrna_eval_loop_pt_v(vrna_fold_compound_t  *fc,
                    int                   i,
                    const short           *pt,
                    int                   verbosity_level);


/**
 * @}
 */


/**
 *  @addtogroup eval_move
 *  @{
 */

/**
 * @brief Calculate energy of a move (closing or opening of a base pair)
 *
 *  If the parameters m1 and m2 are negative, it is deletion (opening)
 *  of a base pair, otherwise it is insertion (opening).
 *
 *  @see              vrna_eval_move_pt()
 *
 *  @param fc         A vrna_fold_compound_t containing the energy parameters and model details
 *  @param structure  secondary structure in dot-bracket notation
 *  @param m1         first coordinate of base pair
 *  @param m2         second coordinate of base pair
 *  @returns          energy change of the move in kcal/mol (#INF / 100. upon any error)
 */
float
vrna_eval_move(vrna_fold_compound_t *fc,
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
 *
 *  @param fc         A vrna_fold_compound_t containing the energy parameters and model details
 *  @param pt         the pair table of the secondary structure
 *  @param m1         first coordinate of base pair
 *  @param m2         second coordinate of base pair
 *  @returns          energy change of the move in 10cal/mol
 */
int
vrna_eval_move_pt(vrna_fold_compound_t  *fc,
                  short                 *pt,
                  int                   m1,
                  int                   m2);


int
vrna_eval_move_pt_simple(const char *string,
                         short      *pt,
                         int        m1,
                         int        m2);


int
vrna_eval_move_shift_pt(vrna_fold_compound_t  *fc,
                        vrna_move_t           *m,
                        short                 *structure);


/**
 * @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup eval_deprecated
 *  @{
 */

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
DEPRECATED(float
           energy_of_structure(const char *string,
                               const char *structure,
                               int        verbosity_level),
           "Use vrna_eval_structure_simple() and vrna_eval_structure() instead");

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
DEPRECATED(float
           energy_of_struct_par(const char    *string,
                                const char    *structure,
                                vrna_param_t  *parameters,
                                int           verbosity_level),
           "Use vrna_eval_structure() instead");

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
DEPRECATED(float
           energy_of_circ_structure(const char  *string,
                                    const char  *structure,
                                    int         verbosity_level),
           "Use vrna_eval_circ_structure_simple() and vrna_eval_structure() instead");

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
DEPRECATED(float
           energy_of_circ_struct_par(const char   *string,
                                     const char   *structure,
                                     vrna_param_t *parameters,
                                     int          verbosity_level),
           "Use vrna_eval_structure() instead");


DEPRECATED(float
           energy_of_gquad_structure(const char *string,
                                     const char *structure,
                                     int        verbosity_level),
           "Use vrna_eval_structure_simple() instead");

DEPRECATED(float
           energy_of_gquad_struct_par(const char    *string,
                                      const char    *structure,
                                      vrna_param_t  *parameters,
                                      int           verbosity_level),
           "Use vrna_eval_structure() instead");


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
DEPRECATED(int
           energy_of_structure_pt(const char  *string,
                                  short       *ptable,
                                  short       *s,
                                  short       *s1,
                                  int         verbosity_level),
           "Use vrna_eval_structure_pt_simple() and vrna_eval_structure_pt() instead");

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
DEPRECATED(int
           energy_of_struct_pt_par(const char   *string,
                                   short        *ptable,
                                   short        *s,
                                   short        *s1,
                                   vrna_param_t *parameters,
                                   int          verbosity_level),
           "Use vrna_eval_structure_pt() instead");


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
DEPRECATED(float
           energy_of_move(const char  *string,
                          const char  *structure,
                          int         m1,
                          int         m2),
           "Use vrna_eval_move() instead");


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
DEPRECATED(int
           energy_of_move_pt(short  *pt,
                             short  *s,
                             short  *s1,
                             int    m1,
                             int    m2),
           "Use vrna_eval_move_pt_simple() and vrna_eval_move_pt() instead");

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
DEPRECATED(int
           loop_energy(short  *ptable,
                       short  *s,
                       short  *s1,
                       int    i),
           "Use vrna_eval_loop_pt() instead");

/**
 *  Calculate the free energy of an already folded RNA
 *
 *  @note This function is not entirely threadsafe! Depending on the state of the global
 *        variable @ref eos_debug it prints energy information to stdout or not...\n
 *
 *  @deprecated This function is deprecated and should not be used in future programs!
 *  Use @ref energy_of_structure() instead!
 *
 *  @see              energy_of_structure, energy_of_circ_struct(), energy_of_struct_pt()
 *
 *  @param string     RNA sequence
 *  @param structure  secondary structure in dot-bracket notation
 *  @return          the free energy of the input structure given the input sequence in kcal/mol
 */
DEPRECATED(float
           energy_of_struct(const char  *string,
                            const char  *structure),
           "Use vrna_eval_structure_simple() instead");

/**
 *  Calculate the free energy of an already folded RNA
 *
 *  @note This function is not entirely threadsafe! Depending on the state of the global
 *        variable @ref eos_debug it prints energy information to stdout or not...\n
 *
 *  @deprecated This function is deprecated and should not be used in future programs!
 *  Use @ref energy_of_structure_pt() instead!
 *
 *  @see              make_pair_table(), energy_of_structure()
 *
 *  @param string     RNA sequence
 *  @param ptable     the pair table of the secondary structure
 *  @param s          encoded RNA sequence
 *  @param s1         encoded RNA sequence
 *  @return          the free energy of the input structure given the input sequence in 10kcal/mol
 */
DEPRECATED(int
           energy_of_struct_pt(const char *string,
                               short      *ptable,
                               short      *s,
                               short      *s1),
           "Use vrna_eval_structure_pt_simple() instead");

/**
 *  Calculate the free energy of an already folded  circular RNA
 *
 *  @note This function is not entirely threadsafe! Depending on the state of the global
 *        variable @ref eos_debug it prints energy information to stdout or not...\n
 *
 *  @deprecated This function is deprecated and should not be used in future programs
 *  Use @ref energy_of_circ_structure() instead!
 *
 *  @see              energy_of_circ_structure(), energy_of_struct(), energy_of_struct_pt()
 *
 *  @param string     RNA sequence
 *  @param structure  secondary structure in dot-bracket notation
 *  @return          the free energy of the input structure given the input sequence in kcal/mol
 */
DEPRECATED(float
           energy_of_circ_struct(const char *string,
                                 const char *structure),
           "Use vrna_eval_circ_structure_simple() and vrna_eval_structure() instead");

#endif

/**
 * @}
 */

#endif
