#ifndef __VIENNA_RNA_PACKAGE_FOLD_H__
#define __VIENNA_RNA_PACKAGE_FOLD_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  \addtogroup mfe_fold
 *  \ingroup folding_routines
 *  \brief This section covers all functions and variables related to the calculation
 *  of minimum free energy (MFE) structures.
 *
 *  The library provides a fast dynamic programming minimum free energy
 *  folding algorithm as described in \cite zuker:1981.
 *  All relevant parts that directly implement the "Zuker & Stiegler" algorithm for single
 *  sequences are described in this section.
 *
 *  Folding of circular RNA sequences is handled as a post-processing step of the forward
 *  recursions. See \cite hofacker:2006 for further details.
 *
 *  Nevertheless, the RNAlib also
 *  provides interfaces for the prediction of consensus MFE structures of sequence alignments,
 *  MFE structure for two hybridized sequences, local optimal structures and many more. For
 *  those more specialized variants of MFE folding routines, please consult the appropriate
 *  subsections (Modules) as listed above.
 *  
 *  \file fold.h
 *  \brief MFE calculations and energy evaluations for single RNA sequences
 * 
 *  This file includes (almost) all function declarations within the RNAlib that are related to
 *  MFE folding...
 */

/**
 *  \defgroup eval Energy evaluation
 *  @{
 *    \brief This module contains all functions and variables related to energy evaluation
 *    of sequence/structure pairs.
 *
 *
 *  @}
 */

/**
 *  \defgroup mfe_fold Calculating Minimum Free Energy Structures
 *  @{
 *    \brief This module contains all functions and variables related to the calculation
 *    of global minimum free energy structures for single sequences.
 *
 *    The library provides a fast dynamic programming minimum free energy
 *    folding algorithm as described by \ref zuker_81 "Zuker & Stiegler (1981)".
 *  @}
 */

/** \brief if nonzero use logarithmic ML energy in energy_of_struct  */
extern  int logML;

/** \brief do ML decomposition uniquely (for subopt)  */
extern  int uniq_ML;

/** \brief set to first pos of second seq for cofolding  */
extern  int cut_point;

/**
 *  \brief verbose info from energy_of_struct
 *  \ingroup eval
 */
extern  int eos_debug;


/**
 *  \brief Compute minimum free energy and an appropriate secondary
 *  structure of an RNA sequence
 * 
 *  The first parameter given, the RNA sequence, must be \a uppercase and should only contain
 *  an alphabet \f$\Sigma\f$ that is understood by the RNAlib\n
 *  (e.g. \f$ \Sigma = \{A,U,C,G\} \f$)\n
 *
 *  The second parameter, \a structure, must always point to an allocated
 *  block of memory with a size of at least \f$\mathrm{strlen}(\mathrm{sequence})+1\f$
 *
 *  If the third parameter is NULL, global model detail settings are assumed for the folding
 *  recursions. Otherwise, the provided parameters are used.
 *
 *  The fourth parameter indicates whether a secondary structure constraint in enhanced dot-bracket
 *  notation is passed through the structure parameter or not. If so, the characters " | x < > " are
 *  recognized to mark bases that are paired, unpaired, paired upstream, or downstream, respectively.
 *  Matching brackets " ( ) " denote base pairs, dots "." are used for unconstrained bases.
 *
 *  To indicate that the RNA sequence is circular and thus has to be post-processed, set the last
 *  parameter to non-zero
 *
 *  After a successful call of fold_par(), a backtracked secondary structure (in dot-bracket notation)
 *  that exhibits the minimum of free energy will be written to the memory \a structure is pointing to.
 *  The function returns the minimum of free energy for any fold of the sequence given.
 *
 *  \note OpenMP: Passing NULL to the 'parameters' argument involves access to several global model
 *        detail variables and thus is not to be considered threadsafe
 *
 *  \ingroup mfe_fold
 *
 *  \see fold(), circfold(), #model_detailsT, set_energy_model(), get_scaled_parameters()
 *
 *  \param sequence       RNA sequence
 *  \param structure      A pointer to the character array where the
 *                        secondary structure in dot-bracket notation will be written to
 *  \param parameters     A data structure containing the prescaled energy contributions
 *                        and the model details. (NULL may be passed, see OpenMP notes above)
 *  \param is_constrained Switch to indicate that a structure contraint is passed via the structure argument (0==off)
 *  \param is_circular    Switch to (de-)activate postprocessing steps in case RNA sequence is circular (0==off)
 *
 *  \return the minimum free energy (MFE) in kcal/mol
 */
float fold_par( const char *sequence,
                char *structure,
                paramT *parameters,
                int is_constrained,
                int is_circular);

/**
 *  \brief Compute minimum free energy and an appropriate secondary structure of an RNA sequence
 *
 *  This function essentially does the same thing as fold_par(). However, it takes its model details,
 *  i.e. #temperature, #dangles, #tetra_loop, #noGU, #no_closingGU, #fold_constrained, #noLonelyPairs
 *  from the current global settings within the library
 *
 *  Use fold_par() for a completely threadsafe variant
 *
 *  \ingroup mfe_fold
 *
 *  \see fold_par(), circfold()
 *
 *  \param sequence RNA sequence
 *  \param structure A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  \return the minimum free energy (MFE) in kcal/mol
 */
float fold( const char *sequence,
            char *structure);

/**
 *  \brief Compute minimum free energy and an appropriate secondary structure of a circular RNA sequence
 *
 *  This function essentially does the same thing as fold_par(). However, it takes its model details,
 *  i.e. #temperature, #dangles, #tetra_loop, #noGU, #no_closingGU, #fold_constrained, #noLonelyPairs
 *  from the current global settings within the library
 *
 *  Use fold_par() for a completely threadsafe variant
 *
 *  \ingroup mfe_fold
 *
 *  \see fold_par(), circfold()
 *
 *  \param sequence RNA sequence
 *  \param structure A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  \return the minimum free energy (MFE) in kcal/mol
 */
float circfold( const char *sequence,
                char *structure);


/**
 *  \addtogroup eval Energy evaluation
 *  \ingroup folding_routines
 *  @{
 *    \brief This module contains all functions and variables related to energy evaluation
 *    of sequence/structure pairs.
 *  @}
 */

/**
 *  \brief Calculate the free energy of an already folded RNA using global model detail settings
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  \note OpenMP: This function relies on several global model settings variables and thus is
 *        not to be considered threadsafe. See energy_of_struct_par() for a completely threadsafe
 *        implementation.
 *
 *  \ingroup eval
 *
 *  \see energy_of_struct_par(), energy_of_circ_structure()
 *
 *  \param string     RNA sequence
 *  \param structure  secondary structure in dot-bracket notation
 *  \param verbosity_level a flag to turn verbose output on/off
 *  \return          the free energy of the input structure given the input sequence in kcal/mol
 */
float energy_of_structure(const char *string,
                          const char *structure,
                          int verbosity_level);

/**
 *  \brief Calculate the free energy of an already folded RNA
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  \ingroup eval
 *
 *  \see energy_of_circ_structure(), energy_of_structure_pt(), get_scaled_parameters()
 *
 *  \param string           RNA sequence in uppercase letters
 *  \param structure        Secondary structure in dot-bracket notation
 *  \param parameters       A data structure containing the prescaled energy contributions and the model details.
 *  \param verbosity_level  A flag to turn verbose output on/off
 *  \return                The free energy of the input structure given the input sequence in kcal/mol
 */
float energy_of_struct_par( const char *string,
                            const char *structure,
                            paramT *parameters,
                            int verbosity_level);

/**
 *  \brief Calculate the free energy of an already folded  circular RNA
 *
 *  \note OpenMP: This function relies on several global model settings variables and thus is
 *        not to be considered threadsafe. See energy_of_circ_struct_par() for a completely threadsafe
 *        implementation.
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  \ingroup eval
 *
 *  \see energy_of_circ_struct_par(), energy_of_struct_par()
 *
 *  \param string           RNA sequence
 *  \param structure        Secondary structure in dot-bracket notation
 *  \param verbosity_level  A flag to turn verbose output on/off
 *  \return                The free energy of the input structure given the input sequence in kcal/mol
 */
float energy_of_circ_structure( const char *string,
                                const char *structure,
                                int verbosity_level);

/**
 *  \brief Calculate the free energy of an already folded circular RNA
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  \ingroup eval
 *
 *  \see energy_of_struct_par(), get_scaled_parameters()
 *
 *  \param string           RNA sequence
 *  \param structure        Secondary structure in dot-bracket notation
 *  \param parameters       A data structure containing the prescaled energy contributions and the model details.
 *  \param verbosity_level  A flag to turn verbose output on/off
 *  \return                The free energy of the input structure given the input sequence in kcal/mol
 */
float energy_of_circ_struct_par(const char *string,
                                const char *structure,
                                paramT *parameters,
                                int verbosity_level);


float energy_of_gquad_structure(const char *string,
                                const char *structure,
                                int verbosity_level);

float energy_of_gquad_struct_par( const char *string,
                                  const char *structure,
                                  paramT *parameters,
                                  int verbosity_level);

/**
 *  \brief Calculate the free energy of an already folded RNA
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  \note OpenMP: This function relies on several global model settings variables and thus is
 *        not to be considered threadsafe. See energy_of_struct_pt_par() for a completely threadsafe
 *        implementation.
 *
 *  \ingroup eval
 *
 *  \see make_pair_table(), energy_of_struct_pt_par()
 *
 *  \param string     RNA sequence
 *  \param ptable     the pair table of the secondary structure
 *  \param s          encoded RNA sequence
 *  \param s1         encoded RNA sequence
 *  \param verbosity_level a flag to turn verbose output on/off
 *  \return          the free energy of the input structure given the input sequence in 10kcal/mol
 */
int energy_of_structure_pt( const char *string,
                            short *ptable,
                            short *s,
                            short *s1,
                            int verbosity_level);

/**
 *  \brief Calculate the free energy of an already folded RNA
 *
 *  If verbosity level is set to a value >0, energies of structure elements are printed to stdout
 *
 *  \ingroup eval
 *
 *  \see make_pair_table(), energy_of_struct_par(), get_scaled_parameters()
 *
 *  \param string           RNA sequence in uppercase letters
 *  \param ptable           The pair table of the secondary structure
 *  \param s                Encoded RNA sequence
 *  \param s1               Encoded RNA sequence
 *  \param parameters       A data structure containing the prescaled energy contributions and the model details.
 *  \param verbosity_level  A flag to turn verbose output on/off
 *  \return                The free energy of the input structure given the input sequence in 10kcal/mol
 */
int energy_of_struct_pt_par(const char *string,
                            short *ptable,
                            short *s,
                            short *s1,
                            paramT *parameters,
                            int verbosity_level);

/**
 *  \brief Free arrays for mfe folding
 *
 *  \ingroup mfe_fold
 *
 */
void  free_arrays(void);


/**
 *  \brief Create a dot-backet/parenthesis structure from backtracking stack
 * 
 *  \note This function is threadsafe
 */
void  parenthesis_structure(char *structure,
                            bondT *bp,
                            int length);

/**
 *  \brief Create a dot-backet/parenthesis structure from backtracking stack
 *  obtained by zuker suboptimal calculation in cofold.c
 * 
 *  \note This function is threadsafe
 */
void parenthesis_zuker( char *structure,
                        bondT *bp,
                        int length);

void letter_structure(char *structure,
                      bondT *bp,
                      int length);


/**
 *  \brief Recalculate energy parameters
 *
 *  \ingroup mfe_fold
 */
void  update_fold_params(void);

/**
 *
 *  \ingroup mfe_fold
 * 
 */
void update_fold_params_par(paramT *parameters);

/**
 *
 *  \ingroup mfe_fold
 * 
 */
char  *backtrack_fold_from_pair(char *sequence,
                                int i,
                                int j);

/** 
 * \brief Calculate energy of a move (closing or opening of a base pair)
 *
 *  If the parameters m1 and m2 are negative, it is deletion (opening)
 *  of a base pair, otherwise it is insertion (opening).
 *
 *  \see              make_pair_table(), energy_of_move()
 *  \param string     RNA sequence
 *  \param structure  secondary structure in dot-bracket notation
 *  \param m1         first coordinate of base pair
 *  \param m2         second coordinate of base pair
 *  \returns          energy change of the move in kcal/mol
 */
float energy_of_move( const char *string,
                      const char *structure,
                      int m1,
                      int m2);


/**
 * 
 * \brief Calculate energy of a move (closing or opening of a base pair)
 *
 *  If the parameters m1 and m2 are negative, it is deletion (opening)
 *  of a base pair, otherwise it is insertion (opening).
 *
 *  \see              make_pair_table(), energy_of_move()
 *  \param pt         the pair table of the secondary structure
 *  \param s          encoded RNA sequence
 *  \param s1         encoded RNA sequence
 *  \param m1         first coordinate of base pair
 *  \param m2         second coordinate of base pair
 *  \returns          energy change of the move in 10cal/mol
 */
int energy_of_move_pt(short *pt,
                   short *s,
                   short *s1,
                   int m1,
                   int m2);

/**
 * \brief Calculate energy of a loop
 *
 *  \param ptable     the pair table of the secondary structure
 *  \param s          encoded RNA sequence
 *  \param s1         encoded RNA sequence
 *  \param i          position of covering base pair
 *  \returns          free energy of the loop in 10cal/mol
 */
int   loop_energy(short *ptable,
                  short *s,
                  short *s1,
                  int i);

/**
 *
 *  \ingroup mfe_fold
 * 
 */
void export_fold_arrays(int **f5_p,
                        int **c_p,
                        int **fML_p,
                        int **fM1_p,
                        int **indx_p,
                        char **ptype_p);

/**
 *
 *  \ingroup mfe_fold
 * 
 */
void export_fold_arrays_par(int **f5_p,
                            int **c_p,
                            int **fML_p,
                            int **fM1_p,
                            int **indx_p,
                            char **ptype_p,
                            paramT **P_p);

/**
 *
 *  \ingroup mfe_fold
 * 
 */
void export_circfold_arrays(int *Fc_p,
                            int *FcH_p,
                            int *FcI_p,
                            int *FcM_p,
                            int **fM2_p,
                            int **f5_p,
                            int **c_p,
                            int **fML_p,
                            int **fM1_p,
                            int **indx_p,
                            char **ptype_p);

/**
 *
 *  \ingroup mfe_fold
 * 
 */
void export_circfold_arrays_par(int *Fc_p,
                                int *FcH_p,
                                int *FcI_p,
                                int *FcM_p,
                                int **fM2_p,
                                int **f5_p,
                                int **c_p,
                                int **fML_p,
                                int **fM1_p,
                                int **indx_p,
                                char **ptype_p,
                                paramT **P_p);


/**
 *  \brief Create a plist from a dot-bracket string
 * 
 *  The dot-bracket string is parsed and for each base pair an
 *  entry in the plist is created. The probability of each pair in
 *  the list is set by a function parameter.
 * 
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 * 
 *  This function is threadsafe
 * 
 *  \param pl     A pointer to the plist that is to be created
 *  \param struc  The secondary structure in dot-bracket notation
 *  \param pr     The probability for each base pair
 */
void assign_plist_from_db(plist **pl,
                          const char *struc,
                          float pr);

/* finally moved the loop energy function declarations to this header...  */
/* BUT: The functions only exist for backward compatibility reasons!      */
/* You better include "loop_energies.h" and call the functions:           */
/* E_Hairpin() and E_IntLoop() which are (almost) threadsafe as they get  */
/* a pointer to the energy parameter datastructure as additional argument */

/**
 *  \deprecated {This function is deprecated and will be removed soon.
 *  Use \ref E_IntLoop() instead!}
 */
DEPRECATED(int LoopEnergy(int n1,
                          int n2,
                          int type,
                          int type_2,
                          int si1,
                          int sj1,
                          int sp1,
                          int sq1));

/**
 *  \deprecated {This function is deprecated and will be removed soon.
 *  Use \ref E_Hairpin() instead!}
 */
DEPRECATED(int HairpinE(int size,
                        int type,
                        int si1,
                        int sj1,
                        const char *string));

/**
 *  Allocate arrays for folding\n
 *  \deprecated {This function is deprecated and will be removed soon!}
 * 
 */
DEPRECATED(void initialize_fold(int length));

/**
 *  Calculate the free energy of an already folded RNA
 * 
 *  \note This function is not entirely threadsafe! Depending on the state of the global
 *  variable \ref eos_debug it prints energy information to stdout or not...\n
 * 
 *  \deprecated This function is deprecated and should not be used in future programs!
 *  Use \ref energy_of_structure() instead!
 * 
 *  \see              energy_of_structure, energy_of_circ_struct(), energy_of_struct_pt()
 *  \param string     RNA sequence
 *  \param structure  secondary structure in dot-bracket notation
 *  \return          the free energy of the input structure given the input sequence in kcal/mol
 */
DEPRECATED(float energy_of_struct(const char *string,
                                  const char *structure));

/**
 *  Calculate the free energy of an already folded RNA
 * 
 *  \note This function is not entirely threadsafe! Depending on the state of the global
 *  variable \ref eos_debug it prints energy information to stdout or not...\n
 * 
 *  \deprecated This function is deprecated and should not be used in future programs!
 *  Use \ref energy_of_structure_pt() instead!
 * 
 *  \see              make_pair_table(), energy_of_structure()
 *  \param string     RNA sequence
 *  \param ptable     the pair table of the secondary structure
 *  \param s          encoded RNA sequence
 *  \param s1         encoded RNA sequence
 *  \return          the free energy of the input structure given the input sequence in 10kcal/mol
 */
DEPRECATED(int energy_of_struct_pt( const char *string,
                                    short *ptable,
                                    short *s,
                                    short *s1));

/**
 *  Calculate the free energy of an already folded  circular RNA
 * 
 *  \note This function is not entirely threadsafe! Depending on the state of the global
 *  variable \ref eos_debug it prints energy information to stdout or not...\n
 * 
 *  \deprecated This function is deprecated and should not be used in future programs
 *  Use \ref energy_of_circ_structure() instead!
 * 
 *  \see              energy_of_circ_structure(), energy_of_struct(), energy_of_struct_pt()
 *  \param string     RNA sequence
 *  \param structure  secondary structure in dot-bracket notation
 *  \return          the free energy of the input structure given the input sequence in kcal/mol
 */
DEPRECATED(float energy_of_circ_struct( const char *string,
                                        const char *structure));

#endif
