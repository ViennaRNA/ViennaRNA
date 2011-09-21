/* functions from part_func.c */
/* calculate partition function and base pair probabilities */
#ifndef __VIENNA_RNA_PACKAGE_PART_FUNC_CO_H__
#define __VIENNA_RNA_PACKAGE_PART_FUNC_CO_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  \file part_func_co.h
 * 
 *  \brief Partition function for two RNA sequences
 * 
 *  As for folding one RNA molecule, this computes the partition function
 *  of all possible structures and the base pair probabilities. Uses the
 *  same global #pf_scale variable to avoid overflows.
 * 
 *  To simplify the implementation the partition function computation is done
 *  internally in a null model that does not include the duplex initiation
 *  energy, i.e. the entropic penalty for producing a dimer from two
 *  monomers). The resulting free energies and pair probabilities are initially
 *  relative to that null model. In a second step the free energies can be
 *  corrected to include the dimerization penalty, and the pair probabilities
 *  can be divided into the conditional pair probabilities given that a re
 *  dimer is formed or not formed.
 * 
 *  After computing the partition functions of all possible dimeres one
 *  can compute the probabilities of base pairs, the concentrations out of
 *  start concentrations and sofar and soaway.
 * 
 *  Dimer formation is inherently concentration dependent. Given the free
 *  energies of the monomers A and B and dimers AB, AA, and BB one can compute
 *  the equilibrium concentrations, given input concentrations of A and B, see
 *  e.g. Dimitrov & Zuker (2004)
 */

/**
 *  \brief Toggles no intrabp in 2nd mol
 */
extern int    mirnatog;

/**
 *  \brief Free energies of the two monomers
 */
extern double F_monomer[2];

/**
 *  \brief Calculate partition function and base pair probabilities
 * 
 *  This is the cofold partition function folding. The second molecule starts
 *  at the #cut_point nucleotide.
 *  
 *  \param  sequence  Concatenated RNA sequences
 *  \param  structure Will hold the structure or constraints
 *  \return           cofoldF structure containing a set of energies needed for
 *                    concentration computations.
 */
cofoldF co_pf_fold( char *sequence,
                    char *structure);

/**
 *  \brief Free the memory occupied by co_pf_fold()
 */
void    free_co_pf_arrays(void);

/**
 *  \brief Recalculate energy parameters
 */
void    update_co_pf_params(int length);

/**
 *  \brief Compute Boltzmann probabilities of dimerization without homodimers
 * 
 *  Given the pair probabilities and free energies (in the null model) for a
 *  dimer AB and the two constituent monomers A and B, compute the conditional pair
 *  probabilities given that a dimer AB actually forms.
 *  Null model pair probabilities are given as a list as produced by
 *  assign_plist_from_pr(), the dimer probabilities 'prAB' are modified in place.
 * 
 *  \param FAB      free energy of dimer AB
 *  \param FEA      free energy of monomer A
 *  \param FEB      free energy of monomer B
 *  \param prAB     pair probabilities for dimer
 *  \param prA      pair probabilities monomer
 *  \param prB      pair probabilities monomer
 *  \param Alength  Length of molecule A
 */
void    compute_probabilities(double FAB,
                              double FEA,
                              double FEB,
                              struct plist  *prAB,
                              struct plist  *prA,
                              struct plist  *prB,
                              int Alength);

/**
 *  \brief Given two start monomer concentrations a and b, compute the
 *  concentrations in thermodynamic equilibrium of all dimers and the monomers.
 * 
 *  This function takes an array  'startconc' of input concentrations with alternating
 *  entries for the initial concentrations of molecules A and B (terminated by
 *  two zeroes), then computes the resulting equilibrium concentrations
 *  from the free energies for the dimers. Dimer free energies should be the
 *  dimer-only free energies, i.e. the FcAB entries from the #cofoldF struct.
 * 
 *  \param FEAB       Free energy of AB dimer (FcAB entry)
 *  \param FEAA       Free energy of AA dimer (FcAB entry)
 *  \param FEBB       Free energy of BB dimer (FcAB entry)
 *  \param FEA        Free energy of monomer A
 *  \param FEB        Free energy of monomer B
 *  \param startconc  List of start concentrations [a0],[b0],[a1],[b1],...,[an][bn],[0],[0]
 *  \return ConcEnt array containing the equilibrium energies and start concentrations
 */
ConcEnt *get_concentrations(double FEAB,
                            double FEAA,
                            double FEBB,
                            double FEA,
                            double FEB,
                            double *startconc);

/*
#################################################
# DEPRECATED FUNCTIONS                          #
#################################################
*/

/**
 *  DO NOT USE THIS FUNCTION ANYMORE
 *  \deprecated{ This function is deprecated and will be removed soon!}
 *  use \ref assign_plist_from_pr() instead!
 */
DEPRECATED(plist  *get_plist( struct plist *pl,
                              int length,
                              double cut_off));
/**
 *  DO NOT USE THIS FUNCTION ANYMORE
 *  \deprecated{ This function is deprecated and will be removed soon!}
 */
DEPRECATED(void   init_co_pf_fold(int length));

#endif
