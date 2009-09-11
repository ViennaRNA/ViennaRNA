#ifndef __VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H__
#define __VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H__

#include "params.h"
/**
*** Compute the Energy of an interior loop
*** @see scale_parameters() 
*** @see paramT
*** \warning Not totally thread safe (yet)!
*** \param  n1      The size of the 'left'-loop (number of unpaired nucleotides)
*** \param  n2      The size of the 'right'-loop (number of unpaired nucleotides)
*** \param  type    The pair type of the base pair closing the interior loop
*** \param  type_2  The pair type of the enclosed base pair
*** \param  si1     The 5'-mismatching nucleotide of the closing pair
*** \param  sj1     The 3'-mismatching nucleotide of the closing pair
*** \param  sp1     The 3'-mismatching nucleotide of the enclosed pair
*** \param  sq1     The 5'-mismatching nucleotide of the enclosed pair
*** \param  P       The datastructure containing scaled energy parameters
*** \return The Free energy of the Interior-loop in dcal/mol
**/
int   E_IntLoop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P);
/**
*** Compute the Energy of a hairpin loop
*** @see scale_parameters() 
*** @see paramT
*** \warning Not totally thread safe (yet)!
*** \param  size  The size of the loop (number of unpaired nucleotides)
*** \param  type  The pair type of the base pair closing the hairpin
*** \param  si1   The 5'-mismatching nucleotide
*** \param  sj1   The 3'-mismatching nucleotide
*** \param  string  The sequence of the loop
*** \param  P     The datastructure containing scaled energy parameters
*** \return The Free energy of the Haurpin-loop in dcal/mol
**/
int   E_Hairpin(int size, int type, int si1, int sj1, const char *string, paramT *P);

/**
*** compute Boltzmann weight of a hairpin loop, multiply by scale[u+2]
*** \return The Boltzmann weight of the Hairpin-loop
**/
double  exp_E_Hairpin(int u, int type, short si1, short sj1, const char *string, pf_paramT *P);
/**
*** compute Boltzmann weight of interior loop, multiply by scale[u1+u2+2] for scaling
*** \return The Boltzmann weight of the Interior-loop
**/
double  exp_E_IntLoop(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1, pf_paramT *P);

#endif
