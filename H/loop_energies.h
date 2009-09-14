#ifndef __VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H__
#define __VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H__

#include "params.h"
/**
*** <H2>Compute the Energy of an interior-loop</H2>
*** This function computes the free energy \f$\Delta G\f$ of an interior-loop with the
*** following structure: <BR>
*** <PRE>
***       3'  5'
***       |   |
***       U - V
***   a_n       b_1
***    .        .
***    .        .
***    .        .
***   a_1       b_m
***       X - Y
***       |   |
***       5'  3'
*** </PRE>
*** This general structure depicts an interior-loop that is closed by the base pair (X,Y).
*** The enclosed base pair is (V,U) which leaves the unpaired bases a_1-a_n and b_1-b_n
*** that constitute the loop. In this example, the length of the interior-loop is \f$(n+m)\f$
*** where n or m may be 0 resulting in a bulge-loop or base pair stack.
*** The mismatching nucleotides for the closing pair (X,Y) are:<BR>
*** 5'-mismatch: a_1<BR>
*** 3'-mismatch: b_m<BR>
*** and for the enclosed base pair (V,U):<BR>
*** 5'-mismatch: b_1<BR>
*** 3'-mismatch: a_n<BR> 
*** \note Base pairs are always denoted in 5'->3' direction. Thus the enclosed base pair
*** must be 'turned arround' when evaluating the free energy of the interior-loop 
*** @see scale_parameters() 
*** @see paramT
*** \warning Not thread safe! A threadsafe implementation will replace this function in a future release!
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
*** <H2>Compute the Energy of a hairpin-loop</H2>
*** To evaluate the free energy of a hairpin-loop, several parameters have to be known.
*** A general hairpin-loop has this structure:<BR>
*** <PRE>
***       a3 a4
***     a2     a5
***     a1     a6
***       X - Y
***       |   |
***       5'  3'
*** </PRE>
*** where X-Y marks the closing pair [e.g. a <B>(G,C)</B> pair]. The length of this loop is 6 as there are
*** six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is
*** a1 while the 3' mismatch is a6. The nucleotide sequence of this loop is &quot;a1.a2.a3.a4.a5.a6&quot; <BR>
*** \note The parameter sequence should contain the sequence of the loop in capital letters of the nucleic acid
*** alphabet if the loop size is below 7. This is useful for unusually stable tri-, tetra- and hexa-loops
*** which are treated differently (based on experimental data) if they are tabulated.
*** @see scale_parameters() 
*** @see paramT
*** \warning Not thread safe! A threadsafe implementation will replace this function in a future release!
*** \param  size  The size of the loop (number of unpaired nucleotides)
*** \param  type  The pair type of the base pair closing the hairpin
*** \param  si1   The 5'-mismatching nucleotide
*** \param  sj1   The 3'-mismatching nucleotide
*** \param  string  The sequence of the loop
*** \param  P     The datastructure containing scaled energy parameters
*** \return The Free energy of the Hairpin-loop in dcal/mol
**/
int   E_Hairpin(int size, int type, int si1, int sj1, const char *string, paramT *P);

/**
*** <H2>compute Boltzmann weight \f$e^{-\Delta G/kT} \f$ of a hairpin loop</H2>
*** multiply by scale[u+2]
*** @see get_scaled_pf_parameters() 
*** @see pf_paramT
*** @see E_Hairpin()
*** \warning Not thread safe! A threadsafe implementation will replace this function in a future release!
*** \param  u       The size of the loop (number of unpaired nucleotides)
*** \param  type    The pair type of the base pair closing the hairpin
*** \param  si1     The 5'-mismatching nucleotide
*** \param  sj1     The 3'-mismatching nucleotide
*** \param  string  The sequence of the loop
*** \param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
*** \return The Boltzmann weight of the Hairpin-loop
**/
double  exp_E_Hairpin(int u, int type, short si1, short sj1, const char *string, pf_paramT *P);
/**
*** <H2>compute Boltzmann weight \f$e^{-\Delta G/kT} \f$ of interior loop</H2>
*** multiply by scale[u1+u2+2] for scaling
*** @see get_scaled_pf_parameters() 
*** @see pf_paramT
*** @see E_IntLoop()
*** \warning Not thread safe! A threadsafe implementation will replace this function in a future release!
*** \param  u1      The size of the 'left'-loop (number of unpaired nucleotides)
*** \param  u2      The size of the 'right'-loop (number of unpaired nucleotides)
*** \param  type    The pair type of the base pair closing the interior loop
*** \param  type2   The pair type of the enclosed base pair
*** \param  si1     The 5'-mismatching nucleotide of the closing pair
*** \param  sj1     The 3'-mismatching nucleotide of the closing pair
*** \param  sp1     The 3'-mismatching nucleotide of the enclosed pair
*** \param  sq1     The 5'-mismatching nucleotide of the enclosed pair
*** \param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
*** \return The Boltzmann weight of the Interior-loop
**/
double  exp_E_IntLoop(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1, pf_paramT *P);

#endif
