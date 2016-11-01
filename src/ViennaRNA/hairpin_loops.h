#ifndef VIENNA_RNA_PACKAGE_HAIRPIN_LOOPS_H
#define VIENNA_RNA_PACKAGE_HAIRPIN_LOOPS_H

#include <math.h>
#include <string.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/**
 *
 *  @file     hairpin_loops.h
 *  @ingroup  loops
 *  @brief    Energy evaluation of hairpin loops for MFE and partition function calculations
 */

/**
 *
 *  @{
 *  @ingroup   loops
 */


/**
 *  @brief Compute the Energy of a hairpin-loop
 *
 *  To evaluate the free energy of a hairpin-loop, several parameters have to be known.
 *  A general hairpin-loop has this structure:<BR>
 *  <PRE>
 *        a3 a4
 *      a2     a5
 *      a1     a6
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  where X-Y marks the closing pair [e.g. a <B>(G,C)</B> pair]. The length of this loop is 6 as there are
 *  six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is
 *  a1 while the 3' mismatch is a6. The nucleotide sequence of this loop is &quot;a1.a2.a3.a4.a5.a6&quot; <BR>
 *  @note The parameter sequence should contain the sequence of the loop in capital letters of the nucleic acid
 *  alphabet if the loop size is below 7. This is useful for unusually stable tri-, tetra- and hexa-loops
 *  which are treated differently (based on experimental data) if they are tabulated.
 *  @see scale_parameters()
 *  @see vrna_param_t
 *  @warning Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *  Energy evaluation may change due to updates in global variable "tetra_loop"
 * 
 *  @param  size  The size of the loop (number of unpaired nucleotides)
 *  @param  type  The pair type of the base pair closing the hairpin
 *  @param  si1   The 5'-mismatching nucleotide
 *  @param  sj1   The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop
 *  @param  P     The datastructure containing scaled energy parameters
 *  @return The Free energy of the Hairpin-loop in dcal/mol
 */
PRIVATE INLINE int
E_Hairpin(int size,
              int type,
              int si1,
              int sj1,
              const char *string,
              vrna_param_t *P);

/**
 *  @brief Compute Boltzmann weight @f$e^{-\Delta G/kT} @f$ of a hairpin loop
 *
 *  multiply by scale[u+2]
 *  @see get_scaled_pf_parameters()
 *  @see vrna_exp_param_t
 *  @see E_Hairpin()
 *  @warning Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *  Energy evaluation may change due to updates in global variable "tetra_loop"
 * 
 *  @param  u       The size of the loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the hairpin
 *  @param  si1     The 5'-mismatching nucleotide
 *  @param  sj1     The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop
 *  @param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
 *  @return The Boltzmann weight of the Hairpin-loop
 */
PRIVATE INLINE FLT_OR_DBL
exp_E_Hairpin(  int u,
                int type,
                short si1,
                short sj1,
                const char *string,
                vrna_exp_param_t *P);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE INLINE int
E_Hairpin(int size,
          int type,
          int si1,
          int sj1,
          const char *string,
          vrna_param_t *P){

  int energy;

  if(size <= 30)
    energy = P->hairpin[size];
  else
    energy = P->hairpin[30] + (int)(P->lxc*log((size)/30.));

  if(size < 3) return energy; /* should only be the case when folding alignments */

  if(P->model_details.special_hp){
    if(size == 4){ /* check for tetraloop bonus */
      char tl[7]={0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl)))
        return (P->Tetraloop_E[(ts - P->Tetraloops)/7]);
    }
    else if(size == 6){
      char tl[9]={0}, *ts;
      strncpy(tl, string, 8);
      if ((ts=strstr(P->Hexaloops, tl)))
        return (energy = P->Hexaloop_E[(ts - P->Hexaloops)/9]);
    }
    else if(size == 3){
      char tl[6]={0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 5);
      if ((ts=strstr(P->Triloops, tl))) {
        return (P->Triloop_E[(ts - P->Triloops)/6]);
      }
      return (energy + (type>2 ? P->TerminalAU : 0));
    }
  }
  energy += P->mismatchH[type][si1][sj1];

  return energy;
}

/**
 *  @brief  Evaluate the free energy of a hairpin loop
 *          and consider hard constraints if they apply
 *
 *  This function evaluates the free energy of a hairpin loop
 *
 *  In case the base pair is not allowed due to a constraint
 *  conflict, this function returns #INF.
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @param vc   The #vrna_fold_compound_t that stores all relevant model settings
 *  @param i    The 5' nucleotide of the base pair (3' to evaluate the pair as exterior hairpin loop)
 *  @param j    The 3' nucleotide of the base pair (5' to evaluate the pair as exterior hairpin loop)
 *  @returns    The free energy of the hairpin loop in 10cal/mol
 */
int
vrna_E_hp_loop( vrna_fold_compound_t *vc,
                int i,
                int j);

/**
 *  @brief  Evaluate the free energy of an exterior hairpin loop
 *          and consider possible hard constraints
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 */
int
vrna_E_ext_hp_loop( vrna_fold_compound_t *vc,
                    int i,
                    int j);

/**
 *  @brief Evaluate free energy of an exterior hairpin loop
 *
 *  @ingroup loops
 *
 */
int
vrna_eval_ext_hp_loop(vrna_fold_compound_t *vc,
                      int i,
                      int j);

/**
 *  @brief Evaluate free energy of a hairpin loop
 *
 *  @ingroup loops
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @param  vc  The #vrna_fold_compound_t for the particular energy evaluation
 *  @param  i   5'-position of the base pair
 *  @param  j   3'-position of the base pair
 *  @returns    Free energy of the hairpin loop closed by @f$ (i,j) @f$ in deka-kal/mol
 */
int
vrna_eval_hp_loop(vrna_fold_compound_t *vc,
                  int i,
                  int j);

/*
*************************************
* Partition function variants below *
*************************************
*/

PRIVATE INLINE FLT_OR_DBL
exp_E_Hairpin(int u,
              int type,
              short si1,
              short sj1,
              const char *string,
              vrna_exp_param_t *P){

  double q, kT;
  kT = P->kT;   /* kT in cal/mol  */

  if(u <= 30)
    q = P->exphairpin[u];
  else
    q = P->exphairpin[30] * exp( -(P->lxc*log( u/30.))*10./kT);

  if(u < 3) return (FLT_OR_DBL)q; /* should only be the case when folding alignments */

  if(P->model_details.special_hp){
    if(u==4){
      char tl[7]={0,0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl))){
        if(type != 7)
          return (FLT_OR_DBL)(P->exptetra[(ts-P->Tetraloops)/7]);
        else
          q *= P->exptetra[(ts-P->Tetraloops)/7];
      }
    }
    else if(u==6){
      char tl[9]={0,0,0,0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 8);
      if ((ts=strstr(P->Hexaloops, tl)))
        return  (FLT_OR_DBL)(P->exphex[(ts-P->Hexaloops)/9]);
    }
    else if(u==3){
      char tl[6]={0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 5);
      if ((ts=strstr(P->Triloops, tl)))
        return (FLT_OR_DBL)(P->exptri[(ts-P->Triloops)/6]);
      if (type>2)
        return (FLT_OR_DBL)(q * P->expTermAU);
      else
        return (FLT_OR_DBL)q;
    }
  }
  q *= P->expmismatchH[type][si1][sj1];

  return (FLT_OR_DBL)q;
}


/**
 *  @brief High-Level function for hairpin loop energy evaluation (partition function variant)
 *
 *  @see vrna_E_hp_loop() for it's free energy counterpart
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
*/
FLT_OR_DBL
vrna_exp_E_hp_loop( vrna_fold_compound_t *vc,
                    int i,
                    int j);

/**
 *  @brief Backtrack a hairpin loop closed by @f$ (i,j) @f$
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 */
int
vrna_BT_hp_loop(vrna_fold_compound_t *vc,
                int i,
                int j,
                int en,
                vrna_bp_stack_t *bp_stack,
                int   *stack_count);

/**
 * @}
 */


#endif
