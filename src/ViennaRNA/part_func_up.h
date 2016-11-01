#ifndef VIENNA_RNA_PACKAGE_PART_FUNC_UP_H
#define VIENNA_RNA_PACKAGE_PART_FUNC_UP_H

#include <ViennaRNA/data_structures.h>

#define   RNA_UP_MODE_1   1U
#define   RNA_UP_MODE_2   2U
#define   RNA_UP_MODE_3   4U

/**
 *  @file     part_func_up.h
 *  @ingroup  pf_fold cofold pf_cofold
 *  @brief    Implementations for accessibility and RNA-RNA interaction as a stepwise process
 */

/**
 *  @addtogroup up_cofold
 *  @brief      RNA-RNA interaction as a stepwise process
 *
 * 
 *  In this approach to cofolding the interaction between two RNA molecules is
 *  seen as a stepwise process. In a first step, the target molecule has to
 *  adopt a structure in which a binding site is accessible. In a second step,
 *  the ligand molecule will hybridize with a region accessible to an
 *  interaction. Consequently the algorithm is designed as a two step process:
 *  The first step is the calculation of the probability
 *  that a region within the target is unpaired, or equivalently, the
 *  calculation of the free energy needed to expose a region. In the second step
 *  we compute the free energy of an interaction for every possible binding site.
 *  @{
 *  @ingroup  up_cofold
 */

/**
 *  @brief Calculate the partition function over all unpaired regions
 *  of a maximal length.
 * 
 *  You have to call function pf_fold() providing the same sequence before calling
 *  pf_unstru(). If you want to calculate unpaired regions for a constrained structure, set
 *  variable 'structure' in function 'pf_fold()' to the constrain string.
 *  It returns a #pu_contrib struct containing four arrays of dimension
 *  [i = 1 to length(sequence)][j = 0 to u-1] containing all possible contributions
 *  to the probabilities of unpaired regions of maximum length u.
 *  Each array in #pu_contrib contains one of the contributions to the
 *  total probability of being unpaired: The probability of being unpaired
 *  within an exterior loop is in array #pu_contrib->E, the probability
 *  of being unpaired within a hairpin loop is in array #pu_contrib->H,
 *  the probability of being unpaired within an interior loop is in array
 *  #pu_contrib->I and probability of being unpaired within a multi-loop
 *  is in array #pu_contrib->M. The total probability of being unpaired
 *  is the sum of the four arrays of #pu_contrib.
 * 
 *  This function frees everything allocated automatically. To
 *  free the output structure call free_pu_contrib().
 * 
 *  @param sequence
 *  @param max_w
 *  @return
 */
pu_contrib *pf_unstru(char *sequence,
                      int max_w);

/**
 *  @brief Calculates the probability of a local interaction between two sequences.
 * 
 *  The function considers the probability that the
 *  region of interaction is unpaired within 's1' and 's2'. The
 *  longer sequence has to be given as 's1'. The shorter sequence has to
 *  be given as 's2'. Function pf_unstru() has to be called
 *  for 's1' and 's2', where the probabilities of  being unpaired
 *  have to be given in 'p_c' and 'p_c2', respectively. If you do
 *  not want to include the probabilities of  being unpaired for 's2' set
 *  'p_c2' to NULL. If variable 'cstruc' is not NULL,
 *  constrained folding is done: The available constrains for intermolecular
 *  interaction are: '.' (no constrain), 'x' (the base has no intermolecular
 *  interaction) and '|' (the corresponding base has to be paired
 *  intermolecularily).\n
 *  The parameter 'w' determines the maximal length of the interaction. The
 *  parameters 'incr5' and 'incr3' allows inclusion of
 *  unpaired residues left ('incr5') and right ('incr3') of the region
 *  of interaction in 's1'. If the 'incr' options are used, function
 *  pf_unstru() has to be called with
 *  w=w+incr5+incr3 for the longer sequence 's1'.
 * 
 *  It returns a structure of type #interact which
 *  contains the probability of the best local interaction including residue i
 *  in Pi and the minimum free energy in Gi, where i is the position in sequence
 *  's1'. The member Gikjl of structure #interact is
 *  the best interaction between region [k,i] k<i in longer sequence
 *  's1' and region [j,l] j<l in 's2'. Gikjl_wo is Gikjl without the
 *  probability of beeing unpaired.\n
 *  Use free_interact() to free the returned structure, all
 *  other stuff is freed inside pf_interact().
 * 
 *  @param  s1
 *  @param  s2
 *  @param  p_c
 *  @param  p_c2
 *  @param  max_w
 *  @param  cstruc
 *  @param  incr3
 *  @param  incr5
 *  @return
 */
interact *pf_interact(const char *s1,
                      const char *s2,
                      pu_contrib *p_c,
                      pu_contrib *p_c2,
                      int max_w,
                      char *cstruc,
                      int incr3,
                      int incr5);

/**
 *  @brief Frees the output of function pf_interact().
 */
void free_interact(interact *pin);

/**
 *  @brief
 */
int Up_plot(pu_contrib *p_c,
            pu_contrib *p_c_sh,
            interact *pint,
            char *ofile,
            int **unpaired_values,
            char *select_contrib,
            char *head,
            unsigned int mode);

/**
 *  @brief
 */
pu_contrib  *get_pu_contrib_struct( unsigned int n,
                                    unsigned int w);

/**
 *  @brief Frees the output of function pf_unstru().
 */
void        free_pu_contrib_struct(pu_contrib *pu);

void
free_pu_contrib(pu_contrib *pu);

/**
 * @}
 */

#endif
