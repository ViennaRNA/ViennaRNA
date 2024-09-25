#ifndef VIENNA_RNA_PACKAGE_LOOPS_HAIRPIN_H
#define VIENNA_RNA_PACKAGE_LOOPS_HAIRPIN_H

#include <math.h>
#include <string.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/params/salt.h>
#include <ViennaRNA/eval/basic.h>


#ifdef VRNA_WARN_DEPRECATED
# if defined(DEPRECATED)
#   undef DEPRECATED
# endif
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

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/**
 *
 *  @file     ViennaRNA/loops/hairpin.h
 *  @ingroup  eval, eval_loops, eval_loops_hp
 *  @brief    Energy evaluation of hairpin loops for MFE and partition function calculations
 */

/**
 *
 *  @addtogroup   eval_loops_hp
 *  @{
 */


/**
 *  @name Basic free energy interface
 *  @{
 */


/**
 *  @brief Retrieve the energy of a hairpin-loop
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
 *
 *  @warning  This function @b only evaluates the free energy of a hairpin loop according to the current
 *            Turner energy parameter set! No additional hard- or soft constraints are applied. See vrna_eval_hp_loop()
 *            for a function that also takes into account any user-supplied constraints!
 *
 *  @note     Whenever one of the mismatch base encodings @p si1 or @p sj1 is negative, terminal mismatch energies are not applied!
 *
 *  @note     The parameter @p sequence is a 0-terminated string of size @p size + 2 that contain the nucleic acid
 *            sequence of the loop in upper-case letters. This parameter is only required for loops of @p size below
 *            7, since it is used for look-up of unusually stable tri-, tetra- and hexa-loops, such as GNRA tetra loops.
 *            Those may have additional sequence-dependent tabulated free energies available.
 *
 *  @see      #vrna_param_t, vrna_eval_hp_loop()
 *
 *  @param  size      The size of the loop (number of unpaired nucleotides)
 *  @param  type      The pair type of the base pair closing the hairpin
 *  @param  si1       The 5'-mismatching nucleotide
 *  @param  sj1       The 3'-mismatching nucleotide
 *  @param  sequence  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P         The datastructure containing scaled energy parameters
 *  @return           The Free energy of the Hairpin-loop in dcal/mol
 */
int
vrna_E_hairpin(unsigned int size,
               unsigned int type,
               int          si1,
               int          sj1,
               const char   *sequence,
               vrna_param_t *P);


/**
 *  @brief Evaluate free energy of a hairpin loop
 *
 *  This function evaluates the free energy of a hairpin loop closed by a base pair
 *  (i,j). By default (@p options = #VRNA_EVAL_DEFAULT), @emph all user-supplied
 *  constraints will be taken into consideration. This means that any hard constraints
 *  that prohibit the formation of this loop will result in an energy contribution
 *  of @b #INF. On the other hand, if, given the set of constraints, the loop is
 *  allowed then its free energy is evaluated according to the Nearest Neighbor
 *  energy parameter set. On top of that, any user-supplied soft-constraints will
 *  be added, if applicable.
 *
 *  The @p options argument allows for (de-)activating certain aspects of the evaluation,
 *  e.g. hard constraints, soft constraints, etc.
 *
 *  @note   If sequence position @p i is larger than @p j, the function assumes a
 *          hairpin loop formed by a circular RNA, where the unpaired loop sequence
 *          spans the n,1-junction.
 *
 *  @note   By default, all user-supplied hard- and soft constraints will be taken
 *          into account! Use the #VRNA_EVAL_LOOP_NO_HC and #VRNA_EVAL_LOOP_NO_SC
 *          bit flags as input for @p options to change the default behavior if necessary.
 *
 *  @note   This function is polymorphic! The provided #vrna_fold_compound_t may be
 *          of type #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @see    vrna_E_hairpin(), vrna_exp_eval_hairpin(),
 *          #VRNA_EVAL_LOOP_NO_HC, #VRNA_EVAL_LOOP_NO_SC, #VRNA_EVAL_LOOP_NO_CONSTRAINTS
 *
 *  @param  fc        The #vrna_fold_compound_t for the particular energy evaluation
 *  @param  i         5'-position of the base pair
 *  @param  j         3'-position of the base pair
 *  @param  options   A bit-field that specifies which aspects (not) to consider during evaluation
 *  @returns          Free energy of the hairpin loop closed by @f$ (i,j) @f$ in deka-kal/mol or #INF if the loop is forbidden
 */
int
vrna_eval_hairpin(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j,
                  unsigned int          options);


/* End basic interface */
/** @} */


/**
 *  @name Boltzmann weight (partition function) interface
 *  @{
 */


/**
 *  @brief Compute Boltzmann weight @f$e^{-\Delta G/kT} @f$ of a hairpin loop
 *
 *  This is the partition function variant of vrna_E_hp() that returns the Boltzmann
 *  weight @f$e^{-\Delta E/kT} @f$ instead of the energy @f$ E @f$.
 *
 *  @note   Whenever one of the mismatch base encodings @p si1 or @p sj1 is negative, terminal mismatch energies are not applied!
 *
 *  @note   Do not forget to scale this Bolztmann factor properly, e.g. by multiplying with scale[u+2]
 *
 *  @see vrna_exp_eval_hp_loop(), #vrna_exp_param_t, vrna_E_hp()
 *
 *  @param  size      The size of the loop (number of unpaired nucleotides)
 *  @param  type      The pair type of the base pair closing the hairpin
 *  @param  si1       The 5'-mismatching nucleotide
 *  @param  sj1       The 3'-mismatching nucleotide
 *  @param  sequence  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P         The datastructure containing scaled Boltzmann weights of the energy parameters
 *  @return The Boltzmann weight of the Hairpin-loop
 */
FLT_OR_DBL
vrna_exp_E_hairpin(unsigned int     size,
                   unsigned int     type,
                   int              si1,
                   int              sj1,
                   const char       *sequence,
                   vrna_exp_param_t *P);


/**
 *  @brief High-Level function for hairpin loop energy evaluation (partition function variant)
 *
 *  This is the partition function variant of vrna_eval_hp_loop() that returns the Boltzmann
 *  weight @f$e^{-\Delta E/kT} @f$ instead of the energy @f$ E @f$. On top of all constraints
 *  application, this function already scales the Boltzmann factor, i.e. it multiplies the
 *  result with scale[u + 2]
 *
 *  The @p options argument allows for (de-)activating certain aspects of the evaluation,
 *  e.g. hard constraints, soft constraints, etc.
 *
 *  @note   If sequence position @p i is larger than @p j, the function assumes a
 *          hairpin loop formed by a circular RNA, where the unpaired loop sequence
 *          spans the n,1-junction.
 *
 *  @note   By default, all user-supplied hard- and soft constraints will be taken
 *          into account! Use the #VRNA_EVAL_LOOP_NO_HC and #VRNA_EVAL_LOOP_NO_SC
 *          bit flags to change the default behavior if necessary.
 *
 *  @note   This function is polymorphic! The provided #vrna_fold_compound_t may be
 *          of type #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @see    vrna_eval_hairpin(), vrna_exp_E_hairpin(),
 *          #VRNA_EVAL_LOOP_NO_HC, #VRNA_EVAL_LOOP_NO_SC, #VRNA_EVAL_LOOP_NO_CONSTRAINTS
 *
 *  @param  fc        The #vrna_fold_compound_t for the particular energy evaluation
 *  @param  i         5'-position of the base pair
 *  @param  j         3'-position of the base pair
 *  @param  options   A bit-field that specifies which aspects (not) to consider during evaluation
 *  @returns          Boltzmann factor of the free energy of the hairpin loop closed by @f$ (i,j) @f$ or 0. if the loop is forbidden
 */
FLT_OR_DBL
vrna_exp_eval_hairpin(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      unsigned int          options);


/* End partition function interface */
/** @} */


/**
 * @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup   eval_deprecated
 *  @{
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
 *
 *  @note The parameter sequence should contain the sequence of the loop in capital letters of the nucleic acid
 *        alphabet if the loop size is below 7. This is useful for unusually stable tri-, tetra- and hexa-loops
 *        which are treated differently (based on experimental data) if they are tabulated.
 *
 *  @see scale_parameters(), vrna_param_t
 *
 *  @warning  Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *            Energy evaluation may change due to updates in global variable "tetra_loop"
 *
 *  @param  size  The size of the loop (number of unpaired nucleotides)
 *  @param  type  The pair type of the base pair closing the hairpin
 *  @param  si1   The 5'-mismatching nucleotide
 *  @param  sj1   The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P     The datastructure containing scaled energy parameters
 *  @return The Free energy of the Hairpin-loop in dcal/mol
 */
DEPRECATED(PRIVATE INLINE int
           E_Hairpin(int          size,
                     int          type,
                     int          si1,
                     int          sj1,
                     const char   *string,
                     vrna_param_t *P),
           "Use vrna_E_hairpin() instead!");


DEPRECATED(int
           vrna_E_hp_loop(vrna_fold_compound_t  *fc,
                          int                   i,
                          int                   j),
           "Use vrna_eval_hairpin() instead!");


DEPRECATED(int
           vrna_E_ext_hp_loop(vrna_fold_compound_t  *fc,
                              int                   i,
                              int                   j),
           "Use vrna_eval_hairpin() instead!");


DEPRECATED(int
           vrna_eval_hp_loop(vrna_fold_compound_t *fc,
                             int                  i,
                             int                  j),
           "Use vrna_eval_hairpin() instead!");


DEPRECATED(int
           vrna_eval_ext_hp_loop(vrna_fold_compound_t *fc,
                                 int                  i,
                                 int                  j),
           "Use vrna_eval_hairpin() instead!");


/**
 *  @brief Compute Boltzmann weight @f$e^{-\Delta G/kT} @f$ of a hairpin loop
 *
 *  @note multiply by scale[u+2]
 *
 *  @see get_scaled_pf_parameters(), vrna_exp_param_t, E_Hairpin()
 *
 *  @warning  Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *            Energy evaluation may change due to updates in global variable "tetra_loop"
 *
 *  @param  u       The size of the loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the hairpin
 *  @param  si1     The 5'-mismatching nucleotide
 *  @param  sj1     The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
 *  @return The Boltzmann weight of the Hairpin-loop
 */
DEPRECATED(PRIVATE INLINE FLT_OR_DBL
           exp_E_Hairpin(int              u,
                         int              type,
                         short            si1,
                         short            sj1,
                         const char       *string,
                         vrna_exp_param_t *P),
           "Use vrna_exp_E_hairpin() instead!");


DEPRECATED(FLT_OR_DBL
           vrna_exp_E_hp_loop(vrna_fold_compound_t  *fc,
                              int                   i,
                              int                   j),
           "Use vrna_exp_eval_hairpin() instead!");


PRIVATE INLINE int
E_Hairpin(int           size,
          int           type,
          int           si1,
          int           sj1,
          const char    *string,
          vrna_param_t  *P)
{
  int energy, salt_correction;

  salt_correction = 0;

  if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
    if (size <= MAXLOOP)
      salt_correction = P->SaltLoop[size + 1];
    else
      salt_correction = vrna_salt_loop_int(size + 1,
                                           P->model_details.salt,
                                           P->temperature + K0,
                                           P->model_details.backbone_length);
  }

  if (size <= 30)
    energy = P->hairpin[size];
  else
    energy = P->hairpin[30] + (int)(P->lxc * log((size) / 30.));

  energy += salt_correction;

  if (size < 3)
    return energy;            /* should only be the case when folding alignments */

  if ((string) && (P->model_details.special_hp)) {
    if (size == 4) {
      /* check for tetraloop bonus */
      char tl[7] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 6);
      tl[6] = '\0';
      if ((ts = strstr(P->Tetraloops, tl)))
        return P->Tetraloop_E[(ts - P->Tetraloops) / 7] + salt_correction;
    } else if (size == 6) {
      char tl[9] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 8);
      tl[8] = '\0';
      if ((ts = strstr(P->Hexaloops, tl)))
        return P->Hexaloop_E[(ts - P->Hexaloops) / 9] + salt_correction;
    } else if (size == 3) {
      char tl[6] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 5);
      tl[5] = '\0';
      if ((ts = strstr(P->Triloops, tl)))
        return P->Triloop_E[(ts - P->Triloops) / 6] + salt_correction;

      return energy + (type > 2 ? P->TerminalAU : 0);
    }
  }

  energy += P->mismatchH[type][si1][sj1];

  return energy;
}


PRIVATE INLINE FLT_OR_DBL
exp_E_Hairpin(int               u,
              int               type,
              short             si1,
              short             sj1,
              const char        *string,
              vrna_exp_param_t  *P)
{
  double q, kT, salt_correction;

  kT              = P->kT; /* kT in cal/mol  */
  salt_correction = 1.;

  if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
    if (u <= MAXLOOP)
      salt_correction = P->expSaltLoop[u + 1];
    else
      salt_correction =
        exp(-vrna_salt_loop_int(u + 1, P->model_details.salt, P->temperature + K0,
                                P->model_details.backbone_length) * 10. / kT);
  }

  if (u <= 30)
    q = P->exphairpin[u];
  else
    q = P->exphairpin[30] * exp(-(P->lxc * log(u / 30.)) * 10. / kT);

  q *= salt_correction;

  if (u < 3)
    return (FLT_OR_DBL)q;         /* should only be the case when folding alignments */

  if ((string) && (P->model_details.special_hp)) {
    if (u == 4) {
      char tl[7] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 6);
      tl[6] = '\0';
      if ((ts = strstr(P->Tetraloops, tl))) {
        if (type != 7)
          return (FLT_OR_DBL)(P->exptetra[(ts - P->Tetraloops) / 7] * salt_correction);
        else
          q *= P->exptetra[(ts - P->Tetraloops) / 7];
      }
    } else if (u == 6) {
      char tl[9] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 8);
      tl[8] = '\0';
      if ((ts = strstr(P->Hexaloops, tl)))
        return (FLT_OR_DBL)(P->exphex[(ts - P->Hexaloops) / 9] * salt_correction);
    } else if (u == 3) {
      char tl[6] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 5);
      tl[5] = '\0';
      if ((ts = strstr(P->Triloops, tl)))
        return (FLT_OR_DBL)(P->exptri[(ts - P->Triloops) / 6] * salt_correction);

      if (type > 2)
        return (FLT_OR_DBL)(q * P->expTermAU);
      else
        return (FLT_OR_DBL)q;
    }
  }

  q *= P->expmismatchH[type][si1][sj1];

  return (FLT_OR_DBL)q;
}


/** @} */

#endif

#endif
