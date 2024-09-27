#ifndef VIENNA_RNA_PACKAGE_EVAL_INTERNAL_H
#define VIENNA_RNA_PACKAGE_EVAL_INTERNAL_H

#include <math.h>

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/params/default.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/constraints/soft.h>
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
 *  @file     ViennaRNA/eval/internal.h
 *  @ingroup  eval, eval_loops, eval_loops_int
 *  @brief    Energy evaluation of internal loops
 */

/**
 *  @addtogroup   eval_loops_int
 *  @{
 */


/**
 *  @name Basic free energy interface
 *  @{
 */


/**
 *  @brief  Compute the Energy of an internal loop
 *
 *  This function computes the free energy @f$ E @f$ of an internal-loop with the
 *  following structure: <BR>
 *  <PRE>
 *        3'  5'
 *        |   |
 *        U - V
 *    a_n       b_1
 *     .        .
 *     .        .
 *     .        .
 *    a_1       b_m
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  This general structure depicts an internal-loop that is closed by the base pair (X,Y).
 *  The enclosed base pair is (V,U) which leaves the unpaired bases a_1-a_n and b_1-b_n
 *  that constitute the loop. In this example, the length of the internal-loop is @f$(n+m)@f$
 *  where n or m may be 0 resulting in a bulge-loop or base pair stack.
 *  The mismatching nucleotides for the closing pair (X,Y) are:<BR>
 *  5'-mismatch: a_1<BR>
 *  3'-mismatch: b_m<BR>
 *  and for the enclosed base pair (V,U):<BR>
 *  5'-mismatch: b_1<BR>
 *  3'-mismatch: a_n<BR>
 *
 *  @note Base pairs are always denoted in 5'->3' direction. Thus the enclosed base pair
 *        must be 'turned arround' when evaluating the free energy of the internal-loop<br>
 *        This function is threadsafe
 *
 *  @see vrna_exp_E_internal()
 *
 *  @param  n1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  n2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the internal loop
 *  @param  type_2  The pair type of the enclosed base pair
 *  @param  si1     The 5'-mismatching nucleotide of the closing pair
 *  @param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  @param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  @param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled energy parameters
 *  @return The Free energy of the internal loop in dcal/mol
 */
int
vrna_E_internal(unsigned int  n1,
                unsigned int  n2,
                unsigned int  type,
                unsigned int  type_2,
                int           si1,
                int           sj1,
                int           sp1,
                int           sq1,
                vrna_param_t  *P);


/**
 *  @brief Evaluate the free energy contribution of an internal loop with delimiting
 *  base pairs @f$(i,j)@f$ and @f$(k,l)@f$
 *
 *  @note This function is polymorphic, i.e. it accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE as well as #VRNA_FC_TYPE_COMPARATIVE
 */
int
vrna_eval_internal(vrna_fold_compound_t *fc,
                   unsigned int         i,
                   unsigned int         j,
                   unsigned int         k,
                   unsigned int         l,
                   unsigned int         options);


int
vrna_eval_stack(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          j,
                unsigned int          options);


/* End basic interface */
/**@}*/


/**
 *  @name Boltzmann weight (partition function) interface
 *  @{
 */


FLT_OR_DBL
vrna_exp_E_internal(unsigned int      n1,
                    unsigned int      n2,
                    unsigned int      type,
                    unsigned int      type_2,
                    int               si1,
                    int               sj1,
                    int               sp1,
                    int               sq1,
                    vrna_exp_param_t  *P);


FLT_OR_DBL
vrna_exp_eval_internal(vrna_fold_compound_t *fc,
                       unsigned int         i,
                       unsigned int         j,
                       unsigned int         k,
                       unsigned int         l,
                       unsigned int         options);


/* End partition function interface */
/**@}*/

/**
 * @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup   eval_deprecated
 *  @{
 */


#ifdef ON_SAME_STRAND
#undef ON_SAME_STRAND
#endif

#define ON_SAME_STRAND(I, J, C)  (((I) >= (C)) || ((J) < (C)))

/**
 *  <H2>Compute the Energy of an internal-loop</H2>
 *  This function computes the free energy @f$\Delta G@f$ of an internal-loop with the
 *  following structure: <BR>
 *  <PRE>
 *        3'  5'
 *        |   |
 *        U - V
 *    a_n       b_1
 *     .        .
 *     .        .
 *     .        .
 *    a_1       b_m
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  This general structure depicts an internal-loop that is closed by the base pair (X,Y).
 *  The enclosed base pair is (V,U) which leaves the unpaired bases a_1-a_n and b_1-b_n
 *  that constitute the loop. In this example, the length of the internal-loop is @f$(n+m)@f$
 *  where n or m may be 0 resulting in a bulge-loop or base pair stack.
 *  The mismatching nucleotides for the closing pair (X,Y) are:<BR>
 *  5'-mismatch: a_1<BR>
 *  3'-mismatch: b_m<BR>
 *  and for the enclosed base pair (V,U):<BR>
 *  5'-mismatch: b_1<BR>
 *  3'-mismatch: a_n<BR>
 *
 *  @note Base pairs are always denoted in 5'->3' direction. Thus the enclosed base pair
 *        must be 'turned arround' when evaluating the free energy of the internal-loop<br>
 *        This function is threadsafe
 *
 *  @see scale_parameters(), vrna_param_t
 *
 *  @param  n1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  n2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the internal loop
 *  @param  type_2  The pair type of the enclosed base pair
 *  @param  si1     The 5'-mismatching nucleotide of the closing pair
 *  @param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  @param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  @param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled energy parameters
 *  @return The Free energy of the Interior-loop in dcal/mol
 */
DEPRECATED(PRIVATE INLINE int
           E_IntLoop(int          n1,
                     int          n2,
                     int          type,
                     int          type_2,
                     int          si1,
                     int          sj1,
                     int          sp1,
                     int          sq1,
                     vrna_param_t *P),
           "Use vrna_E_internal() instead!");


/**
 *  <H2>Compute Boltzmann weight @f$e^{-\Delta G/kT} @f$ of internal loop</H2>
 *  multiply by scale[u1+u2+2] for scaling
 *
 *  @note This function is threadsafe
 *
 *  @see get_scaled_pf_parameters(), vrna_exp_param_t, E_IntLoop()
 *
 *  @param  u1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  u2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the internal loop
 *  @param  type2   The pair type of the enclosed base pair
 *  @param  si1     The 5'-mismatching nucleotide of the closing pair
 *  @param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  @param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  @param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
 *  @return The Boltzmann weight of the Interior-loop
 */
DEPRECATED(PRIVATE INLINE FLT_OR_DBL
           exp_E_IntLoop(int              u1,
                         int              u2,
                         int              type,
                         int              type2,
                         short            si1,
                         short            sj1,
                         short            sp1,
                         short            sq1,
                         vrna_exp_param_t *P),
           "Use vrna_exp_E_internal() instead!");


DEPRECATED(PRIVATE INLINE int
           E_IntLoop_Co(int           type,
                        int           type_2,
                        int           i,
                        int           j,
                        int           p,
                        int           q,
                        int           cutpoint,
                        short         si1,
                        short         sj1,
                        short         sp1,
                        short         sq1,
                        int           dangles,
                        vrna_param_t  *P),
           "This function is obsolete");


PRIVATE INLINE int
__E_IntLoop_Co(int          type,
               int          type_2,
               int          i,
               int          j,
               int          p,
               int          q,
               int          cutpoint,
               short        si1,
               short        sj1,
               short        sp1,
               short        sq1,
               int          dangles,
               vrna_param_t *P);
/*
 *  ugly but fast internal loop evaluation
 *
 *  Avoid including this function in your own code. It only serves
 *  as a fast inline block internally re-used throughout the RNAlib. It
 *  evalutes the free energy of internal loops in single sequences or sequence
 *  hybrids. Soft constraints are also applied if available.
 *
 *  NOTE: do not include into doxygen reference manual!
 */
DEPRECATED(PRIVATE INLINE int
           ubf_eval_int_loop(int            i,
                             int            j,
                             int            p,
                             int            q,
                             int            i1,
                             int            j1,
                             int            p1,
                             int            q1,
                             short          si,
                             short          sj,
                             short          sp,
                             short          sq,
                             unsigned char  type,
                             unsigned char  type_2,
                             int            *rtype,
                             int            ij,
                             int            cp,
                             vrna_param_t   *P,
                             vrna_sc_t      *sc),
           "This function is obsolete");


DEPRECATED(PRIVATE INLINE int
           ubf_eval_int_loop2(int           i,
                              int           j,
                              int           p,
                              int           q,
                              int           i1,
                              int           j1,
                              int           p1,
                              int           q1,
                              short         si,
                              short         sj,
                              short         sp,
                              short         sq,
                              unsigned char type,
                              unsigned char type_2,
                              int           *rtype,
                              int           ij,
                              unsigned int  *sn,
                              unsigned int  *ss,
                              vrna_param_t  *P,
                              vrna_sc_t     *sc),
           "This function is obsolete");


/*
 *  ugly but fast exterior internal loop evaluation
 *
 *  Avoid including this function in your own code. It only serves
 *  as a fast inline block internally re-used throughout the RNAlib. It
 *  evalutes the free energy of internal loops in single sequences or sequence
 *  hybrids. Soft constraints are also applied if available.
 *
 *  NOTE: do not include into doxygen reference manual!
 */
DEPRECATED(PRIVATE INLINE int
           ubf_eval_ext_int_loop(int            i,
                                 int            j,
                                 int            p,
                                 int            q,
                                 int            i1,
                                 int            j1,
                                 int            p1,
                                 int            q1,
                                 short          si,
                                 short          sj,
                                 short          sp,
                                 short          sq,
                                 unsigned char  type,
                                 unsigned char  type_2,
                                 int            length,
                                 vrna_param_t   *P,
                                 vrna_sc_t      *sc),
           "This function is obsolete");


DEPRECATED(int
           vrna_eval_int_loop(vrna_fold_compound_t  *fc,
                              int                   i,
                              int                   j,
                              int                   k,
                              int                   l),
           "Use vrna_eval_internal() instead!");

DEPRECATED(int
           vrna_E_stack(vrna_fold_compound_t  *fc,
                        int                   i,
                        int                   j),
           "Use vrna_eval_stack() instead!");


DEPRECATED(FLT_OR_DBL
           vrna_exp_E_interior_loop(vrna_fold_compound_t *fc,
                                    int                  i,
                                    int                  j,
                                    int                  k,
                                    int                  l),
           "Use vrna_exp_eval_internal() instead!");


PRIVATE INLINE int
ubf_eval_int_loop(int           i,
                  int           j,
                  int           p,
                  int           q,
                  int           i1,
                  int           j1,
                  int           p1,
                  int           q1,
                  short         si,
                  short         sj,
                  short         sp,
                  short         sq,
                  unsigned char type,
                  unsigned char type_2,
                  int           *rtype,
                  int           ij,
                  int           cp,
                  vrna_param_t  *P,
                  vrna_sc_t     *sc)
{
  int energy, u1, u2;

  u1  = p1 - i;
  u2  = j1 - q;

  if ((cp < 0) || (ON_SAME_STRAND(i, p, cp) && ON_SAME_STRAND(q, j, cp))) {
    /* regular internal loop */
    energy = vrna_E_internal(u1, u2, type, type_2, si, sj, sp, sq, P);
  } else {
    /* internal loop like cofold structure */
    short Si, Sj;
    Si      = ON_SAME_STRAND(i, i1, cp) ? si : -1;
    Sj      = ON_SAME_STRAND(j1, j, cp) ? sj : -1;
    energy  = __E_IntLoop_Co(rtype[type], rtype[type_2],
                             i, j, p, q,
                             cp,
                             Si, Sj,
                             sp, sq,
                             P->model_details.dangles,
                             P);
  }

  /* add soft constraints */
  if (sc) {
    if (sc->energy_up)
      energy += sc->energy_up[i1][u1]
                + sc->energy_up[q1][u2];

    if (sc->energy_bp)
      energy += sc->energy_bp[ij];

    if (sc->energy_stack)
      if (u1 + u2 == 0) {
        int a = sc->energy_stack[i]
                + sc->energy_stack[p]
                + sc->energy_stack[q]
                + sc->energy_stack[j];
        energy += a;
      }

    if (sc->f)
      energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
  }

  return energy;
}


PRIVATE INLINE int
ubf_eval_int_loop2(int            i,
                   int            j,
                   int            p,
                   int            q,
                   int            i1,
                   int            j1,
                   int            p1,
                   int            q1,
                   short          si,
                   short          sj,
                   short          sp,
                   short          sq,
                   unsigned char  type,
                   unsigned char  type_2,
                   int            *rtype,
                   int            ij,
                   unsigned int   *sn,
                   unsigned int   *ss,
                   vrna_param_t   *P,
                   vrna_sc_t      *sc)
{
  int energy, u1, u2;

  u1  = p1 - i;
  u2  = j1 - q;

  if ((sn[i] == sn[p]) && (sn[q] == sn[j])) {
    /* regular internal loop */
    energy = vrna_E_internal(u1, u2, type, type_2, si, sj, sp, sq, P);
  } else {
    /* internal loop like cofold structure */
    short Si, Sj;
    Si      = (sn[i1] == sn[i]) ? si : -1;
    Sj      = (sn[j] == sn[j1]) ? sj : -1;
    energy  = __E_IntLoop_Co(rtype[type], rtype[type_2],
                             i, j, p, q,
                             ss[1],
                             Si, Sj,
                             sp, sq,
                             P->model_details.dangles,
                             P);
  }

  /* add soft constraints */
  if (sc) {
    if (sc->energy_up)
      energy += sc->energy_up[i1][u1]
                + sc->energy_up[q1][u2];

    if (sc->energy_bp)
      energy += sc->energy_bp[ij];

    if (sc->energy_stack)
      if (u1 + u2 == 0) {
        int a = sc->energy_stack[i]
                + sc->energy_stack[p]
                + sc->energy_stack[q]
                + sc->energy_stack[j];
        energy += a;
      }

    if (sc->f)
      energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
  }

  return energy;
}


PRIVATE INLINE int
ubf_eval_ext_int_loop(int           i,
                      int           j,
                      int           p,
                      int           q,
                      int           i1,
                      int           j1,
                      int           p1,
                      int           q1,
                      short         si,
                      short         sj,
                      short         sp,
                      short         sq,
                      unsigned char type,
                      unsigned char type_2,
                      int           length,
                      vrna_param_t  *P,
                      vrna_sc_t     *sc)
{
  int energy, u1, u2, u3;

  u1  = i1;
  u2  = p1 - j;
  u3  = length - q;

  energy = vrna_E_internal(u2, u1 + u3, type, type_2, si, sj, sp, sq, P);

  /* add soft constraints */
  if (sc) {
    if (sc->energy_up) {
      energy += sc->energy_up[j1][u2]
                + ((u3 > 0) ? sc->energy_up[q1][u3] : 0)
                + ((u1 > 0) ? sc->energy_up[1][u1] : 0);
    }

    if (sc->energy_stack)
      if (u1 + u2 + u3 == 0)
        energy += sc->energy_stack[i]
                  + sc->energy_stack[p]
                  + sc->energy_stack[q]
                  + sc->energy_stack[j];

    if (sc->f)
      energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
  }

  return energy;
}


PRIVATE INLINE int
E_IntLoop(int           n1,
          int           n2,
          int           type,
          int           type_2,
          int           si1,
          int           sj1,
          int           sp1,
          int           sq1,
          vrna_param_t  *P)
{
  /* compute energy of degree 2 loop (stack bulge or internal) */
  int nl, ns, u, energy, salt_stack_correction, salt_loop_correction, backbones;

  salt_stack_correction = P->SaltStack;
  salt_loop_correction  = 0;

  if (n1 > n2) {
    nl  = n1;
    ns  = n2;
  } else {
    nl  = n2;
    ns  = n1;
  }

  if (nl == 0) {
    return P->stack[type][type_2] + salt_stack_correction;  /* stack */
  }

  backbones = nl + ns + 2;

  if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
    /* salt correction for loop */
    if (backbones <= MAXLOOP + 1)
      salt_loop_correction = P->SaltLoop[backbones];
    else
      salt_loop_correction =
        vrna_salt_loop_int(backbones, P->model_details.salt, P->temperature + K0,
                           P->model_details.backbone_length);
  }

  if (ns == 0) {
    /* bulge */
    energy = (nl <= MAXLOOP) ? P->bulge[nl] :
             (P->bulge[30] + (int)(P->lxc * log(nl / 30.)));
    if (nl == 1) {
      energy += P->stack[type][type_2];
    } else {
      if (type > 2)
        energy += P->TerminalAU;

      if (type_2 > 2)
        energy += P->TerminalAU;
    }
  } else {
    /* internal loop */
    if (ns == 1) {
      if (nl == 1)                    /* 1x1 loop */
        return P->int11[type][type_2][si1][sj1] + salt_loop_correction;

      if (nl == 2) {
        /* 2x1 loop */
        if (n1 == 1)
          energy = P->int21[type][type_2][si1][sq1][sj1];
        else
          energy = P->int21[type_2][type][sq1][si1][sp1];

        return energy + salt_loop_correction;
      } else {
        /* 1xn loop */
        energy =
          (nl + 1 <=
           MAXLOOP) ? (P->internal_loop[nl + 1]) : (P->internal_loop[30] +
                                                    (int)(P->lxc * log((nl + 1) / 30.)));
        energy  += MIN2(MAX_NINIO, (nl - ns) * P->ninio[2]);
        energy  += P->mismatch1nI[type][si1][sj1] + P->mismatch1nI[type_2][sq1][sp1];
        return energy + salt_loop_correction;
      }
    } else if (ns == 2) {
      if (nl == 2) {
        /* 2x2 loop */
        return P->int22[type][type_2][si1][sp1][sq1][sj1] + salt_loop_correction;
      } else if (nl == 3) {
        /* 2x3 loop */
        energy  = P->internal_loop[5] + P->ninio[2];
        energy  += P->mismatch23I[type][si1][sj1] + P->mismatch23I[type_2][sq1][sp1];
        return energy + salt_loop_correction;
      }
    }

    {
      /* generic internal loop (no else here!)*/
      u       = nl + ns;
      energy  =
        (u <=
         MAXLOOP) ? (P->internal_loop[u]) : (P->internal_loop[30] + (int)(P->lxc * log((u) / 30.)));

      energy += MIN2(MAX_NINIO, (nl - ns) * P->ninio[2]);

      energy += P->mismatchI[type][si1][sj1] + P->mismatchI[type_2][sq1][sp1];
    }
  }

  return energy + salt_loop_correction;
}


PRIVATE INLINE FLT_OR_DBL
exp_E_IntLoop(int               u1,
              int               u2,
              int               type,
              int               type2,
              short             si1,
              short             sj1,
              short             sp1,
              short             sq1,
              vrna_exp_param_t  *P)
{
  int     ul, us, no_close = 0;
  double  z           = 0.;
  int     noGUclosure = P->model_details.noGUclosure;
  int     backbones;
  double  salt_stack_correction = P->expSaltStack;
  double  salt_loop_correction  = 1.;

  if ((noGUclosure) && ((type2 == 3) || (type2 == 4) || (type == 3) || (type == 4)))
    no_close = 1;

  if (u1 > u2) {
    ul  = u1;
    us  = u2;
  } else {
    ul  = u2;
    us  = u1;
  }

  /* salt correction for loop */
  backbones = ul + us + 2;

  if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
    if (backbones <= MAXLOOP + 1)
      salt_loop_correction = P->expSaltLoop[backbones];
    else
      salt_loop_correction =
        exp(-vrna_salt_loop_int(backbones, P->model_details.salt, P->temperature + K0,
                                P->model_details.backbone_length) * 10. / P->kT);
  }

  if (ul == 0) {
    /* stack */
    z = P->expstack[type][type2] * salt_stack_correction;
  } else if (!no_close) {
    if (us == 0) {
      /* bulge */
      z = P->expbulge[ul];
      if (ul == 1) {
        z *= P->expstack[type][type2];
      } else {
        if (type > 2)
          z *= P->expTermAU;

        if (type2 > 2)
          z *= P->expTermAU;
      }

      return (FLT_OR_DBL)(z * salt_loop_correction);
    } else if (us == 1) {
      if (ul == 1)                     /* 1x1 loop */
        return (FLT_OR_DBL)(P->expint11[type][type2][si1][sj1] * salt_loop_correction);

      if (ul == 2) {
        /* 2x1 loop */
        if (u1 == 1)
          return (FLT_OR_DBL)(P->expint21[type][type2][si1][sq1][sj1] * salt_loop_correction);
        else
          return (FLT_OR_DBL)(P->expint21[type2][type][sq1][si1][sp1] * salt_loop_correction);
      } else {
        /* 1xn loop */
        z = P->expinternal[ul + us] * P->expmismatch1nI[type][si1][sj1] *
            P->expmismatch1nI[type2][sq1][sp1];
        return (FLT_OR_DBL)(z * P->expninio[2][ul - us] * salt_loop_correction);
      }
    } else if (us == 2) {
      if (ul == 2) {
        /* 2x2 loop */
        return (FLT_OR_DBL)(P->expint22[type][type2][si1][sp1][sq1][sj1] * salt_loop_correction);
      } else if (ul == 3) {
        /* 2x3 loop */
        z = P->expinternal[5] * P->expmismatch23I[type][si1][sj1] *
            P->expmismatch23I[type2][sq1][sp1];
        return (FLT_OR_DBL)(z * P->expninio[2][1] * salt_loop_correction);
      }
    }

    /* generic internal loop (no else here!)*/
    z = P->expinternal[ul + us] * P->expmismatchI[type][si1][sj1] *
        P->expmismatchI[type2][sq1][sp1];
    return (FLT_OR_DBL)(z * P->expninio[2][ul - us] * salt_loop_correction);
  }

  return (FLT_OR_DBL)z;
}


PRIVATE INLINE int
E_IntLoop_Co(int          type,
             int          type_2,
             int          i,
             int          j,
             int          p,
             int          q,
             int          cutpoint,
             short        si1,
             short        sj1,
             short        sp1,
             short        sq1,
             int          dangles,
             vrna_param_t *P)
{
  return __E_IntLoop_Co(type,
                        type_2,
                        i,
                        j,
                        p,
                        q,
                        cutpoint,
                        si1,
                        sj1,
                        sp1,
                        sq1,
                        dangles,
                        P);
}


PRIVATE INLINE int
__E_IntLoop_Co(int          type,
               int          type_2,
               int          i,
               int          j,
               int          p,
               int          q,
               int          cutpoint,
               short        si1,
               short        sj1,
               short        sp1,
               short        sq1,
               int          dangles,
               vrna_param_t *P)
{
  int e, energy, ci, cj, cp, cq, d3, d5, d5_2, d3_2, tmm, tmm_2;
  int salt_loop_correction, backbones;

  salt_loop_correction = 0;

  backbones = p - i + j - q;
  /* salt correction for loop */
  if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
    if (backbones <= MAXLOOP + 1)
      salt_loop_correction = P->SaltLoop[backbones];
    else
      salt_loop_correction =
        vrna_salt_loop_int(backbones, P->model_details.salt, P->temperature + K0,
                           P->model_details.backbone_length);
  }

  energy = 0;
  if (type > 2)
    energy += P->TerminalAU;

  if (type_2 > 2)
    energy += P->TerminalAU;

  if (!dangles)
    return energy + salt_loop_correction;

  ci  = ON_SAME_STRAND(i, i + 1, cutpoint);
  cj  = ON_SAME_STRAND(j - 1, j, cutpoint);
  cp  = ON_SAME_STRAND(p - 1, p, cutpoint);
  cq  = ON_SAME_STRAND(q, q + 1, cutpoint);

  d3    = ci  ? P->dangle3[type][si1]   : 0;
  d5    = cj  ? P->dangle5[type][sj1]   : 0;
  d5_2  = cp  ? P->dangle5[type_2][sp1] : 0;
  d3_2  = cq  ? P->dangle3[type_2][sq1] : 0;

  tmm   = (cj && ci) ? P->mismatchExt[type][sj1][si1]   : d5 + d3;
  tmm_2 = (cp && cq) ? P->mismatchExt[type_2][sp1][sq1] : d5_2 + d3_2;

  if (dangles == 2)
    return energy + tmm + tmm_2 + salt_loop_correction;

  /* now we may have non-double dangles only */
  if (p - i > 2) {
    if (j - q > 2) {
      /* all degrees of freedom */
      e       = MIN2(tmm, d5);
      e       = MIN2(e, d3);
      energy  += e;
      e       = MIN2(tmm_2, d5_2);
      e       = MIN2(e, d3_2);
      energy  += e;
    } else if (j - q == 2) {
      /* all degrees of freedom in 5' part between i and p */
      e       = MIN2(tmm + d5_2, d3 + d5_2);
      e       = MIN2(e, d5 + d5_2);
      e       = MIN2(e, d3 + tmm_2);
      e       = MIN2(e, d3 + d3_2);
      e       = MIN2(e, tmm_2); /* no dangles on enclosing pair */
      e       = MIN2(e, d5_2);  /* no dangles on enclosing pair */
      e       = MIN2(e, d3_2);  /* no dangles on enclosing pair */
      energy  += e;
    } else {
      /* no unpaired base between q and j */
      energy += d3 + d5_2;
    }
  } else if (p - i == 2) {
    if (j - q > 2) {
      /* all degrees of freedom in 3' part between q and j */
      e       = MIN2(tmm + d3_2, d5 + d3_2);
      e       = MIN2(e, d5 + d3_2);
      e       = MIN2(e, d3 + d3_2);
      e       = MIN2(e, d5 + tmm_2);
      e       = MIN2(e, tmm_2);
      e       = MIN2(e, d5_2);
      e       = MIN2(e, d3_2);
      energy  += e;
    } else if (j - q == 2) {
      /* one possible dangling base between either side */
      e       = MIN2(tmm, tmm_2);
      e       = MIN2(e, d3);
      e       = MIN2(e, d5);
      e       = MIN2(e, d5_2);
      e       = MIN2(e, d3_2);
      e       = MIN2(e, d3 + d3_2);
      e       = MIN2(e, d5 + d5_2);
      energy  += e;
    } else {
      /* one unpaired base between i and p */
      energy += MIN2(d3, d5_2);
    }
  } else {
    /* no unpaired base between i and p */
    if (j - q > 2) {
      /* all degrees of freedom in 3' part between q and j */
      energy += d5 + d3_2;
    } else if (j - q == 2) {
      /* one unpaired base between q and j */
      energy += MIN2(d5, d3_2);
    }
  }

  return energy + salt_loop_correction;
}
/**
 * @}
 */

#endif

#endif
