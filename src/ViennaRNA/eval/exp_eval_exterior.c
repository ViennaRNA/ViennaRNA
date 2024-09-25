#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/params/default.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/eval/exterior.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC FLT_OR_DBL
vrna_exp_E_exterior_stem(unsigned int     type,
                         int              n5d,
                         int              n3d,
                         vrna_exp_param_t *p)
{
  double energy = 1.0;

  if (n5d >= 0 && n3d >= 0)
    energy = p->expmismatchExt[type][n5d][n3d];
  else if (n5d >= 0)
    energy = p->expdangle5[type][n5d];
  else if (n3d >= 0)
    energy = p->expdangle3[type][n3d];

  if (type > 2)
    energy *= p->expTermAU;

  return (FLT_OR_DBL)energy;
}


PUBLIC FLT_OR_DBL
vrna_exp_E_exterior_loop(unsigned int n,
                         vrna_md_t    *md)
{
  vrna_md_t md_tmp;

  if (md == NULL) {
    vrna_md_set_default(&md_tmp);
    md = &md_tmp;
  }

  if ((md->circ) &&
      (md->circ_penalty)) {
    double  e   = (double)vrna_E_exterior_loop(n, md) / 100.;
    double  kT  = md->betaScale * (md->temperature + K0) * GASCONST / 1000.;  /* kT in kcal/mol */

    return (FLT_OR_DBL)exp(-e / kT);                                          /* return in dekacal/mol */
  } else {
    return 1.;
  }
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */


/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC FLT_OR_DBL
exp_E_Stem(int              type,
           int              si1,
           int              sj1,
           int              extLoop,
           vrna_exp_param_t *P)
{
  double  energy  = 1.0;
  double  d5      = (si1 >= 0) ? P->expdangle5[type][si1] : 1.;
  double  d3      = (sj1 >= 0) ? P->expdangle3[type][sj1] : 1.;

  if (si1 >= 0 && sj1 >= 0)
    energy = (extLoop) ? P->expmismatchExt[type][si1][sj1] : P->expmismatchM[type][si1][sj1];
  else
    energy = d5 * d3;

  if (type > 2)
    energy *= P->expTermAU;

  if (!extLoop)
    energy *= P->expMLintern[type];

  return (FLT_OR_DBL)energy;
}


PUBLIC FLT_OR_DBL
exp_E_ExtLoop(int               type,
              int               si1,
              int               sj1,
              vrna_exp_param_t  *P)
{
  return vrna_exp_E_exterior_stem((unsigned int)type, si1, sj1, P);
}


PUBLIC FLT_OR_DBL
vrna_exp_E_ext_stem(unsigned int      type,
                    int               n5d,
                    int               n3d,
                    vrna_exp_param_t  *p)
{
  return vrna_exp_E_exterior_stem(type, n5d, n3d, p);
}


#endif
