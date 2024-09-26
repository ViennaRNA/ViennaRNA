#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/eval/multibranch.h"

#include "ViennaRNA/intern/grammar_dat.h"

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
vrna_exp_E_multibranch_stem(unsigned int      type,
                            int               si1,
                            int               sj1,
                            vrna_exp_param_t  *P)
{
  double energy;

  energy = P->expMLintern[type];

  if (si1 >= 0 && sj1 >= 0)
    energy *= P->expmismatchM[type][si1][sj1];
  else if (si1 >= 0)
    energy *= P->expdangle5[type][si1];
  else if (sj1 >= 0)
    energy *= P->expdangle3[type][sj1];

  if (type > 2)
    energy *= P->expTermAU;

  return (FLT_OR_DBL)energy;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
