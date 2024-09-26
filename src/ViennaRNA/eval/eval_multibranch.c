/*
 * WBL 24 Aug 2018 Add AVX512 based on sources_034_578/modular_decomposition_id3.c
 * WBL 22 Aug 2018 by hand d3c17fd3e04e2419c147a1e097d3c4d2c5a6f11d lines 1355-1357
 */

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
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/multibranch.h"

#include "ViennaRNA/intern/grammar_dat.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/multibranch_hc.inc"
#include "ViennaRNA/constraints/multibranch_sc.inc"

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
PUBLIC int
vrna_E_multibranch_stem(unsigned int  type,
                        int           si1,
                        int           sj1,
                        vrna_param_t  *P)
{
  int energy = INF;

  if (P) {
    energy = P->MLintern[type];

    if (si1 >= 0 && sj1 >= 0)
      energy += P->mismatchM[type][si1][sj1];
    else if (si1 >= 0)
      energy += P->dangle5[type][si1];
    else if (sj1 >= 0)
      energy += P->dangle3[type][sj1];

    if (type > 2)
      energy += P->TerminalAU;
  }

  return energy;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
