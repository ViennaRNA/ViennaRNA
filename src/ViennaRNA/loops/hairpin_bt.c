#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/loops/external.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/loops/hairpin.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "hairpin_hc.inc"

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

/**
 *  @brief Backtrack a hairpin loop closed by @f$ (i,j) @f$
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 */
PUBLIC int
vrna_BT_hp_loop(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   en,
                vrna_bp_stack_t       *bp_stack,
                int                   *stack_count)
{
  int       e, u;
  vrna_sc_t *sc;

  sc = NULL;

  u = j - i - 1;

  if (fc->hc->up_hp[i + 1] < u)
    return 0;

  e = vrna_E_hp_loop(fc, i, j);

  if (e == en) {
    switch (fc->type) {
      case  VRNA_FC_TYPE_SINGLE:
        sc = fc->sc;
        break;

      case  VRNA_FC_TYPE_COMPARATIVE:
        if (fc->scs)
          sc = fc->scs[0];

        break;

      default:
        break;
    }

    if (sc) {
      if (sc->bt) {
        vrna_basepair_t *ptr, *aux_bps;
        aux_bps = sc->bt(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
        for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
          bp_stack[++(*stack_count)].i  = ptr->i;
          bp_stack[(*stack_count)].j    = ptr->j;
        }
        free(aux_bps);
      }
    }

    return 1;
  }

  return 0;
}
